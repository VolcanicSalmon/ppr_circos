import os
import json
import pickle
import numpy as np
import faiss
import torch
import torch.nn as nn
from torch.nn import functional as F
from torch_geometric.data import HeteroData
from torch_geometric.nn import HeteroConv, GATConv, Linear
from torch_geometric.utils import from_networkx
import networkx as nx
import pandas as pd
from typing import List, Dict, Tuple, Optional
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment



# ── Geometry helpers ──────────────────────────────────────────────────────────

def calc_dist(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


def calc_angle(vec1, vec2):
    cos = vec1.dot(vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    return np.arccos(np.clip(cos, -1.0, 1.0))


# ── Feature encoding ──────────────────────────────────────────────────────────

BASE_MAP = {'A': 0, 'U': 1, 'G': 2, 'C': 3}
DICER_MAP = {21: 0, 22: 1, 23: 2, 24: 3}
PAIR_MAP  = {'WC': 0, 'wobble': 1, 'mismatch': 2, 'unpaired': 3,
             'BULt': 4, 'BULq': 5, 'AILt': 6, 'AILq': 7, 'SIL': 8, 'UP5': 9, 'UP3': 10}


def encode_base(b: str) -> torch.Tensor:
    t = torch.zeros(4)
    t[BASE_MAP.get(b.upper(), 0)] = 1.0
    return t


def encode_pairing(b1: str, b2: str, annot: str = '') -> torch.Tensor:
    t = torch.zeros(len(PAIR_MAP))
    if annot in PAIR_MAP:
        t[PAIR_MAP[annot]] = 1.0
        return t
    wc     = {('A','U'),('U','A'),('G','C'),('C','G')}
    wobble = {('G','U'),('U','G')}
    pair   = (b1.upper(), b2.upper())
    if pair in wc:
        t[PAIR_MAP['WC']] = 1.0
    elif pair in wobble:
        t[PAIR_MAP['wobble']] = 1.0
    elif b2 != '-':
        t[PAIR_MAP['mismatch']] = 1.0
    else:
        t[PAIR_MAP['unpaired']] = 1.0
    return t


def encode_srna_node(seq: str, mfe: float, stem_len: int,
                     hloop_size: int, size_nt: int) -> torch.Tensor:
    # one-hot bases (padded to 23) + scalar features
    base_feats = torch.zeros(23, 4)
    for i, b in enumerate(seq[:23]):
        base_feats[i] = encode_base(b)
    base_feats = base_feats.flatten()                    # 92
    size_oh = torch.zeros(4)                             # 21/22/23/24 one-hot
    idx = size_nt - 21
    if 0 <= idx < 4:
        size_oh[idx] = 1.0
    scalars = torch.tensor([mfe, float(stem_len), float(hloop_size),
                            float(size_nt)])             # length as continuous feature
    return torch.cat([base_feats, size_oh, scalars])     # dim = 100


def encode_mrna_node(dicer: int, total_reads: float,
                     pt_score: float, ax_score: float,
                     pt_ratio: float, pt_num: int, pt_abun: float,
                     deg_cent: float, bet_cent: float) -> torch.Tensor:
    dicer_oh = torch.zeros(4)
    dicer_oh[DICER_MAP.get(dicer, 3)] = 1.0
    scalars = torch.tensor([
        np.log1p(total_reads),
        pt_score  if not np.isnan(pt_score)        else 0.0,
        ax_score  if not np.isnan(ax_score)        else 0.0,
        float(not np.isnan(pt_score)),
        float(not np.isnan(ax_score)),
        pt_ratio  if not np.isnan(pt_ratio)        else 0.0,
        float(pt_num) if not np.isnan(float(pt_num)) else 0.0,
        np.log1p(pt_abun) if not np.isnan(pt_abun) else 0.0,
        deg_cent, bet_cent,
    ])
    return torch.cat([dicer_oh, scalars])  # dim = 14


def encode_edge(allen_score: float, mfe_site: float,
                mfe_ratio: float, paired_frac: float,
                pairing_annots: List[str]) -> torch.Tensor:
    """
    Edge features for one sRNA-mRNA interaction from GSTAr.
    pairing_annots: list of annotation strings per paired block (WC/wobble/BULt etc.)
    """
    pair_oh = torch.zeros(len(PAIR_MAP))
    for annot in pairing_annots:
        idx = PAIR_MAP.get(annot)
        if idx is not None:
            pair_oh[idx] = 1.0
    scalars = torch.tensor([
        float(allen_score),
        float(mfe_site),
        float(mfe_ratio),
        float(paired_frac),
    ])
    return torch.cat([pair_oh, scalars])  # dim = 11


# ── dot-bracket → base-pair graph ────────────────────────────────────────────

def dotbracket_to_graph(sequence: str, structure: str) -> Tuple[torch.Tensor, torch.Tensor]:
    edges = [(i, i+1) for i in range(len(sequence)-1)]
    stack = []
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            edges += [(i, j), (j, i)]
    x = torch.stack([encode_base(b) for b in sequence])
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    return x, edge_index


# ── RNA structure GNN (per-RNA embedding) ────────────────────────────────────

class RNAStructGNN(nn.Module):
    def __init__(self, in_dim: int = 4, hidden_dim: int = 64,
                 out_dim: int = 32, heads: int = 4):
        super().__init__()
        self.conv1 = GATConv(in_dim, hidden_dim, heads=heads, concat=True)
        self.conv2 = GATConv(hidden_dim * heads, out_dim, heads=1, concat=False)

    def forward(self, x, edge_index):
        x = F.gelu(self.conv1(x, edge_index))
        return self.conv2(x, edge_index)   # [n_nodes, out_dim]

    def embed(self, x, edge_index):
        node_emb = self.forward(x, edge_index)
        return node_emb.mean(dim=0)        # graph-level mean pool → [out_dim]


# ── Heterogeneous interaction graph ──────────────────────────────────────────

class SRNAMRNAHeteroGNN(nn.Module):
    """
    Heterogeneous GNN over sRNA-mRNA interaction graph.
    Node types  : 'srna', 'mrna'
    Edge type   : ('srna', 'targets', 'mrna')
    Edge features: pairing type one-hot (7-dim) from GSTAr alignment
    """
    def __init__(self, srna_in: int, mrna_in: int,
                 hidden: int = 64, out: int = 32, heads: int = 4, dropout: float = 0.3):
        super().__init__()
        proj_dim = hidden * heads
        self.drop       = dropout
        self.srna_proj  = Linear(srna_in, proj_dim)
        self.mrna_proj  = Linear(mrna_in, proj_dim)
        self.convs = nn.ModuleList([
            HeteroConv({
                ('srna', 'targets', 'mrna'): GATConv(
                    (proj_dim, proj_dim), hidden, heads=heads,
                    concat=True, add_self_loops=False, edge_dim=len(PAIR_MAP) + 4
                ),
            }, aggr='sum')
            for _ in range(2)
        ])
        self.out_srna   = Linear(hidden * heads, out)
        self.out_mrna   = Linear(hidden * heads, out)

    def forward(self, data: HeteroData):
        x_dict = {
            'srna': F.dropout(F.gelu(self.srna_proj(data['srna'].x)), p=self.drop, training=self.training),
            'mrna': F.dropout(F.gelu(self.mrna_proj(data['mrna'].x)), p=self.drop, training=self.training),
        }
        edge_index = data['srna', 'targets', 'mrna'].edge_index
        edge_attr  = data['srna', 'targets', 'mrna'].edge_attr

        for conv in self.convs:
            out = conv(x_dict, {('srna','targets','mrna'): edge_index},
                       edge_attr_dict={('srna','targets','mrna'): edge_attr})
            # merge: srna is never a destination so it won't appear in out;
            # keep it from the previous x_dict so the next layer has a source
            x_dict = {**x_dict, **{k: F.dropout(F.gelu(v), p=self.drop, training=self.training)
                                    for k, v in out.items()}}

        return self.out_srna(x_dict['srna']), self.out_mrna(x_dict['mrna'])


# ── Build HeteroData from pipeline outputs ───────────────────────────────────

def build_heterodata(srna_records: List[Dict],
                     mrna_records: List[Dict],
                     edges: List[Tuple[int, int, torch.Tensor]],
                     struct_gnn: Optional[RNAStructGNN] = None) -> HeteroData:
    """
    srna_records: list of dicts with keys:
        seq, structure, mfe, stem_len, hloop_size, size_nt
    mrna_records: list of dicts with keys:
        dicer, total_reads, pt_score, ax_score, pt_ratio,
        pt_num, pt_abun, deg_cent, bet_cent
    edges: list of (srna_idx, mrna_idx, pairing_feat_tensor)
    struct_gnn: optional; if None a zero vector of dim 32 is used in place
                of the structure embedding (useful before RNAfold is available)
    """
    data = HeteroData()

    # sRNA nodes: struct GNN embedding + sequence/structure scalars
    srna_feats = []
    for r in srna_records:
        if struct_gnn is not None:
            x, ei      = dotbracket_to_graph(r['seq'], r['structure'])
            struct_emb = struct_gnn.embed(x, ei)       # [32]
        else:
            struct_emb = torch.zeros(32)
        node_feat  = encode_srna_node(
            r['seq'], r['mfe'], r['stem_len'], r['hloop_size'], r['size_nt']
        )                                               # [99]
        srna_feats.append(torch.cat([struct_emb, node_feat]))  # [132]
    data['srna'].x = torch.stack(srna_feats).float()

    # mRNA nodes
    mrna_feats = []
    for r in mrna_records:
        mrna_feats.append(encode_mrna_node(
            r['dicer'], r['total_reads'],
            r.get('pt_score', float('nan')), r.get('ax_score', float('nan')),
            r.get('pt_ratio', float('nan')), r.get('pt_num', float('nan')),
            r.get('pt_abun', float('nan')),
            r.get('deg_cent', 0.0), r.get('bet_cent', 0.0),
        ))
    data['mrna'].x = torch.stack(mrna_feats).float()

    # edges with pairing features
    src  = torch.tensor([e[0] for e in edges], dtype=torch.long)
    dst  = torch.tensor([e[1] for e in edges], dtype=torch.long)
    eatt = torch.stack([e[2] for e in edges])
    data['srna', 'targets', 'mrna'].edge_index = torch.stack([src, dst])
    data['srna', 'targets', 'mrna'].edge_attr  = eatt

    return data


def _infer_pairing_annots(srna_seq: str, mrna_seq: str) -> List[str]:
    """
    Walk aligned sRNA and mRNA sequences and collect which pairing types
    are present (multi-hot input for encode_edge).
    """
    wc     = {('A','U'),('U','A'),('G','C'),('C','G')}
    wobble = {('G','U'),('U','G')}
    annots: set = set()
    for s, m in zip(srna_seq.upper(), mrna_seq.upper()):
        s = s.replace('T', 'U')
        m = m.replace('T', 'U')
        if s == '-':
            annots.add('BULt')      # gap in sRNA → bulge in target (mRNA)
        elif m == '-':
            annots.add('BULq')      # gap in mRNA → bulge in query (sRNA)
        elif (s, m) in wc:
            annots.add('WC')
        elif (s, m) in wobble:
            annots.add('wobble')
        else:
            annots.add('mismatch')
    return list(annots) if annots else ['unpaired']


def heterodata_from_gstar(
        gstar_file: str,
        di_npy: str = None,
        struct_gnn: Optional[RNAStructGNN] = None,
        allen_max: int = 3,
        min_paired: float = 0.75,
        scaler=None) -> HeteroData:
    """
    Build a HeteroData graph directly from a GSTAr output file.

    sRNA nodes  : one per unique Query ID.
    mRNA nodes  : one per unique Transcript ID. Phase-pipeline scores
                  (pt_score, ax_score, …) are left as NaN placeholders
                  until merged in from map_clusters_to_tx output.
    Edges       : one per filtered GSTAr row.
    data['srna','targets','mrna'].di
                : cross-block DI tensor [l_srna_pad, l_mrna_pad] stored as a
                  graph-level attribute when di_npy is provided; the GNN does
                  not consume it automatically — pass it to whatever module
                  needs positional coevolution signal.

    Parameters
    ----------
    gstar_file  : GSTAr TSV (13/15/16/19-col formats all accepted)
    di_npy      : .npy saved by PottsDCA.fit_evcouplings (full L×L matrix)
    struct_gnn  : RNAStructGNN instance; if None uses zero vectors (dim 32)
    allen_max   : keep rows with AllenScore <= allen_max
    min_paired  : keep rows with paired_frac >= min_paired
    scaler      : fitted EdgeScaler (z-scores AllenScore and MFEsite)
    """
    from util import read_gstar, _parse_paired_frac

    df = read_gstar(gstar_file)
    df['AllenScore'] = pd.to_numeric(df['AllenScore'], errors='coerce')
    df['MFEsite']    = pd.to_numeric(df['MFEsite'],    errors='coerce')
    df['MFEratio']   = pd.to_numeric(df['MFEratio'],   errors='coerce')

    df = df[df['AllenScore'] <= allen_max].copy()
    df['paired_frac'] = df['Paired'].apply(_parse_paired_frac)
    df = df[df['paired_frac'] >= min_paired].copy()

    split = df['Sequence'].str.split('&', n=1, expand=True)
    df['srna_seq'] = (split[0].str.strip().str.upper()
                               .str.replace('T', 'U', regex=False))
    df['mrna_seq'] = (split[1].str.strip().str.upper()
                               .str.replace('T', 'U', regex=False))
    df['srna_len'] = df['srna_seq'].str.replace('-', '', regex=False).str.len()
    df = df[df['srna_len'].isin([21, 22, 23, 24])].copy()

    if df.empty:
        raise ValueError('No interactions passed filters — check allen_max / min_paired')

    if scaler is not None:
        df = scaler.transform(df)

    srna_ids = {s: i for i, s in enumerate(df['Query'].unique())}
    mrna_ids = {t: i for i, t in enumerate(df['Transcript'].unique())}

    # sRNA node records (one per unique Query)
    srna_rows = df.drop_duplicates('Query').set_index('Query')
    srna_records = []
    for sid in srna_ids:
        row = srna_rows.loc[sid]
        seq = row['srna_seq']
        struct_str  = str(row['Structure'])
        srna_struct = struct_str.split('&')[0] if '&' in struct_str else struct_str
        srna_struct = srna_struct.ljust(len(seq), '.')   # pad to seq length
        srna_records.append({
            'seq':        seq,
            'structure':  srna_struct,
            'mfe':        float(row['MFEsite']) if pd.notna(row['MFEsite']) else 0.0,
            'stem_len':   0,    # placeholder until RNAfold integration
            'hloop_size': 0,
            'size_nt':    int(row['srna_len']),
        })

    # mRNA node records — phase scores and ShortStack read counts to be merged
    # from run_phase_pipeline output via merge_phase_scores(); zeros for now
    mrna_records = [
        {'dicer': 24, 'total_reads': 0.0,
         'pt_score': float('nan'), 'ax_score': float('nan'),
         'pt_ratio': float('nan'), 'pt_num':   float('nan'),
         'pt_abun':  float('nan'), 'deg_cent': 0.0, 'bet_cent': 0.0}
        for _ in mrna_ids
    ]

    # edges
    edges = []
    for _, row in df.iterrows():
        si    = srna_ids[row['Query']]
        mi    = mrna_ids[row['Transcript']]
        annots = _infer_pairing_annots(row['srna_seq'], row['mrna_seq'])
        feat  = encode_edge(
            float(row['AllenScore']),
            float(row['MFEsite'])  if pd.notna(row['MFEsite'])  else 0.0,
            float(row['MFEratio']) if pd.notna(row['MFEratio']) else 0.0,
            float(row['paired_frac']),
            annots,
        )
        edges.append((si, mi, feat))

    data = build_heterodata(srna_records, mrna_records, edges, struct_gnn)
    data['mrna'].node_ids = list(mrna_ids.keys())   # transcript ID → node index
    data['srna'].node_ids = list(srna_ids.keys())   # sRNA query ID → node index

    # store DI as graph-level attribute on the edge store
    if di_npy is not None:
        full_di = np.load(di_npy)                       # [L, L]
        l_srna  = int(df['srna_seq'].str.len().max())
        l_mrna  = int(df['mrna_seq'].str.len().max())
        cross   = full_di[:l_srna, l_srna:l_srna + l_mrna]  # [l_srna, l_mrna]
        data['srna', 'targets', 'mrna'].di = torch.tensor(cross, dtype=torch.float32)

    return data


def merge_phase_scores(data: HeteroData,
                       clusters_bed: str) -> HeteroData:
    """
    Fill mRNA node features from clusters_on_tx.bed alone.
    All needed values (PT, AX, DC, MIR, reads) are already encoded
    in the BED name field by map_clusters_to_tx.py.

    Feature vector (dim=9) per transcript:
        log1p(total_reads), log1p(n_clusters),
        max_pt, max_ax, mean_ax,
        frac_phased_pt, frac_phased_ax,
        dominant_dc (21/22/23/24, 0=N/.),
        has_mir (1 if any cluster is .mature or MIR!=N/.)
    """
    from collections import defaultdict, Counter

    def _tag(parts, key):
        for p in parts:
            if p.startswith(key + '='):
                return p.split('=', 1)[1]
        return ''

    tx_clusters: dict = defaultdict(list)
    with open(clusters_bed) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 6:
                continue
            parts = f[3].split('|')
            if len(parts) < 2:
                continue
            tx_id      = parts[0]
            cluster_id = parts[1]
            tags       = parts[2:]
            dc         = _tag(tags, 'DC')
            mir        = _tag(tags, 'MIR')
            pt         = _tag(tags, 'PT')
            ax         = _tag(tags, 'AX')
            reads      = float(f[4]) if f[4].lstrip('-').isdigit() else 0.0
            tx_clusters[tx_id].append((cluster_id, reads, dc, mir, pt, ax))

    node_ids = data['mrna'].node_ids
    mrna_feats = []
    zero_feat = torch.zeros(9, dtype=torch.float32)

    for tx in node_ids:
        entries = tx_clusters.get(tx, [])
        if not entries:
            mrna_feats.append(zero_feat)
            continue

        total_reads = sum(r for _, r, *_ in entries)
        n           = len(entries)
        pt_vals, ax_vals, dc_vals = [], [], []
        has_mir = 0.0

        for cid, _, dc, mir, pt, ax in entries:
            try:    pt_vals.append(float(pt))
            except (ValueError, TypeError): pass
            try:    ax_vals.append(float(ax))
            except (ValueError, TypeError): pass
            dc_vals.append(int(dc) if dc.isdigit() else 0)
            if mir not in ('N', 'NA', '.', '') or '.mature' in cid:
                has_mir = 1.0

        dominant_dc  = Counter(dc_vals).most_common(1)[0][0]
        frac_pt      = len(pt_vals) / n
        frac_ax      = len(ax_vals) / n

        feat = torch.tensor([
            np.log1p(total_reads),
            np.log1p(n),
            max(pt_vals) if pt_vals else 0.0,
            max(ax_vals) if ax_vals else 0.0,
            float(np.mean(ax_vals)) if ax_vals else 0.0,
            frac_pt,
            frac_ax,
            float(dominant_dc),
            has_mir,
        ], dtype=torch.float32)
        mrna_feats.append(feat)

    data['mrna'].x = torch.stack(mrna_feats)
    return data


# ── Network metrics (SRNET) ───────────────────────────────────────────────────

class SRNET:
    def __init__(self, configs: str):
        with open(configs, 'r') as f:
            self.config = json.load(f)

    def build_srna_net(self, targets: List[str], output: str):
        rna_edges = self.get_string_targets(targets)
        G = nx.Graph()
        G.add_edges_from(rna_edges)
        self.calc_net_metrics(G)
        return G

    def calc_net_metrics(self, G: nx.Graph):
        deg_cent  = nx.degree_centrality(G)
        bet_cent  = nx.betweenness_centrality(G)
        clos_cent = nx.closeness_centrality(G)
        try:
            eig_cent = nx.eigenvector_centrality(G, max_iter=500)
        except nx.PowerIterationFailedConvergence:
            eig_cent = {n: 0.0 for n in G.nodes()}
        nx.set_node_attributes(G, deg_cent,  'deg_cent')
        nx.set_node_attributes(G, bet_cent,  'bet_cent')
        nx.set_node_attributes(G, clos_cent, 'clos_cent')
        nx.set_node_attributes(G, eig_cent,  'eig_cent')

    def key_node(self, net: nx.Graph, outfile: str) -> pd.DataFrame:
        if net.number_of_nodes() == 0:
            df = pd.DataFrame(columns=['node','degree','bet_cent','clos_cent','eig_cent'])
            df.to_csv(outfile, index=False)
            return df
        rows = []
        for n, d in net.nodes(data=True):
            rows.append({
                'node':      n,
                'degree':    net.degree(n),
                'bet_cent':  d.get('bet_cent', 0.0),
                'clos_cent': d.get('clos_cent', 0.0),
                'eig_cent':  d.get('eig_cent', 0.0),
            })
        df = pd.DataFrame(rows).sort_values('bet_cent', ascending=False)
        df.to_csv(outfile, index=False)
        return df

    def dotbracket_to_graph(self, sequence: str, structure: str):
        return dotbracket_to_graph(sequence, structure)

    def to_pyg(self, G: nx.Graph,
               node_attrs: List[str] = ['deg_cent','bet_cent','clos_cent','eig_cent']):
        return from_networkx(G, group_node_attrs=node_attrs)


# ── Geometric kernel GP for cross-species prediction ─────────────────────────

class CrossSpeciesKernelGP:
    """
    Uses Matérn geometric kernel on graph space to compute similarity
    between species-level sRNA-mRNA interaction graphs, then applies
    kernel ridge regression to predict edge probabilities in target species.
    """
    def __init__(self, nu: float = 2.5, ridge: float = 1e-3):
        self.nu    = nu
        self.ridge = ridge
        self.alpha = None
        self.X_train = None

    def _build_kernel(self, G: nx.Graph):
        from geometric_kernels.spaces import Graph
        from geometric_kernels.kernels import MaternGeometricKernel
        space  = Graph(nx.to_numpy_array(G))
        return MaternGeometricKernel(space)

    def fit(self, graphs: List[nx.Graph], y: np.ndarray):
        n = len(graphs)
        K = np.zeros((n, n))
        kernels = [self._build_kernel(G) for G in graphs]
        for i in range(n):
            for j in range(i, n):
                # node embeddings as kernel inputs
                Xi = np.arange(graphs[i].number_of_nodes())[:, None]
                Xj = np.arange(graphs[j].number_of_nodes())[:, None]
                K[i, j] = kernels[i].K(Xi, Xj).sum()
                K[j, i] = K[i, j]
        self.alpha   = np.linalg.solve(K + self.ridge * np.eye(n), y)
        self.X_train = graphs
        self.kernels = kernels

    def predict(self, graphs: List[nx.Graph]) -> np.ndarray:
        n_train = len(self.X_train)
        n_test  = len(graphs)
        K_star  = np.zeros((n_test, n_train))
        for i, Gt in enumerate(graphs):
            kt = self._build_kernel(Gt)
            for j, Gs in enumerate(self.X_train):
                Xi = np.arange(Gt.number_of_nodes())[:, None]
                Xj = np.arange(Gs.number_of_nodes())[:, None]
                K_star[i, j] = kt.K(Xi, Xj).sum()
        return K_star @ self.alpha


# ── Potts DCA for sRNA-mRNA coevolution ──────────────────────────────────────

class PottsDCA:
    """
    Direct Coupling Analysis (Potts model) on concatenated sRNA|mRNA MSA.

    Input : ortholog sequences per species — one sRNA and one mRNA per row.
            Species order must match across sRNA and mRNA dicts.

    Output: DI matrix cross-block [L_srna, L_mrna] — direct information
            between each sRNA position and each mRNA position.
            Use as edge features between sRNA codon nodes and mRNA codon nodes.

    Example
    -------
    dca = PottsDCA(l_srna=21, l_mrna=1000)
    dca.build_msa(srna_seqs, mrna_seqs, 'concat_msa.fasta')
    DI = dca.fit('concat_msa.fasta')          # [21, 1000]
    edge_feat = dca.edge_features(DI, srna_pos=5, mrna_pos=42)  # scalar
    """

    RNA_ALPHA = ['A', 'U', 'G', 'C', '-']
    ALPHA_IDX = {c: i for i, c in enumerate(RNA_ALPHA)}

    def __init__(self, l_srna: int, l_mrna: int, q: int = 5,
                 lam: float = 0.01, max_iter: int = 1000):
        self.l_srna   = l_srna
        self.l_mrna   = l_mrna
        self.L        = l_srna + l_mrna
        self.q        = q        # alphabet size (A/U/G/C/-)
        self.lam      = lam      # L2 regularisation
        self.max_iter = max_iter
        self.J        = None     # [L, L, q, q] coupling parameters
        self.h        = None     # [L, q]  local fields
        self.DI       = None     # [L, L]  direct information

    # ── MSA construction ──────────────────────────────────────────────────────

    def build_msa(self, srna_seqs: Dict[str, str],
                  mrna_seqs: Dict[str, str],
                  out_fasta: str):
        """
        Concatenate sRNA and mRNA orthologs per species into a joint MSA.
        srna_seqs / mrna_seqs: {species_id: aligned_sequence}
        Species must be aligned separately before calling (e.g. MAFFT).
        """
        species = sorted(set(srna_seqs) & set(mrna_seqs))
        if not species:
            raise ValueError('No common species between sRNA and mRNA dicts')
        with open(out_fasta, 'w') as f:
            for sp in species:
                s = srna_seqs[sp].upper().replace('T', 'U')
                m = mrna_seqs[sp].upper().replace('T', 'U')
                concat = s + m
                f.write(f'>{sp}\n{concat}\n')
        self.l_srna = len(next(iter(srna_seqs.values())))
        self.l_mrna = len(next(iter(mrna_seqs.values())))
        self.L      = self.l_srna + self.l_mrna

    # ── Encoding ──────────────────────────────────────────────────────────────

    def _encode_msa(self, msa_file: str) -> np.ndarray:
        """Returns [N_seq, L] integer-encoded MSA."""
        seqs = []
        for rec in SeqIO.parse(msa_file, 'fasta'):
            row = [self.ALPHA_IDX.get(c, 4) for c in str(rec.seq).upper()]
            seqs.append(row)
        return np.array(seqs, dtype=np.int32)

    # ── Pseudocount reweighting ───────────────────────────────────────────────

    @staticmethod
    def _reweight(msa: np.ndarray, theta: float = 0.2) -> np.ndarray:
        """Sequence reweighting to correct for phylogenetic bias."""
        N, L = msa.shape
        sim  = np.mean(msa[:, None] == msa[None, :], axis=2)  # [N, N]
        w    = 1.0 / (sim >= (1 - theta)).sum(axis=1)
        return w / w.sum()

    # ── Gaussian (mean-field) DCA — fast approximation ───────────────────────

    def _gaussian_dca(self, msa: np.ndarray,
                      weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Gaussian (mean-field) approximation to plmDCA.
        Returns (J [L,L,q,q], C_inv [Lq, Lq]).
        Fast O(L^2 q^2) — use plmDCA for higher accuracy.
        """
        N, L = msa.shape
        q    = self.q
        Lq   = L * q

        # one-hot encode
        X = np.zeros((N, Lq))
        for n in range(N):
            for i in range(L):
                a = msa[n, i]
                if a < q:
                    X[n, i * q + a] = weights[n]

        # weighted covariance
        mu  = X.sum(axis=0)
        C   = X.T @ X - np.outer(mu, mu)
        C  += self.lam * np.eye(Lq)
        C_inv = np.linalg.inv(C)

        J = np.zeros((L, L, q, q))
        for i in range(L):
            for j in range(L):
                if i != j:
                    J[i, j] = -C_inv[i*q:(i+1)*q, j*q:(j+1)*q]
        return J, C_inv

    # ── Direct Information ────────────────────────────────────────────────────

    def _compute_di(self, J: np.ndarray,
                    freq1: np.ndarray) -> np.ndarray:
        """
        Compute DI[i,j] from coupling matrix J and single-site frequencies.
        freq1: [L, q]
        """
        L, q = freq1.shape
        DI   = np.zeros((L, L))
        for i in range(L):
            for j in range(i+1, L):
                # direct probability via J_ij and marginals
                W = np.exp(J[i, j])
                W = W * freq1[i, :, None] * freq1[j, None, :]
                Z = W.sum()
                if Z == 0:
                    continue
                P  = W / Z
                pi = freq1[i]
                pj = freq1[j]
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.where((pi[:, None] * pj[None, :]) > 0,
                                     P / (pi[:, None] * pj[None, :]), 0)
                    di    = np.where(P > 0, P * np.log(ratio), 0).sum()
                DI[i, j] = di
                DI[j, i] = di
        return DI

    # ── Public API ────────────────────────────────────────────────────────────

    def fit(self, msa_file: str, theta: float = 0.2) -> np.ndarray:
        """
        Fit Potts model and return cross-block DI [l_srna, l_mrna].
        Uses Gaussian mean-field DCA.
        For higher accuracy replace _gaussian_dca with pydca.plmdca.
        """
        msa     = self._encode_msa(msa_file)
        weights = self._reweight(msa, theta)

        freq1 = np.zeros((self.L, self.q))
        for i in range(self.L):
            for a in range(self.q):
                freq1[i, a] = weights[msa[:, i] == a].sum()
        freq1 = np.clip(freq1, 1e-9, None)
        freq1 /= freq1.sum(axis=1, keepdims=True)

        J, _ = self._gaussian_dca(msa, weights)
        self.J  = J
        self.DI = self._compute_di(J, freq1)
        return self.DI[:self.l_srna, self.l_srna:]   # [l_srna, l_mrna]

    def fit_plmdca(self, msa_file: str,
                   python37: str = 'python3.7',
                   conda_env: str = None) -> np.ndarray:
        """
        Higher-accuracy fit using pydca in a separate Python 3.7 env.
        Saves DI matrix to {msa_file}_di.npy then loads it back.
        python37 : path to python3.7 binary or 'python3.7'
        conda_env: conda env name to activate with 'source activate'
        """
        import subprocess, tempfile, textwrap
        out_npy  = msa_file + '_di.npy'
        script   = textwrap.dedent(f"""\
            import numpy as np
            from pydca.plmdca import plmdca
            plm = plmdca.PlmDCA(
                'RNA',
                '{msa_file}',
                num_iter_steps={self.max_iter},
                reg_lambda_pair={self.lam},
            )
            pairs = plm.compute_di()
            L  = {self.L}
            DI = np.zeros((L, L))
            for i, j, di in pairs:
                DI[int(i), int(j)] = di
                DI[int(j), int(i)] = di
            np.save('{out_npy}', DI)
            print('Saved DI to {out_npy}')
        """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py',
                                         delete=False) as tmp:
            tmp.write(script)
            tmp_path = tmp.name
        try:
            if conda_env:
                cmd = f'source activate {conda_env} && {python37} {tmp_path}'
                subprocess.run(cmd, shell=True, check=True,
                               executable='/bin/bash')
            else:
                subprocess.run([python37, tmp_path], check=True)
        finally:
            os.unlink(tmp_path)

        DI      = np.load(out_npy)
        self.DI = DI
        return DI[:self.l_srna, self.l_srna:]

    def fit_evcouplings(self, msa_file: str,
                        theta: float = 0.8,
                        pseudo_count: float = 0.5) -> np.ndarray:
        """
        Mean-field DCA using evcouplings math functions directly.
        Bypasses the Alignment class to avoid header/focus-column issues.
        Returns cross-block DI [l_srna, l_mrna].
        """
        from evcouplings.couplings.mean_field import (
            compute_covariance_matrix, reshape_invC_to_4d,
            direct_information, regularize_frequencies,
            regularize_pair_frequencies,
        )

        msa     = self._encode_msa(msa_file)       # [N, L] int
        N, L    = msa.shape
        q       = self.q

        # sequence weights
        weights = self._reweight(msa, theta=1.0 - theta)

        # single-site frequencies [L, q]
        f_i = np.zeros((L, q))
        for n in range(N):
            for i in range(L):
                a = msa[n, i]
                if a < q:
                    f_i[i, a] += weights[n]

        # pair frequencies [L, L, q, q] via einsum
        X = np.zeros((N, L, q))
        for n in range(N):
            for i in range(L):
                a = msa[n, i]
                if a < q:
                    X[n, i, a] = weights[n]
        f_ij = np.einsum('nia,njb->ijab', X, X)
        for i in range(L):                          # diagonal: f_ij[i,i,a,a] = f_i[i,a]
            for a in range(q):
                f_ij[i, i, a, a] = f_i[i, a]

        reg_f_i  = regularize_frequencies(f_i,  pseudo_count=pseudo_count)
        reg_f_ij = regularize_pair_frequencies(f_ij, pseudo_count=pseudo_count)

        C     = compute_covariance_matrix(reg_f_i, reg_f_ij)
        C_inv = -np.linalg.inv(C)
        J_ij  = reshape_invC_to_4d(C_inv, L, q)
        DI    = direct_information(J_ij, reg_f_i)  # [L, L]

        out_npy = msa_file + '_di.npy'
        np.save(out_npy, DI)
        self.DI = DI
        print(f'  evcouplings DCA done → {out_npy}')
        return DI[:self.l_srna, self.l_srna:]

    def fit_plmdca_direct(self, msa_file: str) -> np.ndarray:
        """Run plmDCA in-process (pydca must be installed in current env)."""
        from pydca.plmdca import plmdca as plmdca_mod
        out_npy = msa_file + '_di.npy'
        plm = plmdca_mod.PlmDCA(
            'RNA',
            msa_file,
            num_iter_steps=self.max_iter,
            reg_lambda_pair=self.lam,
        )
        pairs = plm.compute_di()
        DI = np.zeros((self.L, self.L))
        for i, j, di in pairs:
            DI[int(i), int(j)] = di
            DI[int(j), int(i)] = di
        np.save(out_npy, DI)
        self.DI = DI
        return DI[:self.l_srna, self.l_srna:]

    def edge_features(self, cross_di: np.ndarray) -> torch.Tensor:
        """
        Convert cross-block DI matrix to edge feature tensor.
        cross_di: [l_srna, l_mrna]
        Returns flattened normalised tensor for use as edge_attr.
        """
        t = torch.tensor(cross_di, dtype=torch.float32)
        t = (t - t.min()) / (t.max() - t.min() + 1e-9)
        return t   # [l_srna, l_mrna] — index [i,j] gives DI for sRNA pos i ↔ mRNA pos j

    def top_couplings(self, cross_di: np.ndarray,
                      n: int = 20) -> List[Tuple[int, int, float]]:
        """Return top-n sRNA-mRNA position couplings sorted by DI."""
        flat = [(i, j, cross_di[i, j])
                for i in range(cross_di.shape[0])
                for j in range(cross_di.shape[1])]
        return sorted(flat, key=lambda x: x[2], reverse=True)[:n]


# ── Trainer ───────────────────────────────────────────────────────────────────

class Trainer:
    def __init__(self, model: nn.Module, train_loader, val_loader, config: Dict):
        self.model        = model
        self.train_loader = train_loader
        self.val_loader   = val_loader
        self.config       = config
        self.k_neg        = config.get('k_neg', 10)
        self.neg_temp     = config.get('neg_temperature', 0.5)
        self.optimizer    = torch.optim.Adam(model.parameters(),
                                             lr=config.get('lr', 1e-3))
        self.loss_fn      = nn.BCEWithLogitsLoss()

    def to_device(self, batch, dev):
        return batch.to(dev)

    def _hard_negatives(self, srna_emb: torch.Tensor, mrna_emb: torch.Tensor,
                        edge_index: torch.Tensor,
                        temperature: float = 0.0) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Negative sampling with temperature control.

        temperature=0  : hard — strictly top-k nearest non-positive mRNAs (original behaviour)
        temperature>0  : soft — fetch k_neg*5 candidates, sample with prob ∝ exp(sim/temperature)
                         Higher temperature → closer to uniform random among candidates.
        """
        s = srna_emb.detach().cpu().float().numpy()
        m = mrna_emb.detach().cpu().float().numpy()
        faiss.normalize_L2(s)
        faiss.normalize_L2(m)

        n_mrna, dim = m.shape
        index = faiss.IndexFlatIP(dim)
        index.add(m)

        k_cand = min(self.k_neg * 5 if temperature > 0 else self.k_neg + 5, n_mrna)
        sims, I = index.search(s, k_cand)

        src_list, dst_list = edge_index[0].tolist(), edge_index[1].tolist()
        true_edges = set(zip(src_list, dst_list))
        rng = np.random.default_rng()

        neg_src, neg_dst = [], []
        for i, (sim_row, idx_row) in enumerate(zip(sims, I)):
            valid_sims, valid_idx = [], []
            for sim, j in zip(sim_row, idx_row):
                if j >= 0 and (i, int(j)) not in true_edges:
                    valid_sims.append(sim)
                    valid_idx.append(int(j))

            if not valid_idx:
                continue

            n_pick = min(self.k_neg, len(valid_idx))
            if temperature > 0:
                w = np.exp((np.array(valid_sims) - max(valid_sims)) / temperature)
                w /= w.sum()
                chosen = rng.choice(len(valid_idx), size=n_pick, replace=False, p=w)
            else:
                chosen = range(n_pick)   # top-k (hard)

            for c in chosen:
                neg_src.append(i)
                neg_dst.append(valid_idx[c])

        return (torch.tensor(neg_src, dtype=torch.long),
                torch.tensor(neg_dst, dtype=torch.long))

    def train_step(self, batch, dev):
        self.model.train()
        batch = self.to_device(batch, dev)
        srna_emb, mrna_emb = self.model(batch)

        edge_index = batch['srna', 'targets', 'mrna'].edge_index
        src, dst   = edge_index

        # positive loss
        pos_scores = (srna_emb[src] * mrna_emb[dst]).sum(dim=-1)
        pos_labels = batch['srna', 'targets', 'mrna'].y.float()
        pos_loss   = self.loss_fn(pos_scores, pos_labels)

        # temperature-scaled negative sampling (soft when neg_temp > 0)
        neg_src, neg_dst = self._hard_negatives(srna_emb, mrna_emb, edge_index,
                                                 temperature=self.neg_temp)
        neg_src = neg_src.to(dev)
        neg_dst = neg_dst.to(dev)
        neg_scores = (srna_emb[neg_src] * mrna_emb[neg_dst]).sum(dim=-1)
        neg_labels = torch.zeros(neg_scores.size(0), device=dev)
        neg_loss   = self.loss_fn(neg_scores, neg_labels)

        return pos_loss + neg_loss

    def train_epoch(self, dev):
        for batch in self.train_loader:
            loss = self.train_step(batch, dev)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()

    def valid_epoch(self, dev):
        self.model.eval()
        total = 0.0
        with torch.no_grad():
            for batch in self.val_loader:
                batch = self.to_device(batch, dev)
                loss  = self.train_step(batch, dev)
                total += loss.item()
        return total / len(self.val_loader)

