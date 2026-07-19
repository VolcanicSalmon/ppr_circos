#!/usr/bin/env python3
"""
pval_dist.py — Gibbs-ensemble scoring of predicted slice sites (bscore).

Module of a larger (plotting) object. Reusable components kept as classes:
  * BedReader — reads a pair_slice_to_srna.py slice-site BED and parses the
                `property` (name) field into a trigger -> target relationship.
  * KDETree   — spatial KD-tree (sklearn.neighbors.KDTree) over slice-site
                positions, one per target, for many-vs-many proximity queries.

MODEL (per target i). The trigger candidates associated with target i form the
sampling space. For each candidate j the EVENT is

      a_j = "trigger j is the one that triggered target i"

The events are mutually exclusive; we want P(a_j). Gibbs ensemble of A copies of the
triggered target, a_j copies in microstate j (Sum a_j = A). Multiplicity
W = A!/prod_j a_j!, maximised via Stirling with an energy constraint, gives

      bscore_j = P(a_j) = a_j/A = g_j * exp(-beta*E_j) / Z ,  Z = Sum_i g_i exp(-beta*E_i)

  * PRIOR / degeneracy  g_j = softmax(-beta * z(allen))  — allen is ONLY a prior;
    -beta inverts allen's scale (small allen = better pairing = higher prior).
  * LIKELIHOOD  exp(-beta*E_j) = P(uptick | trigger j fired), the MAIN evidence:
    the read-depth (phasiRNA) uptick, found by scipy.signal find_peaks over the whole
    TARGET TRANSCRIPT REGION (the phasiRNA array spreads kb from the cut, so a fixed
    radius around the slice site is too tight). weight_j = prominence of the peak
    NEAREST trigger j's slice site; peak_dist reports how far.

bscore_j = prior_j * weight_j / Sum_k prior_k * weight_k = P(trigger j -> target i).

COVERAGE: strand is irrelevant — ABSOLUTE read depth |fwd| + |rev| (rev bigWigs may be
negative). Peak thresholds are ABSOLUTE (--min_height/--min_prominence, in reads); no
baseline pseudocount (0 = genuinely no sRNA there).

PERFORMANCE: coverage fetched ONCE per target transcript region, find_peaks once;
targets processed by a thread pool, each worker with its OWN bigWig handles.

The `pval=` in the BED name is the OLD GSTAr/CleaveLand p-value: raw passthrough
(`pred_pval`, may be 'NA'), NEVER used in the bscore.
"""
import argparse
import collections
import json
import math
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
from collections import defaultdict

import numpy as np
import pyBigWig
import pysam
from scipy.signal import find_peaks
from sklearn.neighbors import KDTree


# ---- property-field parsing --------------------------------------------------

def parse_property_field(prop):
    """Parse a pair_slice_to_srna.py name/property field, e.g.

      baldrich_uncoll_Cd33274944_38::RHC10H1G0341.2::chr10_1:7419471-7430291(+)|
      RHC10H1G0341.2|22nt|allen=8|pval=NA|left_of_sRNA=384|RHC10H1G0341.2|
      baldrich_uncoll_Cd17153220_1

    trigger  = first token before the first ::   (the slicer that sets the cut)
    target   = second token
    region   = third token (the target transcript span) -> region_chrom/start/end
    """
    def grab(pat, cast=str, default=None):
        m = re.search(pat, prop)
        try:
            return cast(m.group(1)) if m else default
        except (TypeError, ValueError):
            return default

    tok0 = prop.split('|', 1)[0]              # trigger::target::region(strand)
    parts0 = tok0.split('::')
    trigger = parts0[0]                       # TRIGGER = first token, before the first ::
    target = parts0[1] if len(parts0) > 1 else grab(r'\|([^|]+)\|\d+nt')

    region_chrom = region_start = region_end = None
    if len(parts0) > 2:
        rm = re.match(r'([^:]+):(\d+)-(\d+)', parts0[2])
        if rm:
            region_chrom = rm.group(1)
            region_start = int(rm.group(2))
            region_end = int(rm.group(3))

    nearby = prop.rsplit('|', 1)[-1].strip()     # last field = the closest sRNA read (name_count)
    nm = re.search(r'_(\d+)$', nearby)           # its read count = trailing _N
    return {
        'trigger':      trigger,
        'target':       target,
        'region_chrom': region_chrom,
        'region_start': region_start,
        'region_end':   region_end,
        'size':         grab(r'(\d+)nt', int),
        'allen':        grab(r'allen=([0-9.]+)', float),
        'pred_pval':    grab(r'pval=([^|]+)'),   # raw passthrough string (may be 'NA')
        'left_of_sRNA': grab(r'left_of_sRNA=(-?\d+)', int),   # signed dist to closest sRNA
        'nearby_sRNA':  nearby,
        'nearby_count': int(nm.group(1)) if nm else 1,        # height of the closest sRNA
    }


# ---- prior (allen) -----------------------------------------------------------

def softmax(x):
    """Numerically stable softmax over a 1-D array."""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    x = x - np.nanmax(x)
    e = np.exp(x)
    return e / e.sum()


def allen_prior(allens, beta=1.0):
    """Prior / degeneracy g_j over a target's trigger population.

    allen is z-scored (EdgeScaler convention: (allen-mean)/std, std->1 if 0), then
    softmax(-beta * z). -beta inverts allen's scale so a smaller allen (better
    pairing) gets a larger prior. allen is ONLY a prior, never the main energy.
    """
    a = np.asarray(allens, dtype=float)
    if a.size == 0:
        return a
    mu = np.nanmean(a)
    a = np.where(np.isnan(a), mu, a)                # neutral (z=0) for missing allen
    sd = a.std()
    z = (a - mu) / sd if sd > 0 else np.zeros_like(a)
    return softmax(-beta * z)


# ---- KD-tree over slice positions (kept for the larger object) ---------------

class KDETree(object):
    """Spatial KD-tree (sklearn.neighbors.KDTree) over 1-D genomic positions for a
    (possibly many-vs-many) trigger->target relationship. Typically one per target;
    `meta` holds the parsed slicer record parallel to each position."""

    def __init__(self, positions, meta=None):
        self.positions = np.asarray(list(positions), dtype=float).reshape(-1, 1)
        self.meta = list(meta) if meta is not None else [None] * len(self.positions)
        self.tree = KDTree(self.positions) if len(self.positions) else None

    def within(self, coord, radius):
        """Indices of stored positions within `radius` bp of `coord`."""
        if self.tree is None:
            return []
        return list(self.tree.query_radius(np.array([[coord]], float), r=radius)[0])

    def nearest(self, coord, k=1):
        """(index, distance) of the k nearest stored positions to `coord`."""
        if self.tree is None:
            return []
        k = min(k, len(self.positions))
        dist, idx = self.tree.query(np.array([[coord]], float), k=k)
        return list(zip(idx[0].tolist(), dist[0].tolist()))


# ---- BED reader (kept as a class — component of the bigger object) -----------

class BedReader(object):
    def __init__(self, path):
        self.path = path
        self.nfields = 6
        self.nlines = 0
        self.fields = ['chrom', 'start', 'end', 'property', 'allen', 'strand']
        self.bedinterval = collections.namedtuple('bedinterval', self.fields)
        self.intervals = []
        self.trigger_target = defaultdict(list)   # trigger -> [record, ...]
        self.by_target = defaultdict(list)         # target  -> [record, ...] (sampling space)
        self.kdetrees = {}                         # target  -> KDETree

    def read(self):
        with open(self.path) as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                fld = line.rstrip('\n').split('\t')
                allen = None
                if len(fld) > 4 and fld[4] not in ('.', ''):
                    allen = float(fld[4])
                iv = self.bedinterval(
                    fld[0], int(fld[1]), int(fld[2]),
                    fld[3] if len(fld) > 3 else '',
                    allen,
                    fld[5] if len(fld) > 5 else '+',
                )
                self.intervals.append(iv)
                self.nlines += 1
        return self

    def parse_property(self):
        """Build trigger->records, target->records (the per-target sampling space),
        and a per-target KDETree of slice positions."""
        by_pos = defaultdict(list)
        for iv in self.intervals:
            m = parse_property_field(iv.property)
            m.update(chrom=iv.chrom, site=int(iv.start), strand=iv.strand)
            if iv.allen is not None:            # BED score column is authoritative
                m['allen'] = iv.allen
            self.trigger_target[m['trigger']].append(m)
            self.by_target[m['target']].append(m)
            by_pos[m['target']].append(m['site'])
        self.kdetrees = {t: KDETree(by_pos[t], self.by_target[t])
                         for t in self.by_target}
        return self.trigger_target, self.kdetrees


# ---- coverage / uptick (transcript-region find_peaks) ------------------------

def _coverage(bw, chrom, start, end):
    """0-based half-open ABSOLUTE coverage (|value|, NaN->0), clipped to chrom bounds.
    Returns (values, actual_start) or (None, None)."""
    L = bw.chroms().get(chrom)
    if L is None:
        return None, None
    s, e = max(0, start), min(L, end)
    if e <= s:
        return None, None
    v = np.asarray(bw.values(chrom, s, e), dtype=float)
    v[np.isnan(v)] = 0.0
    return np.abs(v), s                                  # strand-agnostic: |depth|


def combined_coverage(fw, rv, chrom, start, end):
    """Total absolute read depth = |fwd| + |rev| over [start, end). Either handle may
    be None. Returns (values, actual_start) or (None, None)."""
    cf, af = _coverage(fw, chrom, start, end) if fw is not None else (None, None)
    cr, ar = _coverage(rv, chrom, start, end) if rv is not None else (None, None)
    if cf is None and cr is None:
        return None, None
    if cf is None:
        return cr, ar
    if cr is None:
        return cf, af
    if af != ar:                                        # align on the max start
        off = ar - af
        if off > 0:
            cf = cf[off:]
        elif off < 0:
            cr = cr[-off:]
        af = max(af, ar)
    n = min(cf.size, cr.size)
    return cf[:n] + cr[:n], af


def region_upticks(cov, act, records, min_height, min_prominence, pad=421):
    """find_peaks once over the target-transcript coverage `cov` (starting at genomic
    `act`); each site is matched to the NEAREST peak in the region (the transcript is the
    search space).

    weight = an uptick *gradient*, not raw prominence: the angle of the rise from the
    slice site up to that peak MINUS the angle up to the closest sRNA read —
        max(0, atan2(h_peak, d_peak) - atan2(h_srna, d_srna))
    with d_peak = |peak-site|, d_srna = |left_of_sRNA| (already in the name field),
    h_srna = the nearby read's count, and heights normalized by the region baseline /
    distances by `pad` (the vicinity window) so it is scale-free. Rewards peaks that are
    both tall AND close relative to the local sRNA context."""
    out = [None] * len(records)
    if cov is None or cov.size == 0:
        return out
    region_med = float(np.median(cov))
    hscale = max(1.0, region_med)                        # normalize heights by baseline
    peaks, props = find_peaks(cov, height=min_height, prominence=min_prominence)
    prom = props.get('prominences', np.empty(0))
    gp = act + peaks                                     # genomic peak positions

    for i, r in enumerate(records):
        site = r['site']
        if peaks.size == 0:
            out[i] = {'weight': 0.0, 'peak_pos': None, 'peak_dist': None, 'n_peaks': 0,
                      'baseline': region_med, 'has_uptick': False, 'uptick_side': None,
                      'peak_h': 0.0, 'prominence': 0.0}
            continue
        j = int(np.argmin(np.abs(gp - site)))            # nearest peak to slice site
        peak_pos = int(gp[j])
        h_peak = float(cov[int(peaks[j])])               # coverage height AT the peak
        d_peak = max(1.0, float(abs(peak_pos - site)))
        lo = r.get('left_of_sRNA')
        d_srna = max(1.0, float(abs(lo))) if lo is not None else d_peak
        h_srna = float(r.get('nearby_count') or 1)
        a_peak = math.atan2(h_peak / hscale, d_peak / pad)
        a_srna = math.atan2(h_srna / hscale, d_srna / pad)
        grad = abs(a_peak - a_srna)                       # angular separation (test: rescues far peaks)
        out[i] = {'weight': grad, 'peak_pos': peak_pos,
                  'peak_dist': int(peak_pos - site), 'n_peaks': int(peaks.size),
                  'baseline': region_med, 'has_uptick': True,
                  'uptick_side': 'upstream' if peak_pos < site else 'downstream',
                  'peak_h': h_peak, 'prominence': float(prom[j])}
    return out


# ---- ensemble scoring --------------------------------------------------------

# --- Axtell phase score (copied from calc_phase_score.py) ---------------------

def reads_in_cluster(bam, chrom, cl_start, cl_end, strand='.'):
    """{5'-pos: count} for reads whose 5' end is in [cl_start, cl_end]."""
    pos_counts = defaultdict(int)
    try:
        for read in bam.fetch(chrom, cl_start, cl_end + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if strand != '.' and (('-' if read.is_reverse else '+') != strand):
                continue
            pos5 = (read.reference_end - 1) if read.is_reverse else read.reference_start
            if cl_start <= pos5 <= cl_end:
                pos_counts[pos5] += 1
    except (ValueError, KeyError):
        pass
    return dict(pos_counts)


def calc_axtell(pos_counts, cl_start, cl_end, phase_size):
    """Axtell phase score over [cl_start, cl_end]; best across offsets, or None."""
    if not pos_counts:
        return None
    best = None
    for offset in range(phase_size):
        n_max = (cl_end - cl_start) // phase_size + 2
        phased = [cl_start + offset + k * phase_size for k in range(n_max)
                  if cl_start + offset + k * phase_size <= cl_end]
        if not phased:
            continue
        phased_set = set(phased)
        pi = sum(pos_counts.get(p, 0) for p in phased)
        u = sum(v for p, v in pos_counts.items() if p not in phased_set)
        n = sum(1 for p in phased if pos_counts.get(p, 0) > 0)
        if n <= 2:
            continue
        score = (n - 2) * math.log(1 + 10 * pi / (1 + u))
        if best is None or score > best:
            best = score
    return best


def trigger_phasedness(bam_path, alignments, triggers, pad=421, phase_sizes=(21, 22, 23),
                       threads=4):
    """Phasedness of the region each trigger comes from: the Axtell phase score of the
    +/- pad window around the trigger's genomic mapping location, MAX over phase_sizes
    and over the trigger's (multimapped) locations. Returns {trigger: score}.

    UNIQUE locations are scored once (deduped across triggers) and the pysam.fetch pre-
    pass is parallelised (each worker its own BAM handle). NOTE: a BAM on a network
    mount (SMB) makes the random fetches very slow — copy it local first."""
    locs = defaultdict(set)
    uniq = set()
    with open(alignments) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 4 and p[3] in triggers:
                try:
                    loc = (p[0], int(p[1]), int(p[2]))
                except ValueError:
                    continue
                locs[p[3]].add(loc)
                uniq.add(loc)

    uniq = list(uniq)
    sys.stderr.write(f"phasedness: {len(uniq)} unique locations for {len(locs)} triggers\n")
    loc_phase = {}
    q = queue.Queue()
    for loc in uniq:
        q.put(loc)
    lock = threading.Lock()
    done = [0]

    def wk():
        bam = pysam.AlignmentFile(bam_path, 'rb')
        local = {}
        while True:
            try:
                c, s, e = q.get_nowait()
            except queue.Empty:
                break
            ws, we = max(0, s - pad), e + pad
            pc = reads_in_cluster(bam, c, ws, we)
            best = 0.0
            for ps in phase_sizes:
                sc = calc_axtell(pc, ws, we, ps)
                if sc is not None and sc > best:
                    best = sc
            local[(c, s, e)] = best
            with lock:
                done[0] += 1
                if done[0] % 20000 == 0:
                    sys.stderr.write(f"  phasedness: {done[0]}/{len(uniq)} locations\n")
        bam.close()
        with lock:
            loc_phase.update(local)

    ts = [threading.Thread(target=wk) for _ in range(max(1, threads))]
    for t in ts:
        t.start()
    for t in ts:
        t.join()

    return {trig: max((loc_phase.get(l, 0.0) for l in ls), default=0.0)
            for trig, ls in locs.items()}


def score_target(fw, rv, records, min_height, min_prominence, beta, trig_phase):
    """P(a_j) over one target's trigger population. Coverage (|fwd|+|rev|) is fetched
    once over the target transcript region (union of the records' region spans).

        bscore_j = prior_j * weight_j * trigger_phase_j / Z   (== a_j / A in the ensemble)

    trigger_phase_j = Axtell phase score of the +/-421 window at the trigger's mapping
    location (how phased the region the trigger comes from is). If no BAM was given
    (trig_phase empty) the factor is 1.0 -> falls back to the prior*uptick formula.
    """
    prior = allen_prior([r.get('allen') for r in records], beta=beta)   # g_j

    chrom = records[0].get('region_chrom') or records[0]['chrom']
    starts = [r['region_start'] for r in records if r.get('region_start') is not None]
    ends = [r['region_end'] for r in records if r.get('region_end') is not None]
    if starts and ends:
        lo, hi = min(starts), max(ends)
    else:                                               # fallback: span the sites
        s = [r['site'] for r in records]
        lo, hi = min(s) - 3000, max(s) + 3000

    cov, act = combined_coverage(fw, rv, chrom, lo, hi)
    uptick = region_upticks(cov, act, records, min_height, min_prominence)

    covered = [(records[i], float(prior[i]), uptick[i])
               for i in range(len(records)) if uptick[i] is not None]

    def phase_of(r):
        # miRs are never in the alignment BED (they're query-only, no phased origin) -> absent
        # from trig_phase -> keep them with a NEUTRAL 1.0. A mapped-but-unphased phasiRNA IS in
        # trig_phase (with 0.0) and stays downweighted.
        return float(trig_phase.get(r['trigger'], 1.0)) if trig_phase else 1.0

    weights = np.array([g * U['weight'] * phase_of(r) for r, g, U in covered], dtype=float)
    Z = weights.sum()                                   # partition function
    bscore = weights / Z if Z > 0 else np.zeros_like(weights)   # P(a_j) = a_j/A

    out = []
    for (r, g, U), b in zip(covered, bscore):
        out.append({**r, **U, 'prior': g, 'trigger_phase': phase_of(r), 'bscore': float(b)})
    return out


def run(slice_bed, bw_fwd=None, bw_rev=None, min_height=2.0, min_prominence=2.0,
        beta=1.0, num_threads=4, alignments=None, bam=None):
    """Score all targets. |fwd|+|rev| summed; peaks over the target transcript region;
    targets distributed over a thread pool, each worker with its own bigWig handles.
    If --bam (+ --alignments) given, multiply each trigger's score by the phasedness
    (Axtell phase score) of the region it comes from."""
    reader = BedReader(slice_bed).read()
    reader.parse_property()

    trig_phase = {}
    if bam and alignments:
        triggers = set(r['trigger'] for recs in reader.by_target.values() for r in recs)
        trig_phase = trigger_phasedness(bam, alignments, triggers, threads=num_threads)
        nz = sum(1 for v in trig_phase.values() if v > 0)
        sys.stderr.write(f"phasedness computed for {len(trig_phase)} triggers "
                         f"({nz} non-zero, max {max(trig_phase.values(), default=0):.3g})\n")
        if nz == 0:
            sys.stderr.write("WARNING: ALL phasedness = 0 — BAM likely unindexed or chrom names "
                             "don't match the alignment BED. Falling back to prior*uptick "
                             "(phase factor ignored). Fix the BAM and re-run to use phasedness.\n")
            trig_phase = {}

    work = queue.Queue()
    for item in reader.by_target.items():
        work.put(item)

    rows, lock = [], threading.Lock()

    def worker():
        fw = pyBigWig.open(bw_fwd) if bw_fwd else None
        rv = pyBigWig.open(bw_rev) if bw_rev else None
        local = []
        while True:
            try:
                target, records = work.get_nowait()
            except queue.Empty:
                break
            local.extend(score_target(fw, rv, records, min_height, min_prominence, beta, trig_phase))
        if fw:
            fw.close()
        if rv:
            rv.close()
        with lock:
            rows.extend(local)

    n = max(1, num_threads)
    threads = [threading.Thread(target=worker) for _ in range(n)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return rows


def _place(src, dst_dir, link):
    """Copy (or symlink) src into dst_dir; return the placed path."""
    if not src:
        return None
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, os.path.basename(src))
    if os.path.abspath(src) != os.path.abspath(dst):
        if os.path.lexists(dst):          # drop any existing file/symlink first
            os.remove(dst)                # (else copy2 follows an old symlink back to src -> SameFileError)
        if link:
            os.symlink(os.path.abspath(src), dst)
        else:
            shutil.copy(src, dst)     # content+mode only; NOT copy2 (chflags EPERM on SMB sources)
    return dst


def install_dataset(args, bscore_path):
    """Copy inputs into <install>/<name>/{bw,slice,genome}, build origins, and
    register the dataset in <install>/datasets.json (read by the circos app)."""
    root = os.path.abspath(args.install)
    base = os.path.join(root, args.name)
    rel = lambda p: os.path.relpath(p, root) if p else None
    entry = {'name': args.name}

    bw = []
    for b in (args.bw_fwd, args.bw_rev):
        p = _place(b, os.path.join(base, 'bw'), args.link)
        if p:
            bw.append(rel(p))
    entry['bw'] = bw
    entry['bscore'] = rel(bscore_path)

    _place(args.slice_bed, os.path.join(base, 'slice', 'slice_bed'), args.link)

    if args.alignments:                       # build trigger_origins.bed via the sibling script
        obed = os.path.join(base, 'slice', 'slice_bed', 'trigger_origins.bed')
        os.makedirs(os.path.dirname(obed), exist_ok=True)
        here = os.path.dirname(os.path.abspath(__file__))
        subprocess.run([sys.executable, os.path.join(here, 'build_origins.py'),
                        '--bscore', bscore_path, '--alignments', args.alignments,
                        '--out', obed], check=True)
        entry['origins'] = rel(obed)

    if args.slice_fa:
        entry['fasta'] = rel(_place(args.slice_fa, os.path.join(base, 'slice', 'slice_fa'), args.link))
    if args.transcript:
        entry['transcript'] = rel(_place(args.transcript, os.path.join(base, 'genome'), args.link))
    if args.fai:
        entry['fai'] = rel(_place(args.fai, os.path.join(base, 'genome'), args.link))
    if args.bam:                                  # sorted+indexed sRNA BAM -> <name>/bam/ (GW live phase/DC)
        bam_dst = _place(args.bam, os.path.join(base, 'bam'), args.link)
        for idx in (args.bam + '.bai', args.bam + '.csi',
                    os.path.splitext(args.bam)[0] + '.bai', os.path.splitext(args.bam)[0] + '.csi'):
            if os.path.exists(idx):
                _place(idx, os.path.join(base, 'bam'), args.link)
        if bam_dst:
            entry['bam'] = rel(bam_dst)

    manifest = os.path.join(root, 'datasets.json')
    data = {'datasets': []}
    if os.path.exists(manifest):
        try:
            data = json.load(open(manifest))
        except Exception:
            pass
    data.setdefault('datasets', [])
    data['datasets'] = [d for d in data['datasets'] if d.get('name') != args.name]  # replace if re-run
    data['datasets'].append(entry)
    json.dump(data, open(manifest, 'w'), indent=2)
    sys.stderr.write(f"installed '{args.name}' -> {base}\n  registered in {manifest}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--slice_bed', required=True, help='pair_slice_to_srna.py BED')
    ap.add_argument('--bw_fwd', help='forward-strand coverage bigWig')
    ap.add_argument('--bw_rev', help='reverse-strand coverage bigWig')
    ap.add_argument('--min_height', type=float, default=2.0,
                    help='find_peaks absolute height threshold (reads)')
    ap.add_argument('--min_prominence', type=float, default=2.0,
                    help='find_peaks absolute prominence threshold (reads)')
    ap.add_argument('--beta', type=float, default=1.0,
                    help='inverse temperature of the allen prior (small = flat prior)')
    ap.add_argument('--threads', type=int, default=4, help='worker threads')
    ap.add_argument('--bam', help='sorted+indexed sRNA BAM -> multiply score by region phasedness (Axtell)')
    ap.add_argument('--out', default='-', help='output TSV (default stdout)')
    # --- install into the app data tree + register in datasets.json (option B) ---
    ap.add_argument('--install', help='app data root; enables install mode (e.g. bin/data)')
    ap.add_argument('--name', help='dataset name/key under --install (e.g. dm_infec)')
    ap.add_argument('--slice_fa', help='trigger FASTA -> slice/slice_fa/')
    ap.add_argument('--alignments', help='sRNA alignment BED -> build slice/slice_bed/trigger_origins.bed')
    ap.add_argument('--transcript', help='transcript regions BED -> genome/')
    ap.add_argument('--fai', help='genome .fa.fai -> genome/')
    ap.add_argument('--link', action='store_true', help='symlink instead of copy (avoid duplicating big files)')
    args = ap.parse_args()

    if not (args.bw_fwd or args.bw_rev):
        ap.error('provide --bw_fwd and/or --bw_rev')
    if args.install and not args.name:
        ap.error('--install requires --name')

    rows = run(args.slice_bed, bw_fwd=args.bw_fwd, bw_rev=args.bw_rev,
               min_height=args.min_height, min_prominence=args.min_prominence,
               beta=args.beta, num_threads=args.threads,
               alignments=args.alignments, bam=args.bam)

    # in install mode the bscore lands in <install>/<name>/slice/
    out_path = args.out
    if args.install:
        slice_dir = os.path.join(os.path.abspath(args.install), args.name, 'slice')
        os.makedirs(slice_dir, exist_ok=True)
        out_path = os.path.join(slice_dir, f'{args.name}_bscore.tsv')

    cols = ['chrom', 'site', 'strand', 'target', 'trigger', 'size', 'allen',
            'pred_pval', 'left_of_sRNA', 'region_start', 'region_end', 'baseline',
            'n_peaks', 'has_uptick', 'uptick_side', 'peak_pos', 'peak_dist',
            'prior', 'weight', 'trigger_phase', 'bscore']
    if out_path != '-':
        d = os.path.dirname(out_path)
        if d:
            os.makedirs(d, exist_ok=True)
    out = sys.stdout if out_path == '-' else open(out_path, 'w')
    out.write('\t'.join(cols) + '\n')
    for r in rows:
        out.write('\t'.join(
            '' if r.get(c) is None else
            (f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]))
            for c in cols) + '\n')
    if out is not sys.stdout:
        out.close()
    sys.stderr.write(f"Done — {len(rows)} rows -> {out_path}\n")

    if args.install:
        install_dataset(args, out_path)


if __name__ == '__main__':
    main()
