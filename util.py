import argparse
import pickle
import numpy as np
import pandas as pd
import os

# Continuous edge scalars to z-score; ratio/fraction columns left as-is
EDGE_SCALE_COLS = ['AllenScore', 'MFEsite']


class EdgeScaler:
    """
    Fits and applies z-score normalisation to continuous edge scalar columns
    (AllenScore, MFEsite) derived from GSTAr output.
    MFEratio and paired_frac are already in [0,1] and are left untouched.
    """

    def __init__(self):
        self.mean_ = {}
        self.std_  = {}

    def fit(self, df: pd.DataFrame) -> 'EdgeScaler':
        for col in EDGE_SCALE_COLS:
            vals = pd.to_numeric(df[col], errors='coerce').dropna()
            self.mean_[col] = float(vals.mean())
            self.std_[col]  = float(vals.std()) if vals.std() > 0 else 1.0
        return self

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        for col in EDGE_SCALE_COLS:
            if col in df.columns:
                df[col] = (pd.to_numeric(df[col], errors='coerce')
                           - self.mean_[col]) / self.std_[col]
        return df

    def fit_transform(self, df: pd.DataFrame) -> pd.DataFrame:
        return self.fit(df).transform(df)

    def save(self, path: str):
        with open(path, 'wb') as f:
            pickle.dump({'mean': self.mean_, 'std': self.std_}, f)

    @classmethod
    def load(cls, path: str) -> 'EdgeScaler':
        with open(path, 'rb') as f:
            d = pickle.load(f)
        s = cls()
        s.mean_ = d['mean']
        s.std_  = d['std']
        return s


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('--gstar',      required=True, help='GSTAr output file')
    p.add_argument('--out_fasta',  required=True, help='Output concatenated MSA FASTA')
    p.add_argument('--allen_max',  type=int,   default=3,    help='Max AllenScore to keep (default: 3)')
    p.add_argument('--min_paired', type=float, default=0.75, help='Min paired fraction (default: 0.75)')
    p.add_argument('--use_plmdca',      action='store_true', help='Run plmDCA via py37 subprocess')
    p.add_argument('--use_gaussian_dca',action='store_true', help='Run Gaussian mean-field DCA (pure numpy, no pydca needed)')
    p.add_argument('--conda_env',    default='py37',           help='Conda env with pydca (default: py37)')
    p.add_argument('--python37',     default='python3.7',      help='Python 3.7 binary (default: python3.7)')
    p.add_argument('--scaler_out',   default=None,             help='Save fitted EdgeScaler to this path (.pkl)')
    return p.parse_args()

# Original GSTAr columns (16 cols)
GSTAR_COLS = [
    'Query', 'miRNA', 'sde', 'Transcript', 'gde',
    'AllenScore', 'TStart', 'TStop', 'TSlice',
    'MFEperfect', 'MFEsite', 'MFEratio',
    'Paired', 'Unpaired', 'Structure', 'Sequence'
]

# Extended version with cross-species columns (19 cols)
GSTAR_COLS_EXT = [
    'Query', 'miRNA', 'sde', 'Transcript', 'gde',
    'from_org', 'to_org', 'mode',
    'AllenScore', 'TStart', 'TStop', 'TSlice',
    'MFEperfect', 'MFEsite', 'MFEratio',
    'Paired', 'Unpaired', 'Structure', 'Sequence'
]


# 13-col compact format (no miRNA/sde/gde): Query and Transcript run directly
GSTAR_COLS_13 = [
    'Query', 'Transcript', 'TStart', 'TStop', 'TSlice',
    'MFEperfect', 'MFEsite', 'MFEratio', 'AllenScore',
    'Paired', 'Unpaired', 'Structure', 'Sequence'
]


def read_gstar(tsv_file: str) -> pd.DataFrame:
    """
    Read GSTAr output. Handles:
      13-col compact  (no miRNA/sde/gde — Query IS the sRNA id)
      15-col          (miRNA column missing from 16-col)
      16-col standard
      19-col extended (from_org/to_org/mode columns present)
    Header lines (starting with 'Query') are skipped automatically.
    Always returns a DataFrame with the standard GSTAR_COLS columns;
    miRNA/sde/gde are filled with NA for compact files.
    """
    rows    = []
    compact = False
    with open(tsv_file) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip() or line.startswith('#'):
                continue
            fields = [f.strip('"') for f in line.split('\t')]
            n = len(fields)
            # skip header
            if fields[0] == 'Query':
                compact = (n == 13)
                continue
            if n == 13:   # compact: Query Transcript TStart TStop TSlice MFEperfect MFEsite MFEratio AllenScore Paired Unpaired Structure Sequence
                # reorder to match GSTAR_COLS: Query miRNA sde Transcript gde AllenScore TStart TStop TSlice MFEperfect MFEsite MFEratio Paired Unpaired Structure Sequence
                fields = [fields[0], '', '', fields[1], '',
                          fields[8], fields[2], fields[3], fields[4],
                          fields[5], fields[6], fields[7],
                          fields[9], fields[10], fields[11], fields[12]]
            elif n == 19:  # extended with from_org/to_org/mode
                fields = fields[:5] + fields[8:]
            elif n == 15:  # miRNA column missing
                fields.insert(1, '')
            elif n != 16:
                continue
            rows.append(fields)

    df = pd.DataFrame(rows, columns=GSTAR_COLS)
    df['miRNA'] = df['miRNA'].replace('', pd.NA)
    df['Query'] = df['Query'].replace('', pd.NA).ffill()

    return df[GSTAR_COLS]


def build_msa_from_gstar(tsv_file: str, out_fasta: str,
                          allen_max: int = 3,
                          min_paired: float = 0.75) -> tuple:
    """
    Parse GSTAr output and build concatenated sRNA|mRNA_window MSA for DCA.

    Filters:
        allen_max   : keep AllenScore <= allen_max (lower = stronger hybridization energy)
        min_paired  : minimum fraction of paired positions

    Returns (l_srna, l_mrna) — lengths after padding.
    """
    df = read_gstar(tsv_file)

    # parse Sequence column → sRNA and mRNA window
    split       = df['Sequence'].str.split('&', n=1, expand=True)
    df['srna']  = split[0].str.strip()
    df['mrna_win'] = split[1].str.strip()

    df['AllenScore'] = pd.to_numeric(df['AllenScore'], errors='coerce')
    df = df[df['AllenScore'] <= allen_max].copy()
    df['paired_frac'] = df['Paired'].apply(_parse_paired_frac)
    df = df[df['paired_frac'] >= min_paired].copy()

    # filter to canonical sRNA sizes (non-gap length)
    df['srna_len'] = df['srna'].str.replace('-', '', regex=False).str.len()
    size_counts = df['srna_len'].value_counts().sort_index()
    print(f'  sRNA size distribution before filter: {size_counts.to_dict()}')
    df = df[df['srna_len'].isin([21, 22, 23, 24])].copy()

    if df.empty:
        raise ValueError('No interactions passed filters — check allen_max / min_paired')

    # pad to max length per side
    max_s = df['srna'].str.len().max()
    max_m = df['mrna_win'].str.len().max()

    with open(out_fasta, 'w') as f:
        for _, row in df.iterrows():
            s   = row['srna'].upper().replace('T', 'U').ljust(max_s, '-')
            m   = row['mrna_win'].upper().replace('T', 'U').ljust(max_m, '-')
            rid = f"{row['Query']}|{row['Transcript']}"
            f.write(f'>{rid}\n{s}{m}\n')

    print(f'  {len(df)} interactions → {out_fasta}  '
          f'(sRNA padded to {max_s}, mRNA window padded to {max_m})')
    return max_s, max_m


def _parse_paired_frac(paired_str: str) -> float:
    """
    Parse Paired field (e.g. '1-24,594-571') → fraction of positions paired.
    Returns proportion of sRNA length that is in a paired block.
    """
    if not isinstance(paired_str, str) or paired_str in ('NA', '.', ''):
        return 0.0
    total = 0
    for block in paired_str.split(';'):
        parts = block.split(',')
        for part in parts:
            if '-' in part:
                a, b = part.split('-')
                try:
                    total += abs(int(b) - int(a)) + 1
                except ValueError:
                    pass
    # divide by 2 since Paired lists both sRNA and mRNA coordinate spans
    return (total / 2) / 24   # normalise by typical sRNA length


def main():
    args = get_args()
    df = read_gstar(args.gstar)
    scaler = EdgeScaler().fit(df)
    if args.scaler_out:
        scaler.save(args.scaler_out)
        print(f'  EdgeScaler saved → {args.scaler_out}')

    l_srna, l_mrna = build_msa_from_gstar(
        args.gstar, args.out_fasta,
        allen_max=args.allen_max,
        min_paired=args.min_paired
    )
    if args.use_plmdca or args.use_gaussian_dca:
        from tasnet import PottsDCA
        dca = PottsDCA(l_srna=l_srna, l_mrna=l_mrna)
        if args.use_gaussian_dca:
            cross_di = dca.fit_evcouplings(args.out_fasta)
        else:
            cross_di = dca.fit_plmdca(args.out_fasta,
                                      python37=args.python37,
                                      conda_env=args.conda_env)
        print('Top couplings:')
        for i, j, di in dca.top_couplings(cross_di):
            print(f'  sRNA pos {i+1} <-> mRNA pos {j+1}  DI={di:.4f}')


if __name__ == '__main__':
    main()
