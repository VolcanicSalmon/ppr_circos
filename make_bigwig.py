#!/usr/bin/env python3
"""
make_bigwig.py — genome-wide stranded coverage bigWig from an alignment BED, written
directly with pyBigWig (no bedGraphToBigWig / UCSC tools needed).

This is a trackplot APP helper (like build_origins.py) — it is NOT the phasepipe
bed_to_bigwig.py and must not be confused with it. Count-weighted: the read count is
parsed from the name (ShortStack `_Cd{id}_{count}` or `_x{count}`); use --no_count for 1/read.

usage:
  python3 make_bigwig.py --bed baldrich_uncoll_condensed.bed \
     --chrom_sizes rh_chrom.sizes --prefix rh_bw --outdir bw
  -> bw/rh_bw_fwd.bw, bw/rh_bw_rev.bw
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pyBigWig

_X = re.compile(r'_x(\d+)(?:[^0-9]|$)')
_TAILNUM = re.compile(r'_(\d+)$')


def count_from_name(name):
    """Read count = last underscore-delimited field (e.g. ..._Cd18575949_36 -> 36);
    ShortStack `_xN` tried first. Matches bed_to_bigwig.py::count_from_name."""
    m = _X.findall(name)
    if m:
        return int(m[-1])
    m = _TAILNUM.search(name)
    return int(m.group(1)) if m else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--chrom_sizes', required=True, help='chrom<TAB>size')
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--outdir', default='.')
    ap.add_argument('--no_count', action='store_true', help='count 1 per read (ignore name count)')
    args = ap.parse_args()

    sizes = {}
    with open(args.chrom_sizes) as f:
        for l in f:
            if l.strip():
                a = l.split('\t')
                sizes[a[0]] = int(a[1])

    # sweep-line events: strand -> chrom -> [positions[], deltas[]]
    ev = {'+': defaultdict(lambda: [[], []]), '-': defaultdict(lambda: [[], []])}
    n = 0
    with open(args.bed) as f:
        for line in f:
            c = line.split('\t')
            if len(c) < 6:
                continue
            try:
                s, e = int(c[1]), int(c[2])
            except ValueError:
                continue
            strand = c[5].rstrip('\n')
            if strand not in ('+', '-'):
                strand = '+'
            cnt = 1 if args.no_count else count_from_name(c[3])
            L = ev[strand][c[0]]
            L[0].append(s); L[1].append(cnt)
            L[0].append(e); L[1].append(-cnt)
            n += 1
            if n % 5_000_000 == 0:
                sys.stderr.write(f"  read {n:,} bed lines\n")
    sys.stderr.write(f"read {n:,} bed lines\n")

    os.makedirs(args.outdir, exist_ok=True)
    for strand, suffix in (('+', 'fwd'), ('-', 'rev')):
        out = os.path.join(args.outdir, f"{args.prefix}_{suffix}.bw")
        bw = pyBigWig.open(out, 'w')
        chroms = sorted(c for c in ev[strand] if c in sizes)
        bw.addHeader([(c, sizes[c]) for c in chroms])
        total_iv = 0
        for c in chroms:
            pos = np.asarray(ev[strand][c][0], dtype=np.int64)
            dl = np.asarray(ev[strand][c][1], dtype=np.float64)
            order = np.argsort(pos, kind='mergesort')
            pos, dl = pos[order], dl[order]
            csum = np.cumsum(dl)                                   # coverage after each event
            upos = np.unique(pos)
            last = np.searchsorted(pos, upos, side='right') - 1    # last event index per unique pos
            cov = csum[last]                                       # coverage in [upos[i], upos[i+1])
            starts, ends, vals = upos[:-1], upos[1:], cov[:-1]
            mask = vals > 0
            if mask.any():
                bw.addEntries([c] * int(mask.sum()), starts[mask].tolist(),
                              ends=ends[mask].tolist(), values=vals[mask].astype(float).tolist())
                total_iv += int(mask.sum())
        bw.close()
        sys.stderr.write(f"[{suffix}] -> {out}  ({len(chroms)} chroms, {total_iv:,} intervals)\n")


if __name__ == '__main__':
    main()
