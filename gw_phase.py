#!/usr/bin/env python3
"""
gw_phase.py — Axtell phase score + dominant DicerCall (DC) for a single genomic window.

A trackplot GW helper (NOT a phasepipe script). `reads_in_cluster` and `calc_axtell`
are copied VERBATIM from calc_phase_score.py so the score matches the pipeline. DC is the
dominant read length (5'-anchored) in the window.

The app (ofApp / GW view) shells out to this per window and parses one stdout line:
    phase=<float|NA>  DC=<int|NA>  n=<reads>

usage:
    gw_phase.py --bam sorted.bam --chrom chr11_1 --start 9160528 --end 9166528
                [--phase_sizes 21 22 23 24] [--strand + | - | .]
"""
import argparse
import math
import sys
from collections import defaultdict, Counter

import pysam


# ── copied verbatim from calc_phase_score.py ──────────────────────────────────
def reads_in_cluster(bam, chrom, cl_start, cl_end, strand='.'):
    """
    Returns (pos_counts, pos_len_counts) for reads whose 5' end is in [cl_start, cl_end].
    pos_counts:     {pos: total_count}          — all reads
    pos_len_counts: {pos: {read_length: count}} — by size
    """
    pos_counts     = defaultdict(int)
    pos_len_counts = defaultdict(lambda: defaultdict(int))
    try:
        for read in bam.fetch(chrom, cl_start, cl_end + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if strand != '.' and (('-' if read.is_reverse else '+') != strand):
                continue
            pos5 = (read.reference_end - 1) if read.is_reverse else read.reference_start
            if cl_start <= pos5 <= cl_end:
                pos_counts[pos5] += 1
                pos_len_counts[pos5][read.query_length] += 1
    except (ValueError, KeyError):
        pass
    return dict(pos_counts), {p: dict(lc) for p, lc in pos_len_counts.items()}


def calc_axtell(pos_counts, cl_start, cl_end, phase_size):
    """
    Returns best score across all offsets, or None if n <= 2 for all offsets.
    Uses all reads regardless of size.
    """
    if not pos_counts:
        return None
    best = None
    for offset in range(phase_size):
        n_max  = (cl_end - cl_start) // phase_size + 2
        phased = [cl_start + offset + k * phase_size
                  for k in range(n_max)
                  if cl_start + offset + k * phase_size <= cl_end]
        if not phased:
            continue
        phased_set = set(phased)
        pi  = sum(pos_counts.get(p, 0) for p in phased)
        u   = sum(v for p, v in pos_counts.items() if p not in phased_set)
        n   = sum(1 for p in phased if pos_counts.get(p, 0) > 0)
        if n <= 2:
            continue
        score = (n - 2) * math.log(1 + 10 * pi / (1 + u))
        if best is None or score > best:
            best = score
    return best
# ──────────────────────────────────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam',   required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--end',   type=int, required=True)
    ap.add_argument('--phase_sizes', nargs='+', type=int, default=[21, 22, 23, 24])
    ap.add_argument('--strand', default='.')
    a = ap.parse_args()

    try:
        bam = pysam.AlignmentFile(a.bam, 'rb')
    except Exception:
        print("phase=NA DC=NA n=0"); return

    s, e = min(a.start, a.end), max(a.start, a.end)
    pos_counts, pos_len_counts = reads_in_cluster(bam, a.chrom, s, e, a.strand)
    n = sum(pos_counts.values())
    if n == 0:
        print("phase=0 DC=NA n=0"); return

    # dominant DicerCall = most common 5'-anchored read length
    lenhist = Counter()
    for _p, lc in pos_len_counts.items():
        for L, c in lc.items():
            lenhist[L] += c
    dc = lenhist.most_common(1)[0][0] if lenhist else 0

    # phase = Axtell at the DC spacing if DC is a phasing size, else best over sizes
    sizes = [dc] if dc in a.phase_sizes else a.phase_sizes
    best = None
    for ps in sizes:
        v = calc_axtell(pos_counts, s, e, ps)
        if v is not None and (best is None or v > best):
            best = v
    phase = best if best is not None else 0.0
    print(f"phase={phase:.3f} DC={dc} n={n}")


if __name__ == '__main__':
    main()
