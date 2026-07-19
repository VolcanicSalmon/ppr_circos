"""
Calculate phasing scores for ShortStack clusters using two methods:

  PhaseTank (Guo et al.):
      phased_score = phased_ratio * phased_num * ln(phased_abun)
      Uses only reads of exactly phase_size nt.

  Axtell (ShortStack ≤3):
      PS = ln[(1 + 10 * ΣPi / (1 + ΣU))^(n-2)],  n > 2
      Uses all reads regardless of size.

Reads are fetched directly from a sorted, indexed BAM file using pysam.

Usage:
    python3 calc_phase_score.py \
        --results  Results.txt \
        --bam      merged_alignments.bam \
        --output   phasing_scores.tsv \
        --phase_sizes 21 22 23
"""

import math
import sys
import argparse
from collections import defaultdict

import pysam
import pandas as pd


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('--results',       required=True)
    p.add_argument('--bam',           required=True)
    p.add_argument('--output',        required=True)
    p.add_argument('--phase_sizes',   nargs='+', type=int, default=[21, 24])
    p.add_argument('--drift',         type=int, default=2)
    p.add_argument('--min_phased_num',type=int, default=4)
    p.add_argument('--min_reads',     type=int, default=0)
    return p.parse_args()


# ── Read fetch ────────────────────────────────────────────────────────────────

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


# ── PhaseTank score ───────────────────────────────────────────────────────────

def calc_phasetank(pos_len_counts, cl_start, phase_size,
                   phased_drift=2, min_phased_num=4):
    """
    Returns (score, phased_ratio, phased_num, phased_abun) or None.
    Uses only reads of exactly phase_size nt.
    """
    sized = {p: lc[phase_size] for p, lc in pos_len_counts.items()
             if lc.get(phase_size, 0) > 0}
    if not sized:
        return None

    total = sum(sized.values())
    bin_abun, bin_pos = defaultdict(int), defaultdict(list)
    for pos, cnt in sized.items():
        b = (pos - cl_start) % phase_size
        bin_abun[b] += cnt
        bin_pos[b].append(pos)

    sorted_bins = sorted(bin_abun.items(), key=lambda x: x[1], reverse=True)
    max_bin, most = sorted_bins[0]
    if len(sorted_bins) > 1:
        more_bin, more = sorted_bins[1]
        fuzzy = abs(max_bin - more_bin)
        if fuzzy <= phased_drift or fuzzy >= phase_size - phased_drift:
            phased_ratio = (most + more) / total
        else:
            phased_ratio = most / total
    else:
        phased_ratio = most / total

    island   = 5 * phase_size
    pos_list = sorted(bin_pos[max_bin])
    best, cur = [], []
    for pos in pos_list:
        if not cur or pos - cur[0] <= island:
            cur.append(pos)
        else:
            if len(cur) > len(best): best = cur
            cur = [pos]
    if len(cur) > len(best): best = cur

    if len(best) < min_phased_num:
        return None

    phased_abun = sum(sized.get(p, 0) for p in best)
    if phased_abun <= 0:
        return None

    score = phased_ratio * len(best) * math.log(phased_abun)
    return score, phased_ratio, len(best), phased_abun


# ── Axtell score ──────────────────────────────────────────────────────────────

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


# ── Helpers ───────────────────────────────────────────────────────────────────

def phase_sizes_for(dicer_call, defaults):
    if str(dicer_call).isdigit():
        size = int(dicer_call)
        return [size] if size in defaults else []
    return defaults


def fmt(v, decimals=4):
    return round(v, decimals) if v is not None else 'NA'


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = get_args()

    results = pd.read_csv(args.results, sep='\t', header=0)
    results.columns = [c.strip() for c in results.columns]
    print(f"Clusters: {len(results)}  from {args.results}", file=sys.stderr)

    print(f"Opening BAM: {args.bam}", file=sys.stderr)
    bam = pysam.AlignmentFile(args.bam, 'rb')

    rows = []
    for _, row in results.iterrows():
        chrom     = str(row['Chrom']).strip()
        cl_start  = int(row['Start'])
        cl_end    = int(row['End'])
        dicer     = str(row['DicerCall']).strip()
        name      = str(row['Name']).strip()
        cl_strand = str(row['Strand']).strip() if 'Strand' in results.columns else '.'
        total_cl  = int(row['Reads']) if 'Reads' in results.columns else 0

        if args.min_reads > 0 and total_cl < args.min_reads:
            continue

        pos_counts, pos_len_counts = reads_in_cluster(
            bam, chrom, cl_start, cl_end, cl_strand)

        sizes = phase_sizes_for(dicer, args.phase_sizes)

        # PhaseTank — best across target sizes
        pt_best, pt_size, pt_ratio, pt_num, pt_abun = None, None, None, None, None
        for p in sizes:
            r = calc_phasetank(pos_len_counts, cl_start, p,
                               args.drift, args.min_phased_num)
            if r is not None and (pt_best is None or r[0] > pt_best):
                pt_best, pt_size = r[0], p
                _, pt_ratio, pt_num, pt_abun = r

        # Axtell — best across target sizes
        ax_best, ax_size = None, None
        for p in sizes:
            r = calc_axtell(pos_counts, cl_start, cl_end, p)
            if r is not None and (ax_best is None or r > ax_best):
                ax_best, ax_size = r, p

        rows.append({
            'Name':          name,
            'Chrom':         chrom,
            'Start':         cl_start,
            'End':           cl_end,
            'DicerCall':     dicer,
            'TotalReads':    total_cl,
            'PT_BestSize':   pt_size  if pt_best is not None else 'NA',
            'PT_PhaseScore': fmt(pt_best),
            'PT_Ratio':      fmt(pt_ratio),
            'PT_Num':        pt_num   if pt_num  is not None else 'NA',
            'PT_Abun':       pt_abun  if pt_abun is not None else 'NA',
            'AX_BestSize':   ax_size  if ax_best is not None else 'NA',
            'AX_PhaseScore': fmt(ax_best),
        })

    bam.close()

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep='\t', index=False)

    pt_scored = (df['PT_PhaseScore'] != 'NA').sum()
    ax_scored = (df['AX_PhaseScore'] != 'NA').sum()
    print(f"Done — PhaseTank: {pt_scored}/{len(df)} scored, "
          f"Axtell: {ax_scored}/{len(df)} scored → {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
