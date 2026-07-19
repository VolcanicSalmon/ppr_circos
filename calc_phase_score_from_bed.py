"""
Calculate PhaseTank + Axtell phase scores for arbitrary regions in a BED file,
using the sRNA read positions and abundances from a condensed BAM.

Same scoring formulas as calc_phase_score.py, but the region set is defined
by the input BED instead of ShortStack Results.txt.

Region name in the BED should be a unique identifier — used as 'Name' in output.
Reads are weighted by the XW:i: tag in the BAM (ShortStack condensed reads)
so a Cd_1000 alignment counts as 1000 reads.

Usage:
    python3 calc_phase_score_from_bed.py \\
        --bed        regions.bed \\
        --bam        merged_alignments.bam \\
        --output     region_phase_scores.tsv \\
        --phase_sizes 21 22 23 24 \\
        [--min_reads 20]
"""

import math
import sys
import argparse
from collections import defaultdict

import pysam
import pandas as pd


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('--bed',           required=True,
                   help='BED with regions to score (chrom start end name score strand)')
    p.add_argument('--bam',           required=True,
                   help='Sorted, indexed condensed BAM (with XW:i: tags)')
    p.add_argument('--output',        required=True)
    p.add_argument('--phase_sizes',   nargs='+', type=int, default=[21, 22, 23, 24])
    p.add_argument('--drift',         type=int, default=2)
    p.add_argument('--min_phased_num',type=int, default=4)
    p.add_argument('--min_reads',     type=int, default=0,
                   help='Skip regions with fewer than this many weighted reads (default: 0)')
    p.add_argument('--no_normalize', action='store_true',
                   help='Use raw Axtell score (scales with window length). '
                        'Default: length-normalized to remove window-size bias.')
    return p.parse_args()


def get_read_weight(read):
    """Return the XW:i: tag value or 1 if absent."""
    try:
        return read.get_tag('XW')
    except KeyError:
        return 1


def reads_in_region(bam, chrom, r_start, r_end, strand='.'):
    """
    Returns (pos_counts, pos_len_counts, total_reads) for reads whose 5' end
    is in [r_start, r_end], weighted by XW:i: tag.
    """
    pos_counts     = defaultdict(int)
    pos_len_counts = defaultdict(lambda: defaultdict(int))
    total = 0
    try:
        for read in bam.fetch(chrom, r_start, r_end + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if strand != '.' and (('-' if read.is_reverse else '+') != strand):
                continue
            pos5 = (read.reference_end - 1) if read.is_reverse else read.reference_start
            if r_start <= pos5 <= r_end:
                w = get_read_weight(read)
                pos_counts[pos5] += w
                pos_len_counts[pos5][read.query_length] += w
                total += w
    except (ValueError, KeyError):
        pass
    return dict(pos_counts), {p: dict(lc) for p, lc in pos_len_counts.items()}, total


def calc_phasetank(pos_len_counts, r_start, phase_size,
                   phased_drift=2, min_phased_num=4):
    sized = {p: lc[phase_size] for p, lc in pos_len_counts.items()
             if lc.get(phase_size, 0) > 0}
    if not sized:
        return None

    total = sum(sized.values())
    bin_abun, bin_pos = defaultdict(int), defaultdict(list)
    for pos, cnt in sized.items():
        b = (pos - r_start) % phase_size
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


def calc_axtell(pos_counts, r_start, r_end, phase_size, normalize=True):
    """Axtell phase score.

    When normalize=True, divides (n-2) by (n_possible-2) so the score is
    length-independent — a length-500bp and length-2500bp window with the
    same *fraction* of phased positions occupied and the same phased/unphased
    read ratio produce the same score.

    When normalize=False, uses the raw Axtell formula (scales linearly with
    window length via the (n-2) factor). Use for compatibility with the
    ShortStack-cluster path.
    """
    if not pos_counts:
        return None
    best = None
    for offset in range(phase_size):
        n_max  = (r_end - r_start) // phase_size + 2
        phased = [r_start + offset + k * phase_size
                  for k in range(n_max)
                  if r_start + offset + k * phase_size <= r_end]
        if not phased:
            continue
        phased_set = set(phased)
        pi  = sum(pos_counts.get(p, 0) for p in phased)
        u   = sum(v for p, v in pos_counts.items() if p not in phased_set)
        n   = sum(1 for p in phased if pos_counts.get(p, 0) > 0)
        if n <= 2:
            continue
        if normalize:
            denom = max(len(phased) - 2, 1)
            score = ((n - 2) / denom) * math.log(1 + 10 * pi / (1 + u))
        else:
            score = (n - 2) * math.log(1 + 10 * pi / (1 + u))
        if best is None or score > best:
            best = score
    return best


def fmt(v, decimals=4):
    return round(v, decimals) if v is not None else 'NA'


def main():
    args = get_args()

    print(f"Opening BAM: {args.bam}", file=sys.stderr)
    bam = pysam.AlignmentFile(args.bam, 'rb')

    rows = []
    n_regions = 0
    with open(args.bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            n_regions += 1

            chrom   = cols[0].strip()
            r_start = int(cols[1])
            r_end   = int(cols[2]) - 1     # BED is half-open; make inclusive
            name    = cols[3].strip() if len(cols) > 3 else f'{chrom}:{r_start}-{r_end+1}'
            strand  = cols[5].strip() if len(cols) > 5 else '.'
            if strand not in ('+', '-'):
                strand = '.'

            pos_counts, pos_len_counts, total_r = reads_in_region(
                bam, chrom, r_start, r_end, strand)

            if args.min_reads > 0 and total_r < args.min_reads:
                continue

            pt_best, pt_size, pt_ratio, pt_num, pt_abun = None, None, None, None, None
            for p in args.phase_sizes:
                r = calc_phasetank(pos_len_counts, r_start, p,
                                   args.drift, args.min_phased_num)
                if r is not None and (pt_best is None or r[0] > pt_best):
                    pt_best, pt_size = r[0], p
                    _, pt_ratio, pt_num, pt_abun = r

            ax_best, ax_size = None, None
            for p in args.phase_sizes:
                r = calc_axtell(pos_counts, r_start, r_end, p, normalize=not args.no_normalize)
                if r is not None and (ax_best is None or r > ax_best):
                    ax_best, ax_size = r, p

            rows.append({
                'Name':          name,
                'Chrom':         chrom,
                'Start':         r_start,
                'End':           r_end + 1,
                'Strand':        strand,
                'TotalReads':    total_r,
                'PT_BestSize':   pt_size if pt_best is not None else 'NA',
                'PT_PhaseScore': fmt(pt_best),
                'PT_Ratio':      fmt(pt_ratio),
                'PT_Num':        pt_num  if pt_num  is not None else 'NA',
                'PT_Abun':       pt_abun if pt_abun is not None else 'NA',
                'AX_BestSize':   ax_size if ax_best is not None else 'NA',
                'AX_PhaseScore': fmt(ax_best),
            })

    bam.close()

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep='\t', index=False)

    pt_scored = (df['PT_PhaseScore'] != 'NA').sum()
    ax_scored = (df['AX_PhaseScore'] != 'NA').sum()
    print(f"Done — {n_regions} regions read, {len(df)} scored | "
          f"PhaseTank: {pt_scored}, Axtell: {ax_scored} → {args.output}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
