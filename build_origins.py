#!/usr/bin/env python3
"""
build_origins.py — make trigger_origins.bed from a bscore.tsv + the sRNA alignment BED.

Companion to pval_dist.py. The distinct triggers in the bscore's `trigger` column are
looked up in the alignment BED whose name field (col4) IS the trigger id (RH:
baldrich_uncoll_condensed.bed; DM: merged_alignments.bed — NOT the *_on_cluster_flank
subset). Matching rows are written as:  chrom  start  end  trigger  strand

This replaces the previous hand-run awk. The circos app reads the output for the chords.

usage:
  python3 build_origins.py \
     --bscore     rh/slice/rh_bscore.tsv \
     --alignments rh/slice/slice_bed/baldrich_uncoll_condensed.bed \
     --out        rh/slice/slice_bed/trigger_origins.bed
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bscore', required=True, help='pval_dist.py bscore TSV')
    ap.add_argument('--alignments', required=True,
                    help='sRNA alignment BED; name (col4) = trigger id')
    ap.add_argument('--out', required=True, help='output trigger_origins.bed')
    args = ap.parse_args()

    # 1) distinct triggers from the bscore 'trigger' column
    triggers = set()
    with open(args.bscore) as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            tc = header.index('trigger')
        except ValueError:
            sys.exit(f"ERROR: no 'trigger' column in {args.bscore}")
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) > tc:
                triggers.add(p[tc])

    # 2) stream the (large) alignment BED once, keep rows whose name is a trigger
    n = 0
    with open(args.alignments) as fin, open(args.out, 'w') as fout:
        for line in fin:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 6 and p[3] in triggers:
                fout.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[5]}\n")
                n += 1

    sys.stderr.write(f"{len(triggers)} triggers -> {n} origin rows -> {args.out}\n")


if __name__ == '__main__':
    main()
