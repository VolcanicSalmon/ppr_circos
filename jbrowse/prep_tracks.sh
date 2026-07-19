#!/bin/bash
# Prepare the local BED tracks referenced by config.json into jbrowse/tracks/.
# Small BEDs -> copied (BedAdapter). Large BEDs -> sort+bgzip+tabix (BedTabixAdapter).
# sRNA-read tracks are derived from the dataset BAM via bedtools bamtobed.
# Output names match the config URIs. Genomes/GFFs come from the remote hosts via serve.py's proxy.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p tracks
T=/Users/vef25hok/Downloads/ppr_srna_local/trackplot
tabix_bed () { LC_ALL=C sort -k1,1 -k2,2n "$1" | bgzip > "tracks/$2.sorted.bed.gz"; tabix -p bed "tracks/$2.sorted.bed.gz"; }

echo "[copy small BEDs -> BedAdapter]"
cp "$T/bin/data/rh/genome/rh_ppregion_transcript_regions.bed" tracks/RH_PPR_transcript_regions.bed
cp "$T/bin/data/nb/genome/nb_ppregion_transcript_regions.bed" tracks/nb_ppregion_transcript_regions.bed

cp "$T/bin/data/spim/genome/sp_ppregion_transcript_regions.bed" tracks/sp_ppregion_transcript_regions.bed
cp "$T/bin/data/aras/genome/ara_ppregion_transcript_regions.bed" tracks/ara_ppregion_transcript_regions.bed

echo "[copy .fai -> pair with the remote .fa streamed via proxy]"
cp "$T/bin/data/nb/genome/NbeHZ1_genome_1.0.fa.fai"  tracks/NbeHZ1_genome_1.0.fa.fai
cp "$T/bin/data/spim/genome/LA2093_genome_v1.5.fa.fai" tracks/LA2093_genome_v1.5.fa.fai

echo "[slice-site BEDs -> tabix]"
tabix_bed "$T/bin/data/rh/slice/slice_bed/mir_and_baldrich_uncoll_condensedmin20_21_23.fa_gstar_421ft2.bed" rh_slice_sites
tabix_bed "$T/bin/data/nb/slice/slice_bed/mir_and_merged_alignmentsmin6_21_23_gstar_421ft.bed"            nb_slice_sites
tabix_bed "$T/bin/data/spim/slice/slice_bed/mir_and_merged_alignmentsmin6_21_23_gstar_421ft.bed"          spim_slice_sites
tabix_bed "$T/bin/data/aras/slice/slice_bed/mir_and_baldrich_uncoll_condensedmin30_21_23_GSTAr_221ft.bed" aras_slice_sites

echo "[sRNA reads: bamtobed the dataset BAM -> tabix]"
bedtools bamtobed -i "$T/testdata/rh/baldrich_uncoll_condensed.bam" > /tmp/rh_reads.bed 2>/dev/null || \
  cp "$T/testdata/rh/baldrich_uncoll_condensed.bed" /tmp/rh_reads.bed
tabix_bed /tmp/rh_reads.bed baldrich_uncoll_condensed
bedtools bamtobed -i "$T/bin/data/nb/bam/merged_alignments.bam" > /tmp/nb_reads.bed
tabix_bed /tmp/nb_reads.bed nb_merged_alignments
bedtools bamtobed -i "$T/bin/data/spim/bam/merged_alignments.bam" > /tmp/spim_reads.bed
tabix_bed /tmp/spim_reads.bed spim_merged_alignments
bedtools bamtobed -i "$T/bin/data/aras/bam/baldrich_uncoll_condensed.bam" > /tmp/aras_reads.bed
tabix_bed /tmp/aras_reads.bed aras_reads

echo "done -> tracks/"; ls -la tracks/
