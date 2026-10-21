#!/usr/bin/env bash
# -A: keep anomalous (not properly paired) read pairs, as the custom reader does
# Fast version: uses the BAM index to jump to each site (-r) instead of scanning the whole BAM.
# Same filters as the custom reader: MAPQ >= 20; unmapped, secondary, QC-fail and duplicate reads excluded;
# no base-quality filter, BAQ off, overlap correction off, no depth cap.
# Usage: bash run_samtools_check_fast.sh <GHB-732 BAM> <GHB-538 BAM>
set -euo pipefail
BAM732="$1"; BAM538="$2"
for pair in "732:$BAM732" "538:$BAM538"; do
  v="${pair%%:*}"; bam="${pair#*:}"
  : > mpileupA_${v}.txt
  while IFS=$'\t' read -r chr pos; do
    samtools mpileup -A -B -x -Q 0 -q 20 -d 0 --ff UNMAP,SECONDARY,QCFAIL,DUP \
      -r "${chr}:${pos}-${pos}" "$bam" 2>/dev/null >> mpileupA_${v}.txt
  done < sites30.txt
  echo "done: mpileupA_${v}.txt ($(wc -l < mpileupA_${v}.txt) sites)"
done
samtools --version | head -1 > samtools_version.txt
