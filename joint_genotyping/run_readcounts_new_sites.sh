#!/usr/bin/env bash
# Read counts (step 05, same custom reader as the published analysis) at the coding-allele sites
# that are new in the joint-genotyping table. The last 3 lines of each site list are already-published
# sites, included as a consistency control.   Run:  sudo bash run_readcounts_new_sites.sh
set -euo pipefail
cd "$(dirname "$0")"
B732="/mnt/d/732 result/02_align_SO_14026_A1/SO_14026_A1.dedup.bam"
B538="/mnt/d/538 result/538_project/work/aln/538.markdup.bam"
python3 05_bamcheck.py "$B732" "$B732.bai"                                  new_sites_732.tsv self732_new.out
python3 05_bamcheck.py "$B538" "/mnt/d/538 result/538_project/work/aln/538.markdup.bai" new_sites_732.tsv x538_new.out
python3 05_bamcheck.py "$B538" "/mnt/d/538 result/538_project/work/aln/538.markdup.bai" new_sites_538.tsv self538_new.out
python3 05_bamcheck.py "$B732" "$B732.bai"                                  new_sites_538.tsv x732_new.out
wc -l *_new.out
echo "== READ COUNTS DONE"
