#!/bin/bash
# Annotate the four gene-region VCF subsets with SnpEff 5.3 (database: PearlMillet,
# built from the Tift 23D2B1-P1-P5 reference and gene models).
# Usage: bash 02_run_snpeff.sh /path/to/snpEff        (directory containing snpEff.jar and snpEff.config)
set -euo pipefail
SNPEFF="$1"
for s in u732.filt u732.raw u538.filt u538.raw; do
  [ -f vcf/$s.vcf ] || gunzip -k vcf/$s.vcf.gz
  java -Xmx4g -jar "$SNPEFF/snpEff.jar" -c "$SNPEFF/snpEff.config" -noStats -noLog PearlMillet vcf/$s.vcf > vcf/$s.ann.vcf
done
