#!/usr/bin/env bash
# Joint genotyping of GHB-732 and GHB-538 over the 1,370 candidate gene models (mRNA +/- 100 bp).
# Run inside WSL Ubuntu:   sudo bash run_joint_genotyping.sh   (sudo needed: files on D: are not readable by the normal user)
# Steps: (0) install Java + GATK 4.6.2.0 (same version as the original GHB-732 calls)
#        (1) HaplotypeCaller -ERC GVCF for each BAM on the gene intervals
#        (2) CombineGVCFs + GenotypeGVCFs (joint genotyping)
# Work is done in ~/joint_gt (fast Linux disk); final files are copied back to this folder on D:.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RUNUSER="${SUDO_USER:-$USER}"
WORK="$(getent passwd "$RUNUSER" | cut -d: -f6)/joint_gt"; mkdir -p "$WORK"; cd "$WORK"
REF="/mnt/d/732 result/02_align_SO_14026_A1/reference.fa"
BAM732="/mnt/d/732 result/02_align_SO_14026_A1/SO_14026_A1.dedup.bam"
BAM538="/mnt/d/538 result/538_project/work/aln/538.markdup.bam"
cp "$HERE/genes_1370.intervals" .

echo "== step 0: Java + GATK"
if ! command -v java >/dev/null; then apt-get update && apt-get install -y openjdk-17-jre-headless unzip wget; fi
if [ ! -x gatk-4.6.2.0/gatk ]; then
  wget -q --show-progress https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
  unzip -q gatk-4.6.2.0.zip
fi
GATK="python3 $WORK/gatk-4.6.2.0/gatk"   # wrapper is a python script; call python3 explicitly
$GATK --version | head -3

echo "== step 1: per-sample GVCFs on gene intervals"
for pair in "GHB732:$BAM732" "GHB538:$BAM538"; do
  s="${pair%%:*}"; bam="${pair#*:}"
  if [ ! -s "$s.g.vcf.gz.tbi" ]; then
    $GATK --java-options "-Xmx6g" HaplotypeCaller -R "$REF" -I "$bam" -L genes_1370.intervals -ip 100 \
      -ERC GVCF --native-pair-hmm-threads 4 -O "$s.g.vcf.gz" 2> "$s.hc.log"
  fi
  echo "   $s done"
done

echo "== step 2: joint genotyping"
$GATK --java-options "-Xmx6g" CombineGVCFs -R "$REF" -V GHB732.g.vcf.gz -V GHB538.g.vcf.gz -L genes_1370.intervals -ip 100 -O combined.g.vcf.gz 2> combine.log
$GATK --java-options "-Xmx6g" GenotypeGVCFs -R "$REF" -V combined.g.vcf.gz -L genes_1370.intervals -ip 100 -O joint_raw.vcf.gz 2> genotype.log

chown -R "$RUNUSER":"$RUNUSER" "$WORK" || true
cp joint_raw.vcf.gz joint_raw.vcf.gz.tbi GHB732.g.vcf.gz GHB538.g.vcf.gz *.log "$HERE/"
echo "== ALL DONE: $(zcat joint_raw.vcf.gz | grep -vc '^#') joint variant records"
