#!/usr/bin/env bash
set -euo pipefail

# Record-level SNP/INDEL overlap using exact CHROM:POS:REF:ALT identity
# Edit these paths before running.

VCF732="/media/abhinavsingh/EXTERNAL_USB/732 result/03_variants_SO_14026_A1/SO_14026_A1.QD18.SNPbase_INDELrelax.vcf.gz"
VCF538="/media/abhinavsingh/EXTERNAL_USB/538 result/538_project/out/vcf/538.QD18.SNPbase_INDELrelax.vcf.gz"

OUTDIR="variant_record_overlap_output"
mkdir -p "$OUTDIR"

bcftools view -v snps -H "$VCF732" | awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$5}' | sort -u > "$OUTDIR/732.SNP.record.keys.tsv"
bcftools view -v snps -H "$VCF538" | awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$5}' | sort -u > "$OUTDIR/538.SNP.record.keys.tsv"

bcftools view -v indels -H "$VCF732" | awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$5}' | sort -u > "$OUTDIR/732.INDEL.record.keys.tsv"
bcftools view -v indels -H "$VCF538" | awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$5}' | sort -u > "$OUTDIR/538.INDEL.record.keys.tsv"

SNP_COMMON=$(comm -12 "$OUTDIR/732.SNP.record.keys.tsv" "$OUTDIR/538.SNP.record.keys.tsv" | wc -l)
SNP_732_UNIQUE=$(comm -23 "$OUTDIR/732.SNP.record.keys.tsv" "$OUTDIR/538.SNP.record.keys.tsv" | wc -l)
SNP_538_UNIQUE=$(comm -13 "$OUTDIR/732.SNP.record.keys.tsv" "$OUTDIR/538.SNP.record.keys.tsv" | wc -l)

INDEL_COMMON=$(comm -12 "$OUTDIR/732.INDEL.record.keys.tsv" "$OUTDIR/538.INDEL.record.keys.tsv" | wc -l)
INDEL_732_UNIQUE=$(comm -23 "$OUTDIR/732.INDEL.record.keys.tsv" "$OUTDIR/538.INDEL.record.keys.tsv" | wc -l)
INDEL_538_UNIQUE=$(comm -13 "$OUTDIR/732.INDEL.record.keys.tsv" "$OUTDIR/538.INDEL.record.keys.tsv" | wc -l)

cat > "$OUTDIR/variant_record_overlap_summary.tsv" <<SUMMARY
Comparison_layer	Common	Unique_to_GHB732	Unique_to_GHB538	GHB732_total_check	GHB538_total_check
SNP_records	$SNP_COMMON	$SNP_732_UNIQUE	$SNP_538_UNIQUE	$((SNP_COMMON + SNP_732_UNIQUE))	$((SNP_COMMON + SNP_538_UNIQUE))
INDEL_records	$INDEL_COMMON	$INDEL_732_UNIQUE	$INDEL_538_UNIQUE	$((INDEL_COMMON + INDEL_732_UNIQUE))	$((INDEL_COMMON + INDEL_538_UNIQUE))
SUMMARY

cat "$OUTDIR/variant_record_overlap_summary.tsv"
