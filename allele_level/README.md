# Allele-level comparison and read-level validation (manuscript version 2)

Scripts and processed data for the allele-level comparison of GHB-732 (low rancidity, LR) and
GHB-538 (high rancidity, HR) across the combined set of 1,370 rancidity-related candidate genes
(Methods 5.6–5.8 of the manuscript; Additional files 3 and 5).

Run all steps from this directory. Python 3.8+ is required; `openpyxl` and `matplotlib` are needed
only for steps 12–15. SnpEff 5.3 (PearlMillet database built from the Tift 23D2B1-P1-P5 reference),
samtools ≥ 1.19 and the two final BAM files are needed only to re-run steps 02, 05 and the cross-check.

## Steps

| Step | Script | Input | Output |
|---|---|---|---|
| 01 | `01_extract_gene_regions.sh` | filtered and raw genome-wide VCFs, `uni/universe_regions.txt` | `vcf/u{732,538}.{filt,raw}.vcf` (provided, gzipped) |
| 02 | `02_run_snpeff.sh` | `vcf/*.vcf` | `vcf/*.ann.vcf` |
| 03 | `03_allele_table.py` | annotated VCFs, `vcf/u538.gvcf_blocks_at_732_calls.tsv.gz` | `uni/Coding_alleles.tsv` |
| 04 | `04_make_sites.py` | `uni/Coding_alleles.tsv` | `bam/sites_{732,538}.tsv` |
| 05 | `05_bamcheck.py` | BAM + BAI, site list | `bam/self732.out`, `bam/x538.out`, `bam/self538.out`, `bam/x732.out` |
| 06 | `06_validate_reads.py` | allele table + read counts | `uni/coding_alleles_validated.tsv` |
| 07 | `07_gene_classes.py` | validated alleles, `uni/Gene_classes.tsv`, `uni/af3_assign.tsv` | `uni/final_gene_classes.tsv` |
| 08 | `08_select_panel.py` | `uni/panel184.json`, gene classes, validated alleles | `uni/ranking2.json` |
| 09–11 | primer scripts (09 mask, 10 Primer-BLAST screen, 10b PER52 local candidates, 11 PER52 uniqueness) | `uni/cds/`, `vcf/*.filt.vcf.gz`, Primer-BLAST output in `primers/` | masked CDS and specificity screens |
| 12–15 | table and figure scripts | outputs above | `out/Additional_file_{3,5}.xlsx`, `figs/` |
| 16 | `16_threshold_sensitivity.py` | annotated VCFs and all inputs above | `uni/threshold_sensitivity.tsv` (re-runs 03 and 06–08 with alternative paralog and depth thresholds) |

Step 05 usage: `python3 05_bamcheck.py <bam> <bai> <sites> <out>`; run each site list against both
varieties' BAM files (own variety → `self*.out`, other variety → `x*.out`).

## Thresholds
- Reliable genotype: DP ≥ 10, GQ ≥ 20; allele fraction 0.25–0.75 (heterozygous) or ≥ 0.85 (homozygous).
- Paralog-like cluster: ≥ 4 heterozygous coding calls (SnpEff HIGH/MODERATE/LOW) with allele fraction < 0.35 within 100 bp in one gene.
- Read level: MAPQ ≥ 20; unmapped, secondary, QC-failed and duplicate reads excluded. Variety-specific when ≥ 15% of own informative reads are alternate-like, the other variety has ≥ 10 informative reads with < 5% alternate-like and < 25% non-reference reads, and neither BAM exceeds twice its median depth at the site.
- Complex locus: any paralog-like cluster, or gene median coding-variant depth > 2× the variety median.

## Threshold sensitivity
The thresholds were set a priori. Step 16 repeats the classification and panel selection with one threshold changed at a time (paralog cluster ≥ 3 or ≥ 5 calls, allele fraction < 0.30 or < 0.40, 50- or 200-bp window; complex-locus depth > 1.5× or > 3× the median). MAGL was selected in 9/9 settings, SDP1 in 8/9 and the PER52-like peroxidase in 7/9. Alleles that become reliable only under an alternative setting have no read-level data and are counted as not validated.

## Joint genotyping
The same analysis was repeated on joint calls of both varieties (`../joint_genotyping/`). The numbers of validated
variety-specific alleles were unchanged (418 and 257), SDP1 and MAGL were again selected, and the PER52-like
peroxidase was replaced by Pgl_GLEAN_10011170 (PER66-like) because joint calling placed it in a paralog-like cluster.

## Check of the custom BAM reader
`samtools_crosscheck/` compares read counts from `05_bamcheck.py` with `samtools mpileup` v1.19.2
(`-A -B -x -Q 0 -q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP`) at 30 random SNV sites in both BAM files:
reference and alternate counts were identical at all 60 site–BAM combinations.

## Reproducibility note
Re-running steps 02–04 from the provided VCF subsets reproduces all 11,046 rows of the published
`uni/Coding_alleles.tsv` exactly and the two site lists (1,429 and 1,715 sites), with one difference:
five start-codon (`start_lost`, c.1) alleles in three genes (Pgl_GLEAN_10001124, Pgl_GLEAN_10001963,
Pgl_GLEAN_10026293) are additionally reported. Four fail genotype QC; the fifth (GHB-538,
chr7:118913 T>A) is also called in GHB-732, so it cannot become a variety-specific allele. The
published tables in `uni/` are the versions used in the manuscript; running step 03 overwrites
`uni/Coding_alleles.tsv`.

## Data not included
FASTQ, BAM, genome-wide VCF files, the reference genome and the SnpEff database are not included
(raw reads: IBDC-INDA INRP000666; INSDC PRJEB113453). `vcf/` holds only the gene-region subsets.
