# Data dictionary

Column names use `GHB732_LR` for the low-rancidity variety GHB-732 and `GHB538_HR` for the high-rancidity
variety GHB-538.

## manuscript_data/Table1_variant_discovery.tsv
Raw and hard-filtered SNP/INDEL counts, retention and classified impacted-gene entries (manuscript Table 1).

## manuscript_data/Table2_variant_and_candidate_overlap.tsv
Record-level SNP and INDEL overlap based on exact CHROM:POS:REF:ALT identity of filtered VCF records, and
candidate-gene list membership based on exact gene-ID identity (manuscript Table 2; Figure 2).
Generated with `scripts/01_variant_record_level_overlap.sh`.

## manuscript_data/Table3_gene_level_impact_and_annotation.tsv
Gene-level impact classes of all impacted-gene entries, functional-annotation rate and annotation support of the
final candidates (manuscript Table 3).

## manuscript_data/Table4_candidate_impact_classes.tsv
Impact classes of the final putative rancidity-related candidates (manuscript Table 4).

## manuscript_data/Additional_file_6_flagstat_summary.tsv
SAMtools flagstat-based alignment quality-control summary of the final duplicate-processed BAM files
(Additional file 6).

## allele_level/ and joint_genotyping/
See `allele_level/README.md` and `joint_genotyping/README_joint_genotyping.md` for file-by-file descriptions.
