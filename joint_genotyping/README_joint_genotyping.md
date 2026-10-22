# Joint-genotyping robustness check (manuscript v20, Methods 5.6 and Results 2.5)

Both final BAM files were genotyped jointly over the 1,370 candidate gene models and the allele-level
analysis was repeated with unchanged rules.

| Step | File | What it does |
|---|---|---|
| 1 | `run_joint_genotyping.sh` | GATK 4.6.2.0 HaplotypeCaller `-ERC GVCF` on both BAMs (`-L genes_1370.intervals -ip 100`), CombineGVCFs, GenotypeGVCFs → `joint_raw.vcf.gz` |
| 2 | `scripts/hardfilter.py`, `scripts/J01_filter_regions.py` | Original SNP/INDEL hard filters (reproduce the published filtered gene-region VCFs exactly: 36,938 and 34,886 records) and restriction to the gene models |
| 3 | SnpEff 5.3 (`PearlMillet`) | `java -jar snpEff.jar PearlMillet vcf/joint.regions.vcf > vcf/joint.regions.ann.vcf` |
| 4 | `scripts/J03_allele_table_joint.py` | Allele table; absence in the other variety taken from its joint genotype (called, DP ≥ 10) |
| 5 | `04_make_sites.py`, `05_bamcheck.py`, `run_readcounts_new_sites.sh` | Read counts for the 53 new sites (12 published sites re-counted as controls: identical) |
| 6 | `scripts/06_validate_reads.py`, `07_gene_classes.py`, `08_select_panel.py` | Read-level validation, gene classes, panel selection |
| 7 | `scripts/J16_sensitivity_joint.py` | The nine threshold settings on the joint calls |

Input provided: `joint_raw.vcf.gz` (output of step 1, both varieties, gene-model intervals), and `bam_new_sites/` (site lists and read counts for the 53 new sites; append each `*_new.out` to the matching file in `allele_level/bam/` before step 6), so steps 2–7 can be re-run without the BAM files.

Results: `JOINT_GENOTYPING_SUMMARY.txt`, `Coding_alleles.tsv`, `coding_alleles_validated.tsv`,
`final_gene_classes.tsv`, `threshold_sensitivity_joint.tsv`; summary tables in Additional file 5 (sheet "Joint genotyping").
