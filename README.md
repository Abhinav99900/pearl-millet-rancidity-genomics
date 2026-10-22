# Pearl millet rancidity comparative genomics

Processed tables, candidate-gene lists, figures, supplementary files and analysis scripts supporting the study:

**Comparative whole-genome resequencing identifies candidate genes for flour rancidity in pearl millet**

## Study overview

The study compares two pearl millet varieties with contrasting flour rancidity:

- **GHB-732 (LR):** low-rancidity variety
- **GHB-538 (HR):** high-rancidity variety

Both varieties were analysed with the same reference-guided SNP/INDEL discovery, variant-effect prediction,
functional annotation and rancidity-focused prioritization workflow. The two candidate lists (585 and 1,108
genes) were then merged into one set of 1,370 genes and compared allele by allele, with every putative
difference checked in the sequencing reads of both varieties. The validated differences were used to select
three priority candidates (SDP1, a monoacylglycerol lipase and a PER52-like class III peroxidase), for which
variant-aware RT-qPCR primers were designed. SDP1 and the monoacylglycerol lipase are robust to the analysis
thresholds and to joint genotyping; for the peroxidase, the class III peroxidase family rather than the individual
gene is the candidate. Candidates are hypotheses for functional studies, not confirmed causal loci.

## What is new in version 2 (manuscript revision)

- `allele_level/`: allele-level comparison, read-level validation, gene classes, panel selection, primer design,
  samtools cross-check, and the scripts and processed inputs to reproduce them (see `allele_level/README.md`).
- `manuscript_v2/`: the six Additional files and eight figures of the current manuscript.

The earlier list-level outputs (common/variety-unique "top candidate" selection, e.g. `candidate_lists/*TOP*`,
`manuscript_data/top_candidates_Table6.tsv`, `scripts/select_top_candidates.py`) are kept for transparency but are
superseded: the current validation panel is selected from allele-level, read-validated differences.

## What is new in version 2.1

- `allele_level/16_threshold_sensitivity.py` and `allele_level/uni/threshold_sensitivity.tsv`: classification and
  panel selection repeated under nine alternative paralog-cluster and complex-locus depth thresholds.
- `joint_genotyping/`: robustness check in which both varieties were genotyped jointly (GATK HaplotypeCaller GVCF
  mode, CombineGVCFs, GenotypeGVCFs) over the 1,370 gene models and the allele-level analysis was repeated; validated
  alleles were unchanged (418 and 257; 672 of 675 identical) and 97.8% of genes kept their class
  (see `joint_genotyping/README_joint_genotyping.md`).
- `manuscript_v2/`: updated Additional file 5 (sheets "Threshold sensitivity" and "Joint genotyping") and Figure 8.

## Variant filtering thresholds

The genome-wide VCFs were hard-filtered with `bcftools view -i` (Additional file 1;
`docs/Supplementary_Method_S1_Detailed_Bioinformatics_Workflow.docx`):

- SNPs: QD ≥ 1.8, FS ≤ 80, SOR ≤ 4.0, MQ ≥ 40, MQRankSum ≥ −15 (or missing), ReadPosRankSum ≥ −10 (or missing)
- INDELs: QD ≥ 1.8, FS ≤ 200, SOR ≤ 10, ReadPosRankSum ≥ −20 (or missing)

`joint_genotyping/scripts/hardfilter.py` re-implements these filters and reproduces the filtered gene-region VCFs
exactly (36,938 and 34,886 records).

## Repository contents

- `allele_level/`: allele-level comparison and validation (version 2; threshold sensitivity added in 2.1)
- `joint_genotyping/`: joint-genotyping robustness check (version 2.1)
- `manuscript_v2/`: Additional files 1–6 and Figures 1–8 of the current manuscript (version 2)
- `manuscript_data/`: summary tables (variant counts, overlaps, impact classes, flagstat)
- `candidate_lists/`: final candidate lists of each variety and list-level comparisons
- `figures/`, `figures_pdf/`: figure files of the earlier version
- `supplementary/`: supplementary files of the earlier version
- `heatmap_inputs/`: functional-category by impact-class matrix used for the heatmap
- `scripts/`: variant-level analysis and reproducibility scripts
- `docs/`: workflow notes and data dictionary

## Raw sequencing data

FASTQ, BAM, genome-wide VCF files, reference genomes and annotation databases are not included. Raw reads have
been deposited in IBDC-INDA (accession INRP000666) and registered under INSDC Study/BioProject PRJEB113453;
they will be publicly released on 1 June 2027 or upon publication, whichever is earlier.

## Citation

Please cite the associated manuscript and the archived Zenodo release (concept DOI: 10.5281/zenodo.20135136,
which always resolves to the latest version).
