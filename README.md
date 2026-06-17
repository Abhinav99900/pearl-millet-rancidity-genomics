# Pearl millet rancidity comparative genomics

This repository contains processed summary tables, candidate-gene lists, plotting inputs, supplementary files, figures, and analysis scripts supporting the study:

**Genome-wide comparative variant analysis reveals divergent candidate gene architectures associated with flour rancidity in pearl millet**

## Study overview

This study compares two contrasting pearl millet varieties:

- **GHB-732 (LR):** low-rancid variety
- **GHB-538 (HR):** high-rancid variety

Both varieties were analysed using the same reference-guided SNP/INDEL discovery, variant-effect prediction, functional annotation, and rancidity-focused prioritization workflow.

The candidate genes reported in this repository should be treated as validation-ready candidates, not confirmed causal loci.

## Repository contents

- `manuscript_data/`: summary tables used in the manuscript
- `candidate_lists/`: final, common, and variety-unique candidate-gene lists
- `figures/`: main manuscript figure files in PNG format
- `figures_pdf/`: main manuscript figure files in PDF format
- `supplementary/`: supplementary method, supplementary tables, and supplementary figures
- `heatmap_inputs/`: functional-category by impact-class matrix used for the heatmap
- `scripts/`: analysis and reproducibility scripts
- `docs/`: workflow notes, data dictionary, and candidate-selection documentation

## Raw sequencing data

Large raw and intermediate files such as FASTQ, BAM, raw VCF, full SnpEff VCF, reference genomes, and annotation databases are not included in this GitHub repository.

Raw sequencing reads have been submitted separately to IBDC/INDA and are scheduled for public release on 01-Jun-2027. Accession details will be added after public release or once they are available for citation.

## Key reproducibility files

- `manuscript_data/variant_counts_Table1.tsv`
- `manuscript_data/variant_overlap_Table4.tsv`
- `manuscript_data/impact_class_counts_Table3.tsv`
- `manuscript_data/chromosome_distribution_Table5.tsv`
- `manuscript_data/top_candidates_Table6.tsv`
- `manuscript_data/flagstat_summary_SuppTableS5.tsv`
- `candidate_lists/COMMON_ALL_SCORED_NO_ANNOTATION_SUPPORT.tsv`
- `candidate_lists/GHB538_UNIQUE_ALL_SCORED_NO_ANNOTATION_SUPPORT.tsv`
- `candidate_lists/GHB732_UNIQUE_ALL_SCORED_NO_ANNOTATION_SUPPORT.tsv`
- `candidate_lists/TABLE6_READY_TOP_CANDIDATES_NO_ANNOTATION_SUPPORT.tsv`
- `scripts/01_variant_record_level_overlap.sh`
- `scripts/select_top_candidates.py`

## Citation

If you use this repository, please cite the associated manuscript and the archived Zenodo release.

Zenodo DOI: 10.5281/zenodo.20135136

A newer version DOI will be added after the next Zenodo version is released.
