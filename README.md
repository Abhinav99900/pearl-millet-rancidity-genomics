# Pearl millet rancidity comparative genomics

This repository contains processed summary tables, candidate-gene lists, plotting inputs, and analysis scripts supporting the manuscript:

**Genome-wide comparative variant analysis reveals divergent candidate gene architectures associated with flour rancidity in pearl millet**

## Study overview

This study compares two pearl millet genotypes with contrasting flour-rancidity backgrounds:

- **GHB-732**: low-rancidity comparative background
- **GHB-538**: high-rancidity comparative background

The analysis uses reference-guided SNP/INDEL discovery, SnpEff-based variant-effect prediction, functional annotation, and multi-gate rancidity-focused prioritization to identify shared and genotype-specific putative rancidity-associated candidate genes.

## Repository contents

- `manuscript_data/`: summary tables used in the manuscript
- `candidate_lists/`: final, common, and genotype-unique candidate-gene lists
- `heatmap_inputs/`: functional-category by impact-class matrix used for the mechanistic heatmap
- `scripts/`: analysis and plotting scripts
- `figures/`: generated manuscript figure files or figure inputs
- `docs/`: supplementary workflow notes and data dictionary

## Raw sequencing data

Raw sequencing reads are not stored in this repository. They will be deposited in the NCBI Sequence Read Archive under BioProject accession:

`PRJNAxxxxxx`

## Large files

Large intermediate files such as FASTQ, BAM, raw VCF, full SnpEff VCF, reference genomes, and annotation databases are not included in this GitHub repository.

## Key reproducibility files

- `manuscript_data/variant_counts_Table1.tsv`
- `manuscript_data/variant_overlap_Table4.tsv`
- `manuscript_data/impact_class_counts_Table3.tsv`
- `manuscript_data/chromosome_distribution_Table5.tsv`
- `manuscript_data/top_candidates_Table6.tsv`
- `manuscript_data/flagstat_summary_SuppTableS5.tsv`
- `scripts/01_variant_record_level_overlap.sh`

## Citation

If you use this repository, please cite the associated manuscript and the archived Zenodo release DOI:

`DOI: to be added after Zenodo release`
