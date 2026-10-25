# Pearl millet rancidity comparative genomics

Data, figures, supplementary files and analysis scripts supporting the manuscript:

**Comparative whole-genome resequencing identifies candidate genes for flour rancidity in pearl millet**

## Study overview

The study compares two pearl millet varieties with contrasting flour rancidity:

- **GHB-732 (LR):** low-rancidity variety
- **GHB-538 (HR):** high-rancidity variety

Both varieties were analysed with the same reference-guided SNP/INDEL discovery, variant-effect prediction,
functional annotation and rancidity-focused prioritization workflow. The two candidate lists (585 and 1,108
genes) were merged into one set of 1,370 genes and compared allele by allele, with every putative difference
checked in the sequencing reads of both varieties. This confirmed 418 GHB-732-specific and 257 GHB-538-specific
coding alleles (unchanged with joint genotyping) and set aside 283 complex loci. The validated differences were
used to select three priority candidates (SDP1, a monoacylglycerol lipase and a PER52-like class III peroxidase),
for which variant-aware RT-qPCR primers were designed. Candidates are hypotheses for functional studies, not
confirmed causal loci.

## Workflow (start to end)

The complete workflow, with all commands, is described in
`manuscript_v2/additional_files/Additional_file_1.docx` (Supplementary Method S1):

1. **Read QC and trimming:** FastQC, MultiQC and fastp
2. **Alignment:** BWA-MEM2 to the Tift 23D2B1-P1-P5 reference; SAMtools sorting and indexing; Picard
   MarkDuplicates; SAMtools flagstat (`manuscript_data/Additional_file_6_flagstat_summary.tsv`)
3. **Variant calling and filtering:** GATK HaplotypeCaller (GVCF mode) and GenotypeGVCFs; separate SNP and
   INDEL hard filters with bcftools (see *Variant filtering thresholds*)
4. **Variant-effect prediction:** SnpEff, keeping the strongest predicted impact per variant and per gene
5. **Functional annotation:** Swiss-Prot, TrEMBL, eggNOG, Pfam and CDD
6. **Rancidity-focused prioritization:** 585 candidates in GHB-732 and 1,108 in GHB-538 (Additional file 4)
7. **Allele-level comparison, read-level validation, gene classes, panel selection and primer design:**
   scripts 01–16 in `allele_level/` (see `allele_level/README.md`)
8. **Robustness checks:** threshold sensitivity (`allele_level/16_threshold_sensitivity.py`) and joint
   genotyping (`joint_genotyping/`)

Steps 1–6 were run on the genome-wide data, which are not included (see *Raw sequencing data*); steps 7–8 can
be re-run from the files provided here.

## Repository contents

| Folder | Contents | Manuscript |
|---|---|---|
| `manuscript_v2/figures/` | Figures 1–8 (JPG and PDF) and the graphical abstract | Figures 1–8 |
| `manuscript_v2/additional_files/` | Additional files 1–6 | Additional files 1–6 |
| `manuscript_data/` | Summary tables of variant discovery, overlap, impact classes and alignment QC | Tables 1–4, Additional file 6 |
| `allele_level/` | Allele-level comparison, read-level validation, gene classes, panel selection, primer design, samtools cross-check and threshold sensitivity (see `allele_level/README.md`) | Results 2.3–2.6; Methods 5.6–5.8; Tables 5–8; Figures 5–6; Additional files 3 and 5 |
| `joint_genotyping/` | Joint-genotyping robustness check (see `joint_genotyping/README_joint_genotyping.md`) | Results 2.5; Methods 5.6; Additional file 5 |
| `scripts/` | Record-level SNP/INDEL overlap | Table 2; Figure 2 |
| `docs/` | Data dictionary | – |

The candidate lists of each variety, with impact classes and functional categories, are in Additional file 4;
the combined set of 1,370 genes with allele-level classes is in Additional file 3.

## Variant filtering thresholds

The genome-wide VCFs were hard-filtered with `bcftools view -i` (Methods 5.3; Additional file 1):

- SNPs: QD ≥ 1.8, FS ≤ 80, SOR ≤ 4.0, MQ ≥ 40, MQRankSum ≥ −15 (or missing), ReadPosRankSum ≥ −10 (or missing)
- INDELs: QD ≥ 1.8, FS ≤ 200, SOR ≤ 10, ReadPosRankSum ≥ −20 (or missing)

`joint_genotyping/scripts/hardfilter.py` re-implements these filters and reproduces the filtered gene-region VCFs
exactly (36,938 and 34,886 records).

## Raw sequencing data

FASTQ, BAM, genome-wide VCF files, reference genomes and annotation databases are not included. Raw reads have
been deposited in IBDC-INDA (accession INRP000666) and registered under INSDC Study/BioProject PRJEB113453;
they will be publicly released on 1 June 2027 or upon publication, whichever is earlier.

## Version history

- **2.2:** repository aligned with the final manuscript. Figure 1 (workflow) and Figure 6 (Circos view) updated
  and the graphical abstract added; `manuscript_data/` regenerated with the current table numbering; folders of
  the earlier list-level version (`figures/`, `figures_pdf/`, `supplementary/`, `candidate_lists/`,
  `heatmap_inputs/`, old `docs/` files and `scripts/select_top_candidates.py`) removed. They remain available
  in release 2.1 and earlier Zenodo versions.
- **2.1:** threshold-sensitivity analysis (`allele_level/16_threshold_sensitivity.py`) and joint-genotyping
  robustness check (`joint_genotyping/`); validated alleles were unchanged (418 and 257; 672 of 675 identical)
  and 97.8% of genes kept their class.
- **2.0:** allele-level comparison and read-level validation (`allele_level/`) and the manuscript files
  (`manuscript_v2/`).

## Citation

Please cite the associated manuscript and the archived Zenodo release (concept DOI: 10.5281/zenodo.20135136,
which always resolves to the latest version).
