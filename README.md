# **Single-Cell Multiomics Reveals Somatic Copy Number Alterations and DCLK1-Driven Inflammatory Signaling in Injured Proximal Tubule Cells**

This repository contains the analysis code accompanying our study on somatic copy number alterations and DCLK1-associated inflammatory signaling in injured proximal tubule cells during kidney aging. 

It is intended to document the workflows used in the manuscript and is not distributed as a general-purpose software package.

Lab website: [Parker Wilson Lab](https://www.parkerwilsonlab.com)

---

## Contents
- [Data availability](#data-availability)
- [Software requirements](#software-requirements)
- [Computational environment](#computational-environment)
- [Analysis workflow](#analysis-workflow)
- [Figures](#figures)
- [Citation](#citation)
- [Contact](#contact)

---

## Data availability

The raw and processed sequencing data generated in this study are available through GEO under the following accession numbers: (reviewer tokens are provided in the manuscript or upon request)

- **Aging rat kidney cortex single-cell multiome and DEFND-seq**, including corresponding Cell Ranger count matrices and Cell Ranger ATAC fragment files:  
  [GSE312213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312213)  

- **Bulk long-read RNA-seq of aged rat kidney cortex**, including transcript assemblies and transcript- and gene-level count matrices:  
  [GSE312214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312214)  

- **Bulk RNA-seq of DCLK1 overexpression in HEK293T cells**, including count matrices and DESeq2 objects:  
  [GSE312496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312496)  

- **Bulk RNA-seq of rat kidney cortex following DCLK1-IN-1 treatment**, including count matrices and DESeq2 objects:  
  [GSE313195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313195)  

- **Previously published bulk RNA-seq data from aged rat kidney cortex** were obtained from SRA:  
  [PRJNA516151](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA516151)

---

## Software requirements

The main workflows in this repository were developed using:

- **R >= 4.3**
- **Python >= 3.8**
- **Cell Ranger >= 8.0**

Some pipelines were run in containerized environments to improve reproducibility.

---

## Computational environment

Analyses were performed primarily in **R** and **Python** using containerized workflows.

Commonly used R packages in manuscript scripts include:

- `tidyverse`
- `arrow`
- `Seurat`
- `GenomicRanges`
- `patchwork`
- `fs`
- `here`

Docker images used in this project are available on Docker Hub:

- **Long-read analysis**: `garywang7/isoseq:1.0.6`
- **Multiome analysis**: `garywang7/ragemultiome:1.0.4`

For additional runtime notes, see the workflow-specific README or script headers where applicable.

---

## Analysis workflow

### Single-cell preprocessing and gene expression analysis

1. Align and quantify 10x gene expression and ATAC libraries

   * `RAGE24_cellranger/step1_cellranger_gex_count.sh`
   * `RAGE24_cellranger/step1_cellranger_atac_count.sh`

2. Build aggregated gene expression objects and sample-level metadata

   * `RAGE24_gex_aggr_prep/step1_gex_aggr_prep.R`

3. Annotate cell identities and generate analysis-ready objects

   * `RAGE24_gex_aggr_prep/step2_anno.R`

Supporting notes:

* `RAGE24_cellranger/cellranger_PMACS.md`

### Negative binomial CNA analysis

1. Prepare sample sheets, barcode mappings, and genome/bin metadata

   * `RAGE24_negative_binomial_CNA/step0_CNV_prep.R`

2. Count fragments across genomic bins for each library

   * `RAGE24_negative_binomial_CNA/step1_count_frags.R`

3. Aggregate bin-level counts to chromosome and chromosome-arm features

   * `RAGE24_negative_binomial_CNA/step2_chrom_count.R`

4. Combine libraries and run negative binomial CNA modeling

   * `RAGE24_negative_binomial_CNA/step3_combine_chrom_depth.R`

Additional notes:

* `RAGE24_negative_binomial_CNA/README.md`


### Bulk long-read analysis

1. Run Iso-Seq preprocessing for Kinnex reads

   * `bulk_long_read/isoseq/step1_segment_read.sh` to `step9_classify.sh`

2. Run IsoQuant quantification workflows

   * `bulk_long_read/isoquant/isoquant_step1_bam2fastq.sh`
   * `bulk_long_read/isoquant/isoquant_step2_runIsoquant.sh`
   * `bulk_long_read/isoquant/isoquant_step2_runIsoquant_2.sh`

3. Perform downstream isoform analysis

   * `bulk_long_read/isoquant/analysis/step0.Dclk1.sh`

Additional notes:

* `bulk_long_read/IsoSeq_pipeline.md`


### HEK293T DCLK1 overexpression RNA-seq analysis

1. Align and quantify reads for the DCLK1 overexpression experiment

   * `HEK293T_DCLK1_overexpression/HEK293T_DCLK1_oe only/step0_alignment.R`

2. Align and quantify reads for the DCLK1 overexpression + DCLK1-IN-1 experiment

   * `HEK293T_DCLK1_overexpression/HEK293T_DCLK1_oe_DCLK1-IN-1/step0_alignment.R`

3. Combine and harmonize differential expression outputs across conditions

   * `HEK293T_DCLK1_overexpression/step1_combine_DE.R`

Additional experiment-specific notes are provided in the corresponding subdirectories.

### Cell type deconvolution

1. Run deconvolution and summarize inferred cell type composition

   * `RAGE24_celltype_deconvolution/celltype_deconvolution_step1_deconv.R`

### Human atlas integration and analysis

1. Preprocess kidney atlas ATAC/multiome fragments

   * `human_atlas/step1_snapatac2_kidney_prep.md`

2. Integrate multimodal data with MultiVI

   * `human_atlas/step2_multivi.py`

3. Call peaks and generate atlas-wide peak features

   * `human_atlas/step3_callpeak.py`

4. Perform chromosome-level negative binomial CNA modeling and merge results

   * `human_atlas/step4_negbinom.md`
   * `human_atlas/step5_combine.R`
   * `human_atlas/step6_loy.R`

5. Run differential expression and Seurat conversion steps

   * `human_atlas/step7_diffexp.py`
   * `human_atlas/step8_toseurat.py`
   * `human_atlas/step9_convert.R`
   * `human_atlas/step10_deg.R`

---

## Figures

Scripts used to generate manuscript figures are located in:

```text
figures/
```

This directory includes scripts for main figures and selected supplementary figures.

---

## Citation

If you use code or workflows from this repository, please cite our manuscript.

**Manuscript title:**
*Single-Cell Multiomics Reveals Somatic Copy Number Alterations and DCLK1-Driven Inflammatory Signaling in Injured Proximal Tubule Cells*

A DOI or formal citation will be added here when available.

---

## Contact

For questions regarding the analysis or repository, please contact:

**Dr. Parker Wilson**
`parker.wilson@pennmedicine.upenn.edu`

