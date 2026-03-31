# **Single-Cell Multiomics Reveals Somatic Copy Number Alterations and DCLK1-Driven Inflammatory Signaling in Injured Proximal Tubule Cells**

Welcome to this GitHub repository.
This folder is a publication-oriented code release for the RAGE24 analysis project.

Visit the Wilsoin Lab website:
[https://www.parkerwilsonlab.com](https://www.parkerwilsonlab.com)

# Contents
- [Data Availability](#data-availability)
- [Environment](#environment)
- [Analysis Workflow](#analysis-workflow)
- [Figures](#figures)
- [Citation](#citation)
- [Contact](#contact)

---
## **Data Availability**
The raw and processed sequencing data generated in this study have been deposited in the Gene Expression Omnibus (GEO) under the following accession numbers.

- single-cell multiome and DEFND-seq data of aging rat kidney cortex , and corresponding cellranger count matrices and cellranger-atac fragment files can be found at [GSE312213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312213) [Reviewer token: cfavwmimfjehjcf]
 
- bulk long-read RNA sequencing of aged rat kidney cortex, and corresponding transcript assemblies, transcript and gene count matrices can be found at [GSE312214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312214) [Reviewer token: yhulkcckrzyljwl]

- bulk RNA sequencing of DCLK1 overexpression in HEK293T cells, and corresponding count matrices and DESeq objects can be found at [GSE312496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312496) [Reviewer token: ebmtmkwezvefbgj]

- bulk RNA sequencing of rat kidney cortex following DCLK1-IN-1 treatment, and corresponding count matrices and DESeq objects can be found at [GSE313195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313195) [Reviewer token: exqnwwekllwzrup]. 

- Previously published bulk RNA-seq data from aged rat kidney cortex were obtained from the Sequence Read Archive (SRA; [PRJNA516151](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA516151)).

---
## **Environment**
Analysis scripts are written in R and use a containerized workflow in development (see docker run blocks in manuscript scripts).

Common R packages used in manuscript scripts include:
- `tidyverse`
- `arrow`
- `Seurat`
- `GenomicRanges`
- `patchwork`
- `fs`
- `here`

Currently, we developed two docker containers on [docker hub](https://hub.docker.com/u/garywang7).
- long read: in /bulk_long_read/docker_images. 
    `` docker pull garywang7/isoseq:1.0.6 ``
- multiome analysis: in /docker.
    `` docker pull garywang7/ragemultiome:1.0.4 ``  
---
## **Analysis Workflow**

### **Single-cell preprocessing and gene expression analysis**
1. Align and count 10x libraries for gene expression and ATAC data (`RAGE24_cellranger/step1_cellranger_gex_count.sh`, `RAGE24_cellranger/step1_cellranger_atac_count.sh`; PMACS notes in `RAGE24_cellranger/cellranger_PMACS.md`).
2. Build aggregated gene expression objects and sample-level metadata (`RAGE24_gex_aggr_prep/step1_gex_aggr_prep.R`).
3. Annotate cell identities and generate analysis-ready objects for downstream expression analyses (`RAGE24_gex_aggr_prep/step2_anno.R`).

### **Negative binomial CNA identification**
1. Prepare CNV sample sheets, barcode mappings, and genome/bin metadata (`RAGE24_negative_binomial_CNA/step0_CNV_prep.R`).
2. Count fragments for each library across genomic bins (`RAGE24_negative_binomial_CNA/step1_count_frags.R`).
3. Aggregate bin-level counts to chromosome and chromosome-arm features (`RAGE24_negative_binomial_CNA/step2_chrom_count.R`).
4. Combine libraries and run negative binomial depth-based CNA modeling (`RAGE24_negative_binomial_CNA/step3_combine_chrom_depth.R`).
5. Container/runtime notes are documented in `RAGE24_negative_binomial_CNA/README.md`.

### **Bulk long-read analysis (Iso-Seq and IsoQuant)**
1. Run Iso-Seq preprocessing for Kinnex reads (segment, demultiplex/refine, cluster, align, collapse, and classify) using scripts in `bulk_long_read/isoseq/` (`step1_segment_read.sh` to `step9_classify.sh`).
2. Run IsoQuant quantification workflows (`bulk_long_read/isoquant/isoquant_step1_bam2fastq.sh`, `bulk_long_read/isoquant/isoquant_step2_runIsoquant.sh`, `bulk_long_read/isoquant/isoquant_step2_runIsoquant_2.sh`).
3. Continue downstream isoform analyses with `bulk_long_read/isoquant/analysis/step0.Dclk1.sh`.
4. Full long-read processing notes are in `bulk_long_read/IsoSeq_pipeline.md`.

### **HEK293T DCLK1 overexpression mRNA-seq analysis**
1. Align and quantify reads for the DCLK1 overexpression-only experiment (`HEK293T_DCLK1_overexpression/HEK293T_DCLK1_oe only/step0_alignment.R`).
2. Align and quantify reads for the DCLK1 overexpression + DCLK1-IN-1 experiment (`HEK293T_DCLK1_overexpression/HEK293T_DCLK1_oe_DCLK1-IN-1/step0_alignment.R`).
3. Combine and harmonize differential expression outputs across conditions (`HEK293T_DCLK1_overexpression/step1_combine_DE.R`).
4. Experiment-specific metadata, QC, and reference notes are in the two HEK293T README files under `HEK293T_DCLK1_overexpression/`.

### **Cell type deconvolution analysis**
1. Run deconvolution and summarize inferred cell type composition profiles with `RAGE24_celltype_deconvolution/celltype_deconvolution_step1_deconv.R`.

### **Human atlas integration and gene expression analysis**
1. Preprocess kidney atlas ATAC/multiome fragments and build initial objects (`human_atlas/step1_snapatac2_kidney_prep.md`).
2. Integrate multimodal data with MultiVI and generate latent-space representations (`human_atlas/step2_multivi.py`).
3. Call peaks and assemble atlas-wide peak features (`human_atlas/step3_callpeak.py`).
4. Perform negative binomial-based CNA modeling and merge outputs (`human_atlas/step4_negbinom.md`, `human_atlas/step5_combine.R`, `human_atlas/step6_loy.R`).
5. Run differential expression and R/Seurat conversion steps (`human_atlas/step7_diffexp.py`, `human_atlas/step8_toseurat.py`, `human_atlas/step9_convert.R`, `human_atlas/step10_deg.R`).

---
## **Figures**
Manuscript figure scripts are maintained in:
`/figures`

Current figure scripts include:
- `figure1.R`
- `figure2.R`
- `figure3.R`
- `figure4.R`
- `figure5.R`
- `figure6.R`
- `figure7.R`
- `figure8.R`
- other supplementary figures.
---
## **Citation**
If you use any of the code or workflows in this repository please cite our manuscript.

---
## **Contact**
For questions and comments regarding the analysis and publication, please contact the corresponding author, Dr. Parker Wilson, at parker.wilson@pennmedicine.upenn.edu.
