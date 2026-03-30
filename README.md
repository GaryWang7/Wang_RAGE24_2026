# **Single-Cell Multiomics Reveals Somatic Copy Number Alterations and DCLK1-Driven Inflammatory Signaling in Injured Proximal Tubule Cells**

Welcome to this GitHub repository.
This folder is a publication-oriented code release for the RAGE24 analysis project.

Visit the Wilsoin Lab website:
[https://www.parkerwilsonlab.com](https://www.parkerwilsonlab.com)

# Contents
- [Data Availability](#data-availability)
- [Source and Release Locations](#source-and-release-locations)
- [Environment](#environment)
- [Analysis Workflow](#analysis-workflow)
- [Figures](#figures)
- [Citation](#citation)
- [Contact](#contact)

---
## **Data Availability**


---
## **Environment**
Analysis scripts are written in R and use a containerized workflow in development (see docker run blocks in manuscript scripts).

Common R packages used in manuscript scripts include:
- `tidyverse`
- `data.table`
- `arrow`
- `Seurat`
- `GenomicRanges`
- `patchwork`
- `fs`
- `here`

---
## **Analysis Workflow**
The high-level manuscript workflow is:

1. Prepare aggregated data and annotations from the RAGE24 project directories.
2. Run chromosome-level and bin-level CNV analyses (including negative binomial workflows).
3. Generate manuscript figure panels from curated analysis outputs.
4. Run cross-dataset comparisons (for example human atlas mapping) where applicable.
5. Export publication-quality figures and tables.

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

---
## **Citation**
If you use any of the code or workflows in this repository please cite our manuscript.

---
## **Contact**
For questions and comments regarding the analysis and publication, please contact the corresponding author, Dr. Parker Wilson, at parker.wilson@pennmedicine.upenn.edu.
