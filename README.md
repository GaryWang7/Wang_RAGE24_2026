# **RAGE24: Aging Rat Kidney Single-Cell Multiome Analysis (Publication Code Release)**
__Gary Wang__

Welcome to this GitHub repository.
This folder is a publication-oriented code release for the RAGE24 analysis project.
It is used to export and share analysis scripts from the active development workspace:
`/home/gary/garyw/scripts/RAGE`.

The goal of this repository is reproducible manuscript-level analysis and figure generation.

# Contents
- [Repository Purpose](#repository-purpose)
- [Source and Release Locations](#source-and-release-locations)
- [Data Availability](#data-availability)
- [Environment](#environment)
- [Analysis Workflow](#analysis-workflow)
- [Figures](#figures)
- [Citation](#citation)
- [Contact](#contact)

---
## **Repository Purpose**
This repository is the publishable snapshot of code used in the RAGE24 manuscript analysis.
It is not the primary development location; ongoing updates happen in the source repository and are copied here for release.

---
## **Source and Release Locations**
- **Primary source code (active development):** `/home/gary/garyw/scripts/RAGE`
- **Publication code release (this repository):** `/home/gary/garyw/scripts/Wang_RAGE24_2026`

When preparing a release, scripts are copied from the source directory into this repository and documented here.

---
## **Data Availability**
Processed intermediate files and manuscript outputs are generated from the RAGE24 project data directories referenced in scripts (for example under `here("garyw", "RAGE24", ...)`).

Raw and processed sequencing datasets, metadata tables, and supplementary release artifacts should be linked here at manuscript submission time (for example GEO/Zenodo accession links).

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
`/home/gary/garyw/scripts/RAGE/manuscript`

Current figure scripts include:
- `figure1.R`
- `figure2.R`
- `figure3.R`
- `figure4_preprocess.R`
- `figure4.R`
- `figure5.R`
- `figure6.R`
- `pw_figure6.R`
- `figure7.R`
- `figure8.R`

---
## **Citation**
If you use code from this repository, please cite the associated RAGE24 manuscript (citation details to be added after publication/preprint release).

---
## **Contact**
For questions about this code release, contact Gary Wang.
