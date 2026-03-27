# 01/26/2025

# This script: 
# 1. creates a barcode conversion file between ATAC and GEX barcodes.
# 2. outputs a sample sheet for chasm to run
# 3. generates a reference information file for GCF_036323735.1 (with BN7.2 mitochondrial genome), which we used as cellranger-arc reference.

# Why barcode conversion files? 10X multiome gex and atac barcodes are not the same. However, they provide barcode files, which contains
# barcodes for both ATAC and gex. For example, row 1783 of ATAC barcode and row 1783 of the gex barcode are from the same gel bead.
# Refer to this page https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/bam-atac#bam-barcode-translation for details.

# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --workdir $HOME \
# --name CNV \
# -v /mnt/g:$HOME/g \
# -v /mnt/e:$HOME/e \
# -v /mnt/y:$HOME/y \
# -v /mnt/d/garyw:$HOME/garyw \
# -v $HOME:$HOME \
# -v /mnt/g/:$RHOME/g \
# -v /mnt/e:$RHOME/e \
# -v /mnt/y:$RHOME/y \
# -v /mnt/d/garyw:$RHOME/garyw \
# -v /var/run/docker.sock:/var/run/docker.sock \
# -e PASSWORD=garywang \
# -e DISABLE_AUTH=TRUE \
# -p 8787:8787 \
# garywang7/chasm:2.0.4 

library(here)
library(dplyr)
library(fs)

data_out_dir <- here("garyw", "RAGE24","CNV","resources")
cellranger_dir <- here("garyw","tools","cellranger-arc","cellranger-arc-2.0.2")
csv_dir <- here("garyw","RAGE24","gex_aggr_prep_combine")
ref_dir <- here("garyw","RAGE24","reference",'GRCr8')

dir_create(data_out_dir)
# Creating barcode conversion file
cellranger_atac_bc <- data.table::fread(here(cellranger_dir,"lib","python","atac","barcodes","737K-arc-v1.txt.gz"), header = F)
cellranger_gex_bc <- data.table::fread(here(cellranger_dir,"lib","python","cellranger","barcodes","737K-arc-v1.txt.gz"), header = F)
bc_conversion_df <- data.frame("ATAC_bc" = cellranger_atac_bc$V1,
                               "GEX_bc" = cellranger_gex_bc$V1)
write.csv(bc_conversion_df, file = here(data_out_dir, "barcode_conversion_10Xmultiome.csv"), row.names = FALSE, quote = FALSE)

# Creating sample sheet for CNV. We will only analyze the libraries/samples after filtering. This is basically the sample ID and library ID for samples to analyze.
anno_table <- read.csv(here(csv_dir, "gex_aggr_anno.csv"), header = TRUE)
sample_sheet_ATAC <- anno_table %>%
  dplyr::select(sample_id, library_type) %>%
  distinct()%>%
  mutate(library_id = paste0(sample_id, "-",library_type,"-ATAC"))
write.csv(sample_sheet_ATAC, file = here(data_out_dir, "CNV_samples.csv"), row.names = FALSE, quote = FALSE)

# Create genome information based on GRCr8_mt
# The primary assembly (autochromosomes+mitochondrial reads) seqlengths should be the same for all GRCr8

# Download assembly report for GCF_036323735.1, which was used as cellranger-arc reference
# from Genome Reference Consortium https://www.ncbi.nlm.nih.gov/grc/human and RefSeq
# at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_assembly_report.txt
# One change was made to Ref seq accession number of chrMT, which is NC_001665.2. We supplied GRCr8 with this genome.
ref_info <- read.table(here(ref_dir, "GCF_036323735.1_GRCr8_assembly_report.txt"), 
                       skip = 33, header = TRUE,sep = "\t", comment = "") # Include "#" as a character instead of comment
ref_info <- ref_info %>% 
  dplyr::rename("NCBI.Style.Name" = X..Sequence.Name) %>%
  # In this build, the chromosome names for NCBI is the same as UCSC styles. We will use UCSC style name to create bins
  mutate(UCSC.style.name = NCBI.Style.Name)
write.table(ref_info,
            file = here(ref_dir, "GCF_036323735.1_GRCr8_assembly_info.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)




