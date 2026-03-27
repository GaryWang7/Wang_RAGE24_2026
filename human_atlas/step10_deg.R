# load in multiome only anndata converted to h5seurat and analyze in MAST adjusting for covariates
# SCRATCH1=/mnt/h/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/s:$HOME/data \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g:$HOME/g \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/c/scratch" \
# -v /mnt/h/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.3.2 R

# run in rstudio server with NAS mount
# SCRATCH1=/mnt/h/scratch
# workdir=/home/rstudio
# docker run -it --rm \
# -p 8888:8787 \
# -e DISABLE_AUTH=true \
# -v /mnt/s:$workdir/data \
# -v /mnt/g/reference:$workdir/reference \
# -v /mnt/g:$workdir/g \
# -v $workdir:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/c/scratch" \
# -v /mnt/h/scratch:$workdir/scratch \
# p4rkerw/sctools:R4.3.2
# navigate browser to localhost:8888


library(Seurat)
library(SeuratDisk)
library(data.table)
library(tibble)
library(dplyr)
library(future)

BiocManager::install('MAST')
library(MAST)

srat <- LoadH5Seurat("scratch/adata_mvi_model100/annotated_h5ad/multi_only_adata.h5seurat", meta.data = FALSE, misc = FALSE)
# srat <- LoadH5Seurat("scratch/multi_only_adata.h5seurat", meta.data = FALSE, misc = FALSE)
anno <- fread("scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv")

meta <- data.frame(barcode = rownames(srat@meta.data))
meta <- meta %>% left_join(anno, by = "barcode")
srat@meta.data <- meta
rownames(srat@meta.data) <- meta$barcode

DimPlot(srat, group.by="celltype")

srat <- NormalizeData(srat)
VlnPlot(srat, features = "DCLK1", group.by = "celltype", pt.size=0)

saveRDS(srat, 'scratch/adata_mvi_model100/annotated_h5ad/multi_only.rds')
