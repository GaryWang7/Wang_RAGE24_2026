# load in multiome only anndata converted to h5seurat and analyze in MAST adjusting for covariates
# SCRATCH1=/mnt/c/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/s:$HOME/data \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/c/scratch" \
# -v /mnt/c/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.3.2 R

# run in rstudio server with NAS mount
# SCRATCH1=/mnt/c/scratch
# workdir=/home/rstudio
# docker run -it --rm \
# -p 8888:8787 \
# -e PASSWORD=password \
# -v /mnt/s:$workdir/data \
# -v /mnt/g/reference:$workdir/reference \
# -v /mnt/g:$workdir/g \
# -v $workdir:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/c/scratch" \
# -v /mnt/c/scratch:$workdir/scratch \
# p4rkerw/sctools:R4.3.2
# navigate browser to localhost:8888
# username: rstudio
# password: password

library(Seurat)
library(SeuratDisk)
library(data.table)
library(tibble)
library(dplyr)
library(future)

# BiocManager::install('MAST')
# library(MAST)

# convert anndata to seurat obj
# Convert('scratch/adata_mvi_model100/annotated_h5ad/multi_only_adata.h5ad', dest = "h5seurat", overwrite = TRUE)
Convert('scratch/multi_only_adata.h5ad', dest = "h5seurat", overwrite = TRUE)
# srat <- LoadH5Seurat("scratch/adata_mvi_model100/annotated_h5ad/multi_only_adata.h5seurat", meta.data = FALSE, misc = FALSE)
