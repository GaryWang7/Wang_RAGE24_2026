# subset anndata to multiome only RNA assay and prepare for h5seurat conversion

# docker run -it --rm --gpus all \
# -v /mnt/c/scratch:$HOME/scratch \
# -v /mnt/g/scripts:$HOME/scripts \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/c/scratch/kidney_10k_h5ad_output:$HOME/kidney_10k_h5ad_output \
# -v /mnt/s/:$HOME/data \
# -v $HOME:$HOME \
# --workdir $HOME \
# p4rkerw/scvi-tools_py3.11-cu11-runtime-latest-snapatac2-harmonypy-doubletdetection:2.6 python3

# initialize libraries
import snapatac2 as snap
import numpy as np
import os
import magic
import scanpy as sc
import shutil
import snapatac2 as snap
import anndata as ad
import pandas as pd
import scvi
import scipy
import doubletdetection
import polars as pl
from scipy import sparse
from glob import glob
from tqdm import tqdm
from functools import partial
import scipy

# read filtered an annotated adata
adata_mvi = ad.read_h5ad('scratch/adata.h5ad')
# adata_mvi = ad.read_h5ad('scratch/adata_mvi_model100/annotated_h5ad/adata.h5ad')

# read genotype annotations and update adata
anno = pd.read_csv('scratch/filtered_annotations.csv')
# anno = pd.read_csv('scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv')
anno.index = adata_mvi.obs.index
modality = anno['sample'].str.split('_').str.get(1)
adata_mvi.obs['modality'] = modality
adata_mvi.obs['celltype'] = anno['celltype']

# # keep genes and remove peaks from .var
genes = adata_mvi.var.index[0:36600]
adata_mvi = adata_mvi[:, genes]

# # also remove atac only samples
multi_indices = adata_mvi.obs_names[adata_mvi.obs['modality'].isin(['multi'])]
adata_mvi = adata_mvi[adata_mvi.obs_names.isin(multi_indices)]

# save a multiome only anndata with raw counts (ie not normalized)
# remove _scvi_extra_categorical_covs which throws an error with h5seurat conversion
del adata_mvi.obsm['_scvi_extra_categorical_covs']
# adata_mvi.obs.to_csv('scratch/adata_mvi_model100/annotated_h5ad/metadata_multi_only_adata.csv')
adata_mvi.obs.to_csv('scratch/metadata_multi_only_adata.csv')
# adata_mvi.write_h5ad('scratch/adata_mvi_model100/annotated_h5ad/multi_only_adata.h5ad')
adata_mvi.write_h5ad('scratch/multi_only_adata.h5ad')
