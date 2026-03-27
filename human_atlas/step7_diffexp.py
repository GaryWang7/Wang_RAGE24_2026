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

# read filtered an annotated adata
adata_mvi = ad.read_h5ad('scratch/adata_mvi_model100/annotated_h5ad/adata.h5ad')

# read genotype annotations and update adata
anno = pd.read_csv('scratch/adata_mvi_model100/annotated_h5ad/filtered_chrX_annotations.csv')
anno.index = adata_mvi.obs.index
adata_mvi.obs['genotype'] = anno['has_chrX']
adata_mvi.obs['modality'] = anno['modality']

# # keep genes and remove peaks from .var
genes = adata_mvi.var.index[0:36600]
adata_mvi = adata_mvi[:, genes]

# # also remove atac only samples
multi_indices = adata_mvi.obs_names[adata_mvi.obs['modality'].isin(['multi'])]
adata_mvi = adata_mvi[adata_mvi.obs_names.isin(multi_indices)]

# filter genes 
# sc.pp.filter_genes(adata, min_cells=3)

# # Saving count data
adata_mvi.layers["counts"] = adata_mvi.X.copy()

# # Normalizing to median total counts
sc.pp.normalize_total(adata_mvi)

# # Logarithmize the data
sc.pp.log1p(adata_mvi)

# Obtain cluster-specific differentially expressed genes
cell_type_1 = "PCT"
cell_idx1 = (adata_mvi.obs["celltype"] == cell_type_1) & (adata_mvi.obs["modality"] == "multi") & (adata_mvi.obs["genotype"] == 0) 
cell_idx2 = (adata_mvi.obs["celltype"] == cell_type_1) & (adata_mvi.obs["modality"] == "multi") & (adata_mvi.obs["genotype"] == 1)

adata_mvi.obs["comp"] = ""
adata_mvi.obs.loc[cell_idx1, "comp"] = "REF"
adata_mvi.obs.loc[cell_idx2, "comp"] = "TEST"
sc.tl.rank_genes_groups(adata_mvi, groupby="comp", test="wilcoxon", key="genotype", reference="REF", groups=['REF','TEST'])
deg = sc.get.rank_genes_groups_df(adata_mvi, group="TEST", pval_cutoff=0.05)
os.makedirs('scratch/adata_mvi_model100/annotated_h5ad/diffexp', exist_ok=True)
filename=cell_type_1 + "_has_chrX.csv"
deg.to_csv(os.path.join('scratch/adata_mvi_model100/annotated_h5ad/diffexp', filename))

# sc.tl.rank_genes_groups(adata_mvi, groupby="celltype", test="wilcoxon", key="celltype")
# deg = sc.get.rank_genes_groups_df(adata_mvi, group="PCT")

celltypes=pd.unique(adata_mvi.obs['celltype'])
for celltype in celltypes:
  print(celltype)
  cell_type_1 = celltype
  cell_idx1 = (adata_mvi.obs["celltype"] == cell_type_1) & (adata_mvi.obs["modality"] == "multi") & (adata_mvi.obs["genotype"] == 0) 
  cell_idx2 = (adata_mvi.obs["celltype"] == cell_type_1) & (adata_mvi.obs["modality"] == "multi") & (adata_mvi.obs["genotype"] == 1)
  adata_mvi.obs["comp"] = ""
  adata_mvi.obs.loc[cell_idx1, "comp"] = "REF"
  adata_mvi.obs.loc[cell_idx2, "comp"] = "TEST"
  sc.tl.rank_genes_groups(adata_mvi, groupby="comp", test="wilcoxon", key="genotype", reference="REF", groups=['REF','TEST'])
  deg = sc.get.rank_genes_groups_df(adata_mvi, group="TEST", pval_cutoff=0.05)
  os.makedirs('scratch/adata_mvi_model100/annotated_h5ad/diffexp', exist_ok=True)
  filename=cell_type_1 + "_has_chrX.csv"
  deg.to_csv(os.path.join('scratch/adata_mvi_model100/annotated_h5ad/diffexp', filename))


########################################################
# # differential expression with DL model
# # load model
# model_dir = 'adata_mvi_model100'
# model = scvi.model.MULTIVI.load(os.path.join('scratch', model_dir))
# adata_mvi = model.adata

# # read annotations and update adata
# anno = pd.read_csv(os.path.join('scratch', model_dir, 'annotated_h5ad', 'loy_annotations.csv'))
# anno['celltype'] = anno['celltype'].astype('category')
# anno.index = adata_mvi.obs.index
# adata_mvi.obs['celltype'] = anno['celltype']
# adata_mvi.obs['genotype'] = anno['genotype']
# adata_mvi.obs['modality'] = anno['modality']

# os.makedirs(os.path.join('scratch', model_dir, 'differential'))

# # cell-specific markers for cluster
# de_change = model.differential_expression(groupby = "celltype", group1 = "PODO")
# de_change.to_csv(os.path.join('scratch', 'temp6.csv'))
# adata_mvi.obs[['celltype','genotype']].groupby(['celltype','genotype']).value_counts()

# cell_type_1 = "PCT"
# cell_idx1 = (adata_mvi.obs["celltype"] == cell_type_1) & (adata_mvi.obs["modality"] == "multi") & (adata_mvi.obs["genotype"] == "LOY") 
# cell_idx2 = (adata_mvi.obs["celltype"] == cell_type_1) & (adata_mvi.obs["modality"] == "multi") & (adata_mvi.obs["genotype"] == "XY")

# de_change = model.differential_expression(idx1=cell_idx1, idx2=cell_idx2)
# de_change.to_csv(os.path.join('scratch', 'pct_loy.csv'))

# # all PT loy vs xy deg
# cell_type_1 = "PCT"
# cell_type_2 = "PST"
# cell_idx1 = ((adata_mvi.obs["celltype"] == cell_type_1) | (adata_mvi.obs["celltype"] == cell_type_2)) & (adata_mvi.obs["genotype"] == "LOY") 
# cell_idx2 = ((adata_mvi.obs["celltype"] == cell_type_1) | (adata_mvi.obs["celltype"] == cell_type_2)) & (adata_mvi.obs["genotype"] == "XY") 

# # diff exp with batch correction by sample (batch_key is originally set to is_multiome)
# scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="sample")
# de_change = model.differential_expression(adata = adata_mvi, idx1=cell_idx1, idx2=cell_idx2, batch_correction=True)
# de_change.to_csv(os.path.join('scratch', 'pt_loy2.csv'))

# # cell-specific markers
# cell_type_1 = "PT_VCAM1"
# cell_idx1 = (adata_mvi.obs["celltype"] == cell_type_1) 
# cell_idx2 = (adata_mvi.obs["celltype"] != cell_type_1) 

# # diff exp with batch correction by sample (batch_key is originally set to is_multiome)
# scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="sample")
# de_change = model.differential_expression(adata = adata_mvi, idx1=cell_idx1, idx2=cell_idx2, batch_correction=True)
# de_change.to_csv(os.path.join('scratch', 'ptvcam1_deg.csv'))

# # cell-specific markers
# cell_type_1 = "PT_VCAM1"
# cell_type_2 = "PCT"
# cell_type_3 = "PST"
# cell_idx1 = (adata_mvi.obs["celltype"] == cell_type_1) 
# cell_idx2 = ((adata_mvi.obs["celltype"] == cell_type_2) | (adata_mvi.obs["celltype"] == cell_type_3))

# # diff exp with batch correction by sample (batch_key is originally set to is_multiome)
# scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="sample")
# de_change = model.differential_expression(adata = adata_mvi, idx1=cell_idx1, idx2=cell_idx2, batch_correction=True)
# de_change.to_csv(os.path.join('scratch', 'ptvcam1_vs_pt_deg.csv'))

