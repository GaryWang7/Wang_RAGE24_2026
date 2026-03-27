# this script now implements doublet detection with doubletdetection python package for all rna modalities and exports doublet scores to the final h5ad

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

# pbmc = scvi.data.read_10x_multiome('scratch/pbmc_10k/filtered_feature_bc_matrix')
# pbmc.var_names_make_unique()
# rna_only = read_rna_only(sample = "Control_1", rna_path = 'data/cellranger_rna_counts/version_6.1.2/Control_1/outs/filtered_feature_bc_matrix.h5')
# atac_only = read_atac_only(sample = "Control_1", atac_path = 'scratch/kidney_10k_h5ad_output/Control_1.h5ad')

# sample = '2-15Nx'
# rna_path = 'data/cellranger_multi_counts'
# atac_path = 'scratch/kidney_10k_h5ad_output'
# out_dir = 'scratch/10X_multi_h5ad'
def read_10x_multiome(sample, rna_path, atac_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rna = sc.read_10x_h5(os.path.join(rna_path, sample,"outs","filtered_feature_bc_matrix.h5"))
    # filter rna object
    sc.pp.filter_cells(rna, min_genes=100)
    # remove doublets
    clf = doubletdetection.BoostClassifier(
        n_iters=10,
        clustering_algorithm="louvain",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1,
    )
    doublets = clf.fit(rna.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    rna.obs["doublet"] = doublets
    rna.obs["doublet_score"] = doublet_score
    rna = rna[rna.obs['doublet'] == 0] 
    # read in atac data
    atac = ad.read_h5ad(os.path.join(atac_path, sample + ".h5ad"))
    sample = sample + "_multi"
    # intersect indices
    indices = rna.obs_names[rna.obs_names.isin(atac.obs_names)]
    # subset matrices
    rna = rna[rna.obs_names.isin(indices)]
    rna.var['modality'] = "Gene Expression"
    atac = atac[atac.obs_names.isin(indices)]
    atac.var['modality'] = "Peaks"
    # stack the sparse mat counts
    multi_X = sparse.hstack((rna.X,atac.X))
    # concatenate var
    multi_var_names = rna.var_names.tolist() + atac.var_names.tolist()
    multi_var = pd.DataFrame(multi_var_names, columns = ["ID"])
    multi_var.index = multi_var_names
    # convert indices to df
    multi_index_names = atac.obs.index
    multi_index = pd.DataFrame(multi_index_names, columns = ["index"])
    multi_index.index = multi_index_names
    # merge rna and atac into new anndata
    multi = sc.AnnData(multi_X, multi_index, multi_var)
    multi.var_names_make_unique()
    multi.var['modality'] = rna.var['modality'].tolist() + atac.var['modality'].tolist()
    multi.obs['sample'] = sample
    multi.obs['doubletdetection_score'] = rna.obs['doublet_score']
    multi.obs['scrublet_score'] = atac.obs['doublet_score']
    multi.obs.index=multi.obs.index + '_' + multi.obs['sample']
    multi.X.indptr = multi.X.indptr.astype(np.int64)
    multi.X.indices = multi.X.indices.astype(np.int64)
    multi.write_h5ad(os.path.join(out_dir, sample + '.h5ad'))

def read_kpmp_multiome(sample, rna_path, atac_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rna = sc.read_10x_mtx(os.path.join(rna_path, sample + "_single-nuc_expression_matrix"))
    # filter rna object
    sc.pp.filter_cells(rna, min_genes=100)
    # remove doublets
    clf = doubletdetection.BoostClassifier(
        n_iters=10,
        clustering_algorithm="louvain",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1,
    )
    doublets = clf.fit(rna.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    rna.obs["doublet"] = doublets
    rna.obs["doublet_score"] = doublet_score
    rna = rna[rna.obs['doublet'] == 0] 
    # read in atac data
    atac = ad.read_h5ad(os.path.join(atac_path, sample + "_single-nuc_expression_matrix.h5ad"))
    sample = sample + "_multi"
    # intersect indices
    indices = rna.obs_names[rna.obs_names.isin(atac.obs_names)]
    # subset matrices
    rna = rna[rna.obs_names.isin(indices)]
    rna.var['modality'] = "Gene Expression"
    atac = atac[atac.obs_names.isin(indices)]
    atac.var['modality'] = "Peaks"
    # stack the sparse mat counts
    multi_X = sparse.hstack((rna.X,atac.X))
    # concatenate var
    multi_var_names = rna.var_names.tolist() + atac.var_names.tolist()
    multi_var = pd.DataFrame(multi_var_names, columns = ["ID"])
    multi_var.index = multi_var_names
    # convert indices to df
    multi_index_names = atac.obs.index
    multi_index = pd.DataFrame(multi_index_names, columns = ["index"])
    multi_index.index = multi_index_names
    # merge rna and atac into new anndata
    multi = sc.AnnData(multi_X, multi_index, multi_var)
    multi.var_names_make_unique()
    multi.var['modality'] = rna.var['modality'].tolist() + atac.var['modality'].tolist()
    multi.obs['sample'] = sample
    multi.obs['doubletdetection_score'] = rna.obs['doublet_score']
    multi.obs['scrublet_score'] = atac.obs['doublet_score']
    multi.obs.index=multi.obs.index + '_' + multi.obs['sample']
    multi.X = scipy.sparse.coo_matrix.tocsr(multi.X)
    multi.X.indptr = multi.X.indptr.astype(np.int64)
    multi.X.indices = multi.X.indices.astype(np.int64)
    multi.write_h5ad(os.path.join(out_dir, sample + '.h5ad'))

def read_atac_only(sample, atac_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    atac = ad.read_h5ad(os.path.join(atac_path, sample + ".h5ad"))
    atac.var_names_make_unique()
    sample = sample + "_atac"
    atac.obs['sample'] = sample
    atac.obs.index=atac.obs.index + '_' + atac.obs['sample']
    atac_var_names = atac.var_names.tolist()
    atac_var = pd.DataFrame(atac_var_names, columns = ["ID"])
    atac_var.index = atac_var_names
    atac.var = atac_var
    atac.var['modality'] = "Peaks"
    atac_X = atac.X
    atac_index_names = atac.obs.index
    atac_index = pd.DataFrame(atac_index_names, columns = ["index"])
    atac_index.index = atac_index_names
    atac_var = atac.var
    atac = sc.AnnData(atac_X, atac_index, atac_var)
    atac.obs['sample'] = sample
    atac.obs['scrublet_score'] = atac.obs['doublet_score']
    atac.X.indptr = atac.X.indptr.astype(np.int64)
    atac.X.indices = atac.X.indices.astype(np.int64)
    atac.write_h5ad(os.path.join(out_dir, sample + '.h5ad'))
    # return(atac_only)

def read_atac_only_pad(sample, atac_path, out_dir, multi_var):
    os.makedirs(out_dir, exist_ok=True)
    atac = ad.read_h5ad(os.path.join(atac_path, sample + ".h5ad"))
    atac.var_names_make_unique()
    sample = sample + "_atac"
    num_cells = len(atac.obs.index)
    num_rna_features = len(multi_var.index) - len(atac.var.index)
    rna_features = multi_var.index[0:num_rna_features]
    rna_X = np.zeros((num_cells, num_rna_features), dtype = 'float64')
    multi_X = sparse.hstack((rna_X,atac.X))
    # convert indices to df
    multi_index_names = atac.obs.index
    multi_index = pd.DataFrame(multi_index_names)
    multi_index.index = multi_index_names
    # merge rna and atac into new anndata
    multi = sc.AnnData(multi_X, multi_index, multi_var)
    multi.var_names_make_unique()
    multi.obs['sample'] = sample
    multi.obs['scrublet_score'] = atac.obs['doublet_score']
    multi.obs.index=multi.obs.index + '_' + multi.obs['sample']
    multi.X = scipy.sparse.coo_matrix.tocsr(multi.X)
    multi.X.indptr = multi.X.indptr.astype(np.int64)
    multi.X.indices = multi.X.indices.astype(np.int64)
    multi.write_h5ad(os.path.join(out_dir, sample + '.h5ad'))
    # return(multi)

# process 10x in-house multiomes
os.makedirs('scratch/10x_multi_h5ad', exist_ok=True)
samples = []
for filename in os.listdir('data/cellranger_multi_counts'):
    if os.path.isdir(os.path.join('data/cellranger_multi_counts',filename)):
        samples.append(filename)
        print(filename)
for sample in samples:
    read_10x_multiome(sample, rna_path = 'data/cellranger_multi_counts', atac_path = 'scratch/kidney_10k_h5ad_output', out_dir = 'scratch/10X_multi_h5ad')

# process kpmp multiomes 
os.makedirs('scratch/kpmp_multi_h5ad')
kpmp_csv = pd.read_csv('data/KPMP_Data/single-nucleus-multiomes/atlas_repository_filelist-20240227.csv')
kpmp_samples = kpmp_csv['Internal Package ID']
kpmp_samples = kpmp_samples.unique().tolist()

for sample in kpmp_samples:
    read_kpmp_multiome(sample, rna_path = 'data/KPMP_Data/single-nucleus-multiomes/', atac_path = 'scratch/kidney_10k_h5ad_output', out_dir = 'scratch/kpmp_multi_h5ad')

# concat in-house 10X and kpmp multiomes on disk
adatas = glob(os.path.join('scratch','kpmp_multi_h5ad','*.h5ad'), recursive=False) + glob(os.path.join('scratch','10x_multi_h5ad','*.h5ad'), recursive=False)
ad.experimental.concat_on_disk(adatas, out_file = 'scratch/all_multi.h5ad')

# process atac-only samples
atac_samples = []
for filename in os.listdir('data/cellranger_atac_counts/kidney/version_2.1'):
    if os.path.isdir(os.path.join('data/cellranger_atac_counts/kidney/version_2.1',filename)):
        atac_samples.append(filename)
        print(filename)

# retrieve multiome .var and pad the atac samples with zeros for the rna modality
# this makes it possible to do on-disk concat - which only supports inner join
multi = ad.read_h5ad('scratch/10x_multi_h5ad/1-27Nx_multi.h5ad')
for sample in atac_samples:
    read_atac_only_pad(sample, atac_path = 'scratch/kidney_10k_h5ad_output', out_dir = 'scratch/atac_h5ad_pad', multi_var = multi.var)
anndata_list = glob(os.path.join('scratch','atac_h5ad_pad','*.h5ad'), recursive=False)
ad.experimental.concat_on_disk(anndata_list, out_file = 'scratch/atac_only_pad.h5ad')

# here we will use atac-only h5ad padded with zeros for the rna portion and concat with the multiomes
ad.experimental.concat_on_disk(['scratch/all_multi.h5ad', 'scratch/atac_only_pad.h5ad'], out_file = 'scratch/adata_mvi.h5ad')

# load processed adata into memory
adata_mvi = ad.read_h5ad('scratch/adata_mvi.h5ad')

# convert pointers to 64bit indices o/w there are too many non-zero elements for 32bit indices (scanpy will not convert automatically when concatenating)
adata_mvi.X.indptr = adata_mvi.X.indptr.astype('int64')
adata_mvi.X.indices = adata_mvi.X.indices.astype('int64')

# filter features in less than 1% of cells - incredibly slow... and possibly not necessary
# consider filtering the rna features prior to the concat_on_disk inner join
# atac features can be filtered by snapatac2
# sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))

# quantify number of genes and regions for model input
multi = ad.read_h5ad('scratch/10x_multi_h5ad/1-27Nx_multi.h5ad')
all_genes = multi.var.index[multi.var['modality'] == "Gene Expression"].tolist()
all_peaks = multi.var.index[multi.var['modality'] == "Peaks"].tolist()

num_genes = adata_mvi.var.index.isin(all_genes).sum()
num_regions = adata_mvi.var.index.isin(all_peaks).sum()

# setup model
# note that if rna-only data is added this will need to be modified - it only differentiates between multi and non-multi
multiomes = kpmp_samples + samples
multiomes = [ e + '_multi' for e in multiomes]
adata_mvi.obs['is_multiome'] = adata_mvi.obs['sample'].isin(multiomes)

scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="is_multiome", categorical_covariate_keys = ["sample"])

model = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=num_genes,
    n_regions=num_regions,
)
model.view_anndata_setup()

model_dir = 'adata_mvi_model100'
model.train(max_epochs = 100)
model.save(os.path.join('scratch', model_dir), save_anndata = True)
# model = scvi.model.MULTIVI.load(os.path.join('scratch', model_dir))
# adata_mvi = model.adata
os.makedirs(os.path.join('scratch', model_dir, 'figures'))

MULTIVI_LATENT_KEY = "X_multivi"
adata_mvi.obsm[MULTIVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(adata_mvi, use_rep=MULTIVI_LATENT_KEY)
sc.tl.umap(adata_mvi, min_dist=0.2)
pl=sc.pl.umap(adata_mvi, color="_scvi_batch", return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_batch_umap.png"))

# cluster the data without a legend
snap.tl.leiden(adata_mvi)
pl=sc.pl.umap(adata_mvi, color="leiden", return_fig=True)
pl.gca().get_legend().remove()
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_leiden_umap.png"))

# with labels on clusters
pl=sc.pl.umap(adata_mvi, color="leiden", legend_loc='on data', return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_leiden_labeled_umap.png"))

# visualize markers without imputation
pl=sc.pl.umap(adata_mvi, color="VCAM1", vmax=partial(np.percentile, q=99), return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_vcam1_umap.png"))

pl=sc.pl.umap(adata_mvi, color=["CUBN"], vmax=partial(np.percentile, q=99), return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_CUBN_umap.png"))


marker_genes=["CUBN","HAVCR1","SLC5A1","SLC5A2","VCAM1","PROM1","CFH","SLC12A1","CLDN10","CLDN16","SLC12A3","TRPM6","SCNN1G","TRPV5","CALB1","AQP2","ATP6V0D2","SLC4A1","SLC26A7","SLC26A4","NPHS1","NPHS2","WT1","CLIC5","PECAM1","FLT1","IGFBP5","IGFBP7","PLVAP","EHD3","SLC6A6","SLC14A1","AQP1","NOS1","ITGA8","PDGFRB","MEIS2","PIEZO2","REN","ACTA2","CALD1","PROX1","FLT4","PDPN","PTPRC","CD3E","MS4A1","CD19","SDC1","CD14","CSF1R"]
# impute expression for selected genes - not sure this works very well...
# imputed_expression = model.get_normalized_expression(gene_list = marker_genes)

# impute expression for all genes
# imputed_expression = model.get_normalized_expression()

# export DataFrame to text file
# imputed_expression.to_csv(os.path.join('scratch', model_dir, 'imputed_expression.csv'), index=False)

def draw_gene_umap(gene_sel, threshold):
    pl=sc.pl.umap(adata_mvi, color=[gene_sel], vmax=partial(np.percentile, q=threshold), return_fig=True)
    pl.savefig(os.path.join('scratch', model_dir, 'figures', gene_sel + "_umap.png"))

for gene_sel in marker_genes:
  draw_gene_umap(gene_sel, 95)

# def draw_imputed_gene_umap(gene_sel, threshold):
#     gene_idx = marker_genes.index(gene_sel)
#     adata_mvi.obs["gene_imputed"] = imputed_expression.iloc[:, gene_idx]
#     pl=sc.pl.umap(adata_mvi, color=["gene_imputed"], vmax=partial(np.percentile, q=threshold), return_fig=True)
#     pl.savefig(os.path.join('scratch', model_dir, 'figures', "imputed_" + gene_sel + "_umap.png"))
# draw_imputed_gene_umap("CUBN", 99)

pl = sc.pl.dotplot(adata_mvi, var_names = marker_genes, groupby = 'leiden', return_fig = True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "dotplot.png"))

pl=sc.pl.umap(adata_mvi, color=["scrublet_score"], return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scrublet_umap.png"))

# calculate overall mean scrublet score 
# trim the mean and stdev to exclude potential outlier doublet clusters
from scipy import stats
mean_scrublet=stats.trim_mean(adata_mvi.obs['scrublet_score'], 0.1)
sd_scrublet=stats.mstats.trimmed_std(adata_mvi.obs['scrublet_score'], 0.1)
threshold=mean_scrublet + 2 * sd_scrublet
# calculate mean scrublet score grouped by cluster and compare to threshold
is_doublet = adata_mvi.obs.groupby(['leiden'])['scrublet_score'].mean() > threshold
doublet_clusters_to_filter = is_doublet[is_doublet.values].index.to_list()

# calculate proportion of each cluster by sample and filter clusters with > 50% from a single sample
adata_mvi.obs['leiden'].value_counts()
df=pd.DataFrame({'sample': adata_mvi.obs['sample'], 'leiden': adata_mvi.obs['leiden']})
num_sample_by_cluster = df.groupby(['leiden','sample'])['leiden'].count()
num_sample_by_cluster.columns = ['sample','cluster_sample_count']
num_sample = df.groupby(['leiden']).count()
num_sample.columns = ['total_cluster_count']
num_sample = num_sample.merge(num_sample_by_cluster.rename('cluster_sample_count'), left_index=True, right_index=True)
num_sample['prop_sample_cluster'] = num_sample['cluster_sample_count'] / num_sample['total_cluster_count']
num_sample = num_sample.sort_values(by='prop_sample_cluster', ascending=False)
num_sample = num_sample.reset_index()
clusters_to_filter = num_sample[num_sample['prop_sample_cluster'] > 0.50]['leiden']
clusters_to_filter = clusters_to_filter.tolist()

# identify barcodes to keep 
idkeep=~adata_mvi.obs['leiden'].isin(clusters_to_filter + doublet_clusters_to_filter)
idkeep.value_counts()
index_keep = idkeep[idkeep.values].index.to_list()

map_dict = {0: 'PCT',
            1: 'TAL',
            2: 'PST',
            3: 'TAL',
            4: 'PCT',
            5: 'FILTER',
            6: 'DCT2_PC',
            7: 'DCT1',
            8: 'PT_VCAM1',
            9: 'FILTER',
            10: 'ENDO',
            11: 'FIB_VSMC_MC',
            12: 'ICA',
            13: 'MONO',
            14: 'PCT',
            15: 'FILTER',
            16: 'ICB',
            17: 'TL',
            18: 'FILTER',
            19: 'TCELL',
            20: 'FILTER',
            21: 'PEC',
            22: 'PODO',
            23: 'BCELL',
            24: 'ENDO',
            25: 'FILTER',
            26: 'FILTER'}
clusters = adata_mvi.obs['leiden'].astype('int')
adata_mvi.obs['celltype'] = clusters.map(map_dict).astype('category')

# visualize annotation
pl=sc.pl.umap(adata_mvi, color="celltype", legend_loc='on data', return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_celltype_labeled_umap.png"))

# save the annotations in a csv
os.makedirs(os.path.join('scratch', model_dir, 'annotated_h5ad'))
df=pd.DataFrame({'barcode': adata_mvi.obs.index, 'sample': adata_mvi.obs['sample'], 'celltype': adata_mvi.obs['celltype']})
df.to_csv(os.path.join('scratch', model_dir, 'annotated_h5ad', 'annotations.csv'), index=False)

# subset anndata
adata_mvi = adata_mvi[index_keep, ]
# write filtered annotations
df=pd.DataFrame({'barcode': adata_mvi.obs.index, 'sample': adata_mvi.obs['sample'], 'celltype': adata_mvi.obs['celltype']})
df.to_csv(os.path.join('scratch', model_dir, 'annotated_h5ad', 'filtered_annotations.csv'), index=False)
pl=sc.pl.umap(adata_mvi, color="celltype", legend_loc='on data', return_fig=True)
pl.savefig(os.path.join('scratch', model_dir, 'figures', "scvi_filtered_celltype_labeled_umap.png"))

# save the annotated object
adata_mvi.write_h5ad(os.path.join('scratch', model_dir, 'annotated_h5ad', 'adata.h5ad'))

# filter snap anndataset and update annotations
# make a copy of anndataset so it doesnt get corrupted later...
adataset_path = 'scratch/kidney_10k_h5ad_output'
os.makedirs(os.path.join(adataset_path, 'annotated'))
shutil.copy(os.path.join(adataset_path, 'adataset/data.h5ads'), os.path.join(adataset_path, 'annotated'))
adataset = snap.read_dataset(os.path.join(adataset_path, 'annotated', 'data.h5ads'), adata_files_update = adataset_path)

# annotate cell types in anndataset
# reformat adata_mvi barcodes to match with snapatac2 adataset...
df=pd.DataFrame({'sample': adata_mvi.obs['sample'], 'barcode': adata_mvi.obs_names, 'celltype': adata_mvi.obs['celltype']})
sample = df['sample'].str.replace('_atac', '')
sample = sample.str.replace('_multi', '')
df['sample'] = sample
barcode = df['barcode'].str.split('_', n=1, expand=True)[0]
df['barcode'] = barcode
df['query'] = df['sample'] + "_" + df['barcode']
lookup=pd.DataFrame({'barcode': df['query'], 'celltype': df['celltype'], 'sample': df['sample']})

# organize snapatac2 barcodes
df=pd.DataFrame({'sample': adataset.obs['sample'], 'barcode': adataset.obs_names})
df['new_sample'] = df['sample'].str.replace('_single-nuc_expression_matrix', '')
df['query'] = df['new_sample'] + "_" + df['barcode']
query=pd.DataFrame({'barcode': df['query']})

# merge the multivi annotations with snapatac2 barcodes in a pandas dataframe
anno=pd.merge(query, lookup, how='left', on=['barcode'])
anno['celltype'].value_counts()
adata_mvi.obs['celltype'].value_counts()
anno.to_csv(os.path.join('scratch','anno.csv'))
anno = pd.read_csv(os.path.join('scratch','anno.csv'))
anno["celltype"].fillna("Filtered", inplace = True)

# add annotation to anndataset
celltype = pl.Series("celltype", anno['celltype'])
adataset.obs['celltype'] = celltype.to_numpy()

# remove filtered cells
idkeep=~adataset.obs['celltype'].is_in(['Filtered'])
idkeep.value_counts()

# save new obj to adataset_path/annotated/_dataset.h5dads
adataset.subset(idkeep, out = os.path.join(adataset_path, 'annotated'))
adataset.close()

	










