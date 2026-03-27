# bsub -Is \
# -n 10 \
# -M 128000 \
# -R "rusage [mem=128000] span[hosts=1]" \
# bash

# # run pmacs
# singularity shell \
# --bind $HOME/scratch:$HOME/scratch \
# --bind $HOME/scripts:$HOME/scripts \
# --bind $HOME/data:$HOME/data \
# --bind $HOME/reference:$HOME/reference \
# --bind $HOME:$HOME \
# --cleanenv $HOME/images/snapatac2_2.6.sif

# # # run local launch from /mnt/g/hpap
# singularity shell \
# --bind /mnt/g/scratch:$HOME/scratch \
# --bind /mnt/g/scripts:$HOME/scripts \
# --bind /mnt/g/reference:$HOME/reference \
# --home /mnt/g/hpap \
# --cleanenv /mnt/g/docker/snapatac2_2.6.sif

# # python3

# this is necessary for proper multithreading on HPC in a python script
# it is not necessary in an interactive python session or jupyter notebook
if __name__ == '__main__':
	
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
	from glob import glob
	from tqdm import tqdm
	
  	adataset_path = 'scratch/kidney_10k_h5ad_output'
	
	outputdir=os.path.join(adataset_path, 'peaks')
	os.makedirs(outputdir)
	adataset = snap.read_dataset(os.path.join(adataset_path, 'annotated', '_dataset.h5ads'), adata_files_update = os.path.join(adataset_path, 'annotated', 'anndatas'))
	  
	snap.tl.macs3(adataset, groupby='celltype', n_jobs=10)
	
	peaks = snap.tl.merge_peaks(adataset.uns['macs3'], snap.genome.hg38)
	peaks.write_csv(os.path.join(outputdir, 'peaks.csv'))

	# make cellxpeak matrix
	peak_mat = snap.pp.make_peak_matrix(adataset, use_rep=peaks['Peaks'])

	peak_mat.write_h5ad(os.path.join(adataset_path, 'annotated', 'peak_mat.h5ad'))
	
	adataset.close()

	
