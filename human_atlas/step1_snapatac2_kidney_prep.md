# snapatac2
```
# mount nas
sudo mount -t drvfs S: /mnt/s
# restart docker after mount
```

```
# run local docker - gpu integration does work!
docker run -it --rm --gpus all \
-v /mnt/c/scratch:$HOME/scratch \
-v /mnt/g/scripts:$HOME/scripts \
-v /mnt/g/reference:$HOME/reference \
-v /mnt/s/:$HOME/data \
-v $HOME:$HOME \
--workdir $HOME \
p4rkerw/scvi-tools_py3.11-cu11-runtime-latest-snapatac2:2.6 \
python3
```

```
# initialize libraries
import snapatac2 as snap
import numpy as np
import os
import magic
import scanpy as sc
from glob import glob
from tqdm import tqdm

snap.__version__
```

```
# load fragment files
# output dir will be in scratch volume
output_dir = "scratch/kidney_10k_h5ad_output"
os.makedirs(output_dir, exist_ok=True)

atac_fragment_files = glob('data/cellranger_atac_counts/kidney/version_2.1/*/outs/fragments.tsv.gz', recursive=False)
outputs = []
for fl in atac_fragment_files:
    name = fl.split('/')[-3].split('.tsv.gz')[0]
    outputs.append(f'{output_dir}/{name}.h5ad')
    print(name)

multi_fragment_files = glob('data/cellranger_multi_counts/*/outs/atac_fragments.tsv.gz', recursive=False)
for fl in multi_fragment_files:
    name = fl.split('/')[-3].split('.tsv.gz')[0]
    outputs.append(f'{output_dir}/{name}.h5ad')
    print(name)


kpmp_fragment_files = glob('data/KPMP_Data/single-nucleus-multiomes/*/atac_fragments.tsv.gz', recursive=False)
for fl in kpmp_fragment_files:
    name = fl.split('/')[-2].split('.tsv.gz')[0]
    outputs.append(f'{output_dir}/{name}.h5ad')
    print(name)

all_fragment_files = atac_fragment_files + multi_fragment_files + kpmp_fragment_files
len(all_fragment_files)
len(outputs)
```

```
# preprocess fragment files
# NVME runtime ~3h
adatas = snap.pp.import_data(all_fragment_files, file=outputs, chrom_sizes=snap.genome.hg38,
                             min_num_fragments=3000, n_jobs=8, sorted_by_barcode=False)

```

```
# add tile matrix
snap.pp.add_tile_matrix(adatas, bin_size=10000)
snap.metrics.tsse(adatas, snap.genome.hg38)
# filter cells based on fragment count and tsse
snap.pp.filter_cells(adatas, min_counts=3000, min_tsse=10, max_counts=100000)
```

```
# select features and remove doublets
snap.pp.select_features(adatas, n_jobs=10)
snap.pp.scrublet(adatas, n_jobs=10)
snap.pp.filter_doublets(adatas, n_jobs=10)
```

```
# make anndataset
os.makedirs(os.path.join(output_dir, 'adataset'))
adataset = snap.AnnDataSet(
		adatas=[(f.filename.split('/')[-1].split('.h5ad')[0], f) for f in adatas],
		filename=os.path.join(output_dir, 'adataset', 'data.h5ads')
	)
adataset.close()


```
