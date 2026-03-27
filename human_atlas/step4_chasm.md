# CHASM
# mount drives
```
sudo mount -t drvfs S: /mnt/s
```

# set up container for non-interactive session
```
SCRATCH1=/mnt/c/scratch
github=/mnt/g/github_repository/multivi_kidney
docker run -d --name chasm \
--workdir $HOME \
-v $github:$HOME/multivi_kidney \
-v /mnt/s:$HOME/data \
-v /mnt/g/reference:$HOME/reference \
-v $HOME:$HOME \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/mnt/c/scratch" \
-v /mnt/c/scratch:$HOME/scratch \
p4rkerw/sctools:R4.3.2
```

```
# prepare barcode annotation files for all samples
library(data.table)
library(dplyr)
library(tidyr)
library(here)
library(stringr)
anno <- fread(here("scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv"))
barcodes <- str_split(anno$barcode, pattern="_", simplify=TRUE)[,1]
library_id <- gsub("_multi", "", anno$sample)
library_id <- gsub("_atac", "", library_id)

anno <- data.frame(barcode = barcodes, library_id = library_id, celltype = anno$celltype)
library_ids <- unique(anno$library_id)

lapply(library_ids, function(library_id_sel) {
  anno %>%
    dplyr::filter(library_id %in% library_id_sel) %>%
    fwrite(here("scratch/kidney_10k_h5ad_output/chasm",library_id_sel,"bcanno.csv"))
})


```

# set up pipeline script
```
# pipeline_mvi_kidney.sh
function run_step1_count_frags() {
	library_id_sel=$1
	fragment_file=$2
	snap_peaks_csv=$3
	barcode_csv=$4
	outputdir=$5
	
	Rscript multivi_kidney/chasm/step1_count_frags.R $library_id_sel $fragment_file $snap_peaks_csv $barcode_csv $outputdir
}
export -f run_step1_count_frags

function run_step2_aggregate_counts() {
	library_id_sel=$1
	datadir=$2
	feature_count_mat=$3
	cytoband=$4
	
	Rscript multivi_kidney/chasm/step2_aggregate_counts.R $library_id_sel $datadir $feature_count_mat $cytoband
}
export -f run_step2_aggregate_counts

function run_step3_model_counts() {
	library_id=$1
	datadir=$2
	feature_counts=$3
	
	Rscript multivi_kidney/chasm/step3_model_counts.R $library_id $datadir $feature_counts
}
export -f run_step3_model_counts

###############################################################################
barcode_csv=scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv
snap_peaks_csv=scratch/kidney_10k_h5ad_output/peaks/peaks.csv
barcode_csv=scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv
outputdir=scratch/kidney_10k_h5ad_output/chasm

# prepare atac-only jobs
ls data/cellranger_atac_counts/kidney/version_2.1/*/outs/fragments.tsv.gz > /tmp/fragment_files.txt
cut -f5 -d'/' /tmp/fragment_files.txt > /tmp/library_ids.txt
paste -d " " /tmp/library_ids.txt /tmp/fragment_files.txt > /tmp/run_file.txt
awk -v arg1="$snap_peaks_csv" -v arg2="$barcode_csv" -v arg3="$outputdir" '{print $0, arg1, arg2, arg3}' /tmp/run_file.txt > /tmp/parallel_jobs1.txt

# prepare 10x multiomes jobs
ls data/cellranger_multi_counts/*/outs/atac_fragments.tsv.gz > /tmp/fragment_files.txt
 cut -f3 -d'/' /tmp/fragment_files.txt > /tmp/library_ids.txt
paste -d " " /tmp/library_ids.txt /tmp/fragment_files.txt > /tmp/run_file.txt
awk -v arg1="$snap_peaks_csv" -v arg2="$barcode_csv" -v arg3="$outputdir" '{print $0, arg1, arg2, arg3}' /tmp/run_file.txt > /tmp/parallel_jobs2.txt

# prepare kpmp jobs
ls data/KPMP_Data/single-nucleus-multiomes/*/atac_fragments.tsv.gz > /tmp/fragment_files.txt
cut -f4 -d'/' /tmp/fragment_files.txt | sed 's/_single-nuc_expression_matrix//g' > /tmp/library_ids.txt
paste -d " " /tmp/library_ids.txt /tmp/fragment_files.txt > /tmp/run_file.txt
awk -v arg1="$snap_peaks_csv" -v arg2="$barcode_csv" -v arg3="$outputdir" '{print $0, arg1, arg2, arg3}' /tmp/run_file.txt > /tmp/parallel_jobs3.txt

# count all fragment files
cat /tmp/parallel_jobs1.txt /tmp/parallel_jobs2.txt /tmp/parallel_jobs3.txt> /tmp/all_parallel_jobs.txt
# parallel -j2 --dry-run --colsep ' ' run_step1_count_frags :::: /tmp/all_parallel_jobs.txt
parallel -j2 --bar --colsep ' ' --bar run_step1_count_frags :::: /tmp/all_parallel_jobs.txt

feature_count_mats=(cellx_macs3-aggr_peak.rds)
ls $outputdir > /tmp/library_ids.txt
parallel -j1 --bar run_step2_aggregate_counts :::: /tmp/library_ids.txt ::: $outputdir ::: ${feature_count_mats[@]} ::: "chromosome"
parallel -j1 --bar run_step2_aggregate_counts :::: /tmp/library_ids.txt ::: $outputdir ::: ${feature_count_mats[@]} ::: "chromosome_arm"

feature_counts=(cellx_macs3-aggr_peak_chromosome_arm.csv)
parallel -j3 run_step3_model_counts :::: /tmp/library_ids.txt ::: $outputdir ::: ${feature_counts[@]}

feature_counts=(cellx_macs3-aggr_peak_chromosome.csv)
parallel -j3 run_step3_model_counts :::: /tmp/library_ids.txt ::: $outputdir ::: ${feature_counts[@]}
```

# execute pipeline
```
docker exec -i chasm bash < $github/pipeline_mvi_kidney.sh
```
