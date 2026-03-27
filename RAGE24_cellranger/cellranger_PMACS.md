- [Intro](#intro)
- [Make references](#make-references)
- [Align 10X Data on PMACS with cellranger and cellranger-atac](#align-10x-data-on-pmacs-with-cellranger-and-cellranger-atac)
  - [ATAC/LAND](#atacland)
  - [GEX (RNA)](#gex-rna)
  - [Check if cellranger ran successfully](#check-if-cellranger-ran-successfully)
    - [RNA portion check](#rna-portion-check)
    - [ATAC portion check](#atac-portion-check)
- [Transfer and Inspect Alignment Results](#transfer-and-inspect-alignment-results)
  - [Transferring web summaries](#transferring-web-summaries)
  - [Transferring all cellranger runs for individual libraries](#transferring-all-cellranger-runs-for-individual-libraries)
- [Aggregate the GEX libraries](#aggregate-the-gex-libraries)
  - [Aggregate RNA(GEX) libraries](#aggregate-rnagex-libraries)
  - [Transfer and inspect RNA aggr results](#transfer-and-inspect-rna-aggr-results)
# Intro
This markdown file contains commands to run cellranger with step1_cellranger_atac_count, step1_cellranger_gex_count, step2_cellranger_aggr_gex, as well as commands to transfer the data back and forth between PMACS and local workstation. 
**Make sure to generate ATAC and GEX sample sheets first.**

# Make references
Refer to "make_rat_ref.sh" to make rat references.
# Align 10X Data on PMACS with cellranger and cellranger-atac
```sh
# Check total number of fastq.gz files in Fastq_RAGE24/ except the ones in Fastq_RAGE24/Extrafiles
find /home/garyw7/Fastq_RAGE24/ -path /home/garyw7/Fastq_RAGE24/Extrafiles -prune -o -type f -name "*.fastq.gz" -print | wc -l  #There are 2563 files
```
## ATAC/LAND
Note that if some ATAC libraries do not have I1 while others do, cellranger-atac will return an error message. In this case, we will have to rename all I1 libraries so cellranger-atac cannot find them. 
```sh
# split library meta to select ATAC libraries. The library_ids in the meta file are unique each line.
awk -F '\t' 'NR > 1 {print $1" "$2" "$3}' /home/garyw7/RAGE24/cellranger_meta_ATAC.tsv > /home/garyw7/RAGE24/library_ids_ATAC.txt
# (Optional) if some of the Fastq libraries do not contain I1 reads, we have to rename(remove) all I1 reads so that cellranger-atac does not get confused
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  echo "Processing Sample ID: $sample_id, Library ID: $library_id, Cellranger Samples: $cellranger_sample_name"
  fastq_dir=/home/garyw7/Fastq_RAGE24/$sample_id
  find "$fastq_dir" -type f -name "RAGE24*-ATAC-*_I1_*.fastq.gz" | while read -r file; do
    #echo "Renaming $file"
    # Extract the directory and old file name
    dir=$(dirname "$file")
    filename=$(basename "$file")
    new_filename="TempTag$filename"
    # Rename the files
    #mv "$file" "$dir/$new_filename"
    echo "Renamed $file to $dir/$new_filename"
  done
done < /home/garyw7/RAGE24/library_ids_ATAC.txt

# swap sample_id,library id and cellranger_sample_name. The library_ids should be unique in this file.
script_dir=/project/parkercwlab/gary/scripts/RAGE
mkdir -p $script_dir/cellranger/ATAC
log_dir=/project/parkercwlab/gary/logs/cellranger/2024_10_16/
mkdir -p $log_dir

# Process library_id one by one
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  echo "Processing Sample ID: $sample_id, Library ID: $library_id, Cellranger Samples: $cellranger_sample_name"

  # Swap in sample_id and library_id in the script
  sed "s/<sample_id>/$sample_id/g; s/<library_id>/$library_id/g; s/<cellranger_sample_name>/$cellranger_sample_name/g" $script_dir/step1_cellranger_atac_count_PMACS.sh > $script_dir/cellranger/ATAC/$library_id.step1_cellranger_atac_count.sh

  # Submit batch job with 16 cores and 250GB RAM
  bsub \
  -q normal \
  -e $log_dir/$library_id.log.error \
  -o $log_dir/$library_id.log.out \
  -J $library_id.cellranger-atac \
  -n 20 \
  -M 250000 \
  -R "rusage[mem=250000] span[hosts=1]" \
  sh $script_dir/cellranger/ATAC/$library_id.step1_cellranger_atac_count.sh
done < /home/garyw7/RAGE24/library_ids_ATAC.txt
```
## GEX (RNA)
```sh
# split library meta to select RNA libraries
awk -F '\t' 'NR > 1  {print $1" "$2" "$3}' /home/garyw7/RAGE24/cellranger_meta_RNA.tsv > /home/garyw7/RAGE24/library_ids_RNA.txt

# swap fastq_sample_name and sample_id
script_dir=/project/parkercwlab/gary/scripts/RAGE
mkdir -p $script_dir/cellranger/GEX
log_dir=/project/parkercwlab/gary/logs/cellranger/2024_10_16/

# Process library_id one by one
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  echo "Processing Sample ID: $sample_id, Library ID: $library_id, Cellranger Samples: $cellranger_sample_name"
  
  # Swap in sample_id and library_id in the script
  sed "s/<sample_id>/$sample_id/g; s/<library_id>/$library_id/g; s/<cellranger_sample_name>/$cellranger_sample_name/g" $script_dir/step1_cellranger_gex_count_PMACS.sh > $script_dir/cellranger/GEX/$library_id.step1_cellranger_gex_count.sh

  # Submit batch job with 16 cores and 250GB RAM
  bsub \
  -q normal \
  -e $log_dir/$library_id.log.error \
  -o $log_dir/$library_id.log.out \
  -J $library_id.cellranger-gex \
  -n 20 \
  -M 200000 \
  -R "rusage[mem=200000] span[hosts=1]" \
  sh $script_dir/cellranger/GEX/$library_id.step1_cellranger_gex_count.sh
done < /home/garyw7/RAGE24/library_ids_RNA.txt
```
## Check if cellranger ran successfully
If cellranger runs successfully, there will be an "outs" folder.
### RNA portion check
```sh
failed=/home/garyw7/RAGE24/library_ids_RNA_failed.txt
>"$failed"
# Check for each library. Gather all failed library in a file.
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  libdir=/home/garyw7/RAGE24/cellranger/$sample_id/$library_id
  if [ -d "$libdir/outs" ]; then
    echo "Cellranger pipeline successful for $library_id."
  else
    echo "The 'outs' directory does not exist in $libdir."
    echo "$sample_id $library_id $cellranger_sample_name" >> "$failed" #append to a file
  fi
done < /home/garyw7/RAGE24/library_ids_RNA.txt

### Rerun the cellranger on failed files
script_dir=/project/parkercwlab/gary/scripts/RAGE
mkdir -p $script_dir/cellranger/GEX
log_dir=/project/parkercwlab/gary/logs/cellranger/2024_08_14/
mkdir -p $log_dir
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  echo "Processing Sample ID: $sample_id, Library ID: $library_id, Cellranger Samples: $cellranger_sample_name"
  # Swap in sample_id and library_id in the script
  sed "s/<sample_id>/$sample_id/g; s/<library_id>/$library_id/g; s/<cellranger_sample_name>/$cellranger_sample_name/g" $script_dir/step1_cellranger_gex_count_PMACS.sh > $script_dir/cellranger/GEX/$library_id.step1_cellranger_gex_count.sh

  # Submit batch job with 24 cores and 300GB RAM
  bsub \
  -q normal \
  -e $log_dir/$library_id.log.error \
  -o $log_dir/$library_id.log.out \
  -J $library_id.cellranger-gex \
  -n 24 \
  -M 300000 \
  -R "rusage[mem=300000] span[hosts=1]" \
  sh $script_dir/cellranger/GEX/$library_id.step1_cellranger_gex_count.sh
done < $failed
```

### ATAC portion check
```sh
failed=/home/garyw7/RAGE24/library_ids_ATAC_failed.txt
>"$failed"
# Check for each library. Gather all failed library in a file.
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  libdir=/home/garyw7/RAGE24/cellranger/$sample_id/$library_id
  if [ -d "$libdir/outs" ]; then
    echo "Cellranger pipeline successful for $library_id."
  else
    echo "The 'outs' directory does not exist in $libdir."
    echo "$sample_id $library_id $cellranger_sample_name" >> "$failed" #append to a file
  fi
done < /home/garyw7/RAGE24/library_ids_ATAC.txt

### Rerun the job on failed ones
script_dir=/project/parkercwlab/gary/scripts/RAGE
mkdir -p $script_dir/cellranger/ATAC
log_dir=/project/parkercwlab/gary/logs/cellranger/2024_10_20/
mkdir -p $log_dir
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  echo "Processing Sample ID: $sample_id, Library ID: $library_id, Cellranger Samples: $cellranger_sample_name"
  # Swap in sample_id and library_id in the script
  sed "s/<sample_id>/$sample_id/g; s/<library_id>/$library_id/g; s/<cellranger_sample_name>/$cellranger_sample_name/g" $script_dir/step1_cellranger_atac_count_PMACS.sh > $script_dir/cellranger/ATAC/$library_id.step1_cellranger_atac_count.sh

  # Submit batch job with 16 cores and 400GB RAM
  bsub \
  -q normal \
  -e $log_dir/$library_id.log.error \
  -o $log_dir/$library_id.log.out \
  -J $library_id.cellranger-atac \
  -n 24 \
  -M 250000 \
  -R "rusage[mem=250000] span[hosts=1]" \
  sh $script_dir/cellranger/ATAC/$library_id.step1_cellranger_atac_count.sh
done < $failed
```
# Transfer and Inspect Alignment Results
## Transferring web summaries
```sh
# To workstation
top_dir=garyw7@mercury.pmacs.upenn.edu:/home/garyw7/RAGE24/cellranger/
destination_path=/home/gary/garyw/scratch/cellranger
rsync -av -m --include='**/' --include='**/outs/*.html' --exclude='*' $top_dir/ $destination_path

# To NAS from workstation
top_dir=/home/gary/garyw/RAGE24/cellranger_flowcell_combine
destination_path=/home/gary/g/RAGE24/cellranger_flowcell_combine_scratch
rsync -av -m --no-perms --no-group --no-owner --omit-dir-times --progress --include='**/' --include='**/outs/*.html' --exclude='*' $top_dir/ $destination_path
```
## Transferring all cellranger runs for individual libraries
```sh
# To Workstation
# sync everything in cellranger /outs folders
top_dir=garyw7@mercury.pmacs.upenn.edu:/home/garyw7/RAGE24/cellranger
destination_path=/home/gary/garyw/RAGE24/cellranger_flowcell_combine
rsync -av -m --progress --include='**/' --include='**/outs/**' --exclude='*' $top_dir/ $destination_path/

# To NAS
top_dir=garyw7@mercury.pmacs.upenn.edu:/home/garyw7/RAGE24/cellranger
destination_path=/home/gary/g/RAGE24/cellranger_flowcell_combine
rsync -av -m --no-perms --no-group --no-owner --no-t --omit-dir-times --info=progress2 --include='**/' --include='**/outs/**' --exclude='*' $top_dir/ $destination_path/

# To External Drive
top_dir=garyw7@mercury.pmacs.upenn.edu:/home/garyw7/RAGE24/cellranger
destination_path=/mnt/f/RAGE24/cellranger_flowcell_combine
rsync -av -m --no-perms --no-group --no-owner --no-t --omit-dir-times --info=progress2 --include='**/' --include='**/outs/**' --exclude='*' $top_dir/ $destination_path/

# To NAS, using workstation data
top_dir=/home/gary/garyw/RAGE24/cellranger_flowcell_combine
destination_path=/home/gary/g/RAGE24/cellranger_flowcell_combine
rsync -av -m --no-perms --no-group --no-t --no-owner --omit-dir-times --progress --include='**/' --include='**/outs/**' --exclude='*' $top_dir/ $destination_path/
```
# Aggregate the GEX libraries
We don't perform this step for ATAC libraries since we have DEFND libraries, which do not have too many meaningful peaks but cellranger-aggr will redo peak calling.
## Aggregate RNA(GEX) libraries
```sh
# Generate an aggregation csv
bsub -Is -n 10 -M 100000 bash # Request for memory and CPUs

gex_aggr=/home/garyw7/RAGE24/gex_aggr.csv
echo "sample_id,molecule_h5" > "$gex_aggr"
while IFS=' ' read -r sample_id library_id cellranger_sample_name; do
  libdir=/home/garyw7/RAGE24/cellranger/$sample_id/$library_id
  h5dir=$libdir/outs/molecule_info.h5
  if [ -f "$h5dir" ]; then
    echo "File exists for $library_id."
    echo "$library_id,$h5dir" >> "$gex_aggr" #append to a file
  else
    echo "The h5 file does not exist in $libdir."
  fi
done < /home/garyw7/RAGE24/library_ids_RNA.txt
# Here we can manually remove the pancreas outputs (PS-KYC) from gex_aggr since we only focus on kidney cortex. 

# Run cellranger aggr
# create path variables
cellranger=/home/garyw7/wilson_lab/gary/tools/cellranger-8.0.1/cellranger
outputdir=/home/garyw7/RAGE24

# change directory
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir
change_directory $outputdir

# aggregate the output
$cellranger aggr \
--id=cellranger_aggr_gex_flowcell_combine \
--csv=$gex_aggr \
--normalize="none" \
--nosecondary
```
## Transfer and inspect RNA aggr results
```sh
# To workstation
top_dir=garyw7@mercury.pmacs.upenn.edu:/home/garyw7/RAGE24/cellranger_aggr_gex_flowcell_combine
destination_path=/home/gary/garyw/RAGE24/cellranger_aggr_gex_flowcell_combine
rsync -av -m --include='**/' --include='**/outs/**' --exclude='*' $top_dir/ $destination_path/
# To NAS
top_dir=garyw7@mercury.pmacs.upenn.edu:/home/garyw7/RAGE24/cellranger_aggr_gex_flowcell_combine
destination_path=/home/gary/g/RAGE24/cellranger_aggr_gex_flowcell_combine
rsync -av --no-perms --no-group --no-owner --no-t --omit-dir-times --progress --include='**/' --include='**/outs/**' --exclude='*' $top_dir/ $destination_path/
```

