########### For PMACS #############
# assign library_id
sample_id=<sample_id>
library_id=<library_id>
cellranger_sample_name=<cellranger_sample_name>
# library_id=RAGE24-C03-KYC-LN-01-Multiome-RNA-22FMVVLT4
# sample_id=RAGE24-C03-KYC-LN-01
# create path variables
cellranger=/home/garyw7/wilson_lab/gary/tools/cellranger-8.0.1/cellranger
outputdir=/home/garyw7/RAGE24/cellranger/$sample_id
fastq=/home/garyw7/RAGE24/Fastq_RAGE24/$sample_id
reference=/home/garyw7/wilson_lab/reference/GRCr8/GRCr8_mt_cellranger

# change directory to home
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir
change_directory $outputdir

# count the fastq files
${cellranger} count \
--sample=$cellranger_sample_name \
--id=$library_id \
--transcriptome=$reference \
--fastqs="$fastq" \
--chemistry=ARC-v1 \
--create-bam=true
#--localcores=16 \
#--localmem=250 \ # By default, cellranger will use all available cores and 90% of available memory