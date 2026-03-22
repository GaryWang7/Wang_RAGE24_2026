########### For PMACS #############
# assign library_id
sample_id=<sample_id>
library_id=<library_id>
cellranger_sample_name=<cellranger_sample_name>
# sample_id=RAGE24-C03-KYC-LN-01
# library_id=RAGE24-C03-KYC-LN-01-Multiome-ATAC
# cellranger_sample_name=RAGE24-C03-KYC-LN-01-Multiome-ATAC-HTWY3DMXY,RAGE24-C03-KYC-LN-01-Multiome-ATAC-22FMT7LT4-01,RAGE24-C03-KYC-LN-01-Multiome-ATAC-22FMT7LT4-03,RAGE24-C03-KYC-LN-01-Multiome-ATAC-22FMT7LT4-04,RAGE24-C03-KYC-LN-01-Multiome-ATAC-22FMT7LT4-02
# create path variables
cellranger_atac=/home/garyw7/wilson_lab/gary/tools/cellranger-atac-2.1.0/cellranger-atac
outputdir=/home/garyw7/RAGE24/cellranger/$sample_id
# outputdir=/home/garyw7/RAGE24/cellranger/test/$sample_id # Test output directory
fastq=/home/garyw7/RAGE24/Fastq_RAGE24/$sample_id
reference=/home/garyw7/wilson_lab/reference/GRCr8/GRCr8_mt_cellrangerarc

# change directory to home
change_directory () {
  cd $1
}

# set up output dir
mkdir -p $outputdir
change_directory $outputdir

# count the fastq files
${cellranger_atac} count \
--sample=$cellranger_sample_name \
--id=$library_id \
--reference=$reference \
--fastqs="$fastq" \
--chemistry=ARC-v1
#--localcores=16 \
#--localmem=250 \ # By default, cellranger will use all available cores and 90% of available memory