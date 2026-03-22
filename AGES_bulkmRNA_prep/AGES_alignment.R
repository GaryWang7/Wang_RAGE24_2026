# This code aligns and quantifies the bulk mRNA data from Shavlakadze et al.
# https://doi.org/10.1016/j.celrep.2019.08.043

library(Rsubread)
library(tidyverse)
library(here)
library(fs)

#### 1. Define paths ####
project_dir <- "V:/AgeRat/"
data_dir <- here(project_dir, "Data")
fastq_dir <- here(data_dir, "Fastq_trimmed")
ref_dir <- here(project_dir, "Reference")
bam_dir <- here(data_dir, "Alignments")
result_dir <- here(project_dir, "counts")

dir_create(result_dir)
dir_create(bam_dir)

#### 2. Preprocessing reference genome ####

# Create an index file for the reference genome
# Take around 1 hour, depending on the computer.
buildindex(basename = here(ref_dir,"GRCr8_mt"), 
           reference = here(ref_dir,"GCF_036323735.1_GRCr8_mt.fna")) 

#### 3. Read and align fastqs to the genome ####

# Create file names
reads1 <- dir_ls(path = here(fastq_dir),
                     regexp = '*_trimmed_1.fastq.gz$')
reads2 <- dir_ls(path = here(fastq_dir),
                 regexp = '*_trimmed_2.fastq.gz$')
samples <- read.table(here(data_dir, "Accession", "SRR_Acc_List.txt"), 
                      sep = "\t", header = F) %>% pull("V1")
output_files <- here(bam_dir, paste0(samples,".bam"))

# Align to genome. Takes 5-10 min per sample
Rsubread::align(
  index = here(ref_dir,"GRCr8_mt"),
  readfile1 = reads1, readfile2 = reads2,
  input_format = "gzFASTQ",
  output_format = "BAM",
  output_file = output_files,
  nthreads = 12)

#### 4. Count mapped reads for genomic features ####
bam_files <- as.character(dir_ls(here(bam_dir), regexp = ".bam$"))
fc <- featureCounts(files = bam_files, 
                    annot.ext = here(ref_dir, "cellranger_filtered_GCF_036323735.1_GRCr8_mt.gtf"),
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE,
                    nthreads = 12,
                    verbose = TRUE)
saveRDS(fc, file = here(result_dir, "featureCounts.rds"))

# Save the counts
raw_counts <- fc$counts
colnames(raw_counts) <- stringr::str_remove(colnames(raw_counts), 
                                            pattern = ".bam")
data.table::fwrite(raw_counts, file = here(result_dir, "counts.csv"),
                   row.names = TRUE)

# Save annotations used
gtf <- fc$annotation
data.table::fwrite(gtf, file = here(result_dir, "genes_for_counting.csv"))

# Save assign metrics
metrics <- fc$stat
colnames(metrics) <- stringr::str_remove(colnames(metrics), 
                                         pattern = ".bam")
data.table::fwrite(metrics, file = here(result_dir, "featureCount_metrics.csv"))
