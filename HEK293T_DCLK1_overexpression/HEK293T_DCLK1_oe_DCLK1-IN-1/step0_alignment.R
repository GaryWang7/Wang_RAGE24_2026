library(Rsubread)
library(tidyverse)
library(here)
library(fs)

#### 1. Define paths ####
project_dir <- "D:/garyw/RAGE24/DCLK1_HEK293T_oe_DCLK1IN1"
data_dir <- here(project_dir, "Fastq")
fastq_dir <- here(data_dir, "trimmed")
ref_dir <- here("D:/garyw/RAGE24/reference/GRCh38p14")
bam_dir <- here(project_dir, "alignments")
result_dir <- here(project_dir, "counts")

dir_create(result_dir)
dir_create(bam_dir)

#### 2. Preprocessing reference genome ####
# Create an index file for the reference genome
# Take around 1 hour, depending on the computer.
buildindex(basename = here(ref_dir,"GRCh38p14_GFP"), 
           reference = here(ref_dir,"Homo_sapiens.GRCh38.dna.primary_assembly_modified_GFP.fa")) 

#### 3. Read and align fastqs to the genome ####
# meta data
meta <- read_csv(here(project_dir,"metadata_edited.csv")) %>%
  dplyr::rename("sample" = `Customer ID`)

# Create file names
samples <- read.table(here(project_dir, "metadata_edited.csv"), 
                      sep = ",", header = T) %>% pull("Customer.ID") %>%
  gtools::mixedsort()
reads1 <- sapply(samples, function(samp){
  dir_ls(path = here(fastq_dir),
         regexp = paste0('/',samp,'_R1.fastq.gz$'))
})
reads2 <- sapply(samples, function(samp){
  dir_ls(path = here(fastq_dir),
         regexp = paste0('/',samp,'_R2.fastq.gz$'))
})
output_files <- here(bam_dir, paste0(samples,".bam"))

# Align to genome. Takes 5~10 min per sample
align.stat <- Rsubread::align(
  index = here(ref_dir,"GRCh38p14_GFP"),
  readfile1 = reads1, readfile2 = reads2,
  input_format = "gzFASTQ",
  output_format = "BAM",
  output_file = output_files,
  nthreads = 12)

# save alignment statistics
data.table::fwrite(align.stat, file = here(bam_dir, "alignment_summary.csv"), col.names = TRUE, row.names = TRUE)

#### 4. Count mapped reads for genomic features ####
bam_files <- as.character(dir_ls(here(bam_dir), regexp = ".bam$"))
fc <- featureCounts(files = bam_files, 
                    annot.ext = here(ref_dir, "gencode.v47.primary_assembly.annotation_modified_GFP.gtf"),
                    GTF.attrType.extra = "gene_name",
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

# Save mapping metrics
stat <- fc$stat %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("V0") %>%
  janitor::row_to_names(row_number = 1) %>%
  dplyr::rename("sample" = "Status") %>%
  dplyr::mutate(across(-sample, as.numeric)) %>%
  rowwise() %>%
  mutate(Total = sum(c_across(-sample), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    sample = str_remove(sample, ".bam"),
    assign.percent = 100*Assigned/Total
  )
data.table::fwrite(stat, file = here(result_dir, "featureCount_metrics.csv"))