# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --restart unless-stopped \
# --workdir $HOME \
# --name ragegex \
# -v /mnt/y:$HOME/y \
# -v /mnt/j:$HOME/j \
# -v /mnt/g:$HOME/g \
# -v /mnt/e:$HOME/e \
# -v /mnt/d/garyw:$HOME/garyw \
# -v /mnt/v:$HOME/v \
# -v $HOME:$HOME \
# -v /mnt/y:$RHOME/y \
# -v /mnt/j:$RHOME/j \
# -v /mnt/g/:$RHOME/g \
# -v /mnt/e:$RHOME/e \
# -v /mnt/d/garyw:$RHOME/garyw \
# -v /mnt/v:$RHOME/v \
# -v /var/run/docker.sock:/var/run/docker.sock \
# -e PASSWORD=garywang \
# -e DISABLE_AUTH=TRUE \
# -p 8787:8787 \
# garywang7/ragemultiome:1.0.4

library(tidyverse)
library(here)
library(fs)
library(data.table)
library(GenomicRanges)
library(arrow)
library(patchwork)

future::plan("multicore", workers = 6)
options(future.globals.maxSize= 100*1024^3) # 100G

#### Define directories and parameters ####
bin_size_Mb <- 1

project_dir <- here("garyw", "RAGE24")
data_dir <- here(project_dir, "nb_model")
nb_dir <- here(data_dir,"result_aggr","negative_binomial")

output_dir <- here(project_dir,"nb_model","result_aggr",paste0(bin_size_Mb,"_Mb_bins"))
plot_dir <- here(project_dir,"publication/figures/figure2")
dir_create(here(plot_dir), recurse = TRUE)

age_col <- c("16" = "#b8d8ba", "30" = "#d9dbbc", "56" = "#dbac95", "82" = "#555b6e")
celltype_cols2 <- c(
  # Proximal tubule lineage
  "PT"           = "#789C8A",  # darker, higher contrast vs PT-MT
  "PT-injured"  = "#C96A4A",
  
  # Thin limb
  "TL1"  = "#E39A6F",
  "TL2"        = "#AFC7BC",
  
  # Non-PT cells (grey scale)
  "C-TAL"        = "#CFCFCF",
  "M-TAL"        = "#C4C4C4",
  "DCT"          = "#B9B9B9",
  "CNT"          = "#B0B0B0",
  "CNT-PC"       = "#A8A8A8",
  
  "PC"           = "#A0A0A0",
  "IC-A"         = "#9A9A9A",
  "IC-B"         = "#949494",
  
  "POD"          = "#8E8E8E",
  "PEC"          = "#888888",
  "EC"           = "#828282",
  "VSMC/P"       = "#7C7C7C",
  
  "FIB"          = "#767676",
  "IMM"          = "#707070"
)
#### Prep files ####
# Read result. Use arrow to load the data is much faster than fread
setDTthreads(threads = 12)
res.v2 <- read_feather(file = here(output_dir, "bin_depth_per_cell_Feb_02_2026.arrow"))
setDT(res.v2)
gc()

# Sort the bins
all.bins <- unique(res.v2$bins)
chr <- str_extract(all.bins, "chr[0-9XYM]+")
start <- as.integer(str_extract(all.bins, "(?<=:)\\d+"))       # Extract number after ":"
end <- as.integer(str_extract(all.bins, "(?<=-)\\d+"))         # Extract number after "-"
bins.gr <-  GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end))%>%
  # Sort chromosomes
  GenomeInfoDb::sortSeqlevels() %>%
  # Sort within chromosomes
  sort(ignore.strand=TRUE) 
bins.sorted <- as.character(bins.gr)
chromosome_names <- seqlevels(bins.gr)

# Read annotation
anno <- fread(file = here(output_dir, "cell_annotation_Feb_02_2026.csv"), nThread = 10)

# Factor cell types
anno$celltype <- factor(anno$celltype, levels = 
                          c("PT","PT-injured", "TL1", "TL2","C-TAL","M-TAL","DCT","CNT","CNT-PC","PC","IC-A","IC-B","EC","FIB","POD","PEC","VSMC_P","IMM"))

# factor bin and chromosome
res.v2 <- res.v2 %>%
  mutate(
    bins = factor(bins, levels = bins.sorted),
    chromosome = factor(chromosome, levels = chromosome_names)
  ) 

#### 1. Chromosome coverage ####
# Used bin-size = 1Mb for this part.
dir_create(here(plot_dir, "coverage"))
chrom.bin.number <- data.table("chromosome" = seqnames(bins.gr)@values, 
                               "chrom_bins" = seqnames(bins.gr)@lengths)

# Coverage summary                               
coverage.sum <- res.v2 %>%
  dplyr::select(barcode_atac_aggr, bins, chromosome, read_depth) %>%
  left_join(dplyr::select(anno, barcode_atac_aggr, celltype, library_type)) %>%
  group_by(barcode_atac_aggr, library_type, celltype, chromosome) %>%
  summarise(coverage = sum(read_depth > 0)) %>%
  mutate(coverage.pct = coverage/chrom.bin.number[chromosome == chromosome, chrom_bins])

write_csv(coverage.sum, here(output_dir, "coverage_summary.csv"))

## Color palette
method_col <- c(
  "DEFND" = "#D81B60",  # vivid magenta (warm, not orange)
  "Multiome"    = "#1F4E79"   # deep navy (cold)
)

## Whole genome
cov_stats_whole <- coverage.sum %>%
  group_by(barcode_gex, library_type) %>%
  summarise(chrom.cov = sum(coverage)) %>%
  mutate(coverage.pct = chrom.cov/sum(chrom.bin.number$chrom_bins))
write_csv(cov_stats_whole, here(output_dir, "coverage_whole_genome.csv"))

# check coverage stats
cov_stats_whole %>%
  group_by(library_type) %>%
  summarise(median = median(coverage.pct),
            mean = mean(coverage.pct))

p.cov.whole.genome <- cov_stats_whole %>%
  ggplot(aes(x = coverage.pct, color = library_type, fill = library_type, group = library_type))+
  geom_density(alpha = 0.1)+
  geom_histogram(aes(y = after_stat(density)), position = "identity",alpha = 0.3, bins= 100, 
                 linewidth = 0.3)+
  scale_x_continuous(labels = scales::percent, limits = c(0.45,1))+
  labs(
    title = paste0("Coverage of genome: ", bin_size_Mb, " Mb bins with reads"),
    x = "Coverage (%)"
  )+
  scale_color_manual(values = method_col) +
  scale_fill_manual(values = method_col) +
  ggpubr::theme_pubr()
ggsave(here(plot_dir, "coverage", "whole genome coverage.pdf"), plot = p.cov.whole.genome,
       height = 4, width = 3.5)

## Chromosome X
cov_stats_chrX <- coverage.sum %>%
  filter(chromosome == "chrX") %>%
  group_by(barcode_gex, library_type) %>%
  summarise(chrom.cov = sum(coverage)) %>%
  mutate(coverage.pct = chrom.cov/sum(chrom.bin.number[chromosome == "chrX", ]$chrom_bins))

p.cov.chrX <- cov_stats_chrX %>%
  ggplot(aes(x = coverage.pct, color = library_type, fill = library_type, group = library_type))+
  geom_density(alpha = 0.1)+
  geom_histogram(aes(y = after_stat(density)), position = "identity",alpha = 0.3, bins = 79)+
  scale_color_manual(values = method_col) +
  scale_fill_manual(values = method_col) +
  scale_x_continuous(labels = scales::percent)+
  labs(
    title = paste0("Coverage of chrX: ", bin_size_Mb, " Mb bins with reads"),
    x = "Coverage (%)"
  )+
  ggpubr::theme_pubr()

ggsave(here(plot_dir, "coverage", "chrX coverage.pdf"), plot = p.cov.chrX,height = 3, width = 6)

## Chromosome Y
cov_stats_chrY <- coverage.sum %>%
  filter(chromosome == "chrY") %>%
  group_by(barcode_gex, library_type) %>%
  summarise(chrom.cov = sum(coverage)) %>%
  mutate(coverage.pct = chrom.cov/sum(chrom.bin.number[chromosome == "chrY", ]$chrom_bins))

p.cov.chrY <- cov_stats_chrY %>%
  ggplot(aes(x = coverage.pct, color = library_type, fill = library_type, group = library_type))+
  geom_density(alpha = 0.1)+
  geom_histogram(aes(y = after_stat(density)), position = "identity",alpha = 0.3, bins = 62)+
  scale_color_manual(values = method_col) +
  scale_fill_manual(values = method_col) +
  scale_x_continuous(labels = scales::percent)+
  labs(
    title = paste0("Coverage of chrY: ", bin_size_Mb, " Mb bins with reads"),
    x = "Coverage (%)"
  )+
  ggpubr::theme_pubr()

ggsave(here(plot_dir, "coverage", "chrY coverage.pdf"), plot = p.cov.chrY,
       height = 3, width = 6)

## combine coverage plots of chromosome X and chrY 
p.cov.chrX.chrY <- (p.cov.chrX +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none"))/
  (p.cov.chrY +
     theme(legend.position = "none"))
ggsave(here(plot_dir, "coverage", "chrX chrY coverage.pdf"), plot = p.cov.chrX.chrY,
       height = 5, width = 6)

