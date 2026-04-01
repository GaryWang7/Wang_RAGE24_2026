# ==============================================================================
# Example Docker launch command for the analysis environment:
# 
# RHOME=/home/rstudio
# docker run -it \
#   --cpus=14 \
#   --memory="600g" \
#   --restart unless-stopped \
#   --workdir $HOME \
#   --name ragegex \
#   -v /mnt/y:$HOME/y \
#   -v /mnt/j:$HOME/j \
#   -v /mnt/g:$HOME/g \
#   -v /mnt/e:$HOME/e \
#   -v /mnt/d/garyw:$HOME/garyw \
#   -v /mnt/v:$HOME/v \
#   -v $HOME:$HOME \
#   -v /mnt/y:$RHOME/y \
#   -v /mnt/j:$RHOME/j \
#   -v /mnt/g:$RHOME/g \
#   -v /mnt/e:$RHOME/e \
#   -v /mnt/d/garyw:$RHOME/garyw \
#   -v /mnt/v:$RHOME/v \
#   -v /var/run/docker.sock:/var/run/docker.sock \
#   -e PASSWORD=garywang \
#   -e DISABLE_AUTH=TRUE \
#   -p 8787:8787 \
#   garywang7/ragemultiome:1.0.4
# ==============================================================================

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

#### 2. Coverage plot between PT and PT-injured ####
res.nb <- read_csv(here(nb_dir, "negative_binomial_combined_Feb_09_2026.csv"))
anno.nb <- read_csv(here(nb_dir, "cell_annotation_Feb_02_2026.csv"))

# chromosomes
chroms <- unique(res.nb$chromosome) %>% gtools::mixedsort()
res.nb$chromosome <- factor(res.nb$chromosome, levels = chroms)

# Select chromosome 2, 5, X, Y , 10, 12, split by PT vs PT-injured
p <- res.nb %>%
  left_join(anno.nb[,c("barcode_atac_aggr","celltype")]) %>%
  filter(celltype %in% c("PT","PT-injured"),
         chromosome %in% c("chr2","chr5", "chrX","chrY", "chr10", "chr12")) %>%
  mutate(celltype = factor(celltype, levels = c("PT-injured", "PT"))) %>%
  group_by(chromosome, celltype) %>%
  mutate(
    scale_chrom_log = (chrom_log_norm - chrom_log_mode)/sd(chrom_log_norm)
  ) %>%
  ggplot()+
  geom_density(aes(x = scale_chrom_log, y = after_stat(scaled), fill = celltype, color = celltype), 
               alpha = 0.5) +
  facet_wrap(vars(chromosome))+
  # The color palette is later changed in Inkscape.
  scale_fill_manual(values = celltype_cols2)+
  scale_color_manual(values = celltype_cols2)+
  coord_cartesian(xlim = c(-4,4))+
  ggpubr::theme_pubr()
ggsave(here(plot_dir,"chrom read normalized PT vs PTinjured selected chroms.pdf"), width = 6, height = 4, plot = p)

#### 3. PT cells with any CNA in different ages ####
# Read the files again.
res.nb <- read_csv(here(nb_dir, "negative_binomial_combined_Feb_09_2026.csv"))
anno.nb <- read_csv(here(nb_dir, "cell_annotation_Feb_02_2026.csv"))

# chromosomes
chroms <- unique(res.nb$chromosome) %>% gtools::mixedsort()
res.nb$chromosome <- factor(res.nb$chromosome, levels = chroms)

# Function to  calculate CNA in each chromosome
calc_CNA <- function(res, padj_thre = 1e-4, cna_scale = "absolute", cna_thre = 1){
  res.df <- res %>%
    ungroup() %>%
    mutate(
      # Readjust p value
      chrom_p_val_adj_mode = p.adjust(chrom_p_val_mode, method = "BH"),
      # Identify cna
      any_CNA = case_when(
        chrom_p_val_adj_mode < padj_thre ~ TRUE,
        .default = FALSE
      ),
      any_CNA_type = case_when(
        any_CNA ==TRUE & cna_size > 0 ~ "Gain",
        any_CNA ==TRUE & cna_size < 0 ~ "Loss",
        .default = "None"
      ),
      large_CNA = case_when(
        cna_scale == "absolute" & any_CNA ==TRUE & abs(cna_size) >= cna_thre ~ TRUE, 
        cna_scale == "relative" & any_CNA ==TRUE & abs(rel_cna_size) >= cna_thre ~ TRUE, 
        .default = FALSE
      ),
      large_CNA_type = case_when(
        large_CNA==TRUE & cna_size > 0 ~ "Gain",
        large_CNA==TRUE & cna_size < 0 ~ "Loss",
        .default = "None"
      )
    ) %>%
    dplyr::select(barcode_atac_aggr, contains("CNA"), chromosome, chrom_category,
                  chrom_p_val_adj_mode, frag_lib_size) %>%
    ungroup() %>%
    distinct() %>%
    left_join(anno.nb) %>%
    dplyr::select(barcode_atac_aggr, celltype, rat_id, age_wks,
                  contains("CNA"), contains("chrom"), library_type,
                  frag_lib_size)
}

CNA.res.PT <- res.nb %>%
  left_join(anno.nb[,c("barcode_atac_aggr", "library_id_gex", "celltype")]) %>%
  filter(celltype %in% c("PT","PT-injured")) %>%
  ungroup() %>%
  calc_CNA(padj_thre = 1e-4, cna_scale = "relative", cna_thre = 0.5)

CNA.per.PT <- CNA.res.PT %>%
  # Remove rats with abnormal distribution on chr12
 filter(!rat_id %in% c("J10","E05","N14")) %>%
  group_by(barcode_atac_aggr, celltype, rat_id, age_wks,
           library_type, frag_lib_size) %>%
  summarise(
    have_large_CNA = as.integer(any(large_CNA==TRUE)),
    have_CNA = as.integer(any(any_CNA == TRUE)),
    .groups = "drop"
  )

CNA.per.PT.prop <- CNA.per.PT %>%
  group_by(age_wks) %>%
  summarise(
    large_CNA_prop = mean(have_large_CNA),
    any_CNA_prop = mean(have_CNA)
  )

# Any CNA--Cochran–Armitage trend test of binary outcomes
# install.packages("DescTools")
prop_tab <- with(CNA.per.PT, table(age_wks, have_CNA))
CAT <- DescTools::CochranArmitageTest(prop_tab, alternative = "one.sided")
CAT.df <- data.frame(unlist(CAT)) %>% tibble::rownames_to_column("rownames")
write_csv(CAT.df, here(plot_dir,"CNA_prop", "Proportion age PT any CNA test.csv"))

# Any CNA plot
p <- CNA.per.PT.prop %>%
  mutate(age_wks = factor(age_wks, levels = c(16,30,56,82)),
         CNA_pct = any_CNA_prop*100) %>%
  ggplot(aes(x = age_wks, y = CNA_pct, fill = age_wks)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f%%", CNA_pct)), vjust = -0.3, size = 3.5
  )+
  annotate("text", x = 0.6, y = 10, hjust = 0,
           label = paste0(
             CAT$method, "\n",
             "p = ", signif(CAT$p.value, 4)
           ))+
  xlab("age weeks") +
  ylab("Percentage (%) ") +
  ggtitle("Percentage of proximal tubule cells with any CNA")+
  scale_fill_manual(values = age_col)+
  #coord_cartesian(ylim = c(5, 8)) + 
  ggpubr::theme_pubr()
ggsave(here(plot_dir,"CNA_prop","all PT cells_CNA with age_mode.pdf"), height = 3, width = 2.5, plot = p)

#### 4. CNA in PT vs PT-injured cells ####
# Use data imported in section 3.
# Test chi-square for any CNA
tab <- table(celltype = CNA.per.PT$celltype, CNA = CNA.per.PT$have_CNA)
test <- chisq.test(tab, correct = FALSE)

# Plot any CNA in PT vs PT-injured
p <- CNA.per.PT.prop %>%
  mutate(CNA_pct = any_CNA_prop*100) %>%
  ggplot(aes(x = celltype, y = CNA_pct, fill = celltype)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f%%", CNA_pct)), vjust = -0.3, size = 3.5
  )+
  annotate("text", x = 0.6, y = 10, hjust = 0,
           label = paste0(
             test$method, "\n",
             "p = ", signif(test$p.value, 4)
           ))+
  xlab("celltype") +
  ylab("Percentage (%) ") +
  ggtitle("Percentage of PT vs PT-injured to develop any CNA")+
  scale_fill_manual(values = celltype_cols2)+
  #coord_cartesian(ylim = c(5, 8)) + 
  ggpubr::theme_pubr()
ggsave(here(plot_dir,"CNA_prop","PT vs PT injured CNA_mode.pdf"), height = 3, width = 2.5, plot = p)