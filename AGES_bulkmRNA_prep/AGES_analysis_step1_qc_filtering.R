library(ComplexHeatmap)
library(clusterProfiler)
library(AnnotationDbi)
library(DESeq2)
library(dendextend)
library(airway)
library(tidyverse)
library(readxl)
library(biomaRt)
library(edgeR)
library(here)
library(RColorBrewer)
library("factoextra")
library(ggrepel)
library(enrichplot)
library(fs)
library(vsn)
library(pheatmap)
library(rrcov)

#### 1. Define paths ####
project_dir <- "V:/AgeRat/"
counts_dir <- here(project_dir, "counts")
plots_dir <- here(project_dir, "plots","step1")
analysis_dat_dir <- here(project_dir,"analysis_data")
dir_create(plots_dir)
dir_create(analysis_dat_dir)
#### 2. Read and wrangle data ####

# Raw feature count results
fc <- readRDS(here(counts_dir, "featureCounts.rds"))

# Renamed counts
raw_counts <- data.table::fread(here(counts_dir, "counts.csv")) %>%
  tibble::column_to_rownames("V1")

# Metrics of assignment to gtf file
metrics <- data.table::fread(here(counts_dir, "featureCount_metrics.csv"))

# Metadata
meta <- data.table::fread(here(project_dir, "Info", "metadata.csv")) %>%
  dplyr::filter(tissue == "kidney") %>%
  dplyr::select(Run, AGE, sex) %>%
  dplyr::rename("sample" = "Run") %>%
  mutate(age_months = factor(as.numeric(
    str_remove_all(AGE, pattern = "Mo")))) %>%
  tibble::column_to_rownames("sample") %>%
  arrange(age_months)

# Make sure the COLUMNS of count matrix are of the same order as 
# ROWS in metadata (coldata)
raw_counts <- raw_counts[,rownames(meta)]

#### 3. Perform differential expression analysis ####
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = meta,
                              design = ~ age_months)

# Pre-filtering genes for visualization and speed for DESeq
# Require genes to have at least 5 counts in at least 8 samples (smallest group size)
# There is no absolute reason why 5, but 5-10 should be reasonable choice.
smallestGroupSize <- min(table(meta$age_months))
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

# Relevel the reference 
dds$age_months <- relevel(dds$age_months, ref = "9")

# Run DESeq
dds <- DESeq(dds)
#res <- results(dds, contrast = c("age_months", 27, 9), alpha = 0.1)

#### 4. Exploring and QC the results ####
# For PCA, we need to use rlog or vst transformation to avoid 
# the influence of strongly expressed genes
# See https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# and https://hbctraining.github.io/Intro-to-DGE/lessons/03_DGE_QC_analysis.html
vsd <- vst(dds, blind = FALSE) # Setting blind to false to incorporate already estimated dispersion

# Compare the removal of dependence of variance on the mean
pdf(here(plots_dir,"mean_vs_sd_plots.pdf"))
meanSdPlot(assay(normTransform(dds)))
meanSdPlot(assay(vsd))
dev.off()
# Continue with vsd as the variance is more stable

# classical PCA
pcaData <- plotPCA(vsd, intgroup = "age_months", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_prefilter_plot <- ggplot(pcaData, aes(PC1, PC2, color=age_months)) +
  geom_point(size=4, alpha = 0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggpubr::theme_classic2()+
  ggsci::scale_color_nejm()

ggsave(here(plots_dir,"PCA plot before filtering.png"), PCA_prefilter_plot,
       width = 6)

# PC1 variance
pcaData_summary <- pcaData %>%
  group_by(age_months) %>%
  mutate(PC1_mean = mean(PC1),
         PC1_sd = sd(PC1))

# Plot PC1, as PC1 explains the majority of variance
PCA_prefilter_PC1_plot <- pcaData %>%
  mutate(age_months = factor(age_months, levels = sort(unique(meta$age_months))))%>%
  ggplot(aes(x = age_months, y = PC1, color=age_months)) +
  geom_jitter(size=3, alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  xlab(paste0("Age (months)")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  ggpubr::theme_classic2()+
  ggsci::scale_color_nejm()

ggsave(here(plots_dir,"PC1 boxplot before filtering.png"), PCA_prefilter_PC1_plot,
       height = 8, width = 6)

#### 5. Filtering of data ####
# Based on the unfiltered classical PCA data, 6 month
# old rats are not very similar to each other, and the variance of PC1 is much higher than other age groups. So we will drop the 6 month old rats for the downstream analysis.
# Also, one rat from age group 9 is having an outlier PC1 value--SRR8705717
# Therefore, we drop these samples.
samples_to_drop <- meta %>%
  tibble::rownames_to_column("name") %>%
  filter(name == "SRR8705717"| age_months == 6) %>%
  pull(name)

#### 6. Repeat the deseq analysis above with filtered data ####
# Repeat the above data construction 
counts_filt <- raw_counts %>% dplyr::select(-one_of(samples_to_drop))
meta_filt <- meta %>% filter(!rownames(.) %in% samples_to_drop) %>%
  mutate(age_months = droplevels(age_months))
dds_filt <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = meta_filt,
                              design = ~ age_months)

# Filter based on low count genes
smallestGroupSize <- min(table(meta_filt$age_months))
keep <- rowSums(counts(dds_filt) >= 5) >= smallestGroupSize
dds_filt <- dds_filt[keep,]

# Relevel the reference 
dds_filt$age_months <- relevel(dds_filt$age_months, ref = "9")

# Run DESeq
dds_filt <- DESeq(dds_filt)

#### 7. Re-examine the QC metrics ####
# Gene counts
plotCounts(dds_filt, gene = "Il1b", intgroup = "age_months")
plotCounts(dds_filt, gene = "Dclk1", intgroup = "age_months")

# Normalization
ntd <- normTransform(dds_filt) # This gives log2(n+1) 
vsd <- vst(dds_filt, blind = FALSE) # Setting blind to false to incorporate already estimated dispersion

# Compare the removal of dependence of variance on the mean
pdf(here(plots_dir,"Sample_filtered_mean_vs_sd_plots.pdf"))
meanSdPlot(assay(ntd), ylab = "normTransform sd")
meanSdPlot(assay(vsd), ylab = "vst sd")
dev.off()

# classical PCA
pcaData <- plotPCA(vsd, intgroup = "age_months", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_postfilter_plot <- ggplot(pcaData, aes(PC1, PC2, color=age_months)) +
  geom_point(size=4, alpha = 0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggpubr::theme_classic2()+
  ggsci::scale_color_nejm()

ggsave(here(plots_dir,"PCA plot after filtering.png"), PCA_postfilter_plot,
       width = 6)

# Plot PC1, as PC1 explains the majority of variance
PCA_postfilter_PC1_plot <- pcaData %>%
  mutate(age_months = factor(age_months, levels = sort(unique(meta_filt$age_months))))%>%
  ggplot(aes(x = age_months, y = PC1, color=age_months)) +
  geom_jitter(size=3, alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  xlab(paste0("Age (months)")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  ggpubr::theme_classic2()+
  ggsci::scale_color_nejm()

ggsave(here(plots_dir,"PC1 boxplot after filtering.png"), PCA_postfilter_PC1_plot,
       height = 8, width = 6)

# Save data
saveRDS(object = dds_filt, here(analysis_dat_dir, "DESeq_obj.rds"))
saveRDS(object = dds, here(analysis_dat_dir, "DESeq_obj_unfiltered.rds"))
data.table::fwrite(meta, here(analysis_dat_dir,"meta.csv"), row.names = TRUE)
data.table::fwrite(meta_filt, here(analysis_dat_dir,"meta_filtered.csv"), row.names = TRUE)
data.table::fwrite(assay(vsd), here(analysis_dat_dir, "vst_norm_counts.csv"), row.names = TRUE)
