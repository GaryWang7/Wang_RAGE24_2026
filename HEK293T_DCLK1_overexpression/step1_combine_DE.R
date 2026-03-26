# This scripts combines and evaluates the data from DCLK1 overexpression and DCLK1 overexpresison+DCLK1IN1
library(ComplexHeatmap)
library(AnnotationDbi)
library(DESeq2)
library(dendextend)
library(tidyverse)
library(readxl)
library(here)
library(RColorBrewer)
library("factoextra")
library(ggrepel)
library(enrichplot)
library(fs)
library(pheatmap)
library(rrcov)
library(glmpca)
library(rstatix)

#### 1. Define directories ####
data_dir1 <- "D:/garyw/RAGE24/DCLK1_HEK293T"
data_dir2 <- "D:/garyw/RAGE24/DCLK1_HEK293T_oe_DCLK1IN1"
ref_dir <- here("D:/garyw/RAGE24/reference/GRCh38p14")

# Analysis results are deposited in the newer folder
plots_dir <- here(data_dir2, "plots","step2")
analysis_dat_dir <- here(data_dir2,"analysis_data")
dir_create(plots_dir)
dir_create(analysis_dat_dir)

#### 2. Load count data and metadata ####
# Renamed counts
raw_counts1 <- data.table::fread(here(data_dir1,"counts", "counts.csv")) %>%
  tibble::column_to_rownames("V1")
raw_counts2 <- data.table::fread(here(data_dir2,"counts", "counts.csv")) %>%
  tibble::column_to_rownames("V1")

# Combine counts
raw_counts <- cbind(raw_counts1, raw_counts2)

# Metadata 
meta1 <- data.table::fread(here(data_dir1, "metadata_edited.csv")) %>%
  dplyr::rename("sample" = "Customer ID",
                "AdmeraID" = "Sample") %>%
  mutate(Treatment = factor(Treatment),
         Lane = factor(Lane),
         Time = "May25",
         # In this data, lane is having some batch effect
         Batch = paste0(Time,".", Lane)) %>%
  tibble::column_to_rownames("sample") %>%
  arrange(Treatment)
meta2 <- data.table::fread(here(data_dir2, "metadata_edited.csv")) %>%
  dplyr::rename("sample" = "Customer ID",
                "AdmeraID" = "Sample",
                `% >= Q30 bases` = `% >= Q30bases`) %>%
  mutate(Treatment = factor(Treatment),
         Lane = factor(Lane),
         Time = "Aug25",
         # In this data, the lane is having acceptable batch effect
         Batch = paste0(Time,".", Lane)) %>%
  tibble::column_to_rownames("sample") %>%
  arrange(Treatment)

# Combine meta 
meta <- rbind(meta1,meta2)

## 2.1 Genes mapping both to chrX and chrY
# We have to filter out chrY reads in this data, as HEK293 is originated from a female embryo.
# Rsubread default enables multimapping (each mapped gene will get 1 count). 
# For genes that are on both chrX and chrY (pseudoautosomal regions, PARs), the chrX and chrY part are likely to both get 1 count.
# Therefore, removing the chrY counterparts of PAR genes are safe in this regard. 
# chrY-specific genes may be due to contamination.
# Here we perform a QC on how many reads are on chrY.
gtf <- rtracklayer::readGFF(here(ref_dir, "gencode.v47.primary_assembly.annotation_modified_GFP.gtf"))
chrX.genes <- gtf %>%
  dplyr::filter(seqid == "chrX")
chrY.genes <- gtf %>%
  dplyr::filter(seqid == "chrY")
PAR.X.genes <- chrX.genes %>% # PAR genes on chrX
  dplyr::filter(gene_name %in% chrY.genes$gene_name ) %>%
  distinct(gene_name, gene_id)
PAR.Y.genes <- chrY.genes %>%
  dplyr::filter(gene_name %in% chrX.genes$gene_name) %>%
  distinct(gene_name, gene_id)
chrY.only.genes <- chrY.genes %>%
  dplyr::filter(! gene_name %in% chrX.genes$gene_name)%>%
  distinct(gene_name, gene_id)

# Map the categories of gene_ids
cat_map <- cat_map <- rownames(raw_counts) %>%
  tibble(gene_id = .) %>%
  mutate(category = case_when(
    gene_id %in% PAR.X.genes$gene_id      ~ "PAR.X",
    gene_id %in% PAR.Y.genes$gene_id      ~ "PAR.Y",
    gene_id %in% chrY.only.genes$gene_id  ~ "chrY.only",
    .default ="Other"
  ))

# Sum counts per category × sample
cat_pct <- as.data.frame(raw_counts) %>%
  rownames_to_column("gene_id") %>%
  left_join(cat_map, by = "gene_id") %>%
  pivot_longer(-c(gene_id, category),
               names_to = "sample", values_to = "count") %>%
  group_by(sample, category) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(perc = count / sum(count)*100) %>%
  ungroup()

# Plot count percentage and counts
p <- cat_pct %>%
  dplyr::filter(category %in% c("PAR.X","PAR.Y")) %>%
  ggplot(aes(x = sample,y = perc, fill = category)) +
  geom_col() +
  geom_text(aes(label = round(perc,digits = 3)),
            position = position_stack(vjust = 0.5), size = 1.3, color = "black") +
  labs(title = "PAR genes count percentage", y = "% counts", x = NULL) +
  ggpubr::theme_pubclean() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(here(plots_dir,"PAR genes count percentage.png"), p, width = 7, height = 5)

p1 <- cat_pct %>%
  dplyr::filter(category %in% c("PAR.X","PAR.Y")) %>%
  ggplot(aes(x = sample,y = count, fill = category)) +
  geom_col() +
  labs(title = "PAR genes count", y = "counts", x = NULL) +
  ggpubr::theme_pubclean() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(here(plots_dir,"PAR genes counts.png"), p1, width = 7, height = 5)

## 2.2 Remove chrY genes and align the samples
counts <- raw_counts[!rownames(raw_counts) %in% chrY.genes$gene_id, ]

# Make sure the COLUMNS of count matrix are of the same order as 
# ROWS in metadata (coldata)
counts <- counts[,rownames(meta)]

#### 3. Construct DESeq object ####
# We construct this object to evaluate batch effect here. 
# For DEG analysis (next step), we will perform DESeq2 on two batches separately.
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              # DESeq2 is unable to model variable which is a linear combination of other variables
                              design = ~ Treatment)
# Pre-filtering genes for visualization and speed for DESeq
# Require genes to have at least 10 counts in at least 3 samples (smallest group size)
# There is no absolute reason why 10, but 0-10 should be reasonable choice.
smallestGroupSize <- min(table(meta$Treatment))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Run DESeq2 
dds <- DESeq(dds)

#### 4. QC the results ####
resultsNames(dds)

# Plot normalized gene counts
# ENSG00000133083 is DCLK1
plotCounts(dds, gene = "ENSG00000133083", intgroup = "Treatment")
plotCounts(dds, gene = "GFP", intgroup = "Treatment")

# For PCA, we need to use rlog or vst transformation to avoid 
# the influence of strongly expressed genes
# See https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# and https://hbctraining.github.io/Intro-to-DGE/lessons/03_DGE_QC_analysis.html
ntd <- normTransform(dds) # This gives log2(n+1)
vsd <- vst(dds, blind = FALSE)

# classical PCA
pcaData <- plotPCA(vsd, 
                   intgroup = c("Treatment","Lane"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_prefilter_plot <- ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size=4, alpha = 0.8, aes(color = Treatment, shape = Lane)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggpubr::theme_classic2()+
  viridis::scale_color_viridis(discrete = TRUE, option = "turbo")

ggsave(here(plots_dir,"PCA plot before batch correction.png"), PCA_prefilter_plot,
       width = 6)
pcaData_summary <- pcaData %>%
  group_by(Treatment) %>%
  mutate(PC1_mean = mean(PC1),
         PC1_sd = sd(PC1))

# Plot PC1, as PC1 explains the majority of variance
PCA_prefilter_PC1_plot <- pcaData %>%
  ggplot(aes(x = Treatment, y = PC1, color=Treatment, shape = Lane)) +
  geom_jitter(size=3, alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.4, outliers = FALSE)+
  xlab(paste0("Treatment")) +
  ylab(paste0("PC1: ",percentVar[1],"% variance")) +
  ggpubr::theme_classic2()+
  viridis::scale_color_viridis(discrete = TRUE, option = "turbo")

ggsave(here(plots_dir,"PC1 boxplot before batch correction.png"), PCA_prefilter_PC1_plot,
       height = 8, width = 12)

#### 5. save data ####
data.table::fwrite(raw_counts, here(analysis_dat_dir, "combined_raw_counts.csv"),
                   row.names = TRUE)
data.table::fwrite(counts, here(analysis_dat_dir, "combined_filtered_counts.csv"),
                   row.names = TRUE)
data.table::fwrite(meta %>% tibble::rownames_to_column("sample"), here(analysis_dat_dir,"combined_meta.csv"),
                   row.names = FALSE)
data.table::fwrite(assay(vsd), here(analysis_dat_dir, "combined_vst.csv"),
                   row.names = TRUE)
saveRDS(dds, file = here(analysis_dat_dir,"combined_DESeq_obj.rds"))
