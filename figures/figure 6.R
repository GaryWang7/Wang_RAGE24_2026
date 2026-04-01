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


# ==============================================================================
# Analysis notes
# ==============================================================================
# age_wks is modeled as a numeric covariate in the DESeq2 design.
# This tests for a linear association between expression and age.
#
# In this framework, the log fold change represents the expected change in
# expression for a one-week increase in age.
#
# Compared with a categorical age-group analysis, this is a more stringent test
# for age dependence because it evaluates monotonic change across age rather than
# any difference between groups.
#
# Relevant discussion:
# https://support.bioconductor.org/p/106021/
# https://support.bioconductor.org/p/126713/


library(here)
library(tidyverse)
library(fs)
library(DESeq2)

future::plan("multisession", workers = 14)
options(future.globals.maxSize = 100000 * 1024^2)  # ~100 GB


# ==============================================================================
# Input paths and metadata
# ==============================================================================
data_dir <- here("~/v/RAGE24_LongRead/bulk_long_read", "data", "IsoQuant_ambiguous")
project_dir <- here("garyw", "RAGE24")
plot_dir <- here(project_dir, "publication/figures/figure6")

age.col <- c(
  "16" = "#b8d8ba",
  "30" = "#d9dbbc",
  "56" = "#dbac95",
  "82" = "#555b6e"
)

rat_ages <- readxl::read_xlsx(
  here("~/v/RAGE24_LongRead/bulk_long_read", "meta", "meta_ages.xlsx")
)

meta <- readxl::read_xlsx(
  here("~/v/RAGE24_LongRead/bulk_long_read", "meta", "meta_RAGE24_longread_7_23_2025.xlsx")
) %>%
  left_join(rat_ages)

DESeq2_meta <- meta %>%
  distinct(sample_id, age_wks) %>%
  arrange(sample_id) %>%
  column_to_rownames("sample_id") %>%
  mutate(
    age_wks = as.numeric(age_wks),
    age_months = age_wks / 4
  )

# ==============================================================================
# Gene level quantification 
# ==============================================================================
gene.raw.counts <- read_tsv(here(data_dir, "RAGE24_LongRead", "RAGE24_LongRead.gene_grouped_counts.tsv"))
gene.counts <- gene.raw.counts %>%
  tibble::column_to_rownames("gene_id") 

# Make sure the COLUMNS of count matrix are of the same order as 
# ROWS in metadata (coldata)
gene.counts<- gene.counts[,rownames(DESeq2_meta)]

# Construct DESeq object
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = DESeq2_meta,
                              design = ~ age_wks)

# Pre-filtering genes for visualization and speed for DESeq
# Require genes to have at least 5 counts in at least 6 samples (smallest group size)
# There is no absolute reason why 5, but 5-10 should be reasonable choice.
smallestGroupSize <- min(table(DESeq2_meta$age_wks))
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]

vsd <- vst(dds, blind = FALSE)

# Run DESeq
dds <- DESeq(dds, test="Wald")
res <- results(dds, cooksCutoff = FALSE) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

# Function to plot gene expression
vst_plot <- function(gene){
  df <- assay(vsd)[gene,] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")%>%
    dplyr::rename("counts" = ".") %>%
    left_join(DESeq2_meta |> tibble::rownames_to_column("sample") |> mutate(age_wks = factor(age_wks)))
  
  p <- df %>%
    ggplot(aes(x = age_wks, y = counts, color = age_wks, fill = age_wks))+
    geom_point(alpha = 0.9)+
    geom_boxplot(alpha = 0.7)+
    xlab("Age (weeks)")+
    ylab("VST Normalized counts")+
    ggtitle(paste0("RAGE long read: ", gene))+
    ggpubr::theme_pubr() +
    scale_color_manual(values = age.col) +
    scale_fill_manual(values = age.col)
  mod <- lm(counts ~ age_wks, 
            data= df |> mutate(age_wks = as.numeric(age_wks)))
  # summary(mod)$coefficients["age_months","Pr(>|t|)"]
  return(list("plot" = p,"mod" = mod))
}

p <- vst_plot("Dclk1")$plot
ggsave(here(plot_dir,"Dclk1_total_gene.pdf"), plot = p, width = 3, height = 5)

# ==============================================================================
# GSEA analysis 
# ==============================================================================
# Load msigdb of gene pathways
set.seed(1105)
msigDb <- msigdbr::msigdbr(species = "Rattus norvegicus")

# According to https://stephenturner.github.io/deseq-to-fgsea/ and https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
# "stat" from DESeq or Limma are good ranking metrics
# The default metric is signal-to-noise ratio, so we use "stat" here ,which is log2FC/lfcSE

library(fgsea)
library(BiocParallel)
param <- SnowParam(workers = parallel::detectCores() - 2)
library(pheatmap)
# Heatmap colors
heat.up.col <- colorRampPalette(rev( c(
  "#C96A4A",  # high
  "#E39A6F",
  "#F2C9B1"
)))(200)

heat.down.col <- colorRampPalette(
  c(
    "#789C8A",
    "#8FB6A2",
    "#EEF4F1"
  )
)(200)

# To generate stats for each gene, we have to remove duplicates
gsea_stat <- function(res){
  res <- res %>% 
    tibble::rownames_to_column("gene_id") %>%
    dplyr::distinct(gene_id, stat)%>%
    group_by(gene_id) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%  # one row per gene_name
    ungroup() %>%
    arrange(stat)
}

# Function for running gsea
gsea_wrapper <- function(gsea_res=gsea_stat, gset){
  if(length(gset) == 1){
    if(! gset %in% unique(msigDb$gs_collection)) {
      stop("gene set: ",gset," is not in MsigDb")
    }
    Db <- dplyr::filter(msigDb, gs_collection == gset)
    pathway_list <- split(x = Db$gene_symbol,
                          f = Db$gs_name)
  }
  if(length(gset) == 2){
    if(! gset[[1]] %in% unique(msigDb$gs_collection)) {
      stop("gene set: ",gset[[1]]," is not in MsigDb")
    }
    Db <- dplyr::filter(msigDb, gs_collection == gset[[1]],
                        gs_subcollection %in% gset[[2]])
    pathway_list <- split(x = Db$gene_symbol,
                          f = Db$gs_name)
  }
  rank <- deframe(gsea_res)
  logger::log_info("Performing fgsea...")
  fgseaRes <- fgsea(pathways = pathway_list,
                    stats = rank,
                    minSize = 10,
                    maxSize = 500,
                    nproc = parallel::detectCores() - 2,
                    BPPARAM = param)
  return(list(path_list = pathway_list, gsea = fgseaRes))
}

# Relevel the age_wks for group wise comparison
dds <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = DESeq2_meta |> mutate(age_wks = factor(age_wks)),
                              design = ~ age_wks)
dds <- DESeq(dds, test="Wald")

# Get GSEA results for all ages vs 16 weeks
ages <- c(30,56,82)
GSEA_res_all <- lapply(ages, function(.age){
  logger::log_info("Performing GSEA on comparison: ", .age, " vs 16 weeks")
  .res <- as.data.frame(results(dds, contrast = c("age_wks", .age, 16)))
  .gsea_stat <- gsea_stat(.res)
  .gsea_H <- gsea_wrapper(.gsea_stat, gset = "H")
  .df <- .gsea_H$gsea %>%
    arrange(pathway)
  return(.df)
})
names(GSEA_res_all) <- ages

# Construct heatmap matrix
col.viz <- "NES"
temp_plotdir <- here(plot_dir,"GSEA","heatmap")
dir_create(temp_plotdir)

pathways <- GSEA_res_all$`30`[["pathway"]]
mat <- as.matrix(rbind(GSEA_res_all$`30`[[col.viz]],
      GSEA_res_all$`56`[[col.viz]],
      GSEA_res_all$`82`[[col.viz]]))

colnames(mat) <- str_remove(pathways, pattern = "HALLMARK_")
rownames(mat) <- ages

# Split matrix and filter 
up_keep <- colMeans(mat) > 0 & colMeans(abs(mat)) > 1.4
down_keep <- colMeans(mat) < 0 & colMeans(abs(mat)) > 1.7

mat_up <- mat[,up_keep]
mat_down <- mat[,down_keep]

# Heatmap
pheatmap(t(mat_up), cluster_cols = FALSE, color = heat.up.col, scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         filename = here(temp_plotdir,"Upregulated in aging.pdf"),
         width = 6 ,height  = 5)

pheatmap(t(mat_down), cluster_cols = FALSE, color = heat.down.col, scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_method = "average",
         filename = here(temp_plotdir,"Downregulated in aging.pdf"),
         width = 5.3 ,height  = 5)
         

# ==============================================================================
# Transcript counts and annotation
# ==============================================================================
# Discovered transcript groups from IsoQuant are used for downstream analysis.
tx_counts <- data.table::fread(
  here(data_dir, "RAGE24_LongRead", "RAGE24_LongRead.discovered_transcript_grouped_counts.tsv")
) %>%
  column_to_rownames("gene_id") %>%
  select(rownames(DESeq2_meta)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0L) |> round()))

gtf_model <- rtracklayer::import(
  here(data_dir, "RAGE24_LongRead", "RAGE24_LongRead.transcript_models.gtf")
)
gtf_model <- as.data.frame(gtf_model)

# Retrieve all annotated and novel Dclk1 transcripts.
Dclk1_transcripts <- gtf_model %>%
  filter(gene_id == "Dclk1", type == "transcript") %>%
  distinct(transcript_id)

# Example Dclk1 transcripts identified in the GTF model:
# print(Dclk1_transcripts$transcript_id)
#   XM_017591139.3
#   XM_039103250.2
#   XM_039103254.2
#   XM_039103247.2
#   XM_039103260.2
#   XM_039103258.2
#   XM_039103257.2
#   transcript60602.NC_086020.1.nic
#   transcript60613.NC_086020.1.nic
#   transcript60635.NC_086020.1.nic

# Retrieve total transcript counts for each transcript
rowSums(tx_counts[Dclk1_transcripts$transcript_id,]) %>% as.data.frame()
# XM_017591139.3                   32
# XM_039103250.2                  126
# XM_039103254.2                   33
# XM_039103247.2                   44
# XM_039103260.2                  122
# XM_039103258.2                  124
# XM_039103257.2                  352
# transcript60602.NC_086020.1.nic  16
# transcript60613.NC_086020.1.nic  18
# transcript60635.NC_086020.1.nic 211

# Combine the DCX isoforms which correspond to the same protein sequence.
Dclk1_DCX <- tx_counts[c("XM_039103257.2","XM_039103258.2","XM_039103260.2"),] %>%
  colSums()
tx_counts <- rbind(tx_counts, "Dclk1_DCX" = Dclk1_DCX)

# Combine the all L isoforms. They are not exactly the same but share the same DCX domains and kinase domains
Dclk1_L_all <- tx_counts[c("XM_017591139.3","XM_039103250.2","XM_039103254.2","XM_039103247.2"),] %>%
  colSums()
tx_counts <- rbind(tx_counts, "Dclk1_L_all" = Dclk1_L_all)

# Add an annotation for Dclk1_DCX, Dclk1_L_all
gtf_model <- gtf_model%>%
  tibble::add_row(type = "transcript", gene_id = "Dclk1", gene = "Dclk1", transcript_id = "Dclk1_DCX") %>%
  tibble::add_row(type = "transcript", gene_id = "Dclk1", gene = "Dclk1", transcript_id = "Dclk1_L_all")

# Transcript level quantification ----------------------------------------------------
# DESeq object
tx_dds <- DESeq2::DESeqDataSetFromMatrix(countData = tx_counts,
                                 colData = DESeq2_meta,
                                 #design = ~ age_months,
                                 design = ~ age_wks)

# Filter DESeq object
tx_keep <- rowSums(counts(tx_dds) >= 3) >= 6
tx_dds <- tx_dds[tx_keep,]

# VST normalized gene counts
vsd_tx <- vst(tx_dds, blind = FALSE)

# LRT results 
# Turn off cook's outlier detection to retrieve all p-value
tx_dds.LRT <- DESeq(tx_dds, test="LRT", reduced=~1)
res_tx.LRT <- DESeq2::results(tx_dds.LRT, cooksCutoff=FALSE)

# Wald test results. See notes above for using age_wks as a numeric covariate.
# Turn off cook's outlier detection to retrieve all p-value
tx_dds.wald <- DESeq(tx_dds, test = "Wald")
res_tx.wald <- results(tx_dds.wald, cooksCutoff=FALSE)

# Plot graph with statistics (LRT p-value: '~age_wks' vs '~1'; Wald p-value: test for non-zero slope of age_wks)
plot_vst_tx_pval <- function(gene, vsd = vsd_tx, method = c("LRT")){
  group.colors <- c("16" = "#b8d8ba", "30" = "#d9dbbc", "56" = "#dbac95", "82" = "#555b6e")
  .meta <- meta %>%
    distinct(sample_id, age_wks)
  if(length(gene)!=1){
    stop("Only allow 1 gene input")
  }
  avail_tx <- rownames(vsd_tx)
  tx_gene <- filter(gtf_model, gene_id == !!gene, type == "transcript") %>%
    pull(transcript_id) %>%
    unique()
  tx <- intersect(tx_gene, avail_tx)
  if(length(tx) == 0){
    stop("No abundantly transcribed transcript found!")
  }
  if(length(tx) == 1){
    vst_exp <- assay(vsd)[tx,] %>%
      t() %>%
      as.data.frame()
    rownames(vst_exp) <- tx
  }else{
    vst_exp <- assay(vsd)[tx,] %>%
      as.data.frame()
  }
  res <- switch(method,
                "LRT"  = res_tx.LRT,
                "Wald" = res_tx.wald,
                stop("Unknown method")
  )
  res <- res %>% as.data.frame() %>%
    rownames_to_column("transcript_id") 
  vst <-  vst_exp %>%
    tibble::rownames_to_column("transcript_id") %>%
    filter(transcript_id %in% tx) %>%
    tidyr::pivot_longer(cols = -transcript_id,
                        names_to = "sample_id",
                        values_to = "counts") %>%
    left_join(.meta) %>%
    left_join(res) %>%
    mutate(padj2 = p.adjust(pvalue, method = "fdr")) %>%
    group_by(transcript_id) %>%
    mutate(y_pos = min(counts) - 0.05)
  vst %>%
    mutate(age_wks = factor(age_wks)) %>%
    ggplot(aes(x = age_wks, y = counts, color = age_wks, fill = age_wks))+
    geom_boxplot(alpha = 0.6)+
    geom_point(alpha = 0.9)+
    geom_text(
      # Only report pvalues because we are mostly interested in the DCX isoform and do not want to involve other p-values
      aes(y = y_pos, label = paste0("padj = ", signif(padj2, 3))), x = 2.5, color = "#757471",
      fontface = "italic", family = "sans", inherit.aes = FALSE, alpha = 1, size = 4
    )+
    facet_wrap(vars(transcript_id), scales = "free_y", nrow = 1)+
    scale_fill_manual(values = group.colors)+
    scale_color_manual(values = group.colors) + 
    xlab("Age (months)")+
    ylab("VST Normalized counts of transcript")+
    ggtitle(paste0("Normalized expressions in RAGE bulk long read for ", gene))+
    ggpubr::theme_pubr()
}
p <- plot_vst_tx_pval("Dclk1", method = "Wald")
ggsave(here(plot_dir, "DESeq long read isoform expression.pdf"), plot = p, width = 16, height = 4.5)

# ==============================================================================
# DEX-seq
# ==============================================================================
#BiocManager::install("DEXSeq")
library(DEXSeq)

# Read the cleaned up data.
# There were some redundant entries in the exon count matrix. For the same exon, there may be multiple entries, which should be removed.
# count.ex <- readr::read_tsv(here(data_dir,"RAGE24_LongRead", "RAGE24_LongRead.exon_grouped_counts.linear.tsv"))
# 
# # Clean up some redundant entries in the exon count matrix. For the same exon, there may be multiple entries
# count.ex.clean <- count.ex %>%
#   dplyr::rename("seqid" = `#chr`, gene_id = gene_ids) %>%
#   group_by(seqid, start, end, strand, flags, gene_id, group_id) %>%
#   summarise(
#     include_counts = sum(include_counts, na.rm = TRUE),
#     exclude_counts = sum(exclude_counts, na.rm = TRUE),
#     .groups = "drop"
#   )
# readr::write_tsv(count.ex.clean,
#                  here(data_dir,"RAGE24_LongRead", "RAGE24_LongRead.exon_grouped_counts.linear.cleaned.tsv") )

count.ex.clean <- read_tsv(here(data_dir,"RAGE24_LongRead", "RAGE24_LongRead.exon_grouped_counts.linear.cleaned.tsv"))

# wider form of count matrix
count.ex.wide <- count.ex.clean %>%
  dplyr::select(-exclude_counts) %>%
  # DEXSeq dislikes "-" in sample names
  dplyr::mutate(group_id = str_replace(group_id, pattern = "-", replacement = ".")) %>%
  pivot_wider(
    names_from = group_id,
    values_from = include_counts
  ) %>%
  mutate(
    dplyr::across(
      dplyr::ends_with("KYC.LN"),
      ~ tidyr::replace_na(., 0)
    )
  )

# count matrix only
count.mtx <- count.ex.wide %>%
  dplyr::select(dplyr::ends_with("KYC.LN"))

# Filter the count matrix based on number of samples expressing the exons
row.keep <- rowSums(count.mtx > 3, na.rm = TRUE) >= 6 # At least 6 samples should have more than 3 counts
count.ex.df <- count.ex.wide[row.keep,]

# exon meta 
exon.meta <- dplyr::select(count.ex.df,
                           seqid, start, end, strand, flags, gene_id) %>%
  mutate(exon_id = paste(seqid, start, end, strand, sep = "_"),
         # Use a simpler exon number to trace back
         exon_no = paste("E", dplyr::row_number(), sep = "_"),
         strand = case_when(strand == "+-" ~ "*",.default = strand))

gr_exons <- makeGRangesFromDataFrame(exon.meta,
                                     seqnames.field = "seqid",
                                     start.field = "start",
                                     end.field = "end",
                                     strand.field = "strand",
                                     keep.extra.columns = TRUE)
# metadata
meta_dxd <- mutate(meta,
                   group = case_when(age_wks %in% c(16,30) ~ "Young",
                                     age_wks %in% c(56,82) ~ "Old")) %>%
  distinct(sample_id, age_wks, group) %>%
  # Edit sample id according to above
  mutate(sample_id = str_replace(sample_id, "-", ".")) %>%
  mutate(
    age_wks = factor(age_wks, levels = c(16,30,56,82)),
    group = factor(group, levels = c("Old","Young"))
  )

# Create DEXSeq obj
dxd <- DEXSeqDataSet(
  countData = count.ex.df[,meta_dxd$sample_id],
  sampleData = meta_dxd,
  design = ~ sample_id + exon + age_wks:exon,
  featureID = exon.meta$exon_no,
  groupID = exon.meta$gene_id,
  featureRanges = gr_exons
)

# Normalization
dxd <- estimateSizeFactors(dxd)

# Estimate dispersion (parallelization has not been avaiable for glmGamPoi)
dxd <- estimateDispersions(dxd,  BPPARAM = MulticoreParam(14), quiet = FALSE)

# Test for differential exon usage
dxd <- testForDEU(dxd, BPPARAM = MulticoreParam(14))
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "age_wks", BPPARAM = MulticoreParam(14))

# result
dxr <- DEXSeqResults(dxd)
 
pdf(file = here(plot_dir,"Dclk1_exon_usage_all_genes.pdf"),width = 10, height = 8)
plotDEXSeq(dxr,"Dclk1", fitExpToVar = "age_wks", legend=TRUE, splicing = TRUE, expression = TRUE,
           color = c("16" = "#b8d8ba", "30" = "#d9dbbc", "56" = "#dbac95", "82" = "#555b6e"))
dev.off()