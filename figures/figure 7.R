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
# Preprocess data 
# ==============================================================================

library(DESeq2)
library(here)
library(tidyverse)
library(fs)
library(msigdbr)
library(clusterProfiler)

#### Define directories ####
project_dir <- here("garyw", "RAGE24")
plot_dir <- here(project_dir,"publication/figures/figure7")
data_dir <- here(project_dir, "DCLK1_HEK293T_oe_DCLK1IN1/analysis_data")

#### Load data ####
# The raw counts were filtered to remove PAR Y genes, as HEK293T were from a female donor
raw_counts <- read_csv(here(data_dir,"combined_filtered_counts.csv")) %>%
  dplyr::rename( "GeneID" = "...1") 
meta <- read_csv(here(data_dir, "combined_meta.csv"))

# Gene names
gtf <- data.table::fread(here(project_dir,"DCLK1_HEK293T_oe_DCLK1IN1","counts", "genes_for_counting.csv"))
ID2gene <- dplyr::select(gtf, GeneID, gene_name)

# We create two DESeq objects for the two experiments done at different time to avoid confounding batch effects
meta1 <- meta %>%
  dplyr::filter(
    Time == "May25"
  ) %>%
  tibble::column_to_rownames("sample")
meta2 <- meta %>%
  dplyr::filter(Time == "Aug25") %>%
  tibble::column_to_rownames("sample")

# Make sure counts follow the same order as meta
counts <- tibble::column_to_rownames(raw_counts,"GeneID")
counts1 <- counts[,rownames(meta1)]
counts2 <- counts[,rownames(meta2)]

# Create DESeq object
dds1 <- DESeqDataSetFromMatrix(countData = counts1, colData = meta1,
                               # This experiment was having some batch effect (number of genes mapped) in lanes
                               design = ~ Lane + Treatment)
dds2 <- DESeqDataSetFromMatrix(countData = counts2, colData = meta2,
                               # Lane is not having batch effect here
                               design = ~ Treatment)

# Filter based on counts
keep1 <- rowSums(counts(dds1) >= 10) >= 6 # 6 samples in each group 
dds1 <- dds1[keep1,]
keep2 <- rowSums(counts(dds2) >= 10) >= 3 # 3 samples in each group
dds2 <- dds2[keep2,]

# Run DESeq2
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)

# normalized counts
vsd1 <- vst(dds1, blind = FALSE)
vsd2 <- vst(dds2, blind = FALSE)

#### Plot of individual genes ####
plot_vst <- function(.vsd, .meta, genes){
  gene_ids_df <- ID2gene %>%
    dplyr::filter(GeneID %in% rownames(.vsd)) %>%
    dplyr::filter(gene_name %in% genes)
  genes.omit <- genes[!genes %in% gene_ids_df$gene_name]
  if(length(genes.omit) != 0){
    print(paste0("Following genes are not in the query: ", 
                 paste0(genes.omit, collapse = " ")))
  }
  if(nrow(gene_ids_df) == 1){
    genes_vst_exp <- assay(.vsd)[gene_ids_df$GeneID,]
    genes_vst <- as.data.frame(t(genes_vst_exp))
  }else{
    genes_vst <- assay(.vsd)[gene_ids_df$GeneID,] %>%
      as.data.frame()
  }
  .meta <- tibble::rownames_to_column(.meta, "sample")
  genes_vst <- genes_vst %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(cols = -gene_id,
                        names_to = "sample",
                        values_to = "counts") %>%
    left_join(.meta, by = join_by(sample == sample)) %>%
    left_join(gene_ids_df, by = join_by(gene_id == GeneID))
  p <- genes_vst %>%
    ggplot(aes(x = Treatment, y = counts, 
               color = Treatment, fill = Treatment))+
    geom_point(alpha = 0.8, show.legend = FALSE)+
    geom_boxplot(alpha = 0.7, show.legend = FALSE)+
    viridis::scale_color_viridis(discrete = TRUE, option = "turbo")+
    viridis::scale_fill_viridis(discrete = TRUE, option = "turbo") +
    ggpubr::theme_pubr()
  # ggpubr::stat_compare_means(method = "t.test", 
  #                            label.y.npc = "middle",
  #                            alpha = 0.6)+
  if(length(genes)>1){
    p <- p + facet_wrap(vars(gene_name), scales = "free_y")
  }else{
    p <- p + labs(title = genes) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  p +
    xlab("Conditions")+
    ylab("VST Normalized counts")+
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))
}

plot_vst(vsd1,meta1,"JUNB")

# ==============================================================================
# GSEA functions 
# ==============================================================================
library(fgsea)
library(BiocParallel)
param <- SnowParam(workers = parallel::detectCores() - 2)

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

set.seed(1105)
msigDb <- msigdbr::msigdbr(species = "human")
gsea_wrapper <- function(gsea_res=gsea_stat, gset){
  if(length(gset) == 1){
    if(! gset %in% unique(msigDb$gs_collection)) {
      stop("gene set: ",gset," is not in MsigDb")
    }
    Db <- dplyr::filter(msigDb, gs_collection == gset)
    pathway_list <- split(x = Db$ensembl_gene,
                          f = Db$gs_name)
  }
  if(length(gset) == 2){
    if(! gset[[1]] %in% unique(msigDb$gs_collection)) {
      stop("gene set: ",gset[[1]]," is not in MsigDb")
    }
    Db <- dplyr::filter(msigDb, gs_collection == gset[[1]],
                        gs_subcollection %in% gset[[2]])
    pathway_list <- split(x = Db$ensembl_gene,
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


# ==============================================================================
# GSEA for DCLK1-L and DCLK1-S overexpression
# ==============================================================================
#### GSEA using Hallmark gene sets ####
pathway_up_pal <- list(
  "L" =  c(
    "#B65A3C",  # muted brick (high)
    "#D28763",
    "#EBC1A6",
    "#FBF0E8"   # very low
  ),
  "S" = c(
    "#7A4A73",  # deep plum (high)
    "#9E6A97",
    "#C8A6C3",
    "#F1E6F0"   # very low
  )
)

pathway_down_pal <- list(
  "L" =   c(
    "#7A8F6A",  # olive green (high)
    "#9AA87A",
    "#D2DDB7",
    "#F0F4E8"   # very low
  ),
  "S" = pathway_S_down_pal <- c(
    "#4F6F8F",  # steel blue (high)
    "#6B8FA3",
    "#B7CFDE",
    "#EEF4FA"  # very low
  )
)

# GSEA stats for L and S
GSEA_stats <- lapply(c("DCLK1-L","DCLK1-S"), function(treat){
  res <- as.data.frame(results(dds1, contrast = c("Treatment", treat, "Backbone")))
  gsea_stat <- gsea_stat(res)
})
names(GSEA_stats) <- c("L","S")

# Hallmark results
gsea_H <- lapply(c("L","S"), function(treat){
  gsea_wrapper(GSEA_stats[[treat]], gset = "H")
})
names(gsea_H) <- c("L","S")

# Plot pathway enrichment dotplots (upregulated and downregulated for both L and S)
lapply(c("L","S"), function(isoform){
  temp_plotdir <- here(plot_dir,"GSEA","Hallmark",isoform)
  dir_create(temp_plotdir)
  res <- gsea_H[[isoform]]$gsea %>%
    filter(padj < 0.1) %>%
    dplyr::arrange(desc(NES), padj) %>%
    mutate(change = if_else(NES > 0, true = "UP", false = "DOWN"))
  res$pathway <- res$pathway%>%
    tolower() %>%
    str_remove_all("hallmark_")
  lapply(unique(res$change), function(reg){
    if(reg == "UP"){pal <- pathway_up_pal[[isoform]]}else{pal<-pathway_down_pal[[isoform]]}
    p <- res %>%
      filter(change == reg) %>%
      slice_max(order_by = abs(NES), n = 8) %>%
      ggplot(aes(x = reorder(x = pathway, X = abs(NES)), y = NES)) + 
      geom_point(aes(color = padj, size = abs(NES))) +
      coord_flip()+
      labs(x = "", y="Normalized Enrichment Score", size = "|NES|")+
      scale_color_gradientn(colours = pal,
                            guide = guide_colorbar(reverse = TRUE))+
      scale_size_continuous(range = c(2, 6)) +
      ggpubr::theme_classic2() +
      theme( 
        legend.direction = "vertical",
        legend.position = "right"
        #legend.box = "vertical",
        #legend.title.align = 0 ,
        #legend.key.width = unit(0.3, "in")
      ) 
    ggsave(here(temp_plotdir,paste0(isoform,"-",reg,".pdf")), p, width = 6, height = 3.5)
  })
})

# ==============================================================================
# Venn diagram between L, S, L_DIN 
# ==============================================================================
#devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)

#### Generate list of up- and down-regulated genes
# Generate the list of upregulated and downregulated genes

# Function to get significant genes for a given treatment and p-value threshold
get_sig_genes <- function(treat, pval_thre = 0.01){
  if(str_detect(treat, "DIN")){
    dds <- dds2
    ct <- c("Treatment", treat, "Backbone_DIN")
  }else{
    dds <- dds1
    ct <- c("Treatment", treat, "Backbone")
  }
  genes <- as.data.frame(results(dds, contrast = ct)) %>%
    tibble::rownames_to_column("ensembl_gene") %>%
    filter(pvalue < pval_thre) %>%
    mutate(regulation = ifelse(log2FoldChange > 0, yes = "up", no = "down"),
           Treatment = treat) %>%
    dplyr::select(ensembl_gene, regulation, Treatment)
}

sig.df <- rbind(
  get_sig_genes("DCLK1-L", pval_thre = 0.05),
  get_sig_genes("DCLK1-S", pval_thre = 0.05),
  get_sig_genes("DCLK1-L_DIN", pval_thre = 0.05)
)

## For DCX-dependent genes
# Up
DCX.genes.up <- sig.df %>%
  mutate(
    group = paste0(Treatment, ":",regulation)
  ) %>%
  filter(group %in% c("DCLK1-L:up","DCLK1-S:up","DCLK1-L_DIN:up")) %>%
  split(x = .$ensembl_gene,
        f = .$group)

p <- ggVennDiagram(DCX.genes.up, label = "count", label_alpha = 0)
ggsave(here(plot_dir, "Venn_DCX_dependent up genes.pdf"), width = 5, height = 5, plot = p)

# Upregulated genes that are shared between DCLK1-L and DCLK1-L_DIN, but not in DCLK1-S
genes.up <- setdiff(intersect(DCX.genes.up$`DCLK1-L_DIN:up`,DCX.genes.up$`DCLK1-L:up`), DCX.genes.up$`DCLK1-S:up`)

# Down
DCX.genes.down <- sig.df %>%
  mutate(
    group = paste0(Treatment, ":",regulation)
  ) %>%
  filter(group %in% c("DCLK1-L:down","DCLK1-S:down","DCLK1-L_DIN:down")) %>%
  split(x = .$ensembl_gene,
        f = .$group)

p <- ggVennDiagram(DCX.genes.down, label = "count", label_alpha = 0)
ggsave(here(plot_dir, "Venn_DCX_dependent down genes.pdf"), width = 5, height = 5, plot = p)

# Downregulated genes that are shared between DCLK1-L and DCLK1-L_DIN, but not in DCLK1-S
genes.down <- setdiff(intersect(DCX.genes.down$`DCLK1-L_DIN:down`,DCX.genes.down$`DCLK1-L:down`), DCX.genes.down$`DCLK1-S:down`)

## Hallmark analysis of DCX-dependent genes 
ora_H.up <- ora_wrapper(genes.up, gset = "H")
ora_H.down <- ora_wrapper(genes.down, gset = "H")

# plot Hallmark results
temp_plotdir <- here(plot_dir, "ORA", "Hallmark")
dir_create(temp_plotdir)
lapply(c("up","down"), function(resname){
  if(resname == "up"){
    res <- ora_H.up
    pathway_pal <- pathway_up_pal$L
  }else{
    res <- ora_H.down
    pathway_pal <- pathway_down_pal$L
  }
  res@result$Description <- res@result$Description %>%
    str_remove("HALLMARK_") %>%
    tolower()
  res@result <- filter(res@result, Count >= 2)
  p <- enrichplot::dotplot(res, showCategory = 8)+
    scale_fill_gradientn(colours = pathway_pal,
                         guide=guide_colorbar(reverse=TRUE))+
    scale_size_continuous(range = c(2, 6))+
    ggpubr::theme_classic2()+
    theme(legend.direction = "horizontal",  
          legend.position = "bottom",
          legend.box = "vertical",
          legend.title.align = 0
    ) 
  ggsave(here(temp_plotdir, paste0(resname," DCX dependent genes.pdf")),p, 
         width = 3.5, height = 6.5)
})

# Save results in csv
genes.up.df <- data.frame("ensembl_gene" = genes.up) %>%
  left_join(ID2gene, by = join_by(ensembl_gene==GeneID))
genes.down.df <- data.frame("ensembl_gene" = genes.down) %>%
  left_join(ID2gene, by = join_by(ensembl_gene==GeneID))
write_csv(genes.up.df, here(plot_dir,"DCLK1-L and DCLK1-L_DIN both upregulated genes_191.csv"))
write_csv(genes.down.df, here(plot_dir,"DCLK1-L and DCLK1-L_DIN both downregulated genes_77.csv"))