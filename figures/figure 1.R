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

# install.packages("ggpubr")
# install.packages("svglite")
# install.packages("msigdbr")
# install.packages("logger")
library(svglite)
library(tidyverse)
library(here)
library(future)
library(fs)
library(Seurat)
library(scCustomize)
library(clusterProfiler)
library(fgsea)

#### Directories ####
proj.dir <- here("garyw/RAGE24")
plot.all.dir <- here(proj.dir,"publication/figures")
plot.dir <- here(plot.all.dir,"figure1")
dir_create(plot.dir)

#### Load data ####
srat <- readRDS(here(proj.dir,"gex_aggr_prep_combine","step2_anno.rds"))

#### color schemes ####
# Highlight PT cells only--color scheme 1
celltype_cols1 <- c(
  # Proximal tubule lineage
  "PT"           = "#8FAE9E",  # unchanged: stable PT baseline
  "PT-injured"  = "#C96A4A",  # unchanged
  
  # Thin limb
  "TL1"  = "#E39A6F",  # unchanged
  "TL2"        = "#CFDED7",  # slightly lighter + cooler than PT
  
  # Loop of Henle / distal nephron
  "C-TAL"        = "#6B8FA3",
  "M-TAL"        = "#4F6F8F",
  "DCT"          = "#7A8F6A",
  "CNT"          = "#9AA87A",
  "CNT-PC"       = "#B6A58B",
  
  # Collecting duct
  "PC"           = "#C7B299",
  "IC-A"         = "#8C6D8E",
  "IC-B"         = "#A48DB5",
  
  # Glomerular / vascular
  "POD"          = "#5E5A80",
  "PEC"          = "#7C7A9A",
  "EC"           = "#4C7280",
  "VSMC/P"       = "#3F5C5A",
  
  # Stromal / immune
  "FIB"          = "#8B7D6B",
  "IMM"          = "#7A4A73"
)
# Highlight PT cells only--color scheme 3
celltype_cols3 <- c(
  # Proximal tubule lineage (highlighted)
  "PT-S1"        = "#8FB6A2",  # light sage
  "PT-S2"        = "#789C8A",  # core PT
  "PT-S3"        = "#5F8573",  # deeper PT
  "PT-injured"   = "#C96A4A",  # highlighted injury state
  
  # Thin limb (unhighlighted)
  "TL1"          = "#AFAFAF",
  "TL2"          = "#A9A9A9",
  
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
pt_healthy_pal <- c(
  "low"  = "#E8F1EC",  # very light sage
  "high" = "#789C8A"   # core PT green
)
pt_injury_pal <- c(
  "low"  = "#F6E6DC",  # very light warm beige
  "high" = "#C96A4A"   # PT-injured rust
)

pathway_up_pal <- c(
  "#C96A4A",  # high
  "#E39A6F",
  "#F2C9B1",
  "#FBF1EA"  # very low
)
pathway_down_pal <- c(
  "#789C8A",
  "#8FB6A2",
  "#CFE3DA",
  "#EEF4F1"
)
#### 1. UMAP ####
p <-DimPlot_scCustom(srat, reduction = "umap",group.by = "celltype_refined1_updated",
                       colors_use = celltype_cols3, figure_plot = TRUE, 
                       size = 0.01, label = FALSE)
ggsave(here(plot.dir,"UMAP_PT_highlight_PT_segments.png"), p, dpi = 2000) # Save as PNG with high dpi. 

#### 2. iPT characterization ####
# Injury markers including Havcr, Nfkb1, Creb, Adamts1, Dcdc2 are from papers below:
# https://www.nature.com/articles/s41588-025-02285-0 Analysis of individual patient pathway coordination in a cross-species single-cell kidney atlas
# https://www.nature.com/articles/s41467-025-59997-4 Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype

Idents(srat) <- "celltype"
injury.marker <- c("Adamts1", "Creb5","Dcdc2","Havcr1","Nfkb1","Dclk1") # Add Dclk1
healthy.marker <- c("Cubn","Hnf4a","Lrp2","Pax8","Mme")

p_injury_dot <- DotPlot_scCustom(srat, features = injury.marker,
                                 group.by = "celltype", 
                                 idents = c("PT","PT-injured"),
                                 colors_use = pt_injury_pal,
                                 x_lab_rotate = TRUE) +
  scale_size(range = c(1,10)) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical") +
  theme(legend.title.align = 0)
p_healthy_dot <- DotPlot_scCustom(srat, features = healthy.marker,
                                   group.by = "celltype", 
                                   idents = c("PT","PT-injured"),
                                   colors_use = pt_healthy_pal,
                                   x_lab_rotate = TRUE) +
  scale_size(range = c(4, 12)) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical") +
  theme(legend.title.align = 0)

ggsave(here(plot.dir, "injury_markers_PT_legend.png"), p_injury_dot, dpi = 1000,
       width = 5, height = 3)
ggsave(here(plot.dir, "healthy_markers_PT_legend.png"), p_healthy_dot, dpi = 1000,
       width = 5, height = 3)

#### 3. GSEA for iPT vs healthy PT ####
Idents(srat) <- "celltype"
iPT_markers <- FindMarkers(srat, ident.1 = "PT-injured", group.by = "celltype",
                           min.pct = 0.1)
write_csv(iPT_markers, here(plot.dir,"DEG_iPT.csv"))

# For stats, use avg_log2FC as statistic, according to https://rnabio.org/module-08-scrna/0008/05/01/Gene_set_enrichment/
gsea_stat <- function(res){
  res <- res %>% 
    tibble::rownames_to_column("gene_id") %>%
    dplyr::distinct(gene_id,avg_log2FC)%>%
    group_by(gene_id) %>%
    slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%  # one row per gene_name
    ungroup() %>%
    arrange(avg_log2FC)
}
gsea_stat <- gsea_stat(iPT_markers)

# Load msigdb of gene pathways
msigDb <- msigdbr::msigdbr(species = "Rattus norvegicus")

# function for performing GSEA analysis. The function is written for multiple databases. Here we will use the Hallmark gene sets.
set.seed(1105)
library(BiocParallel)
param <- SnowParam(workers = parallel::detectCores() - 2)

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

# # GO database
# gsea_GO <- gsea_wrapper(gsea_stat, gset = list("C5", c("GO:BP","GO:MF","GO:CC")))

# # Reactome database
# gsea_Reactome <- gsea_wrapper(gsea_stat, gset = list("C2", c("CP:REACTOME")))

# Hallmark database
gsea_H <- gsea_wrapper(gsea_stat, gset = "H")

# Plot hallmark GSEA results
temp_plotdir <- here(plot.dir,"GSEA","Hallmark")
dir_create(temp_plotdir)

res <- gsea_H$gsea %>%
  filter(padj < 0.1) %>%
  dplyr::arrange(desc(NES), padj) %>%
  mutate(change = if_else(NES > 0, true = "UP", false = "DOWN"))
lapply(unique(res$change), function(reg){
  if(reg == "UP"){
        pal <- pathway_up_pal
    }else{
        pal<-pathway_down_pal
    }
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
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title.align = 0 ,
      legend.key.width = unit(0.3, "in")
    ) 
  ggsave(here(temp_plotdir,paste0(reg,".pdf")), p, width = 5, height = 5, dpi = 1000)
})

saveRDS(gsea_H, here(plot.dir, "GSEA","Hallmark", "gsea_hallmark.rds"))

#### 4. Ligand/receptor analysis ####
# follow tutorial in paper https://www.nature.com/articles/s41596-024-01045-4 
save.dir <- here(proj.dir,"gex_aggr_prep_combine")

#devtools::install_github("jinworks/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

## 1. prep input data
# For input data, cellchat requires normalized data. 
# Here we use SCT transformed data following https://github.com/sqjin/CellChat/issues/180
data.input <- GetAssayData(srat, assay = 'SCT', layer = "data")
meta <- as.data.frame(srat[[c("celltype","rat_id","age_wks","library_type",
                              "library_id", "celltype_refined1_updated","sample_id")]]) %>%
  mutate(labels = celltype,
         samples = sample_id)
rownames(meta) <- colnames(srat)
## 2. Create cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

## 3. L-R database ##
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB) 
dplyr::glimpse(CellChatDB$interaction) 
CellChatDB.use <- CellChatDB

## 4. Update cellchat object ##
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 

## 5. Identify overexpressed ligands
future::plan("multicore", workers = 8) # Do not use multisession, or the RAM will not be enough
#For 400 Gb RAM (including memory swap). This cannot exceed local disk limit (~100G).
options(future.globals.maxSize = 400000 * 1024^2) 
cellchat <- identifyOverExpressedGenes(cellchat) 

## 6. Identify overexpressed ligand-receptor interactions ##
cellchat <- identifyOverExpressedInteractions(cellchat)

## 7. smooth data ## (Optional) we did not smooth data data here

## 8. Infer cell-cell interaction at L-R pair level
# Check some signaling genes first
computeAveExpr(cellchat, features = c("Cxcl12","Cxcr4"), type = "triMean") # approximate 25% truncated mean
computeAveExpr(cellchat, features = c("Cxcl12","Cxcr4"), type = "truncatedMean", trim = 0.1) 

computeAveExpr(cellchat, features = c("Spp1","Cd44"), type = "triMean")
computeAveExpr(cellchat, features = c("Spp1","Cd44"), type = "truncatedMean", trim = 0.1) 
# Based on the result, may use 10% truncated mean.
# since PT-injured are a lot less than PT-healthy (PT-S1, PT-S2, PT-S3), we may need to consider cell proportion
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE,
                              population.size = TRUE)

## 9. Filter cell-cell communication ##
cellchat <- filterCommunication(cellchat, min.cells = 50)

## 10. Infer cell-cell communication at a signaling pathway level.
cellchat <- computeCommunProbPathway(cellchat)

# Extract the inferred cellcular communication network at L-R and signaling pathway levels as dataframe
df.LR.net <- subsetCommunication(cellchat,slot = "net")
df.P.net <- subsetCommunication(cellchat,slot = "netP")
write_csv(df.LR.net, file = here(save.dir, "cellchat_network_LR.csv"))
write_csv(df.P.net, file = here(save.dir, "cellchat_network_signal_pathways.csv"))


## 11. Calculate the aggregated cell-cell communication network ## 
# Calculate across all groups
cellchat <- aggregateNet(cellchat)

# save 
saveRDS(cellchat, file = here(save.dir,"cellchat_obj.rds"))

## 12. visualization ##
temp_plotdir <- here(plot.dir,"cellchat")
dir_create(temp_plotdir)
pathways.show.all <- cellchat@netP$pathways

# circle plot for PT-injured-outbound L-R communications
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[2, ] <- mat[2, ] # only visualize PT-injured cells' outgoing signal strength

pdf(here(temp_plotdir,"circle_PT-injured_all cell types.pdf")) # Interaction strength (weights)
netVisual_circle(mat2,color.use = celltype_cols1,
                 vertex.weight = groupSize,
                 weight.scale = T, edge.weight.max = max(mat2), 
                 title.name = "PT-injured")
dev.off()