# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
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

library(here)
library(scCustomize)
library(fs)
library(Seurat)
library(tidyverse)

#### Directories ####
proj.dir <- here("garyw/RAGE24")
plot.all.dir <- here(proj.dir,"publication/figures")
plot.dir <- here(plot.all.dir,"supplements")
dir_create(plot.dir)

#### Load data ####
srat <- readRDS(here(proj.dir,"gex_aggr_prep_combine","step2_anno.rds"))

#### 1. QC plot ####
# DEFND-seq vs Multiome gene counts/UMI
method_col <- c(
  "DEFND" = "#D81B60",  # vivid magenta (warm, not orange)
  "Multiome"    = "#1F4E79"   # deep navy (cold)
)
p <- VlnPlot_scCustom(srat, features = c("nFeature_RNA"), group.by = "library_type",
                      colors_use = method_col,
                 pt.size = 0.001, raster = TRUE, alpha = 0.05)+
  scale_y_log10()
ggsave(here(plot.dir, "Gene counts DEFND vs Multiome.pdf"), width = 6, height = 5, 
       plot = p)
p <- VlnPlot_scCustom(srat, features = c("nCount_RNA"), group.by = "library_type",
                      colors_use = method_col, y.max = 20000,
                      pt.size = 0.001, raster = TRUE, alpha = 0.05)+
  scale_y_log10()
ggsave(here(plot.dir, "UMI DEFND vs Multiome.pdf"), width = 6, height = 5, 
       plot = p)

#### 2. Cell type markers ####
# Cell type marker dotplot
marker.genes <- c( "Cubn", # all PT
                   "Slc7a7", # PT-S1
                   "Slc13a1", # PCT (PT-S1 and S2)
                   "Slc7a13", # PST (PT-S2 and S3)
                   "Creb5","Dcdc2", # PT-injured
                   #"Il34","Adamts1", # PT-injured
                   "Cryab", # TL, TL1
                   "Spp1", # TL2
                   "Slc12a1", # TAL1 and TAL2
                   "Cldn16", #TAL1
                   "Cldn10", #TAL2
                   "Slc12a3", # DCT
                   "Slc8a1","Calb1", # CNT
                   "Aqp2","Scnn1g", # PC # PC also high in Scnn1g
                   "Atp6v0d2", # ICA and ICB
                   "Slc4a1", # ICA
                   "Slc26a4", # ICB
                   "Pecam1", # EC
                   "Fbln5", # FIB "C7",
                   "Nphs1", # POD
                   "Akap12", # PEC
                   "Notch3", # VSMC/P (Vascular cmooth muscle cell/pericyte)
                   "Ptprc") # Immune
dotp <- DotPlot_scCustom(srat, features = marker.genes, cols = c("lightyellow", "red"), 
                         group.by = "celltype_refined1_updated",
                         dot.min = 0.15) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust = 1, face = "italic"))+
  theme(axis.text.y = element_text(size = 12))+
  scale_y_discrete(limits = rev)
ggsave(here(plot.dir,"cell marker dotplot.pdf"), width = 12, height = 9, 
       plot = dotp, dpi = 1200)

### 3. circulating biomarkers (human orthologs) in rat ####
# Analogs of circulating injured PT biomarkers found in human. This is manually 
rat_genes <- c("Tnfrsf11a", "Tnfrsf14", "Eda2r", "Il1r1","Tnfrsf4", "Tnfrsf1a", "Tnfrsf1b",
               "Ltbr", "Fas","Cd27", 
               "Tnfrsf10b", # Ortholog for TNFRSF10A as well
               "Tnfrsf12a","Tnfrsf19","Relt","Tnfrsf21",
               "Havcr1", "Epha2", "Tgfbr2",
               "Cd300lg", "Pilra", "Vsig4","Colec12", "Ephb4",
               "Folr1", "Layn",
               "Scarb2", "Creld2", "Ctsz","Klk11", "Cdh3","Dll1",
               "Efna4","Pgf","Hspg2","Ambp","Fstl3","Il18bp", # No ortholog PI3
               "Wfdc2","Cd99l2","Dsc2","Esam","Nectin4","Tff3"
               )

# calculate DEGs
deg <- FindMarkers(srat, ident.1 = "PT-injured", ident.2 = "PT")
df <- deg %>% 
  tibble::rownames_to_column("gene") %>%
  filter(gene %in% rat_genes)

# Plot
df_plot <- df %>%
  mutate(
    sig = ifelse(p_val_adj < 0.05, "Sig", "NS")
  )
sig_cols <- c(
  "NS"  = "#9A9A9A",  # neutral grey (background)
  "Sig" = "#9E4F5A"   # muted rose / brick (statistical highlight)
)

p <- ggplot(df_plot, aes(x = avg_log2FC, y = pct.1)) +
  geom_point(aes(color = sig, fill = sig), size = 1.5, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = gene, color = sig),
    size = 4,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = sig_cols) +
  scale_fill_manual(values = sig_cols) +
  labs(x = "Log2 fold change of PT-injured relative to non-injured PT",
       y = "Proportion of PT-injured expressing gene") +
  ggpubr::theme_pubr()
ggsave(here(plot.dir, "Circulating biomarkers rat.pdf"), plot = p, width = 6)
