# # Written for RAGE24 project
# # This code annotates the filtered single-cell libraries
# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --workdir $HOME \
# --name gexaggr \
# -v /mnt/g/cellranger_atac_counts:$HOME/cellranger_atac_counts \
# -v /mnt/g:$HOME/g \
# -v /mnt/e:$HOME/e \
# -v /mnt/d/garyw:$HOME/garyw \
# -v /mnt/e/Atlas_References:$HOME/Atlas_References \
# -v $HOME:$HOME \
# -v /mnt/g/cellranger_atac_counts:$RHOME/cellranger_atac_counts \
# -v /mnt/g/:$RHOME/g \
# -v /mnt/e:$RHOME/e \
# -v /mnt/d/garyw:$RHOME/garyw \
# -v /mnt/e/Atlas_References:$RHOME/Atlas_References \
# -v /var/run/docker.sock:/var/run/docker.sock \
# -e PASSWORD=garywang \
# -e DISABLE_AUTH=TRUE \
# -p 8787:8787 \
# garywang7/ragemultiome:1.0.4



library(Seurat)
library(ggplot2)
library(here)
library(dplyr)
library(scCustomize)
library(tibble)
library(stringr)
library(fs)
library(openxlsx)
library(harmony)
library(ComplexHeatmap)
library(future)
plan("multicore", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2) #For 100 Gb RAM. Important to set this maxSize or you will encounter error in SCTransfrom

#### Install packages ####
# If needed.
# BiocManager::install("ComplexHeatmap")

#### Define directories ####
# Result output dir
gex_aggr_prep <- here("garyw","RAGE24","gex_aggr_prep_combine")
plot_output_dir <- here(gex_aggr_prep,"plots", "step2", "scDblFinder_workflow")
dir_create(plot_output_dir)

srat <- readRDS(here(gex_aggr_prep,"step1_prep_scDblFinder.rds"))

#### Part 1. Examine markers ####
Idents(srat) <- "seurat_clusters"
## nCount RNA and nFeatures
ggsave(plot = VlnPlot_scCustom(srat, features = c("nCount_RNA"), pt.size = 0.0001) +
  theme(legend.position = NULL), filename = here(plot_output_dir,"nCountRNA.png"),width = 14, height = 8)
ggsave(plot = VlnPlot_scCustom(srat, features = c("nFeature_RNA"), pt.size = 0.0001) +
         theme(legend.position = NULL), filename = here(plot_output_dir,"nFeatureRNA.png"),width = 14, height = 8)
## PT
# Markers from https://www.pnas.org/doi/full/10.1073/pnas.2005477117  and Humphreys Lab Kidney Interactive Transcriptomics
# Use this paper's marker for S1,2 and 3: https://pubmed.ncbi.nlm.nih.gov/37906287/
# Refer to the paper "Analysis of individual patient pathway coordination in a cross-species single-cell kidney atlas" at https://www.nature.com/articles/s41588-025-02285-0 
PT.subset <- c(1,2,4,5,6,8,15,16,18,20,25,30,35,37)
VlnPlot(srat, features = "Prom1")
Clustered_DotPlot(srat, features = c("Slc22a7", "Slc5a12","Slc13a3","Slc16a9","Havcr1"), 
                  flip = TRUE, cluster_feature = FALSE)
Clustered_DotPlot(srat, features = c("Cubn", "Slc13a1","Slc5a12","Slc13a3","Slc22a6", # General PT markers
                                     "Prodh2","Slc5a2","Slc22a8", # PT-S1
                                     "Slc34a1","Slc22a7", # PT-S2 (S2 share a lot of same marker with PT-S1 and S3)
                                     "Slc5a11","Slc22a24","Slc7a13","Satb2"), # PT-S3
                  flip = TRUE, cluster_feature = FALSE, assay = "SCT") # From KPMP. PT-S1, S2 and S3
DotPlot_scCustom(srat, features = c("Slc13a1","Slc7a13"),
                 idents = PT.subset) # From Abedini et al. https://pubmed.ncbi.nlm.nih.gov/37906287/ 
# Slc13a1 is PCT marker and Slc7a13 for PST marker.
# injured PT
# Seee SISKA conserved markers: https://www.nature.com/articles/s41588-025-02285-0
DotPlot_scCustom(srat, features = c("Vcam1","Havcr1", "Pdgfb","Creb5","Adamts1","Dcdc2"),
                 idents = PT.subset)
# Clus8 These cells show similarity to PT cells but is not Havcr1+ or Pdgfb+
# Consider clus8 to be TL (thin limb) population
Cluster_Highlight_Plot(srat, cluster_name = 8)
DotPlot_scCustom(srat, features = c("Cryab","Tacstd2", "Slc44a5","Klrg2","Col26a1","Boc"),
                 idents = NULL) # TL markers from KPMP paper
# clus 18 might be a doublet population
Cluster_Highlight_Plot(srat, cluster_name = 18, pt.size = 0.3)
clus18 <- subset(srat, idents = 18)
hist(clus18$scDblFinder.score)
table(clus18$library_type)
clus18_ident <- ifelse(srat$seurat_clusters == 18, yes= "Yes", no = "No" )
chisq.test(srat$library_type, clus18_ident) # Do a chi-square to test if DEFND is over-represented in clus18

## FIB
DotPlot_scCustom(srat, features = c("Col1a2","C7","Fbln5","Cald1","Pdgfrb"))

## TAL
VlnPlot(srat, features = c("Slc12a1","Cldn10","Cldn16")) # TAL2 or M-TAL(Cldn10+) and TAL1 or C-TAL(Cldn16+) # See https://www.nature.com/articles/s41581-022-00553-4#Sec6
DotPlot_scCustom(srat, features = c("Slc12a1","Esrrb","Egf","Enox1","Tmem207","Nos1","Robo2","Calcr"))
# clus22 might be a doublet population as it scattered all over. Also, the nCount and nFeature are much lower.
clus22 <- subset(srat, idents = 22)
hist(clus22$scDblFinder.score)
Cluster_Highlight_Plot(srat, cluster_name = 22)

## PEC
# SISKA conserved markers.
VlnPlot_scCustom(srat, features= c("Cldn1","Akap12","Pax8","Aldh1a2"))
DotPlot_scCustom(srat, features = c("Cldn1","Akap12","Ncam1","Aldh1a2"))

## CNT
DotPlot_scCustom(srat, features = c("Slc8a1","Scn2a","Hsd11b2","Calb1","Kitlg","Pcdh7"))# CNT
DotPlot_scCustom(srat, features = c("Ralyl","Tox","Sgpp1","Scnn1g","Kcnip1"))

## PC
DotPlot_scCustom(srat,features = c("Gata3","Aqp2","Aqp3","Scnn1g","Scnn1b","Slc8a1","Hsd11b2","Calb1")) # PC-CNT --19

## POD
DotPlot_scCustom(srat, features = c("Wt1","Nphs1","Nphs2"))
FeaturePlot_scCustom(srat, features=c("Wt1","Nphs1","Nphs2"))


## Immune cells
# Lymphocyte markers human analogs: https://www.sciencedirect.com/science/article/pii/S1535610822005487?via%3Dihub
FeaturePlot_scCustom(srat, features = c("Cd4", "Cd247","Cd3e","Cfh","Cd8a","Cd8b"), raster = FALSE)
DotPlot_scCustom(srat, features = c("Cd8a","Cd8b","Cd4","Il2ra"))# T cells Note CD4 can be expressed in Monocytes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4136363/
DotPlot_scCustom(srat, features = c("Cd74", "RT1-Ba", "RT1-Db1")) # Dendritic cells
DotPlot_scCustom(srat, features = c("Cd19","Ms4a1","Cd79a","Pax5","Cd96","Cd247","Cd3e","Cd3g","Cd8a","Cd8b")) # B cells and T cells
DotPlot_scCustom(srat,features = c("Apoe","Lyz2","Csf1r","C1qa"))# Macrophages (Lyz2) and Monocytes(C1qa)

## IC
Stacked_VlnPlot(srat, features = c("Atp6v0d2","Slc4a1","Slc26a4")) # General IC, IC-A and IC-B respectively
Cluster_Highlight_Plot(srat, cluster_name = 34)
# Cluster 34 is very likely a doublet population.
clus34 <- subset(srat, idents=34)
table(clus34$library_type) # 847 DEFND and 12 Multiome
Cluster_Highlight_Plot(srat, cluster_name = 34)+
  ggtitle(paste0("Number of Cells: ",ncol(clus34)), "= 847 DEFND + 12 Multiome")
hist(clus34$scDblFinder.score)

## Endothelial cells (EC)
# Markers from KPMP paper https://www.nature.com/articles/s41586-023-05769-3#Sec8 Supplement Table 5
DotPlot_scCustom(srat, features = c("Ptprb","Pecam1","Flt1","Plvap","Meis2","Emcn")) # From KPMP marker list

## Vascular smooth muscle cell and pericyte VSMC_P
# See KPMP and also table 1 of this paper https://www.nature.com/articles/s41581-021-00474-8
DotPlot_scCustom(srat, features = c("Notch3","Pdgfrb","Ren","Itga8","Pip5k1b","Piezo2")) # MC mesangial cells
DotPlot_scCustom(srat, features = c("Notch3","Pdgfrb","Itga8","Acta2","Ntrk3","Myh11","Rgs5","Adcy3","Pdgfra")) # VSMC vascular smooth muscle cells

#### Part 2. Annotating clusters ####
# Read annotations
anno_table <- read.xlsx(here(gex_aggr_prep,"RAGE24_single_cell_cluster_annotation_updated.xlsx"))
Idents(srat) <- "seurat_clusters"
srat[["orig.clusters"]] <- Idents(srat)
anno.ids.raw.updated <- anno_table$celltype_raw_updated # raw cell type annotation
anno.ids.refined1.updated <- anno_table$celltype_refined1_updated # Updated the TL population based on refined1

names(anno.ids.raw.updated) <- anno_table$cluster
names(anno.ids.refined1.updated) <- anno_table$cluster

# Rename the clusters
Idents(srat) <- "seurat_clusters"
srat <- RenameIdents(srat, anno.ids.raw.updated)
srat[["celltype_raw_updated"]] <- Idents(srat)

Idents(srat) <- "seurat_clusters"
srat <- RenameIdents(srat, anno.ids.refined1.updated)
srat[["celltype_refined1_updated"]] <- Idents(srat)

# Mark doublets
srat[["is_doublet"]] <- ifelse(srat$celltype_raw=="Doublet", "doublet","singlet")
# Examine the cells in the clusters
table(srat$celltype_raw_updated)
table(srat$celltype_refined1_updated)
table(srat$is_doublet)

# Adjust this for the refine level of annotation
srat$celltype <- srat$celltype_raw_updated

# Inspect doublet population
Idents(srat) <- "celltype"
Cluster_Highlight_Plot(srat, cluster_name = "Doublet")

pdf(here(plot_output_dir,"step2_manual_doublets.pdf"))
DimPlot_scCustom(srat, reduction = "umap", group.by = "celltype", raster = FALSE, label = TRUE,
                 alpha = 0.5, pt.size = 0.4)+
  ggtitle("Cell types before manual removal of doublets")
Cluster_Highlight_Plot(srat, cluster_name = "Doublet")
DimPlot_scCustom(srat, reduction = "umap", split.by = "is_doublet", pt.size = 0.01,
                 alpha = 0.5, raster = FALSE, split_seurat = TRUE, figure_plot = FALSE)+
  theme(legend.position = "none")
VlnPlot_scCustom(srat, features = c("nCount_RNA"), group.by = "is_doublet",
                 pt.size = 0.005, alpha = 0.1, raster=FALSE)
VlnPlot_scCustom(srat, features = c("nFeature_RNA"), group.by = "is_doublet", 
                 pt.size = 0.005, alpha = 0.1, raster=FALSE)
dev.off()

# Remove the doublet cell types
Idents(srat) <- "is_doublet"
srat <- subset(srat, idents = "singlet")
srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA"), verbose = TRUE, vst.flavor = "v2")
srat <- RunPCA(srat, verbose = TRUE)
# ElbowPlot(srat, ndims = 50) # to determine number of dimensions for clustering
srat <- RunHarmony(srat, c("library_id","library_type"), plot_convergence = TRUE, assay.use="SCT") #https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
srat <- RunUMAP(srat, dims = 1:50, verbose = TRUE, reduction = "harmony")

# Reorder the idents for raw and refined1
Idents(srat) <- "celltype"
srat$celltype <- factor(srat$celltype, levels = 
                          c("PT","PT-injured", "TL1", "TL2","C-TAL","M-TAL","DCT","CNT","CNT-PC","PC","IC-A","IC-B","EC","FIB","POD","PEC","VSMC/P","IMM"))
srat$celltype_raw_updated <- factor(srat$celltype_raw_updated, levels = 
                          c("PT","PT-injured", "TL1", "TL2","C-TAL","M-TAL","DCT","CNT","CNT-PC","PC","IC-A","IC-B","EC","FIB","POD","PEC","VSMC/P","IMM"))
srat$celltype_refined1_updated <- factor(srat$celltype_refined1_updated, levels = 
                                   c("PT-S1","PT-S2","PT-S3","PT-injured","TL1","TL2", "C-TAL","M-TAL","DCT","CNT","CNT-PC","PC","IC-A","IC-B","EC","FIB","POD","PEC","VSMC/P","IMM"))

# Plot annotation
pdf(here(plot_output_dir,"step2_anno.pdf"))
DimPlot_scCustom(srat, reduction = "umap", group.by = "celltype_raw_updated", label = TRUE)
VlnPlot(srat, features = c("nCount_RNA"), group.by = "celltype_raw_updated", raster = FALSE, pt.size = 0.001, alpha = 0.01)
VlnPlot(srat, features = c("nFeature_RNA"), group.by = "celltype_raw_updated", raster = FALSE, pt.size = 0.001, alpha = 0.01)
DimPlot_scCustom(srat, reduction = "umap", group.by = "celltype_refined1_updated", 
                 label.size = 3.5,
                 pt.size = 0.1,label = TRUE,
                 repel = TRUE,color_seed = 1105)+
  theme(legend.position = "none")
dev.off()

#### Part 3. Plotting signature gene plot ####
marker.genes <- c( "Cubn", # all PT
                  "Slc7a7", # PT-S1
                  "Slc13a1", # PCT （PT-S1 and S2)
                  "Slc7a13", # PST (PT-S2 and S3)
                  #"Creb5", "Adamts1", # PT-injured
                  "Il34","Adamts1", # PT-injured1
                  #"Spp1", # PT-injured2
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
dotp
ggsave(filename = here(plot_output_dir, "step2_markers.png"), 
       plot = dotp, width = 12, height = 10)

#### Saving the object ####
# Save object
Idents(srat) <- srat$celltype
DefaultAssay(srat) <- "SCT"
saveRDS(srat, file = here(gex_aggr_prep,"step2_anno.rds"))
# Save cell identities
write.csv(srat@meta.data, file = here(gex_aggr_prep, "gex_aggr_anno.csv"), row.names = TRUE)
