# # Written for RAGE24 project
# # This code analyzes the RAGE single-cell libraries 
# Reference: https://www.sc-best-practices.org/ and https://www.nature.com/articles/s41467-022-32972-z 

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
library(dplyr)
library(tidyr)
#library(DoubletFinder)
library(harmony)
library(here)
library(stringr)
library(tibble)
library(future)
library(fs)
library(openxlsx)
library(scDblFinder)
library(SingleCellExperiment)
library(scCustomize)
library(BiocParallel)
set.seed(1105)
plan("multicore", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2) #For 100 Gb RAM. Important to set this maxSize or you will encounter error in SCTransfrom

# #### Install additional libraries ####
# # Run this once.
# BiocManager::install('glmGamPoi') # For SCTransform
# devtools::install_github('immunogenomics/presto') # For FindMarkers

#### Defind directories ####
# Data input dir
aggr_input_dir <- here("garyw","RAGE24","cellranger_aggr_gex_flowcell_combine","outs")
# Result output dir
gex_aggr_prep <- here("garyw","RAGE24","gex_aggr_prep_combine")
plot_dir <- here(gex_aggr_prep,"plots","step1")
dir_create(plot_dir)
plot_output_dir <- here(plot_dir,"scDblFinder_workflow")
dir_create(plot_output_dir)

#### Part 1. Load data and create object ####
counts <- Read10X_h5(here(aggr_input_dir, "count", "filtered_feature_bc_matrix.h5"))
aggcsv <- read.csv(here(aggr_input_dir,"aggregation.csv"))
srat <- CreateSeuratObject(counts = counts, min.cells = 10, min.features = 100, project = "RNA")
# extract GEM groups from individual barcodes using string split and the suffix integer
# use the GEM groups to assign sample origin (experimental data) from the aggregation.csv metadata
gemgroup <- sapply(strsplit(rownames(srat@meta.data), split="-"), "[[", 2) %>% as.numeric()
current.gemgroups <- seq(1, length(unique(aggcsv$sample_id)))
sample_id <- plyr::mapvalues(gemgroup, from = current.gemgroups, to = as.character(aggcsv$sample_id))
srat <- AddMetaData(object=srat, metadata=data.frame(library_id=sample_id, row.names=rownames(srat@meta.data)))

#### Part 2. Add other metadata ####
sample_ID_assay_meta <- read.csv(file = here("garyw","RAGE24","sample_ID_assay_meta_08_19.csv"))
meta <- srat@meta.data
meta <- meta %>%
  rownames_to_column("barcodes")%>%
  left_join(distinct(select(sample_ID_assay_meta, rat_id,library_type,sample_id, library_id, age_wks, Lithium_concentration_mM, Lithium_time_min)))%>%
  column_to_rownames("barcodes")
srat <- AddMetaData(object = srat, metadata = meta)
# Other qc metrics
mt_features <- unique(rtracklayer::import(here("e","Atlas_References","Rat","GRCr8","NC_001665.2_edit.gtf"))$gene)
mt_features <- mt_features[mt_features %in% rownames(srat)]
srat <- PercentageFeatureSet(srat, features = mt_features, col.name = "percent.mt") # Mitochondrial genes. We added BN7.2 mt genes 
srat <- PercentageFeatureSet(srat, pattern = "^Rpl", col.name = "percent.rpl") # Ribosomal genes
srat <- PercentageFeatureSet(srat, pattern = "^Rps", col.name = "percent.rps") # Ribosomal genes
## QC plots
pdf(here(gex_aggr_prep,"plots","step1","step1_qc.pdf"))
VlnPlot(object = srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
VlnPlot(object = srat, features = c("percent.rps", "percent.rpl"), pt.size = 0.01, alpha = 0.4)
VlnPlot(object = srat, features = c("nCount_RNA"), group.by = "library_id",pt.size=0)+
  NoLegend()+
  theme(axis.text.x = element_text(size = 6))
VlnPlot(object = srat, features = c("nFeature_RNA"), group.by = "library_id",pt.size=0)+
  NoLegend()+
  theme(axis.text.x = element_text(size = 6))
VlnPlot(object = srat, features = c("percent.mt"), group.by = "library_id",pt.size=0)+
  NoLegend()+
  theme(axis.text.x = element_text(size = 6))
VlnPlot(object = srat, features = c("percent.rpl"), group.by = "library_id",pt.size=0)+
  NoLegend()+
  theme(axis.text.x = element_text(size = 6))
VlnPlot(object = srat, features = c("percent.rps"), group.by = "library_id",pt.size=0)+
  NoLegend()+
  theme(axis.text.x = element_text(size = 6))
dev.off()

#### Part 3. filter the libraries with low number of cells. Also we only want KYC-LN here. ####
table(srat$library_id)
libraries_ditch <- c("RAGE24-L12-PS-LN-01-DEFND-RNA",
                     "RAGE24-O15-KYC-LN-01-Multiome-RNA",
                     "RAGE24-T20-KYC-LN-01-DEFND-RNA",
                     "RAGE24-T20-KYC-LN-01-Multiome-RNA",
                     "RAGE24-X24-PS-LN-01-DEFND-RNA")
for(lib in libraries_ditch){
  srat <- subset(srat, subset = library_id != lib)
}
ncol(srat) # 138430 cells before QC

#### Part 4. filter the aggregated dataset for low quality cells ####
# Metrics
quantile(srat$percent.mt, probs = seq(0,1,0.01))
quantile(srat$percent.rps,probs = seq(0,1,0.05))
quantile(srat$percent.rpl,probs = seq(0,1,0.05))
quantile(srat$nFeature_RNA, probs = seq(0,1,0.05))
# The percentage of rpl and rps is pretty constant among the libraries, so we do not do heavy filtering here.
# Proximal tubule cells may have a high mt content, so we loosen the mt threshold here.
srat <- subset(srat, subset = nFeature_RNA > 400 
               & nFeature_RNA < 7000
               & percent.mt < 30)
table(srat$library_id)
ncol(srat) # After filtering low quality cells, the number of cells is 131721

##### Part 5. Doublet removal--scDblFinder workflow ####
# This workflow is more stringent than DoubletFinder but takes less time
sce <- as.SingleCellExperiment(srat)
# We do not use "cluster" parameter
sce <- scDblFinder(sce, samples="library_id", removeUnidentifiable = TRUE,
                     BPPARAM=MulticoreParam(10, RNGseed=1105), verbose = TRUE)
table(sce$scDblFinder.class) # 
srat$scDblFinder.class <- sce$scDblFinder.class
srat$scDblFinder.score <- sce$scDblFinder.score
rm(sce) # save some space for SCTransform
rm(counts)
## Plot QC graphs for each library
Idents(srat) <- "library_id"
list.doublet.bc <- lapply(unique(srat$library_id), function(lib_id) {
  sobj <- subset(srat, idents = lib_id)
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  sobj <- ScaleData(sobj)
  sobj <- RunPCA(sobj)
  sobj <- RunUMAP(sobj, dims = 1:50)
  p1 <- hist(sobj$scDblFinder.score, plot = F)
  p1$density <- p1$counts/sum(p1$counts)*100
  p2 <- DimPlot(sobj, group.by = "scDblFinder.class")+
    ggtitle("Doublet class", paste0("Cell identities: ", table(sobj$scDblFinder.class)["singlet"],
                                    " singlets and ", table(sobj$scDblFinder.class)["doublet"], " doublets"))
  p3 <- FeaturePlot(sobj, features = "scDblFinder.score")
  p4 <- VlnPlot(sobj, features = c("nFeature_RNA"), group.by = "scDblFinder.class")
  p5 <- VlnPlot(sobj, features = c("nCount_RNA"), group.by = "scDblFinder.class")
  outFile <- here(plot_output_dir, paste0("step1_doublets.",lib_id,".pdf"))
  pdf(outFile)
  plot(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  dev.off()
})
## visualize doublets before proceeding with preprocessing
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat)
srat <- ScaleData(srat)
srat <- RunPCA(srat, verbose = TRUE)
# ElbowPlot(srat, ndims = 50) # to determine number of dimensions for clustering
srat <- RunHarmony(srat, c("library_id","library_type"), reduction.use = "pca", plot_convergence = TRUE)
srat <- FindNeighbors(srat, dims = 1:50, verbose = TRUE, reduction = "harmony")
srat <- FindClusters(srat, verbose = TRUE, resolution = 0.6)
srat <- RunUMAP(srat, dims = 1:50, verbose = TRUE, reduction = "harmony")
# Visualize
Idents(srat) <- "seurat_clusters"
p1 <- DimPlot_scCustom(srat, reduction = "umap", label = TRUE, raster = FALSE, alpha = 0.5) +
  ggtitle("Clustering with Harmony Including Doublets")
p2 <- DimPlot_scCustom(srat, group.by = "scDblFinder.class", order=TRUE, alpha = 0.5, raster = FALSE) + ggtitle("scDblFinder")
p2.0 <- DimPlot_scCustom(srat, reduction = "umap", split.by = "scDblFinder.class",
                         alpha = 0.5, raster = FALSE, split_seurat = TRUE, figure_plot = FALSE)+
  theme(legend.position = "none")
p2.1 <- VlnPlot_scCustom(srat, features = c("nFeature_RNA"), group.by = "scDblFinder.class", 
                pt.size = 0.005, alpha = 0.1, raster=FALSE)
p2.2 <- VlnPlot_scCustom(srat, features = c("nCount_RNA"), group.by = "scDblFinder.class",
                pt.size = 0.005, alpha = 0.1, raster=FALSE)
pdf(here(plot_output_dir,"step1_doublets.pdf"))
p0 <- hist(srat$scDblFinder.score, plot = F)
p0$density <- p0$counts/sum(p0$counts)*100
plot(p0)
print(list(p1,p2,p2.0,p2.1,p2.2))
dev.off()

# save an object with doublets
saveRDS(srat, here(gex_aggr_prep,"step1_doublets_scDblFinder.rds"))

# subset the object for singlets to do additional processing
table(srat$scDblFinder.class) # 115810 singlets vs 15911 doublets
table(srat$scDblFinder.class)[2]/ncol(srat) # 12.08% doublet rate
Idents(srat) <- "scDblFinder.class"
srat <- subset(srat, idents = "singlet")
saveRDS(srat, here(gex_aggr_prep,"step1_no_doublets_scDblFinder.rds"))
# srat <- readRDS(here(gex_aggr_prep,"step1_no_doublets_scDblFinder.rds"))


# #### Part 5b. Doublet removal DoubletFinder workflow. Not used ####
# We leave it here for reference. 
# FindDoublets <- function(library_id, seurat_aggregate) {
#   seurat_obj <- subset(seurat_aggregate, idents = library_id)
#   seurat_obj <- NormalizeData(seurat_obj)
#   seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
#   seurat_obj <- ScaleData(seurat_obj)
#   seurat_obj <- RunPCA(seurat_obj)
#   # ElbowPlot(seurat_obj)
#   seurat_obj <- RunUMAP(seurat_obj, dims = 1:50)
#   # DimPlot(seurat_obj)
  
#   ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#   sweep.res.list_kidney <- paramSweep(seurat_obj, PCs = 1:50, sct = F, num.cores=10)
#   sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
#   bcmvn_kidney <- find.pK(sweep.stats_kidney)
#   pK <- bcmvn_kidney %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
#     filter(BCmetric == max(BCmetric)) %>%
#     select(pK) 
#   pK <- as.numeric(as.character(pK[[1]]))
#   seurat_doublets <- doubletFinder(seurat_obj, PCs = 1:50, pN = 0.25, pK = pK,
#                                       nExp = round(0.05*length(seurat_obj@active.ident)), 
#                                       reuse.pANN = FALSE, sct = F)
  
#   # create doublet groupings and visualize results
#   DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
#   pANN <- names(seurat_doublets@meta.data) %>% str_subset("pANN")
  
#   p1 <- ggplot(bcmvn_kidney, aes(x=pK, y=BCmetric)) +
#     geom_bar(stat = "identity") + 
#     ggtitle(paste0("pKmax=",pK)) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   p2 <- DimPlot(seurat_doublets, group.by = DF.class)
#   p3 <- FeaturePlot(seurat_doublets, features = pANN)
  
#   outFile <- here(gex_aggr_prep,"plots","step1", paste0("step1_doublets.",library_id,".pdf"))
#   pdf(outFile)
#   print(p1) # need to use print() when drawing pdf in a function call
#   print(p2)
#   print(p3)
#   dev.off()
  
#   # create a df of barcodes and doublet designations
#   df_doublet_barcodes <- as.data.frame(cbind(rownames(seurat_doublets@meta.data), seurat_doublets@meta.data[[DF.class]]))
#   return(df_doublet_barcodes)
# }

# #Use the function to get a list of doublet barcodes
# Idents(srat) <- "library_id"
# list.doublet.bc <- lapply(unique(srat$library_id), function(x) {FindDoublets(x, seurat_aggregate = srat)})
# doublet_id <- list.doublet.bc %>%
#   bind_rows() %>%
#   dplyr::rename("doublet_id" = "V2") %>%
#   tibble::column_to_rownames(var = "V1") # this is the barcode column
# table(doublet_id) # quantify total doublet vs. singlet calls (expect ~6% doublets)

# # add doublet calls to aggregated snRNA object as doublet_id in meta.data slot
# srat <- AddMetaData(srat, doublet_id)
# # enable parallel processing via future package
# plan("multicore", workers = 10) 
# options(future.globals.maxSize = 100000 * 1024^2) # for 100 Gb RAM
# plan()
# # visualize doublets before proceeding with preprocessing
# # Here, choosing SCTransform is for consistency. In practice, this visualization step can be carried out by the standard Seurat workflow to save time and RAM.
# srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA"), verbose = TRUE, vst.flavor="v2")
# srat <- RunPCA(srat, verbose = TRUE)
# # ElbowPlot(srat, ndims = 50) # to determine number of dimensions for clustering
# srat <- RunHarmony(srat, "library_id", plot_convergence = TRUE, assay.use="SCT")
# srat <- FindNeighbors(srat, dims = 1:50, verbose = TRUE, assay = "SCT", reduction = "harmony")
# srat <- FindClusters(srat, verbose = TRUE, resolution = 0.6)
# srat <- RunUMAP(srat, dims = 1:50, verbose = TRUE,assay = "SCT", reduction = "harmony")

# srat@meta.data$doublet_viz <- ifelse(srat@meta.data$doublet_id == "Singlet",0,1)

# Idents(srat) <- "seurat_clusters"
# p1 <- DimPlot(srat, reduction = "umap", label = TRUE) +
#   ggtitle("Clustering with Harmony Including Doublets")
# p2 <- FeaturePlot(srat, features = "doublet_viz", order=TRUE) + ggtitle("DoubletFinder")
# pdf(here(gex_aggr_prep,"plots","step1_doublets.pdf"))
# print(list(p1,p2))
# dev.off()

# # save an object with doublets
# # saveRDS(srat, here(gex_aggr_prep,"step1_doublets.rds"))
# # srat <- readRDS(here(gex_aggr_prep,"step1_no_doublets.rds"))
# # subset the object for singlets to do additional processing
# Idents(srat) <- "doublet_id"
# srat <- subset(srat, idents = "Singlet")
# # saveRDS(srat, here(gex_aggr_prep,"step1_no_doublets.rds"))
# # remove unused metadata columns
# meta <- srat@meta.data %>%
#   dplyr::select(-doublet_id, -doublet_viz, -orig.ident)
# srat@meta.data <- meta

#### Part 6. Post filtering of low-quality cells and doublet-removal ####
# Up to this point, you may need to save the srat object and restart the session to free up some memory
# srat <- readRDS(here(gex_aggr_prep,"step1_no_doublets_scDblFinder.rds"))
srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA"), verbose = TRUE, vst.flavor = "v2")
srat <- RunPCA(srat, verbose = TRUE)
# ElbowPlot(srat, ndims = 50) # to determine number of dimensions for clustering
srat <- RunHarmony(srat, c("library_id","library_type"), plot_convergence = TRUE, assay.use="SCT") #https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
srat <- FindNeighbors(srat, dims = 1:50, verbose = TRUE, reduction = "harmony")
srat <- FindClusters(srat, verbose = TRUE, resolution = 1.5, future.seed=TRUE, method = "igraph")
srat <- RunUMAP(srat, dims = 1:50, verbose = TRUE, reduction = "harmony")

# visualize the clustering
Idents(srat) <- "seurat_clusters"
p2 <- DimPlot(srat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Clustering with Harmony No Doublets") +
  NoLegend()

#### Additional QC ####
p3 <- DimPlot(srat, group.by = "library_type",raster = FALSE)
Idents(srat) <- "age_wks"
p4 <- DimPlot(srat)
Idents(srat) <- "library_type"
p5 <- DimPlot(srat, cells.highlight = colnames(subset(srat, idents = "Multiome")),
              sizes.highlight = 0.1, pt.size = 0.1, alpha = 0.5, raster = FALSE)+
  ggtitle("Multiome")+
  theme(legend.position="none")
p6 <- DimPlot(srat, cells.highlight = colnames(subset(srat, idents = "DEFND")), 
              sizes.highlight = 0.1, pt.size = 0.1, alpha = 0.5, raster = FALSE)+
  ggtitle("DEFND")+
  theme(legend.position="none")
p7 <- p5|p6
p8 <- DimPlot_scCustom(srat, group.by = "library_id", alpha = 0.5)+
  theme(legend.text = element_text(size = 5))
# Print UMAP clustering
print("Drawing UMAP clusters")
pdf(here(plot_output_dir,"step1_clusters.pdf"))
tryCatch(print(p2), error=function(e) NULL)
#tryCatch(print(p1), error=function(e) NULL)
tryCatch(print(p3), error=function(e) NULL)
tryCatch(print(p7), error = function(e) NULL)
tryCatch(print(p4), error=function(e) NULL)
tryCatch(print(p8), error = function(e) NULL)
dev.off()

#### Part 7. Plot marker gene expression ####
# Markers
Idents(srat) <- "seurat_clusters"
marker.genes <- c("Cubn","Havcr1","Slc5a1","Slc5a2","Vcam1", "Prom1", # PT and PTVCAM1+ markers
                  "Creb5","Adamts1","Spp1", # PT injury markers
                  "Cfh","Akap12","Ncam1", # PEC
                  "Aqp1", "Cryab","Proser2","Epha7",# TL (Thin Limb)
                  "Slc12a1", # TAL NKCC2
                  "Cldn10", #MTAL
                  "Cldn16", #CTAL
                  #"S100a2", #ATL
                  "Slc12a3","Trpm6", # DCT1 and DCT2 NCC
                  "Scnn1g","Trpv5", # DCT2/CNT ENaC
                  "Calb1", # CNT
                  "Aqp2", # PC # PC also high in Scnn1g
                  "Atp6v0d2", # ICA and ICB
                  "Slc4a1","Slc26a7", # ICA
                  "Slc26a4", # ICB
                  "Nphs1","Nphs2", # PODO
                  "Pecam1","Flt1", # ENDO 
                  "Igfbp5","Igfbp7", # PTC and AVR
                  "Plvap", # PTC and AVR https://www.nature.com/articles/s41467-019-12872-5
                  "Ehd3", # GEC
                  "Slc6a6","Slc14a1","Aqp1", # EA and DVR
                  "Nos1", # MD
                  "Itga8","Pdgfrb","Meis2","Piezo2","Ren", # MES and JGA
                  "Acta2","Cald1", # FIB
                  "Prox1","Flt4","Pdpn", # Lymphatics
                  "Ptprc","Cd3e","Ms4a1", # Lymphocytes
                  "Fcgr3a","Cd14","Csf1r") # Monocyte / Macrophage
print("Drawing UMAP markers")
pdf(here(plot_output_dir,"step1_markers.pdf")) 
lapply(marker.genes, function(gene) {
  tryCatch({plot1 <- FeaturePlot(srat, features=gene, reduction="umap") + theme(axis.title.x = element_blank(),
                                                                      axis.title.y = element_blank(),
                                                                      axis.ticks = element_blank())},
           warning = function(w) return(NULL))
  tryCatch({plot2 <- VlnPlot(srat, features = gene, pt.size = 0.5) + theme(legend.position="none")},
           warning= function(w) return(NULL))
  tryCatch({print(plot1)}, warning=function(w) return(NULL)) 
  tryCatch({print(plot2)}, warning=function(w) return(NULL)) 
})
dev.off()

dotp <- DotPlot(srat, features = marker.genes, cols = c("lightyellow", "red")) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, size = 8))
ggsave(filename = here(plot_output_dir, "step1_markers.png"), 
       plot = dotp, width = 14, height = 12)

# save the preprocessed rna file
saveRDS(srat, file = here(gex_aggr_prep,"step1_prep_scDblFinder.rds"))

#### Part 8. Find differentially expressed genes ####
# Create a folder to write an xlsx file for DE markers
markers_dir <- here(gex_aggr_prep,"markers")
dir_create(here(markers_dir))

# With SCT
Idents(srat) <- "seurat_clusters"
srat <- PrepSCTFindMarkers(srat, assay = "SCT")
deg.markers <- FindAllMarkers(srat, assay = "SCT", min.pct = 0.2)
markers.list <- split(deg.markers, f = deg.markers$cluster)
write.xlsx(markers.list, file = here(markers_dir,"deg_celltype_markers.xlsx"),
           sheetName = levels(srat@meta.data$seurat_clusters),
           rowNamees = F,overwrite = T)

