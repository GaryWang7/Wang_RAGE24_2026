# Written for RAGE24 deconvloluting the celltypes in AGES data
# In this process, we perform both BayesPrism and MuSiC. BayesPrism has a handy function to filter 
# genes in the reference with low specificity and low expression, which is important for deconvolution. 

# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --workdir $HOME \
# --name deconvolution \
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

# Install the libraries
library("devtools")
install_github("Danko-Lab/BayesPrism/BayesPrism")
BiocManager::install("TOAST")
devtools::install_github('xuranw/MuSiC')
install.packages("ggsci")
BiocManager::install("DESeq2")

# Load libraries
library(Seurat)
library(here)
library(fs)
library(tidyverse)
library(BayesPrism)
library(ggsci)
library(MuSiC)
library(DESeq2)
library(SingleCellExperiment)

# Define directories
gex_aggr_prep <- here("garyw","RAGE24","gex_aggr_prep_combine")
ages_dir <- here("~/v/AgeRat")
plot_outdir <- here(ages_dir,"plots","step4_celltype_deconv")
analysis_outdir <- here(ages_dir, "analysis_data","deconv")
dir_create(plot_outdir)
dir_create(analysis_outdir)

# Load data 
fc <- read.csv(here(ages_dir,"counts","counts.csv")) %>%
  column_to_rownames("X")
fc.meta <- read.csv(here(ages_dir,"analysis_data","meta.csv")) %>%
  dplyr::rename("sample_ID" = "X")
srat <- readRDS(here(gex_aggr_prep,"step2_anno.rds"))
srat.filt <- subset(srat, !str_detect(celltype, "TL1")) # Remove TL1 cell type, which is transitional and contain too many mitochondrial counts.

##### 1. BayesPrism ####
# Follow https://github.com/Danko-Lab/BayesPrism Vignette for deconvoluting
# Prepare input data
# As the single cell data cannot detect mitochondrial reads well, we remove the PT population with a lot of mitochondrial genes, and perform the deconvolution
bk.dat <- t(fc)
sc.dat <- LayerData(srat.filt, assay = "RNA", layer = "counts") %>% 
  as.matrix() %>%
  t()
gc()

# The definition of cell type and cell state can be somewhat arbitrary (similar to the issue of assigning cell types for scRNA-seq) 
# and depends on the question of interest. 
# Usually, a good rule of thumb is as follows. 1) Define cell types as the cluster of cells having a sufficient number of significantly differentially expressed genes than other cell types, e.g., greater than 50 or even 100. For clusters that are too similar in transcription, we recommend treating them as cell states, which will be summed up before the final Gibbs sampling. Therefore, cell states are often suitable for cells that form a continuum on the phenotypic manifold rather than distinct clusters. 2) Define multiple cell states for cell types of significant heterogeneity, such as malignant cells, and of interest to deconvolve their transcription.

cell.type.labels <- srat.filt$celltype
cell.state.labels <- srat.filt$celltype_refined1_updated

# QC cell type/state assignment
# merge those cell types/states with the most similar cell types/states, or remove them (if re-clustering and merging is not appropriate).
dir_create(here(plot_outdir,"BayesPrism"))
plot.cor.phi(input=sc.dat,
            input.labels=cell.state.labels,
            title="cell state correlation",
            cexRow=0.8, cexCol=0.8,
            margins=c(6,6),
            pdf.prefix = here(plot_outdir,"BayesPrism","celltype_state_corr"))
gc()
plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              cexRow=0.8, cexCol=0.8,
              margins=c(6,6),
              pdf.prefix = here(plot_outdir,"BayesPrism","celltype_type_corr"))
gc()

# QC if there are outlier genes
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  # use mouse to approximate our data
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE, #return the data used for plotting. 
  pdf.prefix=here(plot_outdir,"BayesPrism","gene_specificity"))
# sc.stat contains the information of those genes
gc()

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="mm", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE,
  pdf.prefix=here(plot_outdir,"BayesPrism","gene_specificity"))
gc()

# Filter genes with low specificity and low expression
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="mm",
                                  gene.group = c("Mrp","chrM","Rb"), # Filter the Mt ribosomal genes and mitochondrial genes--excluding 
                                  # Remove lowly expressed genes (~10% of smallest groups)
                                  exp.cells=30) # filtering 5498 genes with low expression
gc()

# Some MT genes in our reference are not named in a canonical way. Find those genes. 
# The mitochondrial genome has been extracted previously from NC_001665.2
mt_gtf <- rtracklayer::import(here(ages_dir,"Reference","NC_001665.2_edit.gtf"))
mt_genes <- mt_gtf$gene_id
# keep only columns NOT in mt_genes
sc.dat.filtered <- sc.dat.filtered[,!colnames(sc.dat.filtered) %in% mt_genes]
rm(sc.dat)
gc()

# Construct prism object 
myPrism <- new.prism(
  reference=sc.dat.filtered, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1
)

# Run Prism
bp.res <- run.prism(prism = myPrism, n.cores=16)

# Extract results
theta <- get.fraction(bp = bp.res,
                      which.theta = "final",
                      state.or.type = "type")
theta.cv <- bp.res@posterior.theta_f@theta.cv

theta.df <- as.data.frame(theta) %>%
  rownames_to_column("sample_ID") %>%
  pivot_longer(cols = !sample_ID, names_to = "celltype", values_to = "ratio") %>%
  left_join(fc.meta)

# Save the result
dir_create(here(analysis_outdir, "BayesPrism"))
write_csv(theta.df, file = here(analysis_outdir,"BayesPrism","theta_df.csv"))
saveRDS(bp.res, file = here(analysis_outdir,"BayesPrism","BPres.rds"))

theta.df <- read_csv(here(analysis_outdir,"BayesPrism","theta_df.csv"))
bp.res <- readRDS(here(analysis_outdir,"BayesPrism","BPres.rds"))
 
# Visualization
x_order <- theta.df %>%
  filter(!age_months == 6) %>%
  arrange(age_months) %>%
  pull(sample_ID)

theta.df %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = sample_ID, y = ratio*100, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d()+
  #scale_x_discrete(limits = x_order) +
  #facet_wrap(vars( age_months), scale = "free_y")+
  ylab("Percentage (%)") +
  theme_minimal() + 
  theme(axis.text.x = 
          element_text(angle = 90,hjust = 0,vjust = 0.5))

theta.df %>%
  filter(!age_months == 6,
         celltype %in% c("PEC")) %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = sample_ID, y = ratio*100, fill = age_months)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d()+
  scale_x_discrete(limits = x_order) +
  #facet_wrap(vars( age_months), scale = "free_y")+
  ylab("Percentage (%)") +
  theme_minimal() + 
  theme(axis.text.x = 
          element_text(angle = 90,hjust = 0,vjust = 0.5))

# try to see correlation between Dclk1 expression and cell type composition
bk.genes.keep <- colnames(bk.dat) %in% colnames(sc.dat.filtered)
bk.dat.filt <- bk.dat[,bk.genes.keep]
bk.dds <- DESeqDataSetFromMatrix(countData = t(bk.dat.filt),
                                 colData = fc.meta,
                                 design = ~ age_months)
vsd <- vst(bk.dds, blind = FALSE)

# Save vst data
saveRDS(vsd, file = here(analysis_outdir,"vst_matrix.rds"))
vsd <- readRDS(file = here(analysis_outdir,"vst_matrix.rds"))

Dclk1.exp.df <- assay(vsd)["Dclk1",] %>%
  as.data.frame() %>%
  dplyr::rename("expression" = ".") %>%
  tibble::rownames_to_column("sample_ID")

theta.df %>%
  left_join(Dclk1.exp.df) %>%
  filter(celltype == "C-TAL",
         !age_months == 28) %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = expression, y =  log10(ratio), color = age_months)) +
  geom_point()

#### 2. MuSiC ####
# Follow this tutorial: https://xuranw.github.io/MuSiC/articles/MuSiC.html#estimation-of-cell-type-proportions-with-pre-grouping-of-cell-types
# This part: "Estimation of cell type proportions with pre-grouping of cell types"
# where they showed a mouse kidney sample deconvolution. 

# Preparation of single-cell data
sce <- as.SingleCellExperiment(srat)
sce <- swapAltExp(sce, "RNA")
# Similar to above, we filter out the TL1 population. Also, we filter based on the genes used in BayesPrism
genes.use <- colnames(bp.res@reference.update@phi)
sce_filtered <- sce[genes.use, !str_detect( colData(sce)$celltype, "TL1")]

# Preparation of bulk data
varmeta <- data.frame(labelDescription = colnames(fc.meta), row.names = colnames(fc.meta))
rownames(fc.meta) <- fc.meta$sample_ID
fc <- fc[,fc.meta$sample_ID]
bulk.eset <- Biobase::ExpressionSet(assayData = data.matrix(fc),
                                    phenoData = new("AnnotatedDataFrame", 
                                                    data = fc.meta, 
                                                    varMetadata = varmeta))

# Produce the first step information
sce.basis <- music_basis(sce_filtered, clusters = "celltype", samples = "library_id")  

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
dir_create(here(plot_outdir, "MuSiC"))
pdf(here(plot_outdir, "MuSiC","Dendrogram.pdf"))
par(mfrow = c(1, 2))
d <- dist(t(log(sce.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(sce.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')
dev.off()

# Manually split cell types into clusters
clusters.type = list(
  C1 = "IMM", C2 = "POD",
  C3 = c("EC","FIB"), C4 = "VSMC/P",
  C5 = c("IC-B","IC-A","C-TAL","M-TAL","DCT"),
  C6 = c("PEC","PT"),
  C7 = c("PT-injured","TL2"), # TL2 is similar to PT-injured in terms of injury marker genes
  C8 = c("CNT-PC","PC","CNT")
)
cl.type = as.character(sce_filtered$celltype)
for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
sce_filtered$clustertype <- factor(cl.type, levels = c(names(clusters.type)))

# Gather the markers for each cluster group. For each cluster of celltypes, we extract marker genes with p_val_adj < 0.05
sheets <- readxl::excel_sheets(here(gex_aggr_prep,"celltype_markers.xlsx"))
markers.celltypes <- lapply(sheets, function(sh){
  readxl::read_xlsx(here(gex_aggr_prep,"celltype_markers.xlsx"), sheet = sh)
}) 
names(markers.celltypes) <- sheets
names(markers.celltypes) <- gsub("^VSMC_P$", "VSMC/P", names(markers.celltypes))

cluster.markers <- lapply(names(clusters.type), function(clus){
  celltypes.within <- clusters.type[[clus]]
  markers.within <- lapply(celltypes.within, function(ct){
    cl.markers <- markers.celltypes[[ct]] %>%
      dplyr::filter(p_val_adj < 0.1)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::pull(gene) %>%
    unique()
})
names(cluster.markers) <- names(clusters.type)

# Another way is to find markers again using the seurat function
# all.equal(colnames(srat.filt),colnames(sce_filtered)) # Make sure the cell barcodes are of the same order
srat.filt$clustertype <- sce_filtered$clustertype
Idents(srat.filt) <- "clustertype"
cluster.markers2.df <- FindAllMarkers(srat.filt, assay = "SCT", min.pct = 0.2)
cluster.markers2 <- lapply(unique(cluster.markers2.df$cluster), function(clus){
  clus.marker <- cluster.markers2.df %>% 
    dplyr::filter(cluster == clus, p_val_adj < 0.1) %>%
    pull(gene)
})
names(cluster.markers2) <- unique(cluster.markers2.df$cluster)
  
#Run MuSiC proportion with tree-guided procedure
est.bulk.tree <- music_prop.cluster(
  bulk.mtx = exprs(bulk.eset), sc.sce = sce_filtered,
  group.markers = cluster.markers2, clusters = "celltype", groups = "clustertype",
  samples = "library_id", clusters.type = clusters.type, iter.max = 10000, eps = 0.001
)
theta.music.tree <- est.bulk.tree$Est.prop.weighted.cluster %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_ID") %>%
  pivot_longer(cols = -sample_ID, values_to = "ratio", names_to = "celltype") %>%
  left_join(fc.meta)

# Run MuSiC with normal mode (not sensitive enough. Go with tree mode)
est.bulk.normal <- music_prop(bulk.mtx = exprs(bulk.eset),
                        sc.sce = sce_filtered, clusters = 'celltype',
                        samples = 'library_id')

theta.music.normal <- est.bulk.normal$Est.prop.weighted %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_ID") %>%
  pivot_longer(cols = -sample_ID, values_to = "ratio", names_to = "celltype") %>%
  left_join(fc.meta)

# Visualization
x_order <- theta.music.tree %>%
  #filter(!age_months == 6) %>% # Actually if we include 6-month samples, we can find Dclk1 expression is still correlated with PT-injured population, suggesting that the 6-month samples might have some inflammation going on.
  arrange(age_months) %>%
  pull(sample_ID)

theta.music.tree %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = sample_ID, y = ratio*100, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d()+
  scale_x_discrete(limits = x_order) +
  #facet_wrap(vars( age_months), scale = "free_y")+
  ylab("Percentage (%)") +
  theme_minimal() + 
  theme(axis.text.x = 
          element_text(angle = 90,hjust = 0,vjust = 0.5))

theta.music.normal %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = sample_ID, y = ratio*100, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d()+
  scale_x_discrete(limits = x_order) +
  #facet_wrap(vars( age_months), scale = "free_y")+
  ylab("Percentage (%)") +
  theme_minimal() + 
  theme(axis.text.x = 
          element_text(angle = 90,hjust = 0,vjust = 0.5))

theta.music.tree %>%
  filter(
         celltype %in% c("PT")) %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = sample_ID, y = ratio*100, fill = age_months)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d()+
  scale_x_discrete(limits = x_order) +
  #facet_wrap(vars( age_months), scale = "free_y")+
  ylab("Percentage (%)") +
  theme_minimal() + 
  theme(axis.text.x = 
          element_text(angle = 90,hjust = 0,vjust = 0.5))

theta.music.tree %>%
  left_join(Dclk1.exp.df) %>%
  filter(celltype == "TL2") %>%
  mutate(age_months = factor(age_months)) %>%
  ggplot(aes(x = expression, y = ratio*100, color = age_months)) +
  geom_point()

# Save the data
dir_create(here(analysis_outdir,"MuSiC"))
saveRDS(est.bulk.tree, file = here(analysis_outdir,"MuSiC","estimation_tree.rds"))
write_csv(theta.music.tree, file = here(analysis_outdir,"MuSiC","theta_music.csv"))
