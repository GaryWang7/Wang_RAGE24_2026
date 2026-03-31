# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --restart unless-stopped \
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
# garywang7/ragemultiome:1.0.3

library(readxl)
library(here)
library(Seurat)
library(tidyverse)
library(fs)
library(scCustomize)
library(patchwork)

#### 1. Read data ####
# Directories
data_dir <- here("garyw/RAGE24/human_atlas")
project_dir <- here("garyw", "RAGE24")
plot_dir <- here(project_dir,"publication/figures/figure4")
dir_create(plot_dir)

# Read seurat
srat <- readRDS(here(data_dir,"srat.rds"))

# Read metadata
anno <- read_csv(here(data_dir,"filtered_annotations.csv")) %>%
  dplyr::rename(celltype_raw = celltype) %>%
  tibble::column_to_rownames("barcode")

srat <- AddMetaData(srat, anno)

# Update cell types (with PCT and PST)
srat$celltype2 <- dplyr::case_when(
  as.character(srat$celltype) == "PT_VCAM1" ~ "PT-injured",
  as.character(srat$celltype) == "PT_PROM1" ~ "TL",
  TRUE ~ as.character(srat$celltype)
)
celltype2_levels <- c(
  "PCT","PST","PT-injured",
  "TL",
  "TAL","DCT1","DCT2_PC",
  "ICA","ICB",
  "PODO","PEC","ENDO","FIB_VSMC_MC",
  "MONO","TCELL","BCELL"
)
srat$celltype2 <- factor(srat$celltype2, levels = celltype2_levels)

srat$celltype3 <- dplyr::case_when(
  as.character(srat$celltype2) %in% c("PST","PCT") ~ "PT",
  .default = as.character(srat$celltype2)
)
celltype3_levels <- c(
  "PT","PT-injured",
  "TL",
  "TAL","DCT1","DCT2_PC",
  "ICA","ICB",
  "PODO","PEC","ENDO","FIB_VSMC_MC",
  "MONO","TCELL","BCELL"
)
srat$celltype3 <- factor(srat$celltype3, levels = celltype3_levels)

Idents(srat) <- "celltype3"

##### Plot Dclk1 dotplot
# Plot DCLK1 in all celltypes--dotplot
p1 <- DotPlot_scCustom(srat, features = "DCLK1",
                      group.by = "celltype2",
                      x_lab_rotate = TRUE,
                      flip_axes = TRUE)

# Plot Dclk1 in the umap plot via feature plot
p2 <- FeaturePlot_scCustom(srat, reduction = "umap", features = "DCLK1", pt.size = 1,
                           figure_plot = TRUE, colors_use = viridis_plasma_light_high)
# p3 <- DimPlot_scCustom(srat, reduction = "umap", group.by = "celltype2",
#                        colors_use = celltype2_cols, figure_plot = TRUE,
#                        size = 0.01, label = TRUE, repel = TRUE)
ggsave(here(plot_dir, "DCLK1 dotplot.pdf"), plot = p1, height = 4, width = 8)
ggsave(here(plot_dir, "DCLK1 featureplot.png"), plot = p2, dpi = 2000, width = 6)