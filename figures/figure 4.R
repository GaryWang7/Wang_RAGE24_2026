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

#### 2. Plot human atlas UMAP ####
# color palette
celltype2_cols <- c(
  # PT comparison (highlighted but calm)
  "PCT"        = "#8FB6A2",  # light sage (PT-S1–like)
  "PST"        = "#5F8573",  # deeper sage (PT-S3–like)
  
  # Injury (highlighted away from others)
  "PT-injured" = "#C96A4A",  # rust highlight
  
  # Thin limb (unhighlighted)
  "TL"         = "#AFAFAF",
  
  # Distal nephron / TAL
  "TAL"        = "#C4C4C4",
  "DCT1"       = "#B9B9B9",
  "DCT2_PC"    = "#A8A8A8",
  
  # Collecting duct / intercalated
  "ICA"        = "#9A9A9A",
  "ICB"        = "#949494",
  
  # Glomerular / vascular / stromal
  "PODO"       = "#8E8E8E",
  "PEC"        = "#888888",
  "ENDO"       = "#828282",
  "FIB_VSMC_MC"= "#7C7C7C",
  
  # Immune
  "MONO"       = "#666666",
  "TCELL"      = "#6A6A6A",
  "BCELL"      = "#707070"
)

# with labels--The un-changed UMAP is also used in supplemental figure 7c.
p <- DimPlot_scCustom(srat, reduction = "umap", group.by = "celltype2",
                      colors_use = celltype2_cols, figure_plot = TRUE,
                      size = 0.01, label = TRUE)
ggsave(here(plot_dir,"UMAP_PT_highlight_PT_segments_labels.png"), p, dpi = 2000, width = 7)

#### 3. Plot circulating biomarkers ####
# Load circulating biomarkers
circ.markers <- read_excel(here(data_dir,"Proteins_list_46.xlsx")) %>%
  dplyr::select(1:4) %>%
  janitor::remove_empty()

# Load DEGs of PT-injured vs PT
deg <- read_csv(here(data_dir,"deg","deg.PT_VCAM1_vs_PT.csv")) %>%
  rename("gene" = "...1")

df <- deg %>%
  filter(gene %in% circ.markers$NAME_1)


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
  labs(x = "avg_log2FC", y = "pct.1") +
  ggpubr::theme_pubr()
ggsave(here(plot_dir, "biomarkers.pdf"), plot = p, width = 5)

#### 4. Dotplot for PT and PT-injured markers ####
injury.marker <- c("ADAMTS1", "CREB5","DCDC2","HAVCR1","NFKB1","DCLK1") 
healthy.marker <- c("CUBN","HNF4A","LRP2","PAX8","MME")

# color palettes
pt_healthy_pal <- c(
  "low"  = "#E8F1EC",  # very light sage
  "high" = "#789C8A"   # core PT green
)
pt_injury_pal <- c(
  "low"  = "#F6E6DC",  # very light warm beige
  "high" = "#C96A4A"   # PT-injured rust
)

# Dotplots
Idents(srat) <- "celltype3"
p_injury_dot1 <- DotPlot_scCustom(srat, features = injury.marker,
                                  group.by = "celltype3", 
                                  idents = c("PT","PT-injured"),
                                  colors_use = pt_injury_pal,
                                  x_lab_rotate = TRUE,
                                  flip_axes = FALSE) +
  scale_size(range = c(2,8)) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical") +
  theme(legend.title.align = 0)
p_healthy_dot1 <- DotPlot_scCustom(srat, features = healthy.marker,
                                    group.by = "celltype3", 
                                    idents = c("PT","PT-injured"),
                                    colors_use = pt_healthy_pal,
                                   x_lab_rotate = TRUE,
                                   flip_axes = FALSE) +
  scale_size(range = c(2, 10)) +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical") +
  theme(legend.title.align = 0)
ggsave(here(plot_dir,"injury marker_dots.pdf"), plot = p_injury_dot1, height = 3, width = 4.5)
ggsave(here(plot_dir,"healthy marker_dots.pdf"), plot = p_healthy_dot1, height = 3, width = 4.5)

#### 5. RAGE Dclk1 qPCR data ####
# Directories
age_col <- c("16" = "#b8d8ba", "30" = "#d9dbbc", "56" = "#dbac95", "82" = "#555b6e")
qPCR_dir <- here(project_dir, "qPCR", "2025_02_26_RAGE Nfkb Panels")
qPCR_dat <- read_csv(here(qPCR_dir, "RAGE NFkB panel bulk qPCR.csv"))
plots_dir <- here(plot_dir,"RAGE_qPCR")
dir_create(plots_dir)

# Tidy up the qPCR data
df <- qPCR_dat %>%
  filter(Target %in% c("Tnfrsf1a","Il1r1","Dclk1_L+S_all")) %>%
  mutate(
    Target = case_when(
      Target == "Dclk1_L+S_all" ~ "Dclk1",
      .default = Target
    )
  )

# Function to plot qPCR data with linear model p-value
plot_qPCR <- function(gene){
  p <- filter(df, Target == gene)%>%
    mutate(age_wks = factor(age_wks)) %>%
    ggplot(aes(x = age_wks, y = rel.exp, color = age_wks, fill = age_wks))+
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(alpha = 0.9) +
    labs(x = "Age (weeks)", y = "Relative gene expression", title = gene)+
    scale_fill_manual(values = age_col) +
    scale_color_manual(values = age_col)+
    ggpubr::theme_pubr()
  mod <- lm(rel.exp ~ age_wks, 
            data= df %>% mutate(age_wks = as.numeric(age_wks)) %>% 
              filter(Target == gene))
  return(list("plot" = p,"mod" = mod))
}

p <- plot_qPCR("Dclk1")$plot
ggsave(here(plots_dir,"Dclk1_qPCR.pdf"), width = 3)

#### 6. Dclk1 expression in AGES aging rat bulk mRNA-seq data ####
library(DESeq2)

ages_dir <- here("v/AgeRat")
counts_dir <- here(ages_dir, "counts")
analysis_dat_dir <- here(ages_dir,"analysis_data")
plots_dir <- here(plot_dir, "ages_bulk") # plot_dir is defined above (figures/figure3)

dir_create(plots_dir)
# Color palette
Ages7_mo <- c(
  "6"  = "#DDE6E0",  # very light grey-green (will be filtered later)
  "9"  = "#B7C7C1",  # reference anchor (muted cool sage, distinct from PT greens)
  "12" = "#D3D2B8",  # pale khaki
  "18" = "#C8AFA0",  # muted sand/rose
  "21" = "#B98F7E",  # warm taupe (older)
  "24" = "#7E6E74",  # muted mauve-slate
  "27" = "#4F5563"   # deep slate (oldest)
)

## Read bulk-RNAseq data
# Read count data
fc <- readRDS(here(counts_dir, "featureCounts.rds"))

# Renamed counts
raw_counts <- data.table::fread(here(counts_dir, "counts.csv")) %>%
  tibble::column_to_rownames("V1")

# Metrics of assignment to gtf file
metrics <- data.table::fread(here(counts_dir, "featureCount_metrics.csv"))

# Metadata
meta <- data.table::fread(here(ages_dir, "Info", "metadata.csv")) %>%
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

## Perform differential expression analysis
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

## Filter the data (remove age = 6 months)
# Based on the unfiltered classical PCA data, it seems that the 6 month
# old rats are not very similar to each other.
# Also, one rat from age group 9 is having an outlier PC1 value--SRR8705717 (also having lower mapped reads)
# The number of counts for SRR8705717 is much lower than other members in the group
samples_to_drop <- meta %>%
  tibble::rownames_to_column("name") %>%
  filter(name == "SRR8705717"| age_months == 6) %>%
  pull(name)
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
vsd <- vst(dds_filt, blind = FALSE)

## Dclk1 expression in AGES
meta_filt <- meta_filt %>%
  tibble::rownames_to_column("sample")

vst_plot <- function(gene){
  df <- assay(vsd)[gene,] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")%>%
    dplyr::rename("counts" = ".") %>%
    filter(sample %in% meta_filt$sample) %>%
    left_join(meta_filt)
  
  p <- df %>%
    ggplot(aes(x = age_months, y = counts, color = age_months, fill = age_months))+
    geom_point(alpha = 0.9)+
    geom_boxplot(alpha = 0.7)+
    xlab("Age (months)")+
    ylab("VST Normalized counts")+
    ggtitle(paste0("AGES: ", gene))+
    ggpubr::theme_pubr() +
    scale_color_manual(values = Ages7_mo) +
    scale_fill_manual(values = Ages7_mo)
  mod <- lm(counts ~ age_months, 
            data= df |> mutate(age_months = as.numeric(age_months)))
  return(list("plot" = p,"mod" = mod))
}

p <- vst_plot("Dclk1")$plot
ggsave(here(plots_dir,"AGES_Dclk1.pdf"), plot = p, width = 3.7)

#### 7. Deconvolution of AGES bulk mRNAseq data ####
plots_dir <- here(plot_dir, "deconv") # plot_dir is defined above (figures/figure3)
dir_create(plots_dir)

# Use the filtered samples (no 6month and sample "SRR8705717")
deconv <- read_csv(here(ages_dir,"analysis_data","deconv","MuSiC","theta_music.csv"))

# Function to plot deconvolution results with linear model p-value
plot_deconv <- function(.celltype){
  deconv2 <- filter(deconv, 
                    sample_ID %in% meta_filt$sample,
                    celltype == .celltype) %>%
    mutate(celltype = ifelse(.celltype=="PT-injured1", yes = "PT-injured", no = .celltype))
  
  # Plot box plot with lm correlation
  mod <- lm(ratio*100 ~ age_months, data = deconv2)
  pval <- summary(mod)$coefficients["age_months","Pr(>|t|)"] %>% 
    signif(3)
  adj.R2 <- summary(mod)$adj.r.squared %>% signif(3)
  p <- deconv2 %>%
    mutate(age_months = factor(age_months)) %>%
    ggplot(aes(x = age_months, y = ratio * 100)) +
    geom_boxplot(aes(fill = age_months), alpha = 0.7) +
    geom_point(aes(color = age_months), alpha = 0.9)+
    labs(
      title = .celltype,
      subtitle = bquote(
        "formula: " ~ y ~ "~" ~ x * ", " ~
          p == .(pval) * ", " ~
          {R^2}[adj] == .(adj.R2)
      )
    )+
    scale_color_manual(values = Ages7_mo)+
    scale_fill_manual(values = Ages7_mo)+
    xlab("Age (months)") +
    ylab("Percentage (%)") +
    ggpubr::theme_pubr()
  return(p)
}

p <- plot_deconv("PT-injured")
ggsave(here(plots_dir, "PT-injured deconv.pdf"),plot = p, width = 4, height = 6)

#### 8. Dclk1 expression vs PT-injured deconvolution in AGES data ####
Dclk1.exp.df <- assay(vsd)[c("Dclk1"),] %>%
  as.data.frame() %>%
  dplyr::rename("expression" = ".") %>%
  tibble::rownames_to_column("sample_ID")

# Function to plot Dclk1 expression vs cell type proportion
plot_ctprop_Dclk1 <- function(ct){
  theta <- deconv %>%
    filter(celltype %in% ct,
           sample_ID %in% meta_filt$sample) %>%
    left_join(Dclk1.exp.df) %>%
    mutate(age_months = factor(age_months))
  p <- theta %>%
    ggplot(aes(x = expression, y = ratio * 100)) +
    geom_smooth(method = "lm",
                formula = y ~ x,
                se = TRUE,
                color = "grey50",
                linewidth = 0.9,
                alpha = 0.15)+
    geom_point(aes(color = age_months), alpha = 0.9)+
    ggpmisc::stat_poly_eq(formula = y ~ x, ggpmisc::use_label("eq","p.value.label","adj.R2"))+
    scale_color_manual(values = Ages7_mo) +
    xlab("Normalized Dclk1 expression") +
    ylab("Percentage (%)") +
    labs(title=paste0("Percentage of ",ct),
         x = "Normalized Dclk1 expression",
         y = "Percentage (%)") +
    ggpubr::theme_pubr()
  return(p)
}

p <- plot_ctprop_Dclk1("PT-injured")
ggsave(here(plots_dir, "PT-injured with Dclk1.pdf"), width = 4, height = 5.2, plot = p)