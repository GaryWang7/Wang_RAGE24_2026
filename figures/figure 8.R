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

library(tidyverse)
library(here)
library(fs)
library(DESeq2)

#### Directories ####
project_dir <- here("garyw", "RAGE24")
plot_dir <- here(project_dir,"publication/figures/figure8")
data_dir <- here("v/RIDE25")

#### metadata and parameters ####
meta <- read_csv(here(data_dir,"metadata.csv"))%>%
  mutate(sample = library_name)

treat_cols <- c(
  "VEH" = "#7A8BBF",  # softened steel-blue (background, calm)
  "DIN" = "#9B4DFF"   # vivid violet-purple (highlight)
)

# Use DIN up and down palettes
pathway_up_pal <-c(
  "#A9783A",  # deep copper (high)
  "#C59A63",
  "#E2C9A5",
  "#F7F1E8"   # very low
)

pathway_down_pal <- c(
  "#6E8570",
  "#8FA79A",
  "#C9DDD3",
  "#EEF4F1"
)

# ==============================================================================
# Body weight 
# ==============================================================================
library(lme4)
#install.packages("lmerTest")

weight.df <- meta %>%
  filter(dosage == 6) %>%
  dplyr::select(rat_id, treatment, contains("weight_dose")) %>%
  pivot_longer(cols = contains("weight_dose"), names_to = "time_points", values_to = "weight") %>%
  mutate(time = as.numeric(str_extract(time_points, "\\d+$")),
         treatment = factor(treatment),
         treatment = relevel(treatment, ref = "VEH")) 

# linear model
mod <- lmerTest::lmer(
  weight ~ time * treatment + (1 | rat_id),
  data = weight.df
)

# Extract pvalue for the interaction term
stats <- summary(mod)
p_int <- stats$coefficients["time:treatmentDIN", "Pr(>|t|)"]
p_lab <- paste0("p = ", signif(p_int, 3))

# Plot
p <- ggplot(weight.df, aes(x = time, y = weight, color = treatment)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  annotate(
    "text",
    x = 2.5, y = 350,
    label = p_lab,
    hjust = 0,
    size = 4,
    color = "black"
  ) +
  ggpubr::theme_pubr() +
  theme(legend.direction = "vertical",
 legend.position = "right") +
  scale_color_manual(values = treat_cols) +
  labs(x = "Dosage Time", y = "Weight (g)")

ggsave(here(plot_dir,"body weights.pdf"), width = 7, height = 3, plot = p)

# ==============================================================================
# DESeq2 and GSEA analysis 
# ==============================================================================
set.seed(1105)
msigDb <- msigdbr::msigdbr(species = "Rattus norvegicus")

# According to https://stephenturner.github.io/deseq-to-fgsea/ and https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
# "stat" from DESeq or Limma are good ranking metrics
# The default metric is signal-to-noise ratio, so we use "stat" here ,which is log2FC/lfcSE
library(fgsea)
library(BiocParallel)
param <- SnowParam(workers = parallel::detectCores() - 2)

#### Directories 
counts_dir <- here(data_dir, "counts")
analysis_dat_dir <- here(data_dir,"analysis_data")

gsea_dir <- here(plot_dir,"gsea")
dir_create(gsea_dir)

#### Run DESeq2
dds <- readRDS(here(analysis_dat_dir,"DESeq_obj.rds"))
vsd <- vst(dds, blind = FALSE)

# DESeq2 DE result
dds_sub <- dds[,dds$dosage == 6]
res<- DESeq2::results(dds_sub, contrast = c("treatment", "DIN", "VEH")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

# To generate stats for each gene, we have to remove duplicates
gsea_stat <- function(res){
  res <- res %>% dplyr::distinct(gene_id,stat)%>%
    group_by(gene_id) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%  # one row per gene_name
    ungroup() %>%
    arrange(stat)
}

gsea_stat <- gsea_stat(res)

# Function to run gsea
gsea_wrapper <- function(gsea_res, gset){
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
  fgseaRes <- fgsea(pathways = pathway_list,
                    stats = rank,
                    maxSize = 500,
                    nproc = parallel::detectCores() - 2,
                    BPPARAM = param)
  return(list(path_list = pathway_list, gsea = fgseaRes))
}

#### Hallmark GSEA results
gsea_Hallmark <- gsea_wrapper(gsea_stat, gset = "H")

res <- gsea_Hallmark$gsea %>%
  filter(padj < 0.1) %>%
  dplyr::arrange(desc(NES), padj) %>%
  mutate(change = if_else(NES > 0, true = "UP", false = "DOWN"))
res$pathway <- res$pathway%>%
  tolower() %>%
  str_remove_all("hallmark_") 

# Plot dotplot
temp_plotdir <- here(plot_dir,"gsea","Hallmark")
dir_create(temp_plotdir)

lapply(unique(res$change), function(reg){
  if(reg == "UP"){pal <- pathway_up_pal}else{pal<-pathway_down_pal}
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