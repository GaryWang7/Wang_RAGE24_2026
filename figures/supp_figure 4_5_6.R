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
# garywang7/ragemultiome:1.0.4

library(tidyverse)
library(here)
library(fs)
library(data.table)
library(GenomicRanges)
library(arrow)
library(patchwork)

future::plan("multicore", workers = 6)
options(future.globals.maxSize= 100*1024^3) # 100G

#### Define directories and parameters ####
bin_size_Mb <- 1

project_dir <- here("garyw", "RAGE24")
data_dir <- here(project_dir, "nb_model")
nb_dir <- here(data_dir,"result_aggr","negative_binomial")

output_dir <- here(project_dir,"nb_model","result_aggr",paste0(bin_size_Mb,"_Mb_bins"))
plot_dir <- here(project_dir,"publication/figures/figure2")
dir_create(here(plot_dir), recurse = TRUE)

age_col <- c("16" = "#b8d8ba", "30" = "#d9dbbc", "56" = "#dbac95", "82" = "#555b6e")

method_col <- c(
  "DEFND" = "#D81B60",  # vivid magenta (warm, not orange)
  "Multiome"    = "#1F4E79"   # deep navy (cold)
)

#### chrX and chrY scale log-normed read pairwise ####
# Load data from negative binomial output
res.nb <- read_csv(here(nb_dir, "negative_binomial_combined_copynumber_Feb_09_2026.csv"))
anno.nb <- read_csv(here(nb_dir, "cell_annotation_Feb_02_2026.csv"))

# chromosomes
chroms <- unique(res.nb$chromosome) %>% gtools::mixedsort()
res.nb$chromosome <- factor(res.nb$chromosome, levels = chroms)

# Function to plot pairwise dotplot
plot_pairwise <- function(data, chroms_to_plot, title, xlim = c(-0.25,0.25)){
  value_filter <- function(x){
    # Function to limit the plotting range to be -100 to 100
    return(abs(x) <= 100)
  }
  densityp <- function(data, mapping, ...){ 
    ggplot(data = data, mapping = mapping, ...) + 
      #geom_hex(..., alpha = 1, bins = 300) + 
      #viridis::scale_fill_viridis(option = "turbo")+
      ggrastr::rasterise(geom_point(alpha = 0.02, size = 0.01, aes(fill = library_type, color = library_type)), dpi = 1000)+
      scale_fill_manual(values = method_col)+
      scale_color_manual(values = method_col)+
      coord_cartesian(xlim = c(-3,5), ylim = c(-3,5))
  }
  histp <-function(data, mapping,...){
    ggplot(data = data, mapping = mapping, ...) + 
      geom_density(..., aes(y = after_stat(scaled),
                            color = library_type, fill=library_type),alpha = 0.6) + 
      scale_fill_manual(values = method_col)+
      scale_color_manual(values = method_col)+
      coord_cartesian(xlim = c(-4,4))
  }
  corrp <- function(data, mapping, ...){
    ggplot(data = data, mapping = mapping, ...)+
    ggpubr::stat_cor(
      ...,
      aes(color = library_type),
      method = "pearson")+
      scale_color_manual(values = method_col)
  }
  data %>%
    left_join(anno.nb[,c("library_type","barcode_atac_aggr")]) %>%
    group_by(chromosome, library_type) %>%
    mutate(
      scale_chrom_log = (chrom_log_norm - chrom_log_mode)/sd(chrom_log_norm)
    ) %>%
    ungroup() %>%
    dplyr::select(chromosome, scale_chrom_log, library_type, barcode_atac_aggr)%>%
    pivot_wider(names_from = chromosome, values_from = scale_chrom_log)%>%
    dplyr::select(any_of(chroms_to_plot), barcode_atac_aggr, library_type)%>%
    filter(if_all(where(is.numeric), value_filter))%>%
    GGally::ggpairs(
      columns = chroms_to_plot,
      lower = list(continuous = densityp),
      upper = list(continuous = corrp),
      diag = list(continuous = histp),
      title = title)+
    theme_light()+
    theme(
      text = element_text(size  = 15), 
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 18),
      legend.text = element_text(size=10),
      strip.text = element_text(size = 15),
      strip.background = element_rect(fill = "#7C7C7C"))
}

p <- plot_pairwise2(res.nb, chroms_to_plot = c("chrX","chrY"), 
                    title = "Scaled chromosome coverage in all cells")

ggsave(here(plot_dir,"pairwise_chrX_chrY.pdf"), height = 5.5, width = 7, dpi = 1500, plot = p)

#### DEFNDseq vs Multiome chromosome read coverage ####
# Use data loaded from previous section
# Plot 
p <- res.nb %>%
  left_join(anno.nb[,c("barcode_atac_aggr","library_type")]) %>%
  group_by(chromosome, library_type) %>%
  mutate(
    scale_chrom_log = (chrom_log_norm - chrom_log_mode)/sd(chrom_log_norm)
  ) %>%
  ggplot()+
  geom_density(aes(x = scale_chrom_log, y = after_stat(scaled), fill = library_type, color = library_type), 
               alpha = 0.5) +
  facet_wrap(vars(chromosome))+
  scale_fill_manual(values = method_col)+
  scale_color_manual(values = method_col)+
  coord_cartesian(xlim = c(-4,4))+
  ggpubr::theme_pubr()
ggsave(here(plot_dir,"chrom read normalized multiome vs defnd.pdf"), width = 8, height = 8, plot = p)

#### Proportion of PT cells with chrX and chrY gain ####
# Use data from previous sections
# calculate CNA in each chromosome
calc_CNA <- function(res, padj_thre = 1e-4, cna_scale = "absolute", cna_thre = 1){
  res.df <- res %>%
    ungroup() %>%
    mutate(
      # Readjust p value
      chrom_p_val_adj_mode = p.adjust(chrom_p_val_mode, method = "BH"),
      # Identify cna
      any_CNA = case_when(
        chrom_p_val_adj_mode < padj_thre ~ TRUE,
        .default = FALSE
      ),
      any_CNA_type = case_when(
        any_CNA ==TRUE & cna_size > 0 ~ "Gain",
        any_CNA ==TRUE & cna_size < 0 ~ "Loss",
        .default = "None"
      ),
      large_CNA = case_when(
        cna_scale == "absolute" & any_CNA ==TRUE & abs(cna_size) >= cna_thre ~ TRUE, 
        cna_scale == "relative" & any_CNA ==TRUE & abs(rel_cna_size) >= cna_thre ~ TRUE, 
        .default = FALSE
      ),
      large_CNA_type = case_when(
        large_CNA==TRUE & cna_size > 0 ~ "Gain",
        large_CNA==TRUE & cna_size < 0 ~ "Loss",
        .default = "None"
      )
    ) %>%
    dplyr::select(barcode_atac_aggr, contains("CNA"), chromosome, chrom_category,
                  chrom_p_val_adj_mode, frag_lib_size) %>%
    ungroup() %>%
    distinct() %>%
    left_join(anno.nb) %>%
    dplyr::select(barcode_atac_aggr, celltype, rat_id, age_wks,
                  contains("CNA"), contains("chrom"), library_type,
                  frag_lib_size)
}

library(broom)

# chrX and chrY
prop.df <- res.nb %>%
  left_join(anno.nb[,c("barcode_atac_aggr", "library_id_gex", "celltype")]) %>%
  filter(celltype %in% c("PT","PT-injured")) %>%
  ungroup() %>%
  calc_CNA(padj_thre = 0.01, cna_scale = "absolute", cna_thre = 0.5) %>% # Same threshold for calculating PT cell any/large CNAs
  filter(chromosome %in% c("chrX","chrY"),
          !rat_id %in% c("J10","E05","N14"), 
          library_type == "DEFND") %>%
   mutate(
     have_large_CNA = large_CNA==TRUE&large_CNA_type == "Gain",
     have_CNA = any_CNA == TRUE&any_CNA_type == "Gain",
     .groups = "drop"
   ) %>%
   group_by(age_wks, rat_id, chromosome) %>%
   summarise(
     n_cells = n(),
     n_large_CNA = sum(have_large_CNA, na.rm = TRUE),
     n_any_CNA = sum(have_CNA, na.rm = TRUE),
     large_CNA_prop = n_large_CNA / n_cells,
     any_CNA_prop = n_any_CNA / n_cells,
     .groups = "drop"
   )

# Make sure age is ordered correctly
age_levels <- sort(unique(prop.df$age_wks))

# Use linear regression to model sample-level trend of CNA
trend_df <- prop.df %>%
  mutate(
    age_wks_num = as.numeric(as.character(age_wks)),
    age_wks = factor(age_wks, levels = age_levels)
  ) %>%
  group_by(chromosome) %>%
  group_modify(~{
    fit <- lm(any_CNA_prop ~ age_wks_num, data = .x)
    
    broom::tidy(fit) %>%
      filter(term == "age_wks_num") %>%
      mutate(
        p_label = paste0("trend p = ", signif(p.value, 3)),
        x = tail(levels(.x$age_wks), 1),,
        y = max(.x$any_CNA_prop, na.rm = TRUE) * 1.08
      ) %>%
      select(p.value, estimate, p_label, x, y)
  }) %>%
  ungroup()

p <- prop.df %>%
  mutate(age_wks = factor(age_wks, levels = age_levels)) %>%
  ggplot() +
  geom_point(aes(x = age_wks, y = any_CNA_prop, color = age_wks))+
  geom_boxplot(aes(x = age_wks, y = any_CNA_prop, fill = age_wks), alpha = 0.4)+
  facet_wrap(vars(chromosome), scale = "free_y")+
  geom_text(
    data = trend_df,
    aes(x = x, y = y, label = p_label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 4
  ) +
  labs(x = "Age (weeks)", y = "Proportion Gain CNA", title = "PT and PT-injured cells") + 
  scale_color_manual(values = age_col)+
  scale_fill_manual(values = age_col)+
  ggpubr::theme_pubr()
ggsave(plot = p, filename = here(plot_dir, "CNA_prop","chrX chrY Gain Proportion DEFND.pdf"), width = 8, height = 6)
