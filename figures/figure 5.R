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
library(rstatix)
library(here)
library(fs)
library(ggpubr)

#### Directories ####
project_dir <- here("garyw", "RAGE24")
plot_dir <- here(project_dir,"publication/figures/figure5")
data_dir <- here(project_dir, "qPCR")

#### Read merged qPCR data ####
df_merged <- read_csv(here(data_dir,"gene_expression_merged_publication.csv")) %>%
  select(-plasmid, -ctrl_plasmid)

#### !!!Notes!!! ####
# For simplicity in experiments, we used L and L+S to denote "DCX" and "KD" domains, relatively. 
# This is because most literature focus on long (740aa) and short (433aa) isoforms and they share kinase domain but differ by DCX domain.

#### Function to plot gene relative expression +/- SEM with stats comparison on dCt values using unpaired 
make_expr_plot <- function(df, targets = NULL, comparisons = NULL, bar.alpha = 0.6,
                           nrow = 1, outlier_remove = TRUE) {
  
  # Optional target filtering
  if (!is.null(targets)) {
    df <- df %>%
      filter(target %in% targets)
    df$target <- factor(df$target, levels = targets)
  }
  
  # Report outliers
  if ("outlier_any_gene_in_sample" %in% colnames(df) &&
      any(df$outlier_any_gene_in_sample, na.rm = TRUE)) {
    logger::log_info("Outlier detected:")
    outlier_samples <- df %>%
      filter(outlier_any_gene_in_sample == TRUE) %>%
      pull(sample) %>%
      unique()
    print(outlier_samples)
  }
  
  # Remove outliers for plotting/statistics
  if(outlier_remove == TRUE){
    df_plot <- df %>%
      filter(outlier_any_gene_in_sample == FALSE)
  }else{
    df_plot <- df
  }
  
  # Summary stats for plotted bars/error bars
  stats_df <- df_plot %>%
    group_by(group, target) %>%
    summarize(
      mean_rel_exp = mean(rel_exp, na.rm = TRUE),
      sem_rel_exp  = sd(rel_exp, na.rm = TRUE) / sqrt(sum(is.finite(rel_exp))),
      .groups = "drop"
    )
  
  # Statistical tests on d_ct-->d_ct is nearly normal-distribution. Do not use rel_exp because it is right-skewed.
  stat_test <- NULL
  if (!is.null(comparisons)) {
    
    # comparisons should look like:
    # list(c("Vehicle", "DIN"), c("Vehicle", "PMA"))
    
    stat_test <- df_plot %>%
      group_by(target) %>%
      group_modify(~{
        map_dfr(comparisons, function(comp) {
          .x %>%
            filter(group %in% comp) %>%
            t_test(d_ct ~ group, paired = FALSE, alternative = "two.sided") %>%
            mutate(group1 = comp[1], group2 = comp[2])
        })
      }) %>%
      ungroup()
    
    # add p-value labels and bracket positions
    stat_test <- stat_test %>%
      group_by(target) %>%
      mutate(
        y.position = {
          ymax <- max(df_plot$rel_exp[df_plot$target == first(target)], na.rm = TRUE)
          ymax + seq_len(n()) * 0.12 * ymax
        }
      ) %>%
      ungroup() %>%
      mutate(
        p.label = case_when(
          p < 0.0001 ~ "P < 0.0001",
          TRUE ~ paste0("P = ", signif(p, 2))
        )
      )
  }
  
  # Plot
  p <- ggplot() +
    geom_col(
      data = stats_df,
      aes(x = group, y = mean_rel_exp, fill = group),
      width = 0.7,
      alpha = bar.alpha,
      color = NA
    ) +
    geom_errorbar(
      data = stats_df,
      aes(
        x = group,
        ymin = mean_rel_exp - sem_rel_exp,
        ymax = mean_rel_exp + sem_rel_exp
      ),
      width = 0.3,
      linewidth = 0.5
    ) +
    geom_jitter(
      data = df_plot,
      aes(x = group, y = rel_exp),
      size = 1.3,
      width = 0.08,
      alpha = 0.75
    ) +
    facet_wrap(~target, scales = "free_y", nrow = nrow) +
    labs(
      x = NULL,
      y = "Relative expression"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Add statistical annotations
  if (!is.null(stat_test)) {
    p <- p +
      stat_pvalue_manual(
        stat_test,
        label = "p.label",   # or use "p.signif" for stars
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        bracket.size = 0.4,
        tip.length = 0.01,
        size = 3
      )
  }
  
  return(p)
}

#### 1. HEK293T siRNA + DIN experiments ####
# Colors
group_cols <- c(
  "Vehicle + si-L" = "#31688E",
  "DIN + si-L"     = "#35B779",
  
  "Vehicle + si-L+S"  = "#181a19",
  "DIN + si-L+S"      = "#5DC863",
  
  "Vehicle + siNC"  = "#443A83",
  "DIN + siNC"      = "#782c82"
)

# Select experiment
df <- df_merged %>%
  filter(experiment == "2026_02_27 HEK siRNA DIN 02_20_26")
my_comparisons <- list(c("DIN + siNC", "DIN + si-L"), c("DIN + siNC", "DIN + si-L+S"), 
                       c("Vehicle + siNC", "Vehicle + si-L"), c("Vehicle + siNC", "Vehicle + si-L+S"))

# plots
p <- make_expr_plot(df, targets = c("DCLK1-L","DCLK1-L+S","TNFRSF12A","FOSL2"), 
    comparisons = my_comparisons, bar.alpha = 0.7)+
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols)
ggsave(here(plot_dir, "HEK293T_siRNA.pdf"), p, width = 8, height = 3.5)

#### 2. HEK DIN ####
group_cols <- c(
  "Vehicle" = "#1E2B45",  # very dark mako navy
  "DIN"     = "#4E6F8C"   # blue-teal
)
# Experiment of DIN 20uM for 16 hours
df <- df_merged %>%
  filter(source_file == "gene expression_HEK DIN 20uM.csv")
my_comparisons <- list(c("DIN","Vehicle"))

# EMT markers and DCLK1 response to DIN 20uM 16hours
p <- make_expr_plot(df, targets = c("DCLK1-L", "DCLK1-L+S","TNFRSF12A","FOSL2", "TNFA"), comparisons = my_comparisons, bar.alpha = 0.7)+
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols)
ggsave(here(plot_dir, "HEK293T_DIN20uM.pdf"), p, width = 6, height = 3.5)

#### 3. RPTEC TNFa DIN ####
group_cols <-  c(
  "Vehicle"    = "#1E2B45",  # same baseline
  "DIN"        = "#4E6F8C",  # same DIN color
  "TNFa"       = "#3F8C8C",  # teal
  "DIN + TNFa" = "#67C1B8"   # brightest highlight
)
df <- df_merged %>% 
  filter(experiment == "2026_03_03 RPTEC PMA DIN TNFa 02_06 and 02_16",
         group %in% c("Vehicle", "TNFa", "DIN","DIN + TNFa")) 
my_comparisons <- list(c("DIN + TNFa","TNFa"), c("DIN + TNFa", "Vehicle"), c("DIN + TNFa", "DIN"))

# For marker genes
p <- make_expr_plot(df, targets = c("DCLK1-L", "DCLK1-L+S","TNFRSF12A", "TNFA"), comparisons = my_comparisons, bar.alpha = 0.7)+
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols)
ggsave(here(plot_dir, "RPTEC_TNFa DIN.pdf"), p, width = 6, height = 3.5)

#### 4. HEK PMA/Ionomycin ####
# colors
group_cols <- c(
  "Vehicle"           = "#1E2B45",  # unchanged mako navy
  "PMA+Ionomycin"     = "#B7AD4E",  # cividis olive
  "DIN+PMA+Ionomycin" = "#D7C94B"   # brighter cividis yellow
)
# Had two qPCR runs for this experiment, but the samples were the same.
df1 <- df_merged %>%
  filter(source_file == "gene expression selected HEK PMA Ionomycin stim 02_06_2026 exp.csv")
df2 <- df_merged %>%
  filter(source_file == "gene expression_HEK PMA_Io stimulation 02_06_26.csv") 

# In those experiments, two samples were 
df <- rbind(df1,df2) %>%
  mutate(group = str_replace(group, 
                             pattern = fixed("Vehicle for PMA+Ionomycin"), 
                             replacement = "Vehicle")) %>%
  mutate(group = str_replace(group, 
                             pattern = fixed("vehicle"), 
                             replacement = "Vehicle"))
my_comparisons <- list(c("PMA+Ionomycin","Vehicle"))

# Plot for main figure
p <- make_expr_plot(df, targets = c("DCLK1-L", "DCLK1-L+S","TNFRSF12A","FOSL2", "TNFA"), comparisons = my_comparisons, bar.alpha = 0.7)+
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols)
ggsave(here(plot_dir, "HEK293T_PMA.pdf"), p, width = 6, height = 3.8)

#### 5. RPTEC PMA Ionomycin ####
# colors
group_cols <- c(
  "Vehicle"           = "#1E2B45",  # unchanged mako navy
  "PMA + Ionomycin"     = "#B7AD4E",  # cividis olive
  "DIN + PMA + Ionomycin" = "#ecdd52"   # brighter cividis yellow
)

# select experiment
df <- df_merged %>% 
  filter(experiment == "2026_03_03 RPTEC PMA DIN TNFa 02_06 and 02_16",
         group %in% c("Vehicle", "DIN + PMA + Io", "PMA + Io")) %>%
  mutate(
    group = str_replace(group, pattern = "Io", replacement = "Ionomycin")
  )
my_comparisons <- list(c("Vehicle","PMA + Ionomycin"), c("PMA + Ionomycin", "DIN + PMA + Ionomycin"))

# For marker genes
p <- make_expr_plot(df, targets = c("DCLK1-L", "DCLK1-L+S","TNFRSF12A", "TNFA"), comparisons = my_comparisons, bar.alpha = 0.7)+
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols)
ggsave(here(plot_dir, "RPTEC_PMA DIN.pdf"), p, width = 6, height = 4)

