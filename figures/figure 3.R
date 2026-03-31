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

project_dir <- here("garyw", "RAGE24")
dPCR_dir <- here(project_dir, "dPCR", "summarized")

plot_dir <- here(project_dir,"publication/figures/figure3")
dir_create(plot_dir, "dPCR")


## Read dPCR files
# Information:
# dPCR probes--GRCr8 genome:
# chr1: 243,220,858 - 243,220,896 (240Mb-245Mb)
# chrX: 72,655,947 – 72,655,993
# chrY: 947,723 – 947,779
df2025 <- read_csv(here(dPCR_dir, "RAGE24_dPCR_publication.csv"))

#### dPCR chrX and chrY separate with age ####
# chrX copy number with age
p.chrX.chr1 <- df2025 %>%
  dplyr::mutate_at(c('age_weeks'), as.factor) %>%
  ggplot(aes(age_weeks, copy_number_chrX_vs_chr1)) +
  geom_boxplot(aes(fill=age_weeks),outliers = FALSE)  +
  geom_jitter(alpha = 0.7, position = position_jitter(width=0.1), color = chr_cols["chrX"])+
  #ylim(c(1.0,1.25)) +
  theme_light() +
  labs(x = "Age (weeks)", y = "Copy Number chrX", fill = "Age (weeks)")+
  theme(legend.position = "none")+
  # pairwise p-values
  ggpubr::stat_compare_means(comparisons = age_wk_comparison, method = "wilcox.test", 
                     paired = FALSE, p.adjust.methods = "fdr",
                     aes(label = after_stat(p.format)))+
  # global anova p-value
  ggpubr::stat_compare_means( label.y = 1.24,label.x = 3.0) +
  scale_fill_manual(values = age_col)
ggsave(filename = here(plot_dir, "dPCR","chrX copy number with age.pdf"),
       plot = p.chrX.chr1, width = 3, height = 4)

#chrY copy number with age
p.chrY.chr1 <- df2025 %>%
  dplyr::mutate_at(c('age_weeks'), as.factor) %>%
  ggplot(aes(age_weeks, copy_number_chrY_vs_chr1)) +
  geom_boxplot(aes(fill=age_weeks),outliers = FALSE)  +
  geom_jitter(alpha = 0.7, position = position_jitter(width=0.1), color = chr_cols["chrY"])+
  #ylim(c(1.0,1.25)) +
  theme_light() +
  labs(x = "Age (weeks)", y = "Copy Number chrY", fill = "Age (weeks)")+
  theme(legend.position = "none")+
  # pairwise p-values
  ggpubr::stat_compare_means(comparisons = age_wk_comparison, method = "wilcox.test", 
                             paired = FALSE, p.adjust.methods = "fdr",
                             aes(label = after_stat(p.format)))+
  # global anova p-value
  ggpubr::stat_compare_means( label.y = 1.24,label.x = 3.0) +
  scale_fill_manual(values = age_col)
ggsave(filename = here(plot_dir, "dPCR","chrY copy number with age.pdf"),
       plot = p.chrY.chr1, width = 3, height = 4)

#### chrX and chrY correlation ####
chrXY.cor <- df2025 %>%
  dplyr::mutate_at(c('age_weeks'), as.factor) %>%
  ggplot(aes(x = copy_number_chrX_vs_chr1, y = copy_number_chrY_vs_chr1))+
  # Add a grey dashed trend line based on all points (ignoring the group/color)
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  geom_point(size = 3, alpha = 1, aes(color = age_weeks)) +
  # Calculate and display the R² for all points together using ggpmisc
  ggpmisc::stat_poly_eq(
    formula = y ~ x, label.x = "left", label.y = "top", parse = TRUE,
    ggpmisc::use_label("Eq", "adj.R2",  "P", sep = "*\" ; \"*")
  ) +
  labs(
       x = "chrX Copy Number",
       y = "chrY Copy Number",
       color = "Age (weeks)") +
  scale_color_manual(values = age_col)+
  ggpubr::theme_pubr() +
  theme(legend.position = "right")
ggsave(filename = here(plot_dir, "dPCR","chrX 2025 chrY 2024 copy number Association.pdf"),
       plot = chrXY.cor , width = 6, height = 4)
