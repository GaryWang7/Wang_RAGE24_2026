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

#### 6. Base HEK and RPTEC DCLK1 expression compared to housekeeping ####
# HEK cells
df1 <- df_merged %>%
  filter(source_file == "gene expression_HEK DIN 20uM.csv",
         target %in% c("DCLK1-L","DCLK1-L+S"),
         group == "Vehicle") %>%
  mutate(rel_exp_housekeeping = 2^(-d_ct),
         celltype = "HEK293T")

# RPTEC cells
df2 <- df_merged %>% 
  filter(experiment == "2026_03_03 RPTEC PMA DIN TNFa 02_06 and 02_16",
         group == "Vehicle",
         target %in% c("DCLK1-L","DCLK1-L+S")) %>%
  mutate(rel_exp_housekeeping = 2^(-d_ct),
         celltype = "hTERT-RPTEC")
df2 <- filter(df2,
              sample %in% sample(unique(sample),4))
df <- rbind(df1,df2)

cellline_cols <- c(
  "HEK293T"     = "#2F7DE1",
  "hTERT-RPTEC" = "#00A087"
)

# Summary stats for plotted bars/error bars
stats_df <- df %>%
  group_by(celltype, target) %>%
  summarize(
    mean_rel_exp = mean(rel_exp_housekeeping, na.rm = TRUE),
    sem_rel_exp  = sd(rel_exp_housekeeping, na.rm = TRUE) / sqrt(sum(is.finite(rel_exp_housekeeping))),
    .groups = "drop"
  )

# statistical test
comparisons <- list(c("HEK293T","hTERT-RPTEC"))
ymax_tbl <- df %>%
  group_by(target) %>%
  summarise(
    ymax = max(rel_exp_housekeeping, na.rm = TRUE),
    .groups = "drop"
  )
stat_test <- df %>%
  group_by(target) %>%
  group_modify(~{
    map_dfr(comparisons, function(comp) {
      .x %>%
        filter(celltype %in% comp) %>%
        t_test(d_ct ~ celltype, paired = FALSE, alternative = "two.sided") %>%
        mutate(group1 = comp[1], group2 = comp[2])
    })
  }) %>%
  ungroup() %>%
  left_join(ymax_tbl, by = "target") %>%
  # add p-value labels and bracket positions
  group_by(target) %>%
  mutate(
    y.position = ymax + seq_len(n()) * 0.12 * ymax
  ) %>%
  ungroup() %>%
  mutate(
    p.label = case_when(
      p < 0.001 ~ "P < 0.001",
      TRUE ~ paste0("P = ", signif(p, 2))
    )
  )

# Plot
p <- ggplot() +
  geom_col(
    data = stats_df,
    aes(x = celltype, y = mean_rel_exp, fill = celltype),
    width = 0.7,
    alpha = 0.7,
    color = NA
  ) +
  geom_errorbar(
    data = stats_df,
    aes(
      x = celltype,
      ymin = mean_rel_exp - sem_rel_exp,
      ymax = mean_rel_exp + sem_rel_exp
    ),
    width = 0.3,
    linewidth = 0.5
  ) +
  geom_jitter(
    data = df,
    aes(x = celltype, y = rel_exp_housekeeping),
    size = 1.3,
    width = 0.08,
    alpha = 0.75
  ) +
  facet_wrap(~target, scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = "Relative expression to PPIA"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
    stat_pvalue_manual(
      stat_test,
      label = "p.label",   # or use "p.signif" for stars
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      bracket.size = 0.4,
      tip.length = 0.01,
      size = 3
    ) +
# Add color
  scale_fill_manual(values = cellline_cols)+
  scale_color_manual(values = cellline_cols)

ggsave(here(plot_dir, "HEK vs RPTEC DCLK1.pdf"), plot = p, width = 4, height = 4)
