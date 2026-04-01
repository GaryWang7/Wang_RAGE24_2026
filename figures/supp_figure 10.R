# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
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

library(here)
library(fs)
library(tidyverse)

#### Directories ####
proj.dir <- here("garyw/RAGE24")

bin_size_Mb <- 1
output_dir <- here(project_dir,"nb_model","result_aggr",paste0(bin_size_Mb,"_Mb_bins"))

plot.all.dir <- here(proj.dir,"publication/figures")
plot.dir <- here(plot.all.dir,"supplements")
dir_create(plot.dir)

#### Load data ####
# This data relies on the output of figure 2 (coverage)
anno <- fread(file = here(output_dir, "cell_annotation_Feb_02_2026.csv"), nThread = 10)
# Here we characterize lithium concentration 12.mM to be strong depletion, 5mM and 6 mM to be weak depletion
cov_stats2 <- read_csv(here(output_dir, "coverage_whole_genome.csv")) %>%
  left_join(anno) %>%
  mutate(lithium_level = case_when(
    Lithium_concentration_mM > 6 ~ "strong",
    Lithium_concentration_mM > 0 & Lithium_concentration_mM <= 6  ~ "weak",
    .default = "0"))

#### Plot ####
lithium_cols <- c(
  "strong" = "#E76F51",  # coral (warm, not magenta)
  "weak"   = "#2A9D8F"   # teal-cyan (cool, not navy)
)

p <- cov_stats2 %>%
  filter(library_type == "DEFND") %>%
  ggplot(aes(x = coverage.pct, color = lithium_level, fill = lithium_level, group = lithium_level))+
  geom_density(alpha = 0.1)+
  geom_histogram(aes(y = after_stat(density)), position = "identity",alpha = 0.3, bins= 70, 
                 linewidth = 0.3)+
  scale_x_continuous(labels = scales::percent, limits = c(0.40,1))+
  labs(
    title = paste0("Coverage of genome: ", bin_size_Mb, " Mb bins with reads"),
    x = "Coverage (%)"
  )+
  scale_fill_manual(values = lithium_cols) +
  scale_color_manual(values = lithium_cols) +
  ggpubr::theme_pubr()

ggsave(here(plot_dir, "coverage", "lithium strong and weak_coverage.pdf"), plot = p,
       height = 5, width = 6)

