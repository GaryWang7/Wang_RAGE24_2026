# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --workdir $HOME \
# --name chasm \
# -v /mnt/g:$HOME/g \
# -v /mnt/e:$HOME/e \
# -v /mnt/y:$HOME/y \
# -v /mnt/d/garyw:$HOME/garyw \
# -v $HOME:$HOME \
# -v /mnt/g/:$RHOME/g \
# -v /mnt/e:$RHOME/e \
# -v /mnt/y:$RHOME/y \
# -v /mnt/d/garyw:$RHOME/garyw \
# -v /var/run/docker.sock:/var/run/docker.sock \
# -e PASSWORD=garywang \
# -e DISABLE_AUTH=TRUE \
# -e AUTH_TIMEOUT_MINUTES=0 \
# -e SESSION_TIMEOUT_MINUTES=0 \
# -p 8787:8787 \
# garywang7/chasm:2.0.4

library(here)
library(tidyverse)
library(data.table)
library(stringr)
library(future.apply)
library(fs)
library(arrow)
library(logger)

options(future.globals.maxSize= 100*1024^3) # 100G

#### 1. Define directories and parameters ####
qt_lower <- 0
qt_upper <- 100
bin_size_Mb <- 1

project_dir <- here("garyw", "RAGE24")
data_dir <- here(project_dir, "chasm2", "result")
aggcsv <- read.csv(here(project_dir, "chasm2","resources","chasm_samples.csv"), header = T)

output_dir <- here(project_dir,"chasm2","result_aggr",paste0(bin_size_Mb,"_Mb_bins"))
plot_dir <- here(output_dir, "plots")
dir_create(here(plot_dir), recurse = TRUE)

#### 2. Load data ####

# annotation for all samples
anno <- lapply(seq_len(nrow(aggcsv)), function(i){
  samp.id <- aggcsv$sample_id[i]
  lib.id <- aggcsv$library_id[i]
  data.dir <- here(project_dir, "chasm2", "result", samp.id, lib.id)
  anno.file <- dir_ls(data.dir, regex = paste0("*barcode_anno.*",bin_size_Mb,"_Mb_bins_read_depth_filtered\\.csv$"), type = "file")
  if(length(anno.file)==0){
    stop("No annotation files matching the pattern 'barcode_anno_*_read_depth_filtered.csv' found in ", lib.id)
  }
  if(length(anno.file)>1){
    stop("Multiple annotation files matching the pattern 'barcode_anno_*.csv' found in ", lib.id)
  }
  dat <- fread(anno.file[1])
  return(dat)
}) %>% rbindlist()

# save combined annotation file
write_csv(anno, here(output_dir, "cell_annotation_Feb_02_2026.csv"))

# Bin depth summary
res.bins <- lapply(seq_len(nrow(aggcsv)), function(i){
  samp.id <- aggcsv$sample_id[i]
  lib.id <- aggcsv$library_id[i]
    
  data.dir <- here(project_dir, "chasm2", "result", samp.id, lib.id)
  # Find the newest result file
  res.file <- dir_ls(data.dir, regexp = paste0("bin_depth_per_cell_", bin_size_Mb, "_Mb_bins.csv")) %>%
    tibble(file = .) %>%
    mutate(mtime = file_info(file)$modification_time) %>%
    slice_max(mtime) %>%
    pull(file)
  
  if(length(res.file)==0){
    stop("No annotation files matching the pattern 'bin_depth_per_cell_*.csv' found in ", lib.id)
  }
  if(length(res.file)>1){
    stop("Multiple annotation files matching the pattern 'bin_depth_per_cell_*.csv' found in ", lib.id)
  }
  dat <- fread(res.file[1])
  return(dat)
}) %>% rbindlist()

# save bin read depth summary
fwrite(res.bins, here(output_dir, "bin_depth_per_cell_Feb_02_2026.csv"))
write_feather(res.bins, here(output_dir, "bin_depth_per_cell_Feb_02_2026.arrow"))

# negative binomial states summary
res.nb <- lapply(seq_len(nrow(aggcsv)), function(i){
  samp.id <- aggcsv$sample_id[i]
  lib.id <- aggcsv$library_id[i]
  
  data.dir <- here(project_dir, "chasm2", "result", samp.id, lib.id)
  # Find the newest result file
  res.file <- dir_ls(data.dir, regexp = paste0("chrom_neg_binom_read_depth_", bin_size_Mb, "_Mb_bins.csv")) %>%
    tibble(file = .) %>%
    mutate(mtime = file_info(file)$modification_time) %>%
    slice_max(mtime) %>%
    pull(file)
  if(length(res.file)==0){
    stop("No annotation files matching the pattern 'chrom_neg_binom_read_depth_(bin_size)_Mb_bins.csv' found in ", lib.id)
  }
  if(length(res.file)>1){
    stop("Multiple annotation files matching the pattern 'chrom_neg_binom_read_depth_(bin_size)_Mb_bins.csv' found in ", lib.id)
  }
  dat <- fread(res.file[1])
}) %>% rbindlist()

# Tidy up chromosome results
res.nb <- res.nb %>%
    dplyr::select(barcode_atac_aggr, barcode_gex, chromosome, 
                frag_lib_size, chrom_read_depth) %>%
    left_join(anno[,c("barcode_atac_aggr","library_id_gex")] %>%
    mutate(
    # Normalize by library depth
    chrom_read_depth_lib_norm = chrom_read_depth/frag_lib_size,
    # Log transform the data
    chrom_log_norm = log(chrom_read_depth_lib_norm + 1)
  )

#### 3. Find the mode (expected read depth) for each chromosome ####
# Then we find mode. We filter out the cells with too low or too high library sizes to identify the mode.
frag_libs <- res.nb %>%
  distinct(library_id_gex, frag_lib_size, barcode_atac_aggr) %>%
  group_by(library_id_gex) %>%
  summarise(
    # Use 50% to 90% cell counts to find mode
    lower_frag_lib_thre = quantile(frag_lib_size, prob = 0.5),
    upper_frag_lib_thre = quantile(frag_lib_size, prob = 0.9)
  )
mode_df <- res.nb %>%
left_join(frag_libs) %>%
filter(frag_lib_size >= lower_frag_lib_thre,
         frag_lib_size <= upper_frag_lib_thre) %>%
  group_by(library_id_gex, chromosome) %>%
  mutate(
    # Find the mode of log and read depth normalized chromosome read 
    chrom_log_mode = {
      d <- density(chrom_log_norm)
      d$x[which.max(d$y)]
    }
  ) %>% 
  ungroup() %>%
  distinct(library_id_gex,chromosome, chrom_log_mode)

# Use mode to calculate expected chromosome read depth.
res.nb2 <- res.nb %>%
  left_join(frag_libs) %>%
  left_join(mode_df) %>%
  mutate(
    # Use mode to calculate expected chromosome reads
    expected_chrom_read_depth_mode = frag_lib_size*(exp(chrom_log_mode) - 1)
  )

#### 4. Negative binomial modeling ####
# For each chromosome, assign the cells into different read depth groups so that we can calculate an empirical variance and size factor for negative binomial model.
ncell_df <- res.nb2 %>%
  distinct(library_id_gex, barcode_gex) %>%
  count(library_id_gex, name = "ncell") 

chrom_group <- res.nb2 %>%
  dplyr::select(library_id_gex, barcode_gex,frag_lib_size, chromosome, chrom_read_depth,chrom_log_norm, expected_chrom_read_depth_mode) %>%
  left_join(ncell_df) %>%
  group_by(library_id_gex, chromosome) %>%
  # calculate read depth group based on expected read depth. Each group must have at least 300 cells
  mutate(
    chrom_read_depth_group_mode = 
           paste0(library_id_gex,":",chromosome,":",
                  Hmisc::cut2(expected_chrom_read_depth_mode, m=max(300, ceiling(0.05*ncell))))) %>%
  ungroup()

chrom_phi <- chrom_group %>%
  # Within each read depth group, set a cutoff to be 3% at both head and tail
  group_by(chrom_read_depth_group_mode) %>%
  filter(chrom_read_depth > quantile(chrom_read_depth, 0.03) &
           chrom_read_depth < quantile(chrom_read_depth, 1 - 0.03)) %>%
  mutate(var_read_depth = var(chrom_read_depth),
            mean_read_depth = mean(chrom_read_depth),
         ncell_group = n()) %>%
  mutate(
    phi_mode = mean_read_depth^2/(var_read_depth - mean_read_depth)) %>%
  dplyr::distinct(chrom_read_depth_group_mode, phi_mode)%>%
  ungroup()

chrom_mean_var <- chrom_group %>%
  left_join(chrom_phi) %>%
  ungroup()

# Plot size factor phi vs chromosomes
chr_order <- gtools::mixedsort(unique(chrom_mean_var$chromosome))
p <- chrom_mean_var %>%
  distinct(library_id_gex, chromosome, chrom_read_depth_group_mode, phi_mode) %>%
  mutate(chromosome = factor(chromosome, levels = chr_order)) %>%
  ggplot(aes(x = chromosome, y = phi_mode))+
  geom_jitter(width = 0.15, alpha = 0.5, size = 0.3)+
  geom_boxplot(alpha = 0.4, outliers = FALSE) +
  scale_y_log10() +
  ggtitle("Size factor Phi of different read depth groups for Negative Binomial",
          "Larger Phi will converge to Poisson")+
  ggpubr::theme_pubr()
ggsave(here(plot_dir, "Phi calculation based on MODE.png"), plot = p, width = 12, height = 6, dpi = 900)

# Save the mean, variance and phi for each chromosome read depth group for QC purposes.
write_csv(chrom_mean_var, here(output_dir, "chrom_variance_dispersion_mode_Feb_09_2026.csv"))

# Test NB
res.nb3 <- res.nb2 %>% 
  left_join(dplyr::select(chrom_mean_var,
                          barcode_gex, chromosome, chrom_read_depth_group_mode, phi_mode), by = c("barcode_gex", "chromosome")) %>%
  mutate(
    chrom_p_val_mode = pnbinom(q = chrom_read_depth, mu = expected_chrom_read_depth_mode, size = phi_mode),
    # adjust the tail probabilities
    chrom_p_val_mode = if_else(chrom_read_depth > expected_chrom_read_depth_mode, 1- chrom_p_val_mode, chrom_p_val_mode)
  ) %>%
  group_by(library_id_gex) %>%
  mutate(
    # FDR control
    chrom_p_val_adj_mode = p.adjust(chrom_p_val_mode, method = "BH"),
  ) %>%
  ungroup() %>%
  dplyr::select(
    barcode_atac_aggr, frag_lib_size, chromosome, chrom_read_depth,
    chrom_read_depth_lib_norm, chrom_log_norm, chrom_log_mode,
    expected_chrom_read_depth_mode, chrom_read_depth_group_mode, phi_mode,
    chrom_p_val_mode, chrom_p_val_adj_mode
  )

# Save to file
write_csv(res.nb3, file = here(output_dir,"negative_binomial_combined_Feb_09_2026.csv"))


