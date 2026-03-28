# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --workdir $HOME \
# --name CNV \
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
# -p 8787:8787 \
# garywang7/chasm:2.0.4

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(wavethresh)
  library(ggplot2)
  library(here)
  library(openxlsx)
  library(argparse)
  library(logger)
  library(fs)
})


#### 1. Define arguments and parameters ####
# Create a parser
parser <- ArgumentParser(description = 'Step 2. filtering',
                         formatter_class = 'argparse.ArgumentDefaultsHelpFormatter')

# Add arguments
parser$add_argument('--library_id', type = 'character', 
                    default = 'RAGE24-C03-KYC-LN-01-DEFND-ATAC',
                    help = 'Library id in the format of "RAGE24-C03-KYC-LN-01-DEFND-ATAC"')
parser$add_argument('--out_dir', type = 'character', 
                    default = 'garyw/RAGE24/chasm2/result',
                    help = 'Output path for all chasm results. 
                    Results will be stored in out_dir/sample_id/library_id, where sample_id is extracted from library_id.')
parser$add_argument('--bin_size_Mb', type = 'numeric', 
                    default = 5,
                    help = 'Size of genome bins in megabases used for counting fragments in step 1.')
parser$add_argument('--frag_library_size_thre', type = 'numeric',
                    default = 8000,
                    help = 'Minimum number of fragments required for a cell to be included in the analysis. 
                            Cells with fewer fragments will be excluded to ensure analysis quality. ')

# Parse the arguments
args <- parser$parse_args()

# Assign each argument into the namespace
for (arg_name in names(args)) {
  assign(arg_name, args[[arg_name]])
}


#### 1. Edit parameters and directories ####
# Directories
sample_id <- str_remove_all(library_id,"-ATAC|-Multiome|-DEFND")
data_dir <- here(out_dir, sample_id, library_id) # output path for the current library_id (different from out_dir)
plot_dir <- here(data_dir, "plot")
dir_create(plot_dir)

log_info(paste0("Step 2. Now processing ", library_id,"\n"))

#### 2. Read and preprocess data ####
# Read cell x bin matrix
read_depth <- fread(here(data_dir, paste0("cell_x_", bin_size_Mb, "_Mb_bins.csv")))%>%
  dplyr::rename("barcode_atac" = "barcode")

# extract bin names
bins <- colnames(read_depth) 
bins <- bins[!grepl("barcode|celltype|library_id", bins)]

# Read cell type annotations
anno <- fread(here(data_dir,paste0("barcode_anno_", library_id,".csv")))

# Merge annotation with cell x bin matrix. 
# Note cell x bin matrix uses barcodes from ATAC (without aggregation, ending with "-1")
read_depth <- read_depth %>% 
  left_join(anno[,c("barcode_atac","barcode_gex","barcode_atac_aggr","library_id_gex")], 
            by = "barcode_atac")

#### 3. Wrangle and filter based on read depth ####
frag_lib_size <- read_depth %>% 
  tibble::column_to_rownames("barcode_atac_aggr") %>%
  dplyr::select(all_of(bins)) %>%
  mutate(frag_lib_size = rowSums(.)) %>%
  tibble::rownames_to_column("barcode_atac_aggr") %>%
  dplyr::select(barcode_atac_aggr, frag_lib_size)

anno <- left_join(anno, frag_lib_size, by = "barcode_atac_aggr")
read_depth <- left_join(read_depth, frag_lib_size, by = "barcode_atac_aggr")

# Save updated annotation file. The annotation is not filtered yet but we will filter the read count data.
fwrite(anno, file = here(data_dir,paste0("barcode_anno_", library_id,"_",bin_size_Mb,"_Mb_bins_read_depth_filtered.csv")))

# Plot library size histogram
png(here(plot_dir, paste0("Read depth histogram_",bin_size_Mb,"_Mb_bins.png")), width = 600, height = 600)
hist(read_depth$frag_lib_size, xlim = c(0, 1e5), 
     main = paste0("Read depth--Number of fragments per cell\n Median: ", median(read_depth$frag_lib_size)),
xlab = "Number of fragments", breaks = 1000)
dev.off()

# Filter out cells with low number of reads
read_depth <- dplyr::filter(read_depth, frag_lib_size >= frag_library_size_thre)

# Pivot_longer. Store all bin counts (read depth) into one column.
bin_depth_per_cell <- read_depth %>%
  tidyr::pivot_longer(cols = all_of(bins),
                      names_to = 'bins', names_transform = list(bins = as.factor),
                      values_to = 'read_depth')

# Save bin-level read depth data.
fwrite(bin_depth_per_cell, file = here(data_dir, paste0("bin_depth_per_cell_", bin_size_Mb, "_Mb_bins.csv")), row.names = FALSE)

#### 4 Generate chromosome read depth data for negative binomial ####
log_info("Carrying out negative binomial modeling...")
# Calculate beta_i_c similar to above
chrom_depth_per_cell <- bin_depth_per_cell %>%
  dplyr::select(matches("barcode|celltype|library"), frag_lib_size, read_depth, bins) %>%
  mutate(chrom = str_split_i(bins,":",1),
         chrom = factor(chrom, levels = gtools::mixedsort(unique(chrom)))) %>%
  group_by(across(matches("barcode|celltype|library|lib_size")), chrom) %>%
  summarise(chrom_read_depth = sum(read_depth)) 

# Save whole chromosome results 
# The chromosome-level read depth is not affected by the bin size in this version, but it could change the counts if some low-confidence bins (e.g., bins with too many Ns or too repetitive) are filtered out
fwrite(chrom_depth_per_cell, file=here(data_dir, paste0("chrom_neg_binom_read_depth_",bin_size_Mb,"_Mb_bins.csv"))) 
log_success("Step 2 finished for ", library_id, "!")
