library(data.table)
library(tidyr)
library(ggplot2)
library(here)

compile_atac_frags <- function(library_ids, data_dir) {
  results.ls <- lapply(library_ids, function(library_id) {
    fragment_file <- here(data_dir, library_id, "cellx_macs3-aggr_peak_chromosome.csv")
    print(fragment_file)
    fread(fragment_file)
  }) 
  names(results.ls) <- library_ids
  results.df <- rbindlist(results.ls, idcol = "library_id") %>% as.data.frame()
}
library_ids <- list.dirs(here("C:/scratch","kidney_10k_h5ad_output","chasm"), full.names=FALSE, recursive=FALSE)

df <- compile_atac_frags(library_ids, 
                         data_dir = here("C:/scratch","kidney_10k_h5ad_output","chasm"))


df$total_atac_frags <- rowSums((df[,3:ncol(df)]))
df <- df %>%
  pivot_longer(cols = colnames(df)[grepl("chr", colnames(df))], names_to = "seqnames", values_to = "atac_frags")

df <- df %>%
  dplyr::mutate(log_atac_frags = log10(atac_frags / total_atac_frags + 1))

df <- df %>%
  dplyr::mutate(scaled_atac_frags = scale(log_atac_frags, center = 0))

# sex anno
library(openxlsx)
library(dplyr)
meta <- read.xlsx(here("C:/scratch","kidney_10k_h5ad_output","clinical_metadata","clinical_meta.xlsx"))
df <- df %>% left_join(meta[,c("library_id", "sex", "tissue_category", "age_group")], by = "library_id")

library(stringr)
anno <- fread(here("C:/scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv"))
barcodes <- str_split(anno$barcode, pattern="_", simplify=TRUE)[,1]
library_id <- gsub("_multi", "", anno$sample)
library_id <- gsub("_atac", "", library_id)

anno <- data.frame(barcode = barcodes, library_id = library_id, celltype = anno$celltype, sample=anno$sample)

df <- df %>% left_join(anno, by = c("library_id","barcode"))

toplot <- df %>%
  group_by(seqnames) %>%
  mutate(global_median = median(log_atac_frags[sex == "male"])) %>%
  group_by(library_id, celltype) %>%
  summarize(median=median(log_atac_frags), global_median = global_median, sex = sex) %>%
  distinct()
global_median <- unique(toplot$global_median)
toplot %>%
  ggplot(aes(celltype, median, fill=as.factor(sex))) + 
  geom_boxplot() +
  geom_hline(yintercept = global_median, color="red")

# compute global median for each chromosome
# for chrY only consider male samples to compute the median
df <- df %>%
  group_by(sex, seqnames) %>%
  mutate(global_median_counts = median(log_atac_frags))

# compute median for each celltype
# for chrY only consider male samples to compute the median
df <- df %>%
  group_by(sex, seqnames, celltype) %>%
  mutate(celltype_median_counts = median(log_atac_frags))   

# compute celltype correction factor
# and correct the counts
df <- df %>%
  group_by(sex, seqnames, celltype) %>%
  mutate(correction_factor_counts = global_median_counts / celltype_median_counts) %>%
  mutate(corrected_atac_frags = log_atac_frags * correction_factor_counts)

df %>%
  filter(seqnames == "chrY", sex == "male") %>%
  ggplot(aes(celltype, log_atac_frags)) + 
  geom_boxplot() 

df %>%
  filter(seqnames == "chrY", sex == "male") %>%
  ggplot(aes(celltype, corrected_atac_frags)) + 
  geom_boxplot() 


df <- df %>%
  group_by(sex, seqnames, celltype) %>%
  dplyr::mutate(scaled_atac_frags = scale(corrected_atac_frags, center = 0))


df %>%
  dplyr::filter(seqnames == "chrY", sex == "male") %>%
  ggplot(aes(scaled_atac_frags, fill=sex)) +
  geom_density() 

df %>%
  na.omit() %>%
  dplyr::filter(sex == "male", seqnames == "chrY") %>%
  ggplot(aes(scaled_atac_frags, fill=sex)) +
  geom_density(aes(y = ..scaled..)) +
  xlim(c(-0.05,0.1))

kd <- df %>%
  dplyr::filter(seqnames == "chrY", sex == "male", total_atac_frags > 1000)

density <- density(kd$scaled_atac_frags, n=10000)
# search for the first trough between zero and the maximum x coord
x <- density$x
y <- density$y
# calculate the difference in the y coord for successive steps in density function
diffy <- diff(y)
# find the sign of the change for the difference
signy <- sign(diffy)
# find points where the difference in the sign of change for successive steps in density function
# is negative 2 (ie. the first point was ascending and the second point was descending)
diffpeak <- diff(signy) == 2
# find the first two peak values
local_trough_index <- which(diffpeak)[1]
trough <- x[local_trough_index]
toplot <- df %>%
  na.omit() %>%
  dplyr::filter(sex == "male", seqnames == "chrY") %>%
  dplyr::mutate(fill = ifelse(scaled_atac_frags < trough, "LOY", "Haploid"))
# log-normalize, scale and center on the maximum of the density plot
p2 <- toplot %>% 
  ggplot(aes(scaled_atac_frags, fill=sex)) + 
  geom_density(aes(y = ..scaled..)) +
  xlab(paste0("Scaled coverage: ", "chrY")) +
  ylab("Density") +
  theme(legend.position="none")
# create a shaded density plot
p <-  toplot %>% 
  ggplot(aes(scaled_atac_frags)) + 
  geom_density(aes(y = ..scaled..)) +
  geom_vline(xintercept = trough, linetype = "dotted") + 
  xlab(paste0("Scaled coverage: ", "chrY")) +
  ylab("Density") +
  xlim(c(-0.05, 5))
d = ggplot_build(p)$data[[1]]
p = p + geom_area(data = subset(d, x < trough), aes(x=x,y=y), fill = "#00AFBB", alpha = 0.5)
p = p + geom_area(data = subset(d, x > trough), aes(x=x,y=y), fill = "darkgrey", alpha = 0.5)
res <- toplot %>%
  ungroup() %>%
  summarize(loy = scaled_atac_frags < trough) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::mutate(total = sum(Freq)) %>%
  dplyr::mutate(loy_prop = Freq / total)
colnames(res)[1] <- "LOY"
res <- res[res$LOY == "TRUE",]

p3 <- p +
  annotate(geom="text",
           x = 3,
           y = 0.8,
           label = paste0("LOY: ", round(res$loy_prop, digits=2), "\n", "Cells: ", res$total)) +
  theme_bw() +
  xlim(c(-0.05, 5))
rm(p)


loy <- df %>% 
  dplyr::filter(sex == "male", seqnames == "chrY") %>%
  na.omit() %>%
  group_by(library_id, celltype) %>%
  summarize(loy = scaled_atac_frags < trough) %>%
  table() %>%
  as.data.frame() %>%
  pivot_wider(names_from = loy, values_from = Freq, names_prefix = "LOSS_") %>%
  dplyr::mutate(loy_prop = LOSS_TRUE/(LOSS_TRUE + LOSS_FALSE)) %>%
  dplyr::mutate(total_cells = LOSS_TRUE + LOSS_FALSE) %>%
  dplyr::filter(total_cells > 100)

loy %>%
  ggplot(aes(celltype, loy_prop, color=celltype, fill=celltype)) + 
  geom_boxplot() +
  ylab("Proportion LOY in cell type")

df$genotype = ifelse(df$scaled_atac_frags < trough & df$sex == "male", "LOY", "XY")
df$genotype = ifelse(df$sex == "female" & df$total_atac_frags > 1000, "XX", df$genotype)


# read unfiltered annotations
no_filter <- fread("C:/scratch/adata_mvi_model100/annotated_h5ad/annotations.csv") %>%
  dplyr::select(barcode, sample)

# write loy annotations to file
df <- df %>%
  dplyr::filter(seqnames == "chrY") %>%
  dplyr::mutate(barcode = paste0(barcode,"_",sample)) %>%
  dplyr::select(barcode, genotype) 

no_filter <- no_filter %>%
  dplyr::left_join(df, by = "barcode")
no_filter$modality <- gsub(".*_(.*)", "\\1", no_filter$barcode)

fwrite(no_filter, "C:/scratch/adata_mvi_model100/annotated_h5ad/loy_annotations.csv", col.names = TRUE)




