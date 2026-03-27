# SCRATCH1=/mnt/h/scratch
# docker run -it --rm \
# --workdir $HOME \
# -v /mnt/s:$HOME/data \
# -v /mnt/g/reference:$HOME/reference \
# -v /mnt/g/ckd:$HOME/ckd \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# -v /mnt/h/scratch:$HOME/scratch \
# p4rkerw/sctools:R4.3.2 R

# run in rstudio server with NAS mount
# SCRATCH1=/mnt/h/scratch
# workdir=/home/rstudio
# docker run -it --rm \
# -p 8888:8787 \
# -e PASSWORD=password \
# -v /mnt/s:$workdir/data \
# -v /mnt/g/reference:$workdir/reference \
# -v /mnt/g/ckd:$workdir/ckd \
# -v $workdir:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/h/scratch" \
# -v /mnt/h/scratch:$workdir/scratch \
# p4rkerw/sctools:R4.3.2
# navigate browser to localhost:8888
# username: rstudio
# password: password
##################################################
compile_results <- function(library_ids, result_dir, data_dir) {
  results.ls <- lapply(library_ids, function(library_id) {
    fread(here(data_dir, library_id, result_dir, "chasm_all.csv"))
                      }) 
    names(results.ls) <- library_ids
    results.df <- rbindlist(results.ls, idcol = "library_id") %>% as.data.frame()
  }

plot_prop_cna <- function(results, chrom_sel, celltype_sel, pval, group.by, fill.by, facet.by = NULL) {
  results %>%
    dplyr::count(chrom, celltype, library_id, sig = p_val_adjust < pval, sex) %>%
    group_by(chrom, celltype, library_id) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE") %>%
    dplyr::filter(celltype %in% celltype_sel,
                  chrom %in% chrom_sel) %>%
    ggplot(aes(!! sym(group.by), prop, fill=!! sym(fill.by))) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(facet.by)
  }

plot_residual_cna <- function(results, pval, facet.by = NULL) {
  x <- seq(0, 2000, 5)
  y_sqrtx <- sqrt(x) # one copy gain
  y_sqrtx_0.5 <- sqrt(x) * 0.5 # two copies gain
  y_sqrtx_1.5 <- sqrt(x) * 1.5 # three copies gain
  y_sqrtx_neg0.5 <- sqrt(x) * 0.5 * (-1) # one copy loss
  y_sqrtx_neg <- sqrt(x) * (-1) # two copies loss
  df <- data.frame(x, y_sqrtx, y_sqrtx_0.5, y_sqrtx_1.5, y_sqrtx_neg0.5, y_sqrtx_neg)

  ggplot() +
  geom_point(data = results[results$p_val_adjust > pval,],
             aes(x = expected_depth, y = diff_normed2),
             color='#C9C9C9', alpha=0.5, size=0.1) +
  geom_point(data = results[results$p_val_adjust < pval,],
             aes(x = expected_depth, y = diff_normed2),
             color='red',alpha=0.8, size=0.5) +
  geom_line(data = df, aes(x = x, y = y_sqrtx_0.5), color='#FFD1DA') +
  geom_line(data = df, aes(x = x, y = y_sqrtx),  color='#FBA1B7') +
  geom_line(data = df, aes(x = x, y = y_sqrtx_1.5), color='#DB005B') +
  geom_line(data = df, aes(x = x, y = y_sqrtx_neg0.5), color='#8CABFF') +
  geom_line(data = df, aes(x = x, y = y_sqrtx_neg), color='#1D5D9B') +
  facet_wrap(~ celltype) +
  theme_light() +
  theme(legend.position = "none") +
  ylim(c(-25, 25)) +
  xlim(c(0,2000))
}
##################################################
library(here)
library(dplyr)
library(data.table)
library(ggplot2)
library(openxlsx)

install.packages('sjPlot')
install.packages('ggeffects')
install.packages('emmeans')
library(ggeffects)
library(sjPlot) 
library(emmeans)       

outputdir <- here("scratch","kidney_10k_h5ad_output","chasm_aggregate")
plotdir <- here(outputdir, "plots")
dir.create(here(plotdir), recursive=TRUE)

library_ids <- list.dirs(here("scratch","kidney_10k_h5ad_output","chasm"), full.names=FALSE, recursive=FALSE)
result_dirs = c("cellx_macs3-aggr_peak_chromosome")
# compile results
results.df <- compile_results(library_ids, 
                              result_dir = result_dirs,
                              data_dir = here("scratch","kidney_10k_h5ad_output","chasm"))

# estimate copy number stat based on normalized difference
x <- seq(0, 2000, 5)
y_sqrtx_0.5 <- sqrt(x) * 0.5 # one copy gain
y_sqrtx <- sqrt(x) # two copy gain
y_sqrtx_1.5 <- sqrt(x) * 1.5 # three copies gain
y_sqrtx_neg0.5 <- sqrt(x) * 0.5 * (-1) # one copy loss
y_sqrtx_neg <- sqrt(x) * (-1) # two copies loss
df <- data.frame(x, y_sqrtx, y_sqrtx_0.5, y_sqrtx_1.5, y_sqrtx_neg0.5, y_sqrtx_neg)

results.df <- results.df %>%
  dplyr::mutate(copy_number_state = case_when(p_val_adjust < 0.05 & diff_normed2 > 1.25 * sqrt(expected_depth) & diff_normed2 ~ "three_copy_gain",
                                              p_val_adjust < 0.05 & diff_normed2 > 0.75 * sqrt(expected_depth) & diff_normed2 < 1.25 * sqrt(expected_depth) ~ "two_copy_gain",
                                              p_val_adjust < 0.05 & diff_normed2 > 0.25 & diff_normed2 < 0.75 * sqrt(expected_depth) ~ "one_copy_gain",
                                              p_val_adjust < 0.05 & diff_normed2 < -0.25 & diff_normed2 > -0.75 * sqrt(expected_depth) ~ "one_copy_loss",
                                              p_val_adjust < 0.05 & diff_normed2 < -0.75 * sqrt(expected_depth) ~ "two_copy_loss",
                                              .default = "no_copy_change"))

results.df %>%
  dplyr::filter(library_id %in% "CKD_1") %>%
  ggplot() +
  geom_point(aes(x = expected_depth, y = diff_normed2, fill=copy_number_state, color=copy_number_state), size=0.5) +
  geom_line(data = df, aes(x = x, y = y_sqrtx_0.5, color='y = 0.5 * sqrt(x)')) +
  geom_line(data = df, aes(x = x, y = y_sqrtx, color='y = sqrt(x)')) +
  geom_line(data = df, aes(x = x, y = y_sqrtx_1.5, color='y = 1.5 * sqrt(x)')) +
  geom_line(data = df, aes(x = x, y = y_sqrtx_neg0.5, color='y = -0.5 * sqrt(x)')) +
  geom_line(data = df, aes(x = x, y = y_sqrtx_neg, color='y = -1 * sqrt(x)')) +
  theme_light() 

# join with sex annotation
meta <- read.xlsx(here("scratch","kidney_10k_h5ad_output","clinical_metadata","clinical_meta.xlsx"))
results.df <- results.df %>% left_join(meta[,c("library_id", "sex", "tissue_category", "age_group")], by = "library_id")

# count number of barcodes in each library
results.df <- results.df %>%
  group_by(library_id) %>%
  dplyr::mutate(num_cells_library = n_distinct(barcode))
#########################################################################
# LINEAR MODELING
num_autosomal_cna <- results.df %>%
  dplyr::filter(chrom != "chrX") %>%
  dplyr::count(library_id, barcode, calledCNA, name = "num_auto_cna") %>%
  dplyr::filter(calledCNA == "YES") %>%
  dplyr::distinct(library_id, barcode, num_auto_cna)

model.df <- results.df %>% left_join(num_autosomal_cna, by = c("library_id","barcode"))

model.df <- model.df %>% dplyr::mutate(has_chrX = ifelse(p_val_adjust < 0.05 & chrom == "chrX", 1, 0))

model.df <- model.df %>%
  dplyr::mutate(age_decade = case_when(age_group == "10-19" ~ 2,
                                       age_group == "20-29" ~ 3,
                                       age_group == "30-39" ~ 4,
                                       age_group == "40-49" ~ 5,
                                       age_group == "50-59" ~ 6,
                                       age_group == "60-69" ~ 7,
                                       age_group == "70-79" ~ 8,
                                       age_group == "80-89" ~ 9,
                                       age_group == "90-99" ~ 10,
                                       .default = TRUE))
# write chrX cna annotations to file
library(stringr)                  
chrX_anno <- model.df %>% 
    dplyr::filter(chrom == "chrX") %>%
    dplyr::mutate(bc = paste0(barcode, "_", library_id)) %>%
    dplyr::select(-library_id, -celltype, -barcode)
anno <- fread("scratch/adata_mvi_model100/annotated_h5ad/annotations.csv")
anno$library_id <- gsub("_atac|_multi", "", anno$sample)
anno <- anno %>%
    dplyr::mutate(bc = str_split(barcode, pattern = "_", simplify=TRUE)[,1]) %>%
    dplyr::mutate(bc = paste0(bc, "_", library_id))
chrX_anno <- anno %>% left_join(chrX_anno, by = c("bc"))
chrX_anno <- chrX_anno %>% mutate(modality = str_extract(sample, pattern="atac|multi"))                  
fwrite(chrX_anno, "scratch/adata_mvi_model100/annotated_h5ad/chrX_annotations.csv")

# write filtered annotations to file
chrX_anno <- model.df %>% 
    dplyr::filter(chrom == "chrX") %>%
    dplyr::mutate(bc = paste0(barcode, "_", library_id)) %>%
    dplyr::select(-library_id, -celltype, -barcode)
filtered_anno <- fread("scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv")
filtered_anno$library_id <- gsub("_atac|_multi", "", filtered_anno$sample)
filtered_anno <- filtered_anno %>%
    dplyr::mutate(bc = str_split(barcode, pattern = "_", simplify=TRUE)[,1]) %>%
    dplyr::mutate(bc = paste0(bc, "_", library_id))
chrX_filtered_anno <- filtered_anno %>% left_join(chrX_anno, by = c("bc"))
chrX_filtered_anno <- chrX_filtered_anno %>% mutate(modality = str_extract(sample, pattern="atac|multi"))                  
fwrite(chrX_filtered_anno, "scratch/adata_mvi_model100/annotated_h5ad/filtered_chrX_annotations.csv")    
                  
df <- model.df %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("PT_PROM1") ~ "TL",
                                     celltype %in% c("DCT1","DCT2_PC") ~ "DCT",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype)))


df <- df %>%
  dplyr::filter(chrom == "chrX") %>%
  distinct(library_id, barcode, has_chrX, num_auto_cna, celltype, age_decade, sex, lib_size, tissue_category)
df$num_auto_cna[is.na(df$num_auto_cna)] <- 0

df <- df %>%
  dplyr::mutate(is_control = case_when(tissue_category %in% c("AKI","CKD") ~ 0,
                                     tissue_category %in% c("Control") ~ 1,
                                     TRUE ~ NA))

df$celltype <- factor(df$celltype, levels = c("PT","PT_VCAM1","PEC","TL","TAL","DCT","IC","PODO","ENDO","FIB_VSMC_MC","TCELL","BCELL","MONO"))
df$sex <- as.factor(df$sex)

# check for libraries with too few cells to model as a fixed or random effect
df <- df %>% 
  group_by(library_id) %>% 
  mutate(num_cells_library = n_distinct(barcode)) 


# model.df <- model.df %>%
#   dplyr::group_by(library_id) %>%
#   dplyr::mutate(scaled_cnv = scale(num_CNA))

install.packages(c('ggeffects','sjPlot','emmeans'))
library(lme4)
library(ggeffects)
library(sjPlot)
library(emmeans)

mod1 <- glm(has_chrX ~ celltype + age_decade,
            family = binomial(link = "logit"),
            data = df)
plot_model(mod1)
plot_model(mod1, type = "emm", terms = c("age_decade", "celltype"))

plot1 <- plot_model(mod1, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot1 + theme_bw() + 
  theme(legend.title=element_blank(),
        legend.pos=c(0.2,0.75),
        legend.key.size = unit(1,"line"),
        plot.title = element_text(face="bold"),
        axis.text = element_text(color="black"),
        panel.grid = element_blank()) +
  ylab("Predicted Probability for chrX") +
  xlab("Age (decade)") 


mod2 <- glm(has_chrX ~ celltype + age_decade + num_auto_cna,
            family = binomial(link = "logit"),
            data = df)
pmod2 <- plot_model(mod2)
pmod2

plot_model(mod2, type = "emm", terms = c("age_decade","celltype"))
plot_model(mod2, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod2, type = "emm", terms = c("num_auto_cna","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod2, type = "emm", terms = c("num_auto_cna","celltype"))
plot_model(mod2, type = "emm", terms = c("age_decade","num_auto_cna"))
plot_model(mod2, type = "emm", terms = c("num_auto_cna","age_decade"))

mod3 <- glm(has_chrX ~ celltype + age_decade + sex,
            family = binomial(link = "logit"),
            data = df)
pmod3 <- plot_model(mod3)
pmod3

plot_model(mod3, type = "emm", terms = c("age_decade","sex"))
plot_model(mod3, type = "emm", terms = c("age_decade","celltype"))
plot_model(mod3, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))


mod4 <- glm(has_chrX ~ celltype + age_decade + sex + num_auto_cna,
            family = binomial(link = "logit"),
            data = df[df$num_cells_library > 500,])
pmod4 <- plot_model(mod4)
pmod4

plot_model(mod4, type = "emm", terms = c("age_decade","sex"))
plot_model(mod4, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod4, type = "emm", terms = c("age_decade","sex","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod4, type = "emm", terms = c("sex","num_auto_cna"))

mod5 <- glm(has_chrX ~ celltype + age_decade + sex + is_control,
            family = binomial(link = "logit"),
            data = df)
pmod5 <- plot_model(mod5)
pmod5

plot_model(mod5, type = "emm", terms = c("age_decade","sex"))
plot_model(mod5, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod5, type = "emm", terms = c("age_decade","sex","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod5, type = "emm", terms = c("sex","is_control"))

# with a mixed effect
mod6 <- glmer(has_chrX ~ celltype + scale(age_decade) + sex + scale(num_auto_cna) + (1|library_id),
              family = binomial(link = "logit"),
              data = df)
pmod6 <- plot_model(mod6)
pmod6
plot_model(mod6)
plot_model(mod6, type = "emm", terms = c("age_decade","sex"))
plot_model(mod6, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod6, type = "emm", terms = c("age_decade","sex","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod6, type = "emm", terms = c("age_decade", "num_auto_cna"))


mod7 <- glmer(has_chrX ~ celltype + age_decade + sex + is_control + num_auto_cna + (1|library_id),
              family = binomial(link = "logit"),
              data = df)
pmod7 <- plot_model(mod7)
pmod7
plot_model(mod7, type = "emm", terms = c("age_decade","sex"))
plot_model(mod7, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod7, type = "emm", terms = c("age_decade","sex","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod7, type = "emm", terms = c("age_decade", "num_auto_cna"))
plot_model(mod7, type = "emm", terms = c("age_decade","is_control"))
plot_model(mod6, type = "emm", terms = c("age_decade","is_control"))

mod8 <- glmer(has_chrX ~ celltype + age_decade + sex + is_control + num_auto_cna + (1|library_id),
              family = binomial(link = "logit"),
              glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)),
              data = df)
pmod8 <- plot_model(mod8)
pmod8
plot_model(mod8, type = "emm", terms = c("age_decade","sex"))
plot_model(mod8, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod8, type = "emm", terms = c("age_decade","sex","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod8, type = "emm", terms = c("age_decade", "num_auto_cna"))
plot_model(mod8, type = "emm", terms = c("age_decade","is_control"))
plot_model(mod8, type = "emm", terms = c("age_decade","is_control"))

# only analyze cells from libraries with > 500 cells to allow fitting mixed effect for library_id
mod9 <- glmer(has_chrX ~ celltype + age_decade + sex + is_control + num_auto_cna + (1|library_id),
              family = binomial(link = "logit"),
              glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)),
              data = df[df$num_cells_library > 500,])
pmod9 <- plot_model(mod9)
pmod9
plot_model(mod9, type = "emm", terms = c("age_decade","sex"))
plot_model(mod9, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod9, type = "emm", terms = c("age_decade","sex","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
plot_model(mod9, type = "emm", terms = c("age_decade", "num_auto_cna"))
plot_model(mod9, type = "emm", terms = c("age_decade","is_control"))
plot_model(mod9, type = "emm", terms = c("age_decade","is_control"))




########################################################################


# proportion of celltype with selected CNA
filename <- paste0("selected_chrom_cna_", result_dirs, ".png")
p <- plot_prop_cna(results = results.df, chrom_sel = c("chrX","chr1","chr7"),
               celltype_sel = c("PCT","PST","PT_VCAM1","PT_PROM1","PODO"),
               pval = 0.01,
               group.by = "chrom",
               fill.by = "celltype",
               facet.by = "sex")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
print(p)
dev.off()

results.df %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::filter(p_val_adjust < 0.05) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(mean_diff_normed2 = mean(abs(diff_normed2))) %>%
    dplyr::filter(!(copy_number_state %in% c("one_copy_gain","one_copy_loss"))) %>%
    dplyr::distinct(library_id, chrom, celltype, mean_diff_normed2, copy_number_state, sex, tissue_category) %>%
    ggplot(aes(celltype, mean_diff_normed2, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom) +
      coord_cartesian(ylim=c(0,20))


results.df %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, copy_number_state, sig = p_val_adjust < 0.01, sex, tissue_category) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chrX") %>%
    dplyr::distinct(library_id, chrom, celltype, copy_number_state, prop, sex, tissue_category) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom + copy_number_state)

results.df %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, copy_number_state, sig = p_val_adjust < 0.01, sex, tissue_category) %>%
    group_by(library_id, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chr7") %>%
    dplyr::distinct(library_id, chrom, celltype, copy_number_state, prop, sex, tissue_category) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~copy_number_state) 

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex, tissue_category) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chrX") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex, tissue_category) %>%
    ggplot(aes(celltype, prop)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom + sex) +
      coord_cartesian(ylim = c(0,0.2))

# calculate overall proportion for selected chromosome CNA by tissue category (AKI vs CKD vs control) in a single cell type
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2500) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex, tissue_category, age_group) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chrX", celltype == "PCT") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex, tissue_category, age_group) %>%
    ggplot(aes(age_group, prop, color=sex, fill=sex)) +
      geom_boxplot() +
      geom_point() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom + tissue_category) +
      coord_cartesian(ylim = c(0,0.2))

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chr7") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom) +
      coord_cartesian(ylim = c(0,0.2))

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 1000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 50) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom != "chrY") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom) +
      coord_cartesian(ylim = c(0,0.2))

# quantify number of CNA per cell
results.df %>%
  dplyr::select(library_id, barcode, celltype, sex, p_val_adjust) %>%
  group_by(library_id, barcode) %>%
  dplyr::add_count(name = "num_barcode_cna", sig = p_val_adjust < 0.05) %>%
  dplyr::filter(sig == "TRUE") %>%
  dplyr::distinct(library_id, barcode, celltype, num_barcode_cna, sex) %>%
  ggplot(aes(celltype, num_barcode_cna, fill=celltype)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~sex)

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.05, sex) %>%
    group_by(chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE") %>%
    dplyr::filter(chrom %in% "chrX") %>%
    ggplot(aes(celltype, prop, fill=celltype)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~sex)

results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::filter(chrom != "chrY") %>%
    group_by(library_id) %>%
    dplyr::mutate(num_cells_library_id = n_distinct(barcode)) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_id_celltype = n_distinct(barcode)) %>%
    dplyr::add_count(p_val_adjust < 0.05, name="num_cna_library_id_celltype") %>%
    dplyr::add_count(p_val_adjust >= 0.05, name="num_non_cna_library_id_celltype") %>%
    mutate(prop = num_cna_library_id_celltype / sum(num_cells_library_id_celltype + num_non_cna_library_id_celltype)) %>%
    distinct(celltype, sex, prop) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      coord_cartesian(ylim = c(0,0.02))

results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::filter(chrom != "chrY") %>%
    group_by(library_id) %>%
    dplyr::mutate(num_cells_library_id = n_distinct(barcode)) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_id_celltype = n_distinct(barcode)) %>%
    dplyr::add_count(p_val_adjust < 0.01, name="num_cna_library_id_celltype") %>%
    dplyr::add_count(p_val_adjust >= 0.01, name="num_non_cna_library_id_celltype") %>%
    dplyr::group_by(library_id, celltype, chrom) %>%
    dplyr::mutate(num_cells_library_id_celltype_chrom = n_distinct(barcode)) %>%
    dplyr::add_count(p_val_adjust < 0.01, name="num_cna_library_id_celltype_chrom") %>%
    dplyr::add_count(p_val_adjust >= 0.01, name="num_non_cna_library_id_celltype_chrom") %>%
    mutate(prop = num_cna_library_id_celltype_chrom / sum(num_cells_library_id_celltype_chrom + num_non_cna_library_id_celltype_chrom)) %>%
    distinct(celltype, sex, prop) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      coord_cartesian(ylim = c(0,0.02)) 


p <- plot_prop_cna(results = results.df, 
               chrom_sel = c("chrX"),
               celltype_sel = unique(results.df$celltype),
               pval = 0.01,
               group.by = "celltype",
               fill.by = "celltype")
filename <- paste0("chrX_celltype_cna_", result_dir, ".png")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
print(p)
dev.off()
  
# Visualization
filename <- paste0("combined_normalized_residual_cna_state_plot_annotated_", result_dir, ".png")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
plot_residual_cna(results.df, pval = 0.05, facet.by = "celltype")
dev.off()

chroms <- unique(results.df$chrom)
names(chroms) <- chroms
 p <- plot_prop_cna(results = results.df, 
             chrom_sel = chroms[names(chroms) != "chrY"],
             celltype_sel = c("PCT","PT_VCAM1"),
             pval = 0.01,
             group.by = "chrom",
             fill.by = "celltype")
filename <- paste0("summary_celltype_", result_dir, ".png")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
print(p)
dev.off()

#####################################################################################################
#####################################################################################################
#####################################################################################################
library(here)
library(dplyr)
library(data.table)
library(ggplot2)
library(openxlsx)

outputdir <- here("scratch","kidney_10k_h5ad_output","chasm_aggregate")
plotdir <- here(outputdir, "plots")
dir.create(here(plotdir), recursive=TRUE)

library_ids <- list.dirs(here("scratch","kidney_10k_h5ad_output","chasm"), full.names=FALSE, recursive=FALSE)
result_dirs = c("cellx_macs3-aggr_peak_chromosome_arm")
# compile results
results.df <- compile_results(library_ids, 
                              result_dir = result_dirs,
                              data_dir = here("scratch","kidney_10k_h5ad_output","chasm"))
  
# join with sex annotation
meta <- read.xlsx(here("scratch","kidney_10k_h5ad_output","clinical_metadata","clinical_meta.xlsx"))
results.df <- results.df %>% left_join(meta[,c("library_id", "sex", "tissue_category", "age_group")], by = "library_id")

# count number of barcodes in each library
results.df <- results.df %>%
  group_by(library_id) %>%
  dplyr::mutate(num_cells_library = n_distinct(barcode))

# proportion of celltype with selected CNA
plot_prop_cna(results = results.df, chrom_sel = c("chrXq","chr1q","chr7q"),
               celltype_sel = c("PCT","PST","PT_VCAM1","PT_PROM1","PODO"),
               pval = 0.01,
               group.by = "chrom",
               fill.by = "celltype",
               facet.by = "sex")

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex, tissue_category) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chrXq") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex, tissue_category) %>%
    ggplot(aes(celltype, prop, fill=tissue_category)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom + sex) +
      coord_cartesian(ylim = c(0,0.2))

# calculate overall proportion for selected chromosome CNA by tissue category (AKI vs CKD vs control) in a single cell type
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex, tissue_category, age_group) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chrXp", celltype == "PCT") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex, tissue_category, age_group) %>%
    ggplot(aes(age_group, prop, fill=sex)) +
      geom_boxplot() +
      geom_point() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom + tissue_category) +
      coord_cartesian(ylim = c(0,0.2))

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 2000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 100) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom == "chrXq") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom) +
      coord_cartesian(ylim = c(0,0.2))




# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::group_by(library_id) %>%
    dplyr::mutate(num_cells_library = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library > 1000) %>%
    dplyr::group_by(library_id, celltype) %>%
    dplyr::mutate(num_cells_library_celltype = n_distinct(barcode)) %>%
    dplyr::filter(num_cells_library_celltype > 50) %>%
    dplyr::filter(lib_size > 1000) %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.01, sex) %>%
    group_by(library_id, chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom != "chrY") %>%
    dplyr::distinct(library_id, chrom, celltype, prop, sex) %>%
    ggplot(aes(celltype, prop, fill=sex)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom) +
      coord_cartesian(ylim = c(0,0.2))

# quantify number of CNA per cell
results.df %>%
  dplyr::select(library_id, barcode, celltype, sex, p_val_adjust) %>%
  group_by(library_id, barcode) %>%
  dplyr::add_count(name = "num_barcode_cna", sig = p_val_adjust < 0.05) %>%
  dplyr::filter(sig == "TRUE") %>%
  dplyr::distinct(library_id, barcode, celltype, num_barcode_cna, sex) %>%
  ggplot(aes(celltype, num_barcode_cna, fill=celltype)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~sex)

# calculate overall proportion for selected chromosome CNA
results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.05, sex) %>%
    group_by(chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE") %>%
    dplyr::filter(chrom %in% "chrX") %>%
    ggplot(aes(celltype, prop, fill=celltype)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~sex)

results.df %>%
    na.omit() %>%
    ungroup() %>%
    dplyr::count(chrom, celltype, sig = p_val_adjust < 0.05, sex) %>%
    group_by(chrom, celltype) %>% 
    mutate(prop = n / sum(n)) %>%
    dplyr::filter(sig == "TRUE", chrom != "chrY") %>%
    ggplot(aes(celltype, prop, fill=celltype)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(~chrom + sex)

p <- plot_prop_cna(results = results.df, 
               chrom_sel = c("chrX"),
               celltype_sel = unique(results.df$celltype),
               pval = 0.01,
               group.by = "celltype",
               fill.by = "celltype")
filename <- paste0("chrX_celltype_cna_", result_dir, ".png")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
print(p)
dev.off()
  
# Visualization
filename <- paste0("combined_normalized_residual_cna_state_plot_annotated_", result_dir, ".png")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
plot_residual_cna(results.df, pval = 0.05, facet.by = "celltype")
dev.off()

chroms <- unique(results.df$chrom)
names(chroms) <- chroms
 p <- plot_prop_cna(results = results.df, 
             chrom_sel = chroms[names(chroms) != "chrY"],
             celltype_sel = c("PCT","PT_VCAM1"),
             pval = 0.01,
             group.by = "chrom",
             fill.by = "celltype")
filename <- paste0("summary_celltype_", result_dir, ".png")
png(here(plotdir, filename), width=3000, height=3000,res = 300)
print(p)
dev.off()

#####################################################################################    
#####################################################################################
                 
library_ids <- list.dirs(here("scratch","kidney_10k_h5ad_output","chasm"), full.names=FALSE, recursive=FALSE)
result_dirs = c("cellx_macs3-aggr_peak_chromosome")
# compile results
results.df <- compile_results(library_ids, 
                              result_dir = result_dirs,
                              data_dir = here("scratch","kidney_10k_h5ad_output","chasm"))
                  
results.df$celltype <- ifelse(results.df$celltype == "PT_PROM1","TL",results.df$celltype)

# join with sex annotation
meta <- read.xlsx(here("scratch","kidney_10k_h5ad_output","clinical_metadata","clinical_meta.xlsx"))
results.df <- results.df %>% left_join(meta[,c("library_id", "sex", "tissue_category", "age_group")], by = "library_id")

results.df <- results.df %>%
  dplyr::mutate(copy_number_state = case_when(p_val_adjust < 0.05 & diff_normed2 > 1.25 * sqrt(expected_depth) ~ "three_copy_gain",
                                              p_val_adjust < 0.05 & diff_normed2 > 0.75 * sqrt(expected_depth) & diff_normed2 < 1.25 * sqrt(expected_depth) ~ "two_copy_gain",
                                              p_val_adjust < 0.05 & diff_normed2 > 0.25 & diff_normed2 < 0.75 * sqrt(expected_depth) ~ "one_copy_gain",
                                              p_val_adjust < 0.05 & diff_normed2 < -0.25 & diff_normed2 > -0.75 * sqrt(expected_depth) ~ "one_copy_loss",
                                              p_val_adjust < 0.05 & diff_normed2 < -0.75 * sqrt(expected_depth) ~ "two_copy_loss",
                                              .default = "no_copy_change"))
                  
# count number of barcodes in each library
results.df <- results.df %>%
  group_by(library_id) %>%
  dplyr::mutate(num_cells_library = n_distinct(barcode))

num_autosomal_cna <- results.df %>%
  dplyr::filter(chrom != "chrX") %>%
  dplyr::count(library_id, barcode, calledCNA, name = "num_auto_cna") %>%
  dplyr::filter(calledCNA == "YES") %>%
  dplyr::distinct(library_id, barcode, num_auto_cna)

# num_total_cna <- results.df %>%
#  dplyr::count(library_id, barcode, calledCNA, name = "num_total_cna") %>%
#  dplyr::filter(calledCNA == "YES") %>%
#  dplyr::distinct(library_id, barcode, num_total_cna)
                  
model.df <- results.df %>% left_join(num_autosomal_cna, by = c("library_id","barcode"))
# model.df <- results.df %>% left_join(num_total_cna, by = c("library_id","barcode"))

# model.df$num_total_cna[is.na(model.df$num_total_cna)] <- 0                
model.df$num_auto_cna[is.na(model.df$num_auto_cna)] <- 0                  

model.df$celltype <- factor(model.df$celltype, levels = c("PCT","PST","PT_VCAM1","PEC","TL","TAL","DCT1","DCT2_PC","ICA","ICB","PODO","ENDO","FIB_VSMC_MC","TCELL","BCELL","MONO"))
                  
model.df %>%
  distinct(library_id, barcode, celltype, num_auto_cna) %>%
  na.omit() %>%
  dplyr::group_by(library_id, celltype) %>%
  dplyr::summarize(mean_auto_cna = sum(num_auto_cna) / n_distinct(library_id, barcode)) %>%                
  ggplot(aes(celltype, mean_auto_cna, fill=celltype)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(0,5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of autosomal CNA") +
  xlab("")

model.df <- model.df %>% dplyr::mutate(has_chrX = as.factor(ifelse(p_val_adjust < 0.05 & chrom == "chrX", 1, 0)))

model.df <- model.df %>%
  dplyr::mutate(age_decade = case_when(age_group == "10-19" ~ 2,
                                       age_group == "20-29" ~ 3,
                                       age_group == "30-39" ~ 4,
                                       age_group == "40-49" ~ 5,
                                       age_group == "50-59" ~ 6,
                                       age_group == "60-69" ~ 7,
                                       age_group == "70-79" ~ 8,
                                       age_group == "80-89" ~ 9,
                                       age_group == "90-99" ~ 10,
                                       .default = TRUE))
# write chrX cna annotations to file
library(stringr)                  
chrX_anno <- model.df %>% 
    dplyr::filter(chrom == "chrX") %>%
    dplyr::mutate(bc = paste0(barcode, "_", library_id)) %>%
    dplyr::select(-library_id, -celltype, -barcode)
anno <- fread("scratch/adata_mvi_model100/annotated_h5ad/annotations.csv")
anno$library_id <- gsub("_atac|_multi", "", anno$sample)
anno <- anno %>%
    dplyr::mutate(bc = str_split(barcode, pattern = "_", simplify=TRUE)[,1]) %>%
    dplyr::mutate(bc = paste0(bc, "_", library_id))
chrX_anno <- anno %>% left_join(chrX_anno, by = c("bc"))
chrX_anno <- chrX_anno %>% mutate(modality = str_extract(sample, pattern="atac|multi"))                  
fwrite(chrX_anno, "scratch/adata_mvi_model100/annotated_h5ad/chrX_annotations.csv")                   

# write filtered annotations to file
chrX_anno <- model.df %>% 
    dplyr::filter(chrom == "chrX") %>%
    dplyr::mutate(bc = paste0(barcode, "_", library_id)) %>%
    dplyr::select(-library_id, -celltype, -barcode)
filtered_anno <- fread("scratch/adata_mvi_model100/annotated_h5ad/filtered_annotations.csv")
filtered_anno$library_id <- gsub("_atac|_multi", "", filtered_anno$sample)
filtered_anno <- filtered_anno %>%
    dplyr::mutate(bc = str_split(barcode, pattern = "_", simplify=TRUE)[,1]) %>%
    dplyr::mutate(bc = paste0(bc, "_", library_id))
chrX_filtered_anno <- filtered_anno %>% left_join(chrX_anno, by = c("bc"))
chrX_filtered_anno <- chrX_filtered_anno %>% mutate(modality = str_extract(sample, pattern="atac|multi"))                  
fwrite(chrX_filtered_anno, "scratch/adata_mvi_model100/annotated_h5ad/filtered_chrX_annotations.csv")    

df <- model.df %>%
  dplyr::mutate(celltype = case_when(celltype %in% c("PCT","PST") ~ "PT",
                                     celltype %in% c("DCT1","DCT2_PC") ~ "DCT",
                                     celltype %in% c("ICA","ICB") ~ "IC",
                                     TRUE ~ as.character(celltype)))


df <- df %>%
  dplyr::filter(chrom == "chrX") %>%
  distinct(library_id, barcode, has_chrX, num_auto_cna, celltype, age_decade, sex, lib_size, tissue_category)
df$num_auto_cna[is.na(df$num_auto_cna)] <- 0

df <- df %>%
  dplyr::mutate(is_control = case_when(tissue_category %in% c("AKI","CKD") ~ 0,
                                     tissue_category %in% c("Control") ~ 1,
                                     TRUE ~ NA))
df$is_control <- as.factor(df$is_control)
df$celltype <- factor(df$celltype, levels = c("PT","PT_VCAM1","PEC","TL","TAL","DCT","IC","PODO","ENDO","FIB_VSMC_MC","TCELL","BCELL","MONO"))
df$sex <- as.factor(df$sex)

# check for libraries with too few cells to model as a fixed or random effect
df <- df %>% 
  group_by(library_id) %>% 
  mutate(num_cells_library = n_distinct(barcode)) 
        
mod0 <- lm(num_auto_cna ~ celltype + age_decade, data = df)
plot0 <- plot_model(mod0)
# plot0 <- plot_model(mod0, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1,FIB_VSMC_MC,BCELL,TCELL,MONO]"))
# swap out yaxis labels
newlabels <- str_replace(plot0$data$term,"celltype","")
newlabels <- factor(newlabels, levels=newlabels)

plot0 <- plot0 + theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(plot.title = element_text(face="bold"), axis.text = element_text(color="black")) + 
  ylab("Estimate for number of autosomal CNA relative to PT") + 
  scale_x_discrete(labels=rev(newlabels)) +
  ggtitle("")
plot0
                  
model.df %>%
  filter(chrom == "chrX") %>%
  distinct(library_id, barcode, celltype, has_chrX) %>%
  na.omit() %>%
  dplyr::group_by(library_id, celltype) %>%
  dplyr::summarize(prop_chrX_cna = sum(ifelse(has_chrX == 1, 1, 0)) / n_distinct(library_id, barcode)) %>%                
  ggplot(aes(celltype, prop_chrX_cna, fill=celltype)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(0,0.5)) +
  theme_bw() +
  theme(legend.pos = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Proportion of cells with chrX CNA") +
  xlab("")

#################################################################################################
#################################################################################################
#################################################################################################                  
# chrX models                  
mod1 <- glm(has_chrX ~ celltype + age_decade,
            family = binomial(link = "logit"),
            data = df)
plot_model(mod1, type = "emm", terms = c("age_decade", "celltype"))
plot_model(mod1, type = "emm", terms = c("age_decade","celltype [PT,PT_VCAM1]"))







