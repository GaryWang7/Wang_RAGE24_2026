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
  library(Signac)
  library(Seurat)
  library(hdf5r)
  library(here)
  library(GenomicRanges)
  library(BSgenome)
  library(dplyr)
  library(future)
  library(GenomeInfoDb)  
  library(plyranges)
  library(data.table)
  library(stringr)
  library(future)
  library(future.apply)
  library(fs)
  library(argparse)
  library(logger)
})

#### Define functions ####
# Convert gex barcodes to corresponding atac barcode in the same 10X gel bead
convertGEXtoATAC_bc <- function(gex_bc, aggr_ident = FALSE){
  gex_bc_pure <- str_split_i(gex_bc, pattern = "-", i = 1)
  atac_bc_pure <- barcode_conv_df$ATAC_bc[which(barcode_conv_df$GEX_bc == gex_bc_pure)]
  if(aggr_ident == FALSE){
    # If aggr_ident = FALSE, we return barcodes that always end with -1
    return(paste(atac_bc_pure, 1,sep = "-"))
  }else{
    gex_bc_no <- str_split_i(gex_bc, pattern = "-", i = 2)
    return(paste(atac_bc_pure, gex_bc_no,sep = "-"))
  }
}

#### 1. Define arguments and parameters ####
# Create a parser
parser <- ArgumentParser(description = 'Step 1. Count fragments',
                         formatter_class = 'argparse.ArgumentDefaultsHelpFormatter')

# Add arguments
parser$add_argument('--library_id', type = 'character', 
                    default = 'RAGE24-C03-KYC-LN-01-DEFND-ATAC',
                    help = 'Library id in the format of "RAGE24-C03-KYC-LN-01-DEFND-ATAC"')
parser$add_argument('--cellranger_dir', type = 'character',
                    default = 'garyw/RAGE24/cellranger_flowcell_combine',
                    help = 'The family directory of all Cell Ranger outputs. The fragment file should be found in cellranger_dir/.../library_id/outs/.') 
parser$add_argument('--out_dir', type = 'character', 
                    default = 'garyw/RAGE24/CNV/result',
                    help = 'Output path for all CNV results. If it does not exists, it will be created. 
                    Results will be stored in out_dir/sample_id/library_id, where sample_id is extracted from library_id.')
parser$add_argument('--bin_size_Mb', type = 'numeric', 
                    default = 5,
                    help = 'Size of genome bins in megabases used for counting fragments')
parser$add_argument('--anno_csv', type = 'character',
                    default = 'garyw/RAGE24/gex_aggr_prep_combine/gex_aggr_anno.csv',
                    help = 'Path to annotation csv file containing info from all cells from gex_aggr.')
parser$add_argument('--ref_dir', type = 'character',
                    default = 'garyw/RAGE24/reference/GRCr8/',
                    help = 'Path to folder for storing reference genomes and information, etc.')
parser$add_argument('--ref_info_tsv', type = 'character',
                    default = 'GCF_036323735.1_GRCr8_assembly_info.tsv',
                    help = 'Name of the sequence info tsv file containing chromosome names and lengths. This should be placed in the ref_dir')
parser$add_argument('--genome_fasta', type = 'character',
                    required = FALSE,
                    default = 'GCF_036323735.1_GRCr8_mt.fna',
                    help = 'Name of genome fasta used for cellranger-atac mapping. This should be placed in the ref_dir.')
parser$add_argument('--bin_N_metrics', type = 'logical',
                    default = TRUE,
                    help = 'Whether to calculate the percentage of N and GC of the genome bin sequences. If TRUE, will create a tsv file in the ref_dir. The genome is specified by --genome_fasta.')
parser$add_argument('--barcode_conv_csv', type = 'character', 
                    default = 'garyw/RAGE24/chasm2/resources/barcode_conversion_10Xmultiome.csv',
                    help = 'Path to the barcode conversion file converting barcodes between gex and atac.')

# Parse the arguments
args <- parser$parse_args()

# Assign each argument into the namespace
for (arg_name in names(args)) {
  assign(arg_name, args[[arg_name]])
}

# Edit bin size
bin_size <- bin_size_Mb*1e6

# Edit/create directories
sample_id <- str_remove_all(library_id,"-ATAC|-Multiome|-DEFND")
log_info(paste0("Locating cellranger output folder for ", library_id, "\n"))
input_dir <- dir_ls(
  path = cellranger_dir, regexp = paste0("/.*/",library_id,"/outs", "$"),
  recurse = TRUE, type = "directory"
) %>% here()
output_dir <- here(out_dir, sample_id, library_id) # output path for the current library_id (different from out_dir)
plot_dir <- here(output_dir, "plot")

dir_create(output_dir)
dir_create(plot_dir)

log_info(paste0("Step 1. Now processing ", library_id,"\n"))

# chromosome information for reference genome
ref_info <- read.table(file = here(ref_dir,ref_info_tsv), sep = '\t', header = TRUE)

#### Part 2. Read in/Create genomic windows ####
bin_file <- here(ref_dir, paste0(genome_fasta, "_",bin_size_Mb,"_Mb_bins.rds"))
if(file_exists(bin_file)){
  log_info("Found genome bin file!\n")
  bins.raw.gr <- readRDS(file = bin_file)
}else{
  log_info("Creating reference genome bin file... \n")
  # Create genome bins from chromosome name and lengths
  ref_seqlengths <- ref_info$Sequence.Length
  names(ref_seqlengths) <- ref_info$UCSC.style.name
  bins.raw.gr <- tileGenome(seqlengths =ref_seqlengths, tilewidth = bin_size, cut.last.tile.in.chrom = FALSE)%>%
    as_granges()%>%
    keepStandardChromosomes(pruning.mode = "coarse") %>%
    dropSeqlevels(value = "chrM", pruning.mode = "coarse")
  bins.raw.gr$group_name <- NULL
  bins.raw.gr$region <- as.character(bins.raw.gr)
  saveRDS(bins.raw.gr, bin_file)
}

#### Part 2.1 Generate N and GC content metrics for the reference genome ####
if(bin_N_metrics==TRUE){
  bin_N_GC_metrics_file <- here(ref_dir, paste0("N_GC_stats_",genome_fasta,"_",bin_size_Mb,"_Mb_bins.tsv"))
  if(file_exists(bin_N_GC_metrics_file)){
    log_info("Calculate N and GC content for genome bins: TRUE. \n but metrics file already exists! Skipping...\n")
  }else{
    log_info("Calculate N and GC content for genome bins: TRUE. Reading reference genomes")
    # Read the reference genome. Make sure this is the same one you used for cellranger-arc mkref
    ref_genome <- biomartr::read_genome(file = here(ref_dir, genome_fasta))
    # Edit seqnames in ref_genome. This is applicable to data downloaded from Ensembl data base and GRC/ncbi (refseq)
    names(ref_genome) <- str_split_i(names(ref_genome), " ", 1)
    
    # For rats/some species: we need to convert accessions to chromosome names. For example, "NC_086019.1" to "chr1"
    refseq_accession <- ref_info$RefSeq.Accn
    chr_names <- ref_info$UCSC.style.name
    bins_chr_names <- unique(seqnames(bins.raw.gr))
    if(all(names(ref_genome) %in% refseq_accession)){
      names(ref_genome) <- chr_names[match(names(ref_genome), refseq_accession)]
      log_info("All sequence names in reference genome were converted to UCSC style.")
      }else if(all(bins_chr_names %in% names(ref_genome))){
      log_info("Bin sequences can be extracted from reference genome, 
               though there might be differences in 
               sequence names for unlocalized/unplaced scaffolds.")
    }else{
      stop("Error. The chromosomes do not match RefSeq accessions")
    }
    # Extract sequences for the bins
    bins.seqs <- BSgenome::getSeq(ref_genome, bins.raw.gr)
    # Count N and GC%
    bins.raw.gr$N_counts <- as.integer(letterFrequency(letters = "N", bins.seqs))
    bins.raw.gr$GC_counts <- rowSums(letterFrequency(bins.seqs, letters = c("G", "C")))
    bins.raw.gr$seq_lengths <- width(bins.seqs)
    bins.raw.gr$N_percentage <- round(bins.raw.gr$N_counts/bins.raw.gr$seq_lengths*100, digits = 4)
    bins.raw.gr$GC_percentage <- round(bins.raw.gr$GC_counts/bins.raw.gr$seq_lengths*100, digits = 4)
    # Save a summary as genome N statistics
    write.table(as.data.frame(bins.raw.gr), 
                file = bin_N_GC_metrics_file, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    log_info(paste0("N and GC content metrics saved to ", bin_N_GC_metrics_file), ".\nPlease check histogram of N% and GC% of the bins.")
    png(filename = here(ref_dir, paste0("hist_N_stats_",genome_fasta,"_",bin_size_Mb,"_Mb_bins.png")), width = 1000, height = 1000, res = 200)
      hist(bins.raw.gr$N_percentage, breaks = 50, main = paste0("N% for bins of ", bin_size_Mb, "Mb"),
         ylab = "Counts", xlab = "N%")
    dev.off()
    png(filename = here(ref_dir, paste0("hist_GC_stats_",genome_fasta,"_",bin_size_Mb,"_Mb_bins.png")), width = 1000, height = 1000, res = 200)
      hist(bins.raw.gr$GC_percentage, breaks = 50, main = paste0("GC% for bins of ", bin_size_Mb, "Mb"),
         ylab = "Counts", xlab = "N%")
    dev.off()
  }
}

# In case filtering is needed:
# bins.gr <- bins.raw.gr[which(bins.raw.gr$N_percentage < 50)]
bins.gr <- bins.raw.gr

#### Part 3. Read and process fragment file ####
log_info(paste0("Reading fragment file for ", library_id,"\n"))
fragment_file <- list.files(here(input_dir), pattern = "fragments.tsv.gz$|atac_fragments.tsv.gz$")
fragments <- fread(here(input_dir,fragment_file), showProgress = TRUE, verbose = FALSE, nThread=4)
colnames(fragments) <- c('seqnames','start','end','barcode','pcr')
fragments.gr <- plyranges::as_granges(fragments)
rm(fragments)
invisible(gc())

## For some genomes, we will have to convert RefSeq accession number to chromosome number
refseq_accession <- ref_info$RefSeq.Accn
chr_names <- ref_info$UCSC.style.name
seqs_fragments <- seqlevels(fragments.gr)
seqs_bins <- seqlevels(bins.gr)
if(all(seqs_fragments %in% refseq_accession)){
  log_info("Converting refseq accession numbers of fragment files to UCSC style (chr1, chr2...)")
  seqlevels(fragments.gr) <- chr_names[match(seqlevels(fragments.gr), refseq_accession)]
}else if (all(seqs_fragments %in% chr_names)){
  log_info("All chromosome names are UCSC style.")
}else if (all(seqs_bins %in% seqs_fragments)){
  log_info("Fragment file contains non-USCS style names for unplaced/unscaffold contigs,
           but all chromosomes named as UCSC style.")
}else{
  log_error("Error. The chromosomes do not match RefSeq accessions")
  stop()
}

#### Part 4. Read and process annotation file for the donor ####
log_info(paste0("Processing gex annotation for ", library_id, "\n"))
# We need information of barcode conversion from gex to atac in the file. If that file does not exist, we 
# create a file 
modified.anno_csv <- paste0(path_dir(anno_csv),
                            "/atac_conv_",basename(anno_csv))
if(file_exists(modified.anno_csv)){
  meta.dat <- read.csv(modified.anno_csv, header = TRUE)
}else{
  log_info("Did not find annotation csv file with atac barcode conversion. Creating one. This may take up to 20mins...\n")
  meta.dat <- read.csv(anno_csv, header = TRUE)
  # Read gex-atac barcode conversion file
  # 10X gex and atac barcodes are not the same. However, they provide barcode files, which contains
  # barcodes for both ATAC and gex. For example, row 1783 of ATAC barcode and row 1783 of the gex barcode are from the same gel bead.
  # Refer to this page https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/bam-atac#bam-barcode-translation for details.
  barcode_conv_df <- read.csv(barcode_conv_csv, header = TRUE)
  # select annotations and reformat the barcode field 
  meta.dat <- meta.dat %>%
    dplyr::rename(c(
      library_id_gex = library_id,
      # Some times the first column for barcode_gex is not named. It is replaced with "X"
      barcode_gex = X
    )) %>%
    select(-contains("snn_res"), -orig.clusters) %>%
    mutate(
      # Convert the barcodes to atac. Also the suffix need to map "-1" in the fragment file
      "barcode_atac" = future_mapply(convertGEXtoATAC_bc, barcode_gex, aggr_ident = FALSE),
      # Keep the original gex barcode suffix as a record
      "barcode_atac_aggr" = future_mapply(convertGEXtoATAC_bc, barcode_gex, aggr_ident = TRUE))
  write.csv(meta.dat, file = modified.anno_csv, row.names = FALSE, quote = FALSE)
  log_info(paste0("Created atac barcode converted annotation csv at ", modified.anno_csv, "\n"))
}

# Select only entries related to our library_id
library.id.gex <- str_replace(library_id, pattern = "-LAND|-ATAC", replacement = "-RNA")
meta.dat <- filter(meta.dat, library_id_gex == library.id.gex)

# Check how many annotated barcodes are captured in ATAC fragment file
fragment_bc <- unique(fragments.gr$barcode)
log_info("Check how many annotated cells are found in the fragment file \n")
print(table(meta.dat$barcode_atac %in% fragment_bc))

# Select cells with ATAC fragments
meta.dat <- subset(meta.dat, barcode_atac %in% fragment_bc)

# Save subseted meta data
fwrite(meta.dat, file = here(output_dir,paste0("barcode_anno_", library_id,".csv")),row.names = FALSE)

#### Step 5. Counting fragments in cells ####
# This setp will use GRanges intersect to count fragments within the bins
# Subset fragments to fragments in cells
log_info(paste0("Counting fragments for ", library_id))
fragments.gr <- keepSeqlevels(fragments.gr, seqlevels(bins.gr), pruning.mode = "coarse")
fragments.gr <- fragments.gr[fragments.gr$barcode %in% meta.dat$barcode_atac,]
fragments.gr <- sort(sortSeqlevels(fragments.gr))
invisible(gc())

# Identify the bin for each fragment
fragments.gr <- join_overlap_left(fragments.gr, bins.gr)
invisible(gc())

# Convert to data table and aggregate by bins
dt <- as.data.table(fragments.gr)
invisible(gc())

# Count by bins
counts_by_bin <- dt[, .N, by = .(barcode, region)]
invisible(gc())

# Identify regions with zero counts across cells and assign them a value of zero in final matrix
all_bins <- unique(bins.gr$region)
counted_bins <- unique(counts_by_bin$region)
missing_bins <- all_bins[!(all_bins %in% counted_bins)]

counts_by_bin <- tidyr::pivot_wider(counts_by_bin,
                                    names_from = "region", 
                                    values_from = "N", values_fill = 0)
counts_by_bin[, missing_bins] <- 0
counts_by_bin <- counts_by_bin[, c("barcode", all_bins)]

# Save counts file
fwrite(counts_by_bin, file = here(output_dir, paste0("cell_x_", bin_size_Mb, "_Mb_bins.csv")))
log_success(paste0("Step 1. Finished ", library_id,"! \n"))
