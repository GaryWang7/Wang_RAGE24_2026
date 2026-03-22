# This script generates the reference genome for rat GRCr8 with mitochondrial genome from mRatBN7.2, which will be used for both bulk long read and single cell analysis.
# Download GRCr8 fasta and gtf from 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036323735.1/
# The current RefSeq GRCr8 does not include mitonchondrial genome, 
# the GenBank GRCr8 has mitochondrion (CM070413.1) but is not annotated.
# CM070413.1 is only different from the MT sequence in mRatBN7.2 (NC_001665.2) by one mismatch and two insertions
# so we will use mRatBN7.2 mitochondrial genome by concatenating it to the GRCr8 genome.
# Last update: 10/16/2024
# On PMACS:
bsub -Is \
-n 10 \
-M 200000 \
-R "rusage [mem=200000] span[hosts=1]" \
bash

# Download genome fasta:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz

# Uncompress fasta
gunzip -v -f GCF_036323735.1_GRCr8_genomic.fna.gz
gunzip -v -f chrMT.fna.gz

# Download genome gtf:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz 

# Uncompress gtf
gunzip -v -f GCF_036323735.1_GRCr8_genomic.gtf.gz
gunzip -v -f GCF_015227675.2_mRatBN7.2_genomic.gtf.gz

# Extract NC_001665.2 from mRatBN7.2 gtf 
grep NC_001665.2 GCF_015227675.2_mRatBN7.2_genomic.gtf > NC_001665.2.gtf
cat NC_001665.2.gtf | awk '{a[$3]++}END{for(k in a){print k,a[k]}}' # check indeed there are 37 genes 

# Edit NC_001665.2 to modify "unassigned_transcript", which may contain conflicting transcript_id against GRCr8
# e.g. unassigned_transcript_1825 corresponds to different genes in GRCr8 and NC_001665.2
sed 's/transcript_id "unassigned_transcript_/transcript_id "unassigned_transcript_mt/g' NC_001665.2.gtf > NC_001665.2_edit.gtf

# Append mitochondrial genome NC_001665.2 to GRCr8 genome and gtf
cat GCF_036323735.1_GRCr8_genomic.fna chrMT.fna > GCF_036323735.1_GRCr8_mt.fna
cat GCF_036323735.1_GRCr8_genomic.gtf NC_001665.2_edit.gtf > GCF_036323735.1_GRCr8_mt.gtf

#### For Cellranger ####
# Filter GTF
cellranger mkgtf \
GCF_036323735.1_GRCr8_mt.gtf cellranger_filtered_GCF_036323735.1_GRCr8_mt.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lncRNA \
--attribute=gene_biotype:antisense \
--attribute=gene_biotype:IG_LV_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:IG_V_pseudogene \
--attribute=gene_biotype:IG_D_gene \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_J_pseudogene \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_C_pseudogene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:TR_V_pseudogene \
--attribute=gene_biotype:TR_D_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_J_pseudogene \
--attribute=gene_biotype:TR_C_gene

#Run mkref
cellranger mkref \
--genome=GRCr8_mt_cellranger \
--fasta=GCF_036323735.1_GRCr8_mt.fna \
--genes=cellranger_filtered_GCF_036323735.1_GRCr8_mt.gtf \
--ref-version=1.0.1 \
--nthreads=10 \
--memgb=90 \
--uiport=3000 # optional. Haven't tested UI interface on PMACS yet.

#### For cellranger-arc ####
# This takes a long time, so I would run it locally or as a job for PMACS, not in an interactive session.
cellranger-arc mkgtf \
GCF_036323735.1_GRCr8_mt.gtf cellrangerarc_filtered_GCF_036323735.1_GRCr8_mt.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lncRNA \
--attribute=gene_biotype:antisense \
--attribute=gene_biotype:IG_LV_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:IG_V_pseudogene \
--attribute=gene_biotype:IG_D_gene \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_J_pseudogene \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_C_pseudogene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:TR_V_pseudogene \
--attribute=gene_biotype:TR_D_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_J_pseudogene \
--attribute=gene_biotype:TR_C_gene

cellranger-arc mkref --config=GRCr8_with_mt.config --memgb=70 --nthreads=10