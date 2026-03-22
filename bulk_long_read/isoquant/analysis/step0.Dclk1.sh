# Run the script below to generate reads specific to Dclk1 (bed file and fastq files)
# Run the following commands to get the fasta file:
mamba activate ucsc
proj_dir=/home/gary/v/RAGE24_LongRead/bulk_long_read
data_dir=/home/gary/v/RAGE24_LongRead/bulk_long_read/data/IsoQuant_ambiguous
ref_fa=$proj_dir/reference/pigeon/GRCr8_mt.fa
ref_gtf=$proj_dir/reference/pigeon/pigeon_GRCr8_mt.gtf
Dclk1_gtf=$data_dir/Dclk1/Dclk1_IsoQuant.gtf
Dclk1_bed=$data_dir/Dclk1/Dclk1_IsoQuant.bed
Dclk1_transcript_fa=$data_dir/Dclk1/Dclk1_IsoQuant.fa
#Use samtools to index the genome first
samtools faidx $ref_fa
gtftoGenePred -genePredExt -ignoreGroupsWithoutExons \
$Dclk1_gtf stdout \
| genePredToBed stdin $Dclk1_bed
bedtools getfasta -fi $ref_fa -bed $Dclk1_bed -fo $Dclk1_transcript_fa -s -split -name

# Then use the ORFfinder to predict the ORF from those sequences.