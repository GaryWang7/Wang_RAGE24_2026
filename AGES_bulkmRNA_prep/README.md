# Introduction
This pipline downloads bulk RNA-seq data from SRA database. The data is from this paper: Age-Related Gene Expression Signature in Rats Demonstrate Early, Late, and Linear Transcriptional Changes from Multiple Tissues [link](https://www.sciencedirect.com/science/article/pii/S2211124719310915)
- SRA accession: SRA: PRJNA516151

# Data download
On SRA Run Selector, select accession number PRJNA516151. Select data that contains "kidney" and download accession list. Follow the instruction to download the data.
```sh
mamba activate general
file_list=/home/gary/v/AgeRat/Data/Accession/SRR_Acc_List.txt
fastq_dir=/home/gary/v/AgeRat/Data/Fastq

while read -r line
do
  echo "Donwloading $line"
  fasterq-dump --split-files --verbose --progress --threads 10 --outdir $fastq_dir $line
done < $file_list

# zip all the files with all available cores as fasterq-dump does not support gzip
cd $fastq_dir
zstd --format=gzip -T0 --verbose --rm *.fastq
```

# QC and trim data
Filtering notes:
- Adaptor: The RNA libraries were prepared using the Illumina TruSeq Stranded Total RNASample Preparation protocol with the Ribo-Zero Gold Kit. This means that the TruSeq adapters are used.
- Quality filtering is enabled by default in fastp. The default quality filtering consider reads with average quality score (Phred score) lower than 15, which corresponds to a base call accuracy of 97%, as **unqualified bases**. The default N-base limit is 5, and percents of unqualitified bases is 40%. 
- Length filtering: reads shorter than 15bp is by default removed.
- Poly-G trimming: enabled to detect a minimum of 10-based polyG by default.


```sh
data_dir=/home/gary/v/AgeRat/Data
fastq_trimmed_dir=/home/gary/v/AgeRat/Data/Fastq_trimmed
fastp_report_dir=/home/gary/v/AgeRat/Data/Fastp_report
cd $data_dir

while read -r line
do
  echo "Processing sample $line"
    fastp --in1 $fastq_dir/${line}_1.fastq.gz \
    --in2 $fastq_dir/${line}_2.fastq.gz \
    --out1 $fastq_trimmed_dir/${line}_trimmed_1.fastq.gz \
    --out2 $fastq_trimmed_dir/${line}_trimmed_2.fastq.gz \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --html $fastp_report_dir/${line}.html \
    --json $fastp_report_dir/${line}.json \
    --report_title $line \
    --thread 10
done < $file_list

cd $fastp_report_dir
mkdir multiqc_report
multiqc . --outdir multiqc_report --title AGES_Fastp_report
```
# Reference transcriptome
Here we use the same reference as RAGE24 data, which is GRCr8 with mitochondrial genome added. 

# Option 1: Mapping and read counting with Subread
The steps (alignment and quantification) is performed in "AGES_alignment.R" to generate a count matrix.

# Option 2: Counting exon-level reads with Salmon
There is an issue of multi-mapping if we would like to quantify exon usage. Subread is not good in handling this but Salmon can take care of this. However, here we did not pursue this option as short read sequencing is not ideal for exon usage quantification.

# Resources:
  - Harvard Chan Bioinformatics Core (HBC) traning [link](https://hbctraining.github.io/main/)
    - Count normalization [link](https://hbctraining.github.io/Intro-to-DGE/lessons/02_DGE_count_normalization.html)