# Introduction
This markdown contain works for preprocessing mRNAseq data from DCLK1 over-expression in HEK293T cells. In this experiment, we have 5 groups, each having 6 replicates. We transfected HEK293T cells with DCLK1-L (NM_001330071.2), DCLK1-S (NM_001195416.2), DCLK1-ultrashort (NM_001195430.2) or empty vector plasmid pcDNA3.1(+)-C-eGFP. We also have a group of HEK293T cells without any transfection. After 30 hours of transfection, cells were collected in lysis buffer and sent for bulk mRNAseq. The groups are as follows:
- DCLK1-L: HEK overexpressing DCLK1-L plasmid
- DCLK1-S: HEK overexpressing DCLK1-S plasmid
- DCLK1-US: HEK overexpressing DCLK1-US plasmid
- Backbone: HEK overexpressing empty vector plasmid pcDNA3.1(+)-C-eGFP
- Vehicle: HEK293T cells without any transfection (OptiMEM + lipofectamine only)

## Transfection
Cell were cultured in 24-well plates. 2hr before transfection, the media was replaced with fresh HEK293T complete media without antibiotics. We follow the lipoectamine 3000 protocol for mixing the reagents. For each well, the master mix contains 0.5µg plasmid DNA, 1.5µL Lipofectamine 3000 reagent, 1µL P3000 reagent and 50µL Opti-MEM medium. The mixture was then added to each well in a dropwise manner. The cells were incubated for 48 hours before harvesting.

## RNA collection
In each well, the cells were then harvested in 700 µL Buffer RLT Plus supplied with β-mercaptoethanol (β-ME, at a concentration of 10 µl β-ME per 1 ml of Buffer RLT Plus). 350 µL of the lysate were sent to to extract RNA and construct libraries.

## Library preparation and sequencing
The RNA is extracted from the cells, prepared with NEBNext® Poly(A) mRNA Magnetic Isolation Module to enrich polyA mRNA, and NEBNext® Ultra™ II Directional RNA Library Prep Kit to construct RNA libraries. The libraries are multiplexed with NEBNext® Multiplex Oligos for Illumina® and sequenced on Illumina NovaSeq X platform. The sequencing is done in paired-end mode with 2x150bp read length.

## Notes from Admera:
Isolated RNA sample quality was assessed by High Sensitivity RNA Tapestation (Agilent Technologies Inc., California, USA) and quantified by AccuBlue® Broad Range RNA Quantitation assay (Biotium, California, USA). Paramagnetic beads coupled with oligo d(T)25 are combined with total RNA to isolate poly(A)+ transcripts based on NEBNext® Poly(A) mRNA Magnetic Isolation Module manual (New England BioLabs Inc., Massachusetts, USA). Prior to first strand synthesis, samples are randomly primed (5´ d(N6) 3´ [N=A,C,G,T]) and fragmented based on manufacturer’s recommendations. The first strand is synthesized with the Protoscript II Reverse Transcriptase with a longer extension period, approximately 40 minutes at 42⁰C. All remaining steps for library construction were used according to the NEBNext® Ultra™ II Directional RNA Library Prep Kit for Illumina® (New England BioLabs Inc., Massachusetts, USA). Final libraries quantity was assessed by Qubit 2.0 (ThermoFisher, Massachusetts, USA) and quality was assessed by TapeStation D1000 ScreenTape (Agilent Technologies Inc., California, USA). Final library size was about 430bp with an insert size of about 300bp. Illumina® 8-nt dual-indices were used. Equimolar pooling of libraries was performed based on QC values and sequenced on an Illumina® NovaseqX plus platform(Illumina, California, USA) with a read length configuration of 150 PE for 40M PE reads per sample (20M in each direction).

# QC and trim data
**Filtering notes:**
- Adaptor: The libraries are multiplexed with NEBNext® Multiplex Oligos (Dual Index) for Illumina®. According to the [FAQ from NEB](https://www.neb.com/en-us/faqs/2021/01/15/what-sequences-need-to-be-trimmed-for-nebnext-libraries-that-are-sequenced-on-an-illumina-instrument),
  ```
  The NEBNext libraries for Illumina resemble TruSeq libraries and can be trimmed like TruSeq:

    Adaptor Read1   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    Adaptor Read2   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   
  ```
  This is exactly the same as the [Illumina TruSeq adapters](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/UDIndexes.htm).
- Quality filtering is enabled by default in fastp. The default quality filtering consider reads with average quality score (Phred score) lower than 15, which corresponds to a base call accuracy of 97%, as **unqualified bases**. The default N-base limit is 5, and percents of unqualitified bases is 40%. 
- Length filtering: reads shorter than 15bp is by default removed.
- Poly-G trimming: enabled to detect a minimum of 10-based polyG by default.

Below is the command to trim the data using fastp v0.24.1
```sh
# conda environment 
mamba activate general

# Define the directories
data_dir=/home/gary/garyw/RAGE24/DCLK1_HEK293T/Fastq_overexpression
fastq_dir=${data_dir}/untrimmed_raw
fastq_trimmed_dir=${data_dir}/trimmed
fastp_report_dir=${data_dir}/fastp_report
mkdir -p $fastq_trimmed_dir
mkdir -p $fastp_report_dir

# Get all the raw sample names into one file. Note for metadata, we changed the "customer ID" column for vehicle samples to match the second experiment.
file_list=${data_dir}/raw_sample_list.txt
metadata=${data_dir}/metadata_edited.csv
file_list_modified=${data_dir}/raw_sample_list_modified.txt

cd $fastq_dir
ls -1 *.fastq.gz | awk -F'_' '{print $1"_"$2"_"$3","$1}' | sort -u > $file_list

# Map from Admera to Sample ID
awk -F',' 'NR==FNR {meta[$3] = $2; next} {print $0, meta[$2]}' OFS=',' $metadata $file_list > $file_list_modified

cd $data_dir
while IFS=',' read -r sample_name_raw Admera_name sample_name  
do
  echo "Processing $sample_name_raw, which is $sample_name"
  in_1=${fastq_dir}/${sample_name_raw}_R1_001.fastq.gz
  in_2=${fastq_dir}/${sample_name_raw}_R2_001.fastq.gz
  out_1=${fastq_trimmed_dir}/${sample_name}_R1.fastq.gz
  out_2=${fastq_trimmed_dir}/${sample_name}_R2.fastq.gz
  fastp --in1 $in_1 \
  --in2 $in_2 \
  --out1 $out_1 \
  --out2 $out_2 \
  --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --html $fastp_report_dir/${sample_name}.html \
  --json $fastp_report_dir/${sample_name}.json \
  --report_title $sample_name \
  --thread 12
done < $file_list_modified

# Generate a multiqc report for all the fastp reports
cd $fastp_report_dir
mkdir multiqc_report
multiqc . --outdir multiqc_report --title DCLK1_HEK_Fastp_report
```
# Reference genome and annotation
Here we use GRCh38p14 genome from Ensembl release 113. The annotation is from GENCODE release 47. We made some modifications to the files based on CellRanger's requirements. No filtering is done on fastq and annotation files though. As we are introducing GFP plasmid, we add custom entries to detect the GFP gene. Follow the guide line from [10X](https://www.10xgenomics.com/support/software/cell-ranger/7.2/tutorials/cr-tutorial-mr#runmkref) to add fasta and gtf entries. Find the format of those files in another [10X tutorial](https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/inputs/cr-3p-references).

An analysis without mapping the GFP protein is in the result folder as well.

## a) Download and modify the reference genome and annotation files
```sh
# Download fasta and gtf
fasta_url="http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz"
gtf_in="gencode.v47.primary_assembly.annotation.gtf"

if [ ! -f "$fasta_in" ]; then
    wget -qO- "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    wget -qO- "$gtf_url" | zcat > "$gtf_in"
fi

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"
```

## b) supply the custom entries to the fasta and gtf files
The fasta sequence of eGFP. It is 717 bp. Save it as GFP.fa
```sh
>GFP
ATGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCACGGCATGGACGAGCTGTACAAGTAA
```
Now make a custom GTF for EGFP

```sh
cd /mnt/d/garyw/RAGE24/reference/GRCh38p14
# Note this has been modified from the tutorial to match GENCODE format.
echo -e 'GFP\tunknown\texon\t1\t717\t.\t+\t.\tgene_id "GFP"; transcript_id "GFP"; gene_type "protein_coding"; gene_name "GFP"; transcript_type "protein_coding";' > GFP.gtf
```

We then create a separate fasta and gtf file.
```sh
# copy fasta file
cp Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified Homo_sapiens.GRCh38.dna.primary_assembly_modified_GFP.fa

# Append EGFP.fa to the fasta file
cat GFP.fa >> Homo_sapiens.GRCh38.dna.primary_assembly_modified_GFP.fa

# Confirm the GFP entry was added to the FASTA file. Use the grep ">" command to search for lines with ">" character
grep ">" Homo_sapiens.GRCh38.dna.primary_assembly_modified_GFP.fa

# copy the gtf file and append GFP.gtf to it
cp gencode.v47.primary_assembly.annotation.gtf.modified gencode.v47.primary_assembly.annotation_modified_GFP.gtf
cat GFP.gtf >> gencode.v47.primary_assembly.annotation_modified_GFP.gtf

# Check that gtf has been appended
tail gencode.v47.primary_assembly.annotation_modified_GFP.gtf
```

# Mapping and read counting with Subread
The steps (alignment and quantification) is performed in "step0_alignment.R" to generate a count matrix. 