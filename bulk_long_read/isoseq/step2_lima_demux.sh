input_bam=<input_bam>
output_bam=<output_bam>
num_threads=<num_threads>
isoseq_primer="/home/garyw7/RAGE_longread/reference/IsoSeq_v2_primers_12.fasta"

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

lima $input_bam $isoseq_primer $output_bam \
    --isoseq --peek-guess \
    -j $num_threads --log-level INFO \
    --overwrite-biosample-names