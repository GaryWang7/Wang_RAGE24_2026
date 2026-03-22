input_bam=<input_bam>
output_prefix=<output_prefix>
num_threads=<num_threads>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Run IsoSeq bam2fastq to convert HiFi reads to fastq:
bam2fastq $input_bam \
    --output $output_prefix \
    --num-threads $num_threads \
    --with-biosample-prefix