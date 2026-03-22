input_bam=<input_bam>
output_bam=<output_bam>

sample_id=<sample_id>
new_header=<new_header>

num_threads=<num_threads>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Extract and rewrite header in sam format
samtools view -H "$input_bam" | \
  awk -v sm="$sample_id" '/^@RG/ { sub(/SM:[^[:space:]]+/, "SM:" sm) }1' \
  > "$new_header"

# Rebuild the BAM files with the modified headers
samtools reheader $new_header $input_bam > $output_bam

# Index the modified BAM file
pbindex $output_bam \
    --num-threads $num_threads \
    --log-level INFO
