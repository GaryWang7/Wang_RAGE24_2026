input_bam=<input_bam>
output_bam=<output_bam>
num_threads=<num_threads>
adapter=/home/garyw7/RAGE_longread/reference/mas8_primers.fasta

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

skera split $input_bam $adapter $output_bam -j $num_threads --log-level INFO 