config=<config>
ref_genome=<ref_genome>
ref_gtf=<ref_gtf>
output_dir=<output_dir>
num_threads=<num_threads>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

cd $output_dir

isoquant.py \
    --output $output_dir \
    --high_memory \
    --threads $num_threads \
    --resume