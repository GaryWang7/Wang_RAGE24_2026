fofn=<fofn>
output_bam=<output_bam>
num_threads=<num_threads>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Create a temporary directory for the clustering process
# This is to avoid a bug in temporary directory handling in isoseq cluster2
# instead of creating /home/garyw7/temp.bam, it creates /home/garyw7temp.bam, which leads to an error.
mkdir -p /home/garyw7/temp/
cd /home/garyw7/temp/

# Run IsoSeq cluster
isoseq cluster2 --num-threads ${num_threads} \
    --log-level INFO \
    ${fofn} \
    ${output_bam}