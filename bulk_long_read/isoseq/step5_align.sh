num_threads=<num_threads>
ref_genome=<ref_genome>
clustered_bam=<clustered_bam>
mapped_bam=<mapped_bam>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Create a temporary directory for the clustering process
# This is to avoid a bug in temporary directory handling in isoseq cluster2
# instead of creating /home/garyw7/temp.bam, it creates /home/garyw7temp.bam, which leads to an error.
mkdir -p /home/garyw7/temp/
cd /home/garyw7/temp/

# Run alignment
pbmm2 align --preset ISOSEQ \
    --sort \
    --sort-memory 12G \
    --sort-threads ${num_threads} \
    --num-threads ${num_threads} \
    --log-level INFO \
    ${ref_genome} \
    ${clustered_bam} \
    ${mapped_bam}