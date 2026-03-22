num_threads=<num_threads>

fofn=<fofn>
mapped_bam=<mapped_bam>
collapsed=<collapsed>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Create a temporary directory for the clustering process
# This is to avoid a bug in temporary directory handling in isoseq cluster2
# instead of creating /home/garyw7/temp.bam, it creates /home/garyw7temp.bam, which leads to an error.
mkdir -p /home/garyw7/temp/
cd /home/garyw7/temp/

# Run collapse

isoseq collapse --do-not-collapse-extra-5exons \
    --num-threads ${num_threads} \
    --log-level INFO \
    ${mapped_bam} \
    ${fofn} \
    ${collapsed}