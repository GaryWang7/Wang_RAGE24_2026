# This scripts generates both FLNC (full-length non-chimeric, poly-A trimmed) reads according to IsoSeq
# and FLNC-polyA reads for IsoQuant.
input_bam=<input_bam>
flnc_bam=<flnc_bam>
flnc_polyA_bam=<flnc_polyA_bam>
num_threads=<num_threads>
isoseq_primer="/home/garyw7/RAGE_longread/reference/IsoSeq_v2_primers_12.fasta"

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Run IsoSeq refine. Remove chimeric reads and trim poly-A tails.
# --require-polya. This filters for FL reads that have a poly(A) tail with at least 20 base pairs (--min-polya-length) and removes identified tail.
isoseq refine -v --log-level INFO \
    --require-polya \
    --num-threads ${num_threads} \
    ${input_bam} ${isoseq_primer} ${flnc_bam}

# Run IsoSeq refine again. But this time only removes chimeric reads
# Ideally we would also like to filter for reads containing poly-A tails. However, if we identify poly-A in isoseq refine, the algorithm will trim the poly-A tails. We avoid an additional layer of complexity here.
isoseq refine -v --log-level INFO \
    --num-threads ${num_threads} \
    ${input_bam} ${isoseq_primer} ${flnc_polyA_bam}
