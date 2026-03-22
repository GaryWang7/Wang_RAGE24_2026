num_threads=<num_threads>
ref_genome=<ref_genome>
ref_gtf=<ref_gtf>
collapsed_gff=<collapsed_gff>
fl_count=<fl_count>
classification_dir=<classification_dir>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

# Create a temporary directory for the clustering process
# This is to avoid a bug in temporary directory handling in isoseq cluster2
# instead of creating /home/garyw7/temp.bam, it creates /home/garyw7temp.bam, which leads to an error.
mkdir -p /home/garyw7/temp/
cd /home/garyw7/temp/

# Classify and filtering
pigeon classify \
    $collapsed_gff \
    $ref_gtf \
    $ref_genome \
    --out-dir $classification_dir \
    --num-threads $num_threads \
    --log-level INFO \
    --fl $fl_count

pigeon filter $classification_dir/collapsed_classification.txt \
    --isoforms $collapsed_gff \
    --polya-percent 0.8 \
    --num-threads $num_threads \
    --log-level INFO