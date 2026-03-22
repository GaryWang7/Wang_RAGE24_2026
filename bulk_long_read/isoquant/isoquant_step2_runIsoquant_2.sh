# Specified without fl_data (as 5' end may be degraded) and a transcript counting strategy to include 
# ambiguous reads by splitting between features with equal weights 
# (e.g. 1/2 when a read is assigned to 2 features simultaneously)
config=<config>
ref_genome=<ref_genome>
ref_gtf=<ref_gtf>
output_dir=<output_dir>
num_threads=<num_threads>

# Initialize environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoseq_all

cd $output_dir

isoquant.py -d pacbio_ccs \
    --yaml $config \
    --reference $ref_genome \
    --genedb $ref_gtf \
    --complete_genedb \
    --output $output_dir \
    --report_novel_unspliced false \
    --transcript_quantification with_ambiguous \
    --count_exons \
    --check_canonical \
    --sqanti_output \
    --threads $num_threads