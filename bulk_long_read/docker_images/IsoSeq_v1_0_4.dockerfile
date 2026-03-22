# docker build --tag=garywang7/isoseq:1.0.4 -f IsoSeq_v1_0_4.dockerfile .
# docker push garywang7/isoseq:1.0.4
# singularity pull -F docker://garywang7/isoseq:1.0.4
# To tes: singularity shell isoseq_1.0.4.sif
FROM mambaorg/micromamba:2.3.1

# create the environment with all needed packages in one step and clean up
RUN micromamba create -n isoseq_all
RUN micromamba install -n isoseq_all -c conda-forge \
    libgcc libstdcxx-ng
RUN micromamba install -n isoseq_all -c bioconda \
      isoseq lima pbtk pbmm2 pbskera pbpigeon=1.4.0 \
      # cleans the cache
    && micromamba clean --all --yes
RUN micromamba install -n isoseq_all -c conda-forge -c bioconda python=3.8 isoquant=3.7.1

# make isoseq_all the default active env at container runtime
ENV ENV_NAME=isoseq_all