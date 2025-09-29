# Dockerfile
FROM mambaorg/micromamba:1.5.8

# Create an env with pinned versions (fast & reproducible)
# - We keep bash available for Nextflow's default task shell
# - Add gawk for awk-compatible scripts
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    samtools=1.20 \
    bcftools=1.20 \
    bedtools=2.31.1 \
    minimap2=2.28 \
    multiqc=1.31 \
    gawk \
  && micromamba clean -a -y

# Ensure tools are on PATH for both root and runtime user
ENV PATH=/opt/conda/bin:$PATH

# (Optional) slim down image a bit
RUN rm -rf /opt/conda/pkgs/*

# Workdir and default entry
WORKDIR /data
SHELL ["/bin/bash", "-c"]
# Remove problematic ENTRYPOINT - let Nextflow handle command execution
# ENTRYPOINT ["/bin/bash", "-lc"]
