# Container for the pipeline. All software comes from environment.yml (pinned
# versions). The pipeline's own scripts are NOT copied in: Nextflow supplies
# bin/ at run time, so the image only changes when the environment does.
#
# Build (the image is amd64; on Apple Silicon this runs under emulation):
#   docker build --platform linux/amd64 -t <repo>:<version> .
FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Make the environment's tools available however the container is started
# (Nextflow, docker run, Singularity), not only via the micromamba entrypoint
ENV PATH=/opt/conda/bin:$PATH
