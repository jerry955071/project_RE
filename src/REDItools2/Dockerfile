# Set python:2.7.12 as the base-image
FROM blcdsdockerregistry/reditools2:1.0.0

# Labeling
LABEL version="0.1"

# Args
ARG DEBIAN_FRONTEND=noninteractive

# Install htslib
RUN cd / && \
    wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -vxjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure --prefix=/htslib && \
    make && \
    make install
# Add htslib to path
ENV PATH="/htslib/bin:${PATH}"

# Install samtools
RUN cd / && \
    wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar -vxjf samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    ./configure --prefix=/samtools && \
    make && \
    make install
# Add samtools to path
ENV PATH="/samtools/bin:${PATH}"

# COPY annotate_with_DNA.py /reditools2.0/src/cineca/
COPY extract_coverage.sh /reditools2.0
