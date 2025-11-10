# Dockerfile for SCAR (Single-Cell ATAC Router)
FROM ubuntu:22.04

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Basic utilities and build tools
RUN apt-get update && \
    apt-get install -y \
    git \
    build-essential \
    cmake \
    wget \
    pkg-config \
    ca-certificates \
    autoconf \
    autoheader \
    automake \
    libtool \
    && apt-get clean

# Install compression libraries required by htslib
RUN apt-get update && \
    apt-get install -y \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean

# Show installed versions for debugging
RUN g++ --version && \
    cmake --version

# Clone SCAR repository
WORKDIR /opt
RUN git clone https://github.com/davidebolo1993/scar.git
WORKDIR /opt/scar

# Build SCAR using CMake
RUN mkdir -p build
WORKDIR /opt/scar/build
RUN cmake .. && \
    make -j4

# Add SCAR to PATH
ENV PATH="/opt/scar/build:${PATH}"

# Set working directory for data
WORKDIR /data

# Default entrypoint
ENTRYPOINT ["/opt/scar/build/scar"]

# Usage examples:
# 
# Build:
#   docker build -t scar .
#
# Filter step:
#   docker run --rm -v $PWD:/data scar filter \
#     -i metadata.tsv \
#     -o metadata.filtered.tsv \
#     -c celltypes_LEV1 \
#     -m 10 \
#     -p peak_matrix.tsv \
#     -b celltype_beds/
#
# Split step:
#   docker run --rm -v $PWD:/data scar split \
#     -i metadata.filtered.tsv \
#     -o output_dir/ \
#     -c celltypes_LEV1 \
#     -d

