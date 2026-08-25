# syntax=docker/dockerfile:1

# ============================================================
# Build Transanno
# ============================================================
FROM rust:1.97-slim AS transanno-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    ca-certificates \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 --branch v0.4.5 \
    https://github.com/informationsea/transanno.git \
    /src/transanno \
    && cd /src/transanno \
    && cargo build --release --locked


# ============================================================
# Main pipeline image
# ============================================================
FROM python:3.11-slim

WORKDIR /app

ARG BCFTOOLS_VERSION=1.22

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    bzip2 \
    make \
    gcc \
    g++ \
    default-jre-headless \
    tabix \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*


# ============================================================
# Python tools
# ============================================================
RUN python -m pip install --no-cache-dir --upgrade \
    pip setuptools wheel \
    && python -m pip install --no-cache-dir \
    CrossMap \
    pysam \
    snakemake


# ============================================================
# bcftools 1.22 + liftover plugin
# ============================================================
RUN wget -q \
    "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
    -O /tmp/bcftools.tar.bz2 \
    && tar xjf /tmp/bcftools.tar.bz2 -C /tmp \
    && cd "/tmp/bcftools-${BCFTOOLS_VERSION}" \
    && rm -f plugins/liftover.c \
    && wget -q \
       "https://raw.githubusercontent.com/freeseek/score/master/liftover.c" \
       -O plugins/liftover.c \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make plugins/liftover.so \
    && mkdir -p /usr/local/lib/bcftools_plugins \
    && install -m 755 bcftools /usr/local/bin/bcftools \
    && install -m 755 plugins/liftover.so \
       /usr/local/lib/bcftools_plugins/liftover.so \
    && rm -rf \
       "/tmp/bcftools-${BCFTOOLS_VERSION}" \
       /tmp/bcftools.tar.bz2

ENV BCFTOOLS_PLUGINS=/usr/local/lib/bcftools_plugins


# ============================================================
# Picard 3.4.0
# ============================================================
RUN mkdir -p /usr/local/lib \
    && wget -q \
       "https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar" \
       -O /usr/local/lib/picard.jar \
    && printf '#!/bin/sh\nexec java -jar /usr/local/lib/picard.jar "$@"\n' \
       > /usr/local/bin/picard \
    && chmod +x /usr/local/bin/picard


# ============================================================
# Transanno binary from builder stage
# ============================================================
COPY --from=transanno-builder \
    /src/transanno/target/release/transanno \
    /usr/local/bin/transanno

RUN chmod +x /usr/local/bin/transanno


# ============================================================
# Project files
#
# Docker build context = crossbuild-sv root
# ============================================================
COPY liftover/ ./liftover/
COPY snake/ ./snake/
COPY regions/ ./regions/

RUN chmod +x ./liftover/*.py \
    && test -x ./liftover/merge_tools.py \
    && test -x ./liftover/normalize_sv_coordinates.py


# ============================================================
# Runtime mount points
#
# Reference/ will be mounted at /ref
# pilot input directory will be mounted at /data
# output directory will be mounted at /results
# ============================================================
RUN mkdir -p \
    /ref \
    /data \
    /results \
    /home/crossbuild


# ============================================================
# Non-root user
# ============================================================
RUN groupadd -r crossbuild \
    && useradd -r -g crossbuild -d /home/crossbuild crossbuild \
    && chown -R crossbuild:crossbuild \
       /app \
       /home/crossbuild \
       /results

USER crossbuild

ENV HOME=/home/crossbuild
ENV PATH="/usr/local/bin:${PATH}"


# ============================================================
# Verify required tools during image build
# ============================================================
RUN command -v CrossMap \
    && command -v bcftools \
    && command -v bgzip \
    && command -v tabix \
    && command -v picard \
    && command -v transanno \
    && bcftools plugin -l | grep -qx liftover \
    && echo "All liftover tools and bcftools liftover plugin verified"


CMD ["snakemake", "--help"]
