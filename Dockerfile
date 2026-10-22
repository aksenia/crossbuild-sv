# Build context is the repo root:  docker build -t crossbuild-sv:latest .

# Transanno from its tagged source (no prebuilt binary for every arch).
FROM rust:1.97-slim AS transanno-builder
RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates pkg-config zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch v0.4.5 https://github.com/informationsea/transanno.git /src/transanno \
    && cd /src/transanno && cargo build --release --locked

FROM python:3.11-slim
WORKDIR /app
ARG BCFTOOLS_VERSION=1.22

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget bzip2 make gcc g++ git default-jre-headless tabix \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Versions are pinned: an unpinned CrossMap once rejected most sequence-resolved
# DELs as Fail(REF==ALT). setuptools<81 keeps pkg_resources for CrossMap;
# svdb and truvari build with Cython/numpy without build isolation.
RUN python -m pip install --no-cache-dir --upgrade pip "setuptools==80.10.2" wheel Cython numpy \
    && python -m pip install --no-cache-dir --no-build-isolation \
        CrossMap==0.7.3 pysam snakemake==8.30.0 svdb==2.12.0 truvari==5.3.0

# bcftools with the freeseek liftover plugin
RUN wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
        -O /tmp/bcftools.tar.bz2 \
    && tar xjf /tmp/bcftools.tar.bz2 -C /tmp && cd "/tmp/bcftools-${BCFTOOLS_VERSION}" \
    && wget -q https://raw.githubusercontent.com/freeseek/score/master/liftover.c -O plugins/liftover.c \
    && ./configure --prefix=/usr/local && make -j"$(nproc)" && make plugins/liftover.so \
    && install -m 755 bcftools /usr/local/bin/bcftools \
    && install -D -m 755 plugins/liftover.so /usr/local/lib/bcftools_plugins/liftover.so \
    && rm -rf "/tmp/bcftools-${BCFTOOLS_VERSION}" /tmp/bcftools.tar.bz2
ENV BCFTOOLS_PLUGINS=/usr/local/lib/bcftools_plugins

# Picard 3.4.0
RUN wget -q https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar \
        -O /usr/local/lib/picard.jar \
    && printf '#!/bin/sh\nexec java $JAVA_OPTS -jar /usr/local/lib/picard.jar "$@"\n' > /usr/local/bin/picard \
    && chmod +x /usr/local/bin/picard

COPY --from=transanno-builder /src/transanno/target/release/transanno /usr/local/bin/transanno

COPY liftover/ ./liftover/
COPY snake/ ./snake/
COPY regions/ ./regions/

RUN groupadd -r crossbuild && useradd -r -g crossbuild -d /home/crossbuild crossbuild \
    && mkdir -p /home/crossbuild && chown -R crossbuild:crossbuild /app /home/crossbuild
USER crossbuild
ENV HOME=/home/crossbuild

RUN which CrossMap bcftools picard transanno svdb truvari tabix \
    && bcftools plugin -l | grep -qx liftover

CMD ["snakemake", "--help"]
