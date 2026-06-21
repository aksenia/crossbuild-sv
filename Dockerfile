FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    tar \
    bzip2 \
    unzip \
    make \
    gcc \
    g++ \
    default-jre-headless \
    tabix \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && apt-get install -y --only-upgrade liblcms2-2 \
    && rm -rf /var/lib/apt/lists/*

# Python liftover tools
RUN pip install --no-cache-dir CrossMap pysam snakemake && \
    pip install --no-cache-dir --upgrade setuptools

# bcftools with liftover plugin
RUN wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar xjf bcftools-1.22.tar.bz2 && \
    cd bcftools-1.22 && \
    /bin/rm -f plugins/*.{c,h,mk} && \
    wget -P plugins https://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    /bin/rm -f plugins/pgs.mk || true && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make plugins/liftover.so && \
    mkdir -p /usr/local/bin /usr/local/lib/bcftools_plugins && \
    cp bcftools /usr/local/bin/ && \
    cp plugins/liftover.so /usr/local/lib/bcftools_plugins/ && \
    cd .. && rm -rf bcftools-1.22 bcftools-1.22.tar.bz2

ENV BCFTOOLS_PLUGINS=/usr/local/lib/bcftools_plugins

# Picard 3.4.0
RUN mkdir -p /usr/local/lib \
    && wget -q https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar \
       -O /usr/local/lib/picard.jar \
    && printf '#!/bin/sh\nexec java -jar /usr/local/lib/picard.jar "$@"\n' \
       > /usr/local/bin/picard \
    && chmod +x /usr/local/bin/picard

# Transanno (statically linked, arch-aware)
RUN ARCH=$(uname -m) \
    && wget -q "https://github.com/informationsea/transanno/releases/download/v0.4.5/transanno-${ARCH}-unknown-linux-musl-v0.4.5.zip" \
       -O /tmp/transanno.zip \
    && unzip /tmp/transanno.zip -d /tmp/transanno \
    && install -m 755 "/tmp/transanno/transanno-${ARCH}-unknown-linux-musl-v0.4.5/transanno" /usr/local/bin/transanno \
    && rm -rf /tmp/transanno /tmp/transanno.zip

# Build context must be the repo root:
#   docker build -f sv-preprocess/Dockerfile -t crossbuild-sv:latest .
COPY sv-preprocess/liftover/ ./liftover/
COPY sv-preprocess/snake/ ./snake/
COPY preprocess/regions/ ./regions/

RUN chmod +x ./liftover/*.py

RUN groupadd -r crossbuild && useradd -r -g crossbuild crossbuild && \
    mkdir -p /home/crossbuild && \
    chown -R crossbuild:crossbuild /app /home/crossbuild

USER crossbuild

ENV HOME=/home/crossbuild
ENV PATH="/usr/local/bin:${PATH}"

RUN which CrossMap bcftools picard transanno && echo "All four liftover tools verified"

CMD ["python", "--help"]
