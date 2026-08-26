FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG PLINK2_VERSION=2.0.0-a.6.9
ARG NXF_VER=26.04.6

RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools ca-certificates curl gzip openjdk-21-jre-headless python3 python-is-python3 python3-pip unzip \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL -o /tmp/plink2.zip \
      "https://www.cog-genomics.org/static/bin/plink${PLINK2_VERSION}_linux_avx2.zip" \
    && unzip /tmp/plink2.zip -d /usr/local/bin \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip

ENV NXF_VER=${NXF_VER}
RUN cd /usr/local/bin && curl -fsSL https://get.nextflow.io | bash && chmod +x /usr/local/bin/nextflow

WORKDIR /work
COPY pyproject.toml ./
COPY scripts ./scripts
RUN python3 -m pip install --break-system-packages --no-cache-dir .
COPY main.nf nextflow.config README.md ./
COPY config ./config

ENTRYPOINT ["nextflow", "run", "main.nf", "-profile", "local", "--outdir", "/work/results"]
