# bonefish container: built from this source tree, the binary only on Debian slim.
#   docker build -t bonefish .
#   docker run --rm -v "$PWD":/data -w /data bonefish bonefish -t 16 -o out.bam ref.fa reads.fq.gz
# procps is included because Nextflow needs `ps` inside task containers.
FROM rust:1.98-bookworm AS build
WORKDIR /src
COPY . .
RUN cargo build --release --locked --features jemalloc \
 && strip target/release/bonefish \
 && target/release/bonefish --version

FROM debian:bookworm-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/bonefish /usr/local/bin/bonefish
COPY LICENSE README.md /usr/share/doc/bonefish/
LABEL org.opencontainers.image.title="bonefish" \
      org.opencontainers.image.description="Fast aligner for Oxford Nanopore DNA reads" \
      org.opencontainers.image.source="https://github.com/DylanCosto/bonefish" \
      org.opencontainers.image.licenses="MIT"
CMD ["bonefish", "--help"]
