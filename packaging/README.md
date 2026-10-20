# Packaging

- `bioconda/`: the Bioconda recipe (`meta.yaml`, `build.sh`), as submitted to
  [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) under `recipes/bonefish/`. For a
  new release, update `version` and `sha256` (of
  `https://github.com/DylanCosto/bonefish/archive/refs/tags/vX.Y.Z.tar.gz`); Bioconda's bot usually
  opens that update itself. Local test (Bioconda's CI supplies the `stdlib('c')`
  variant; locally pass it in a config file):
  `printf 'c_stdlib:\n  - sysroot\nc_stdlib_version:\n  - "2.17"\n' > cbc.yaml && conda-build packaging/bioconda -m cbc.yaml -c conda-forge -c bioconda`.
- `../Dockerfile`: container image built from the source tree (`docker build -t bonefish .`).

Both were checked on 16 Oct 2026 (v0.1.6): the conda-built and the container binaries give output
byte-identical to the release binary on the genome-wide HG001 sample (MD5
`9bc7ff3e9f498cc1530b2b1d78c47e10`), and the container's sorted-BAM index matches `samtools index`.
