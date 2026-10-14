# Upstream rammap documentation

These documents come from [rammap](https://github.com/jwanglab/rammap), the Rust port of minimap2
that bonefish was forked from. They are kept because the index format, SIMD kernels and much of
the mapping pipeline are rammap's, but they describe rammap, not bonefish:

- they cover features bonefish removed (spliced/RNA alignment, paired-end short reads, PacBio and
  assembly presets, homopolymer-compressed indexes) and minimap2-identical output, which bonefish
  does not produce;
- paths and commands may say `rammap` / `rammap-core` (now `bonefish` / `bonefish-core`);
- `performance.md` measures rammap 1.0.0 on workloads bonefish does not support.

For bonefish itself see the [README](../../README.md) and [benchmarks](../benchmarks.md).
