> **Note:** inherited from upstream [rammap](https://github.com/jwanglab/rammap) and kept for reference. It describes rammap (all minimap2 presets, minimap2-identical output), not bonefish; paths may still say `rammap`/`rammap-core` (now `bonefish`/`bonefish-core`).

# WebAssembly Support

rammap compiles to WebAssembly with SIMD128 support, running the full alignment
pipeline (index build, seeding, chaining, DP extension, output formatting) in
the browser or via WASI runtimes like wasmtime.

Two build modes:
- **Single-threaded** — stable Rust, works in all browsers, no special headers
- **Multi-threaded** — nightly Rust, rayon via `wasm-bindgen-rayon`, requires `SharedArrayBuffer`

---

## Quick Start

```bash
# Single-threaded build (stable Rust)
bash scripts/build-wasm.sh
python3 -m http.server 8080
# Open http://localhost:8080/web/

# Multi-threaded build (nightly Rust)
bash scripts/build-wasm-threads.sh
python3 scripts/serve.py
# Open http://localhost:8080/web/
```

---

## Architecture

### WASM-Specific Code

| File | Purpose |
|------|---------|
| `src/align/wasm_lib.rs` | wasm-bindgen entry points (`align_wasm_full`, `force_align_wasm`) |
| `src/align/dp/common.rs` | SIMD128 compatibility layer (34 SSE intrinsic → WASM v128 mappings) |
| `src/align/chain_simd.rs` | WASM SIMD128 chaining (4-wide, mirrors NEON) |
| `web/index.html` | Browser demo UI |
| `pkg/worker.js` | Web Worker for background alignment |
| `scripts/build-wasm.sh` | Single-threaded build script |
| `scripts/build-wasm-threads.sh` | Multi-threaded build script |
| `scripts/serve.py` | Dev server with COOP/COEP headers |

### SIMD128 Coverage

All DP kernels (single-affine, dual-affine, splice) compile for WASM SIMD128
via a compatibility layer that maps SSE intrinsic names to `core::arch::wasm32`
equivalents. The SSE/SSE4.1 macros are instantiated for WASM unchanged — the
compat layer handles the translation at the function level.

Chaining also has a native WASM SIMD128 implementation (4-wide `i32x4`/`f32x4`
batch scoring), dispatched automatically on `wasm32`.

| Component | SIMD128 | Notes |
|-----------|---------|-------|
| Single-affine DP | Yes | Via SSE compat layer |
| Dual-affine DP | Yes | Via SSE compat layer |
| Splice-aware DP | Yes | Via SSE compat layer |
| Lightweight i16 SW | Yes | Native WASM intrinsics |
| Chaining | Yes | Native `i32x4`/`f32x4`, 4-wide |
| Index sort | No | Scalar (no SIMD sort for WASM) |

### Threading Model

| Build | Parallelism | Mechanism |
|-------|-------------|-----------|
| Single-threaded | None | Sequential sketch, chain, align |
| Multi-threaded (browser) | rayon thread pool | `wasm-bindgen-rayon` → Web Workers + SharedArrayBuffer |
| WASI (wasmtime) | None | Single-threaded only (WASI threads not yet supported) |

When multi-threaded, rayon's `par_iter` distributes query reads across worker
threads. The index build uses parallel bucket sort. The thread count is
configurable in the demo UI (defaults to 4, max = system threads via `navigator.hardwareConcurrency`).

---

## Reference size ceiling

The browser-WASM demo can handle FASTA references up to roughly **1.5 GB of
bases** (or ~0.5 GB compressed for gzipped input). Above that, the index is
too large to fit in WebAssembly's 4 GB linear-memory cap and the build fails
partway through with a trap.

The steady-state memory needed to index an N-base reference is:

| Component | Size |
|---|---|
| Packed 4-bit reference | N / 2 bytes |
| Bucket entries (`~N/w × 16 B` for w=10 default) | ~1.6 × N bytes |
| Working memory | ~0.2 GB |

For hs1 (T2T human, 3.1 Gb): packed ≈ 1.5 GB + buckets ≈ 2.6 GB ≈ 4.3 GB of
index, before counting any working memory. That exceeds the WASM32 ceiling
even under perfect allocation. Use the native rammap CLI for genome-scale
references.

The demo's drop-zone surfaces a warning when a reference file is dropped that
will likely exceed the ceiling.

---

## API Reference

### `AlignSession` (streaming, recommended)

Class-based API for multi-GB inputs. Bytes flow through streaming FASTA/FASTQ
parsers; the reference is packed-and-dropped as records complete, queries are
aligned one read at a time. Peak WASM memory stays bounded regardless of
total input size — though still subject to the 4 GB linear-memory cap (see
above).

```js
const session = new AlignSession(preset, output_sam, output_cigar);
session.reserve_ref_bases(BigInt(uncompressed_ref_bytes));  // optional hint
for await (const chunk of streamChunks) session.append_ref(chunk);
session.finalize_ref();
let out = '';
for await (const chunk of queryChunks) out += session.append_query(chunk);
out += session.finalize();   // trailing output + "---LOG---\n<log>"
```

### `align_wasm_full(target, query, preset, output_sam, output_cigar) → String`

Single-call entry point. Builds index from `target` (FASTA text), aligns all
sequences in `query` (FASTA or FASTQ text), returns results. Constrained by
V8's ~512 MB max-string length on both inputs — prefer `AlignSession` for
anything beyond ~100 MB.

**Parameters:**
- `target`: FASTA reference text (multi-sequence supported)
- `query`: FASTA or FASTQ query text (auto-detected by `@` vs `>` header)
- `preset`: `"map-ont"`, `"map-hifi"`, `"sr"`, `"splice"`, `"asm20"`
- `output_sam`: `true` for SAM, `false` for PAF
- `output_cigar`: `true` to compute CIGAR strings

**Return value:** String with format `"<output>\n---LOG---\n<log>"`. Split on
`\n---LOG---\n` to separate alignment output from timing/progress messages.

### `align_wasm(target, query, output_sam, is_splice) → String`

Legacy API. Wraps `align_wasm_full` with preset = `"splice"` or `"map-ont"`.
Returns alignment output only (no log section).

### `force_align_wasm(tseq, qseq) → String`

Force-align two sequences via full DP (no seeding/chaining). Returns CIGAR string.
Useful for small pairwise alignments.

### `initThreadPool(num_threads) → Promise` (threaded build only)

Initialize the rayon thread pool. Must be called once before any parallel
alignment. Only available when built with `wasm-threads` feature.

---

## Building

### Prerequisites

```bash
# For all builds
cargo install wasm-pack

# For threaded build only
rustup component add rust-src --toolchain nightly
cargo install wasm-bindgen-cli
```

### Single-Threaded Build

```bash
bash scripts/build-wasm.sh
```

Uses stable Rust. Output in `pkg/` at the repo root. No special browser requirements.

The wasm bindings live in the `rammap-core` library crate (whose crate name is
`rammap`, so the wasm-pack artifact is `pkg/rammap.js`). The script builds that
crate and relocates its `pkg/` to the repo root so the `web/` demo's
`../pkg/rammap.js` import resolves. Equivalent to:
```bash
(cd rammap-core && wasm-pack build --target web --release --no-default-features)
mv rammap-core/pkg pkg
```
(Don't pass `wasm-pack`'s `--out-dir` — it's forwarded to cargo as the
nightly-only `--artifact-dir`.)

### Multi-Threaded Build

```bash
bash scripts/build-wasm-threads.sh
```

Requires nightly Rust. Uses `-Z build-std` to rebuild std with atomics support.
Output in `pkg/`.

The build script:
1. Compiles with `+simd128,+atomics,+bulk-memory` target features
2. Links with `--shared-memory` and TLS exports
3. Runs `wasm-bindgen` to generate JS glue
4. Patches `workerHelpers.js` for browser module loading (fixes a
   `import('../../..')` directory resolution issue)

**RUSTFLAGS used:**
```
-C target-feature=+simd128,+atomics,+bulk-memory
-C link-arg=--shared-memory
-C link-arg=--max-memory=4294967296
-C link-arg=--import-memory
-C link-arg=--export=__wasm_init_tls
-C link-arg=--export=__tls_size
-C link-arg=--export=__tls_align
-C link-arg=--export=__tls_base
```

### Cargo.toml Features

The WASM bindings build from the `rammap-core` library crate. Its features:

| Feature | Effect |
|---------|--------|
| `parallel` (default) | Enables `rayon` for native parallelism |
| `wasm-threads` | Enables `rayon` + `wasm-bindgen-rayon` for browser threading |
| (no features, i.e. `--no-default-features`) | Single-threaded WASM |

The CLI is a separate crate (`rammap`, which depends on `rammap-core`) and owns
the binary plus the allocator opt-ins (`jemalloc`/`mimalloc`); there is no longer
a `cli` feature on the library.

The `web_spin_lock` feature on rayon is always enabled — it replaces futex-based
locks with spin locks on WASM (required because `Atomics.wait` traps on the
browser main thread). It's a no-op on native targets.

### `.cargo/config.toml`

SIMD128 is enabled globally for WASM targets:
```toml
[target.wasm32-wasip1]
rustflags = ["-C", "target-feature=+simd128"]

[target.wasm32-unknown-unknown]
rustflags = ["-C", "target-feature=+simd128"]
```

---

## Running

### Browser Demo

```bash
# Single-threaded (any HTTP server works)
python3 -m http.server 8080

# Multi-threaded (needs COOP/COEP headers for SharedArrayBuffer)
python3 scripts/serve.py
```

Open `http://localhost:8080/web/`. The demo page provides:
- Drag-and-drop file input for reference and queries
- Preset selector (map-ont, map-hifi, sr, splice, asm20)
- Output format (PAF/SAM) with CIGAR toggle
- Thread count selector (1 to `navigator.hardwareConcurrency`)
- Real-time log panel (streams progress from Web Worker)
- Output panel with Raw/Table view toggle
- Table view parses PAF/SAM fields with sticky headers, CIGAR expansion on hover

The alignment runs in a Web Worker (`pkg/worker.js`) to keep the UI responsive.
When built with `wasm-threads`, rayon's internal workers provide additional
parallelism for the alignment loop.

**Browser requirements (threaded build):**
- `SharedArrayBuffer` support (Chrome 68+, Firefox 79+, Safari 15.2+)
- Server must send headers:
  - `Cross-Origin-Opener-Policy: same-origin`
  - `Cross-Origin-Embedder-Policy: require-corp`

---

## Testing

### Unit Tests (wasm-bindgen-test)

```bash
wasm-pack test --node -- --lib
```

Tests in `wasm_lib.rs`:
- `test_force_align_exact_match` — 8bp exact match → 8M CIGAR
- `test_force_align_with_mismatch` — single mismatch
- `test_force_align_with_insertion` — 1bp insertion
- `test_force_align_with_deletion` — 1bp deletion
- `test_align_wasm_basic` — full pipeline (48bp target, 40bp query)
- `test_align_wasm_longer_sequence` — 400bp target, 200bp query
- `test_align_wasm_sam_output` — SAM format output

### Concordance Testing

WASM alignment output should be identical to native for the same input. The
`wasm-pack test --node` unit tests above exercise this on small cases; for
larger inputs, diff the browser demo's output (`web/`) against the native
`rammap` CLI on the same reference and reads.

---

## Performance

Benchmarked on chr20 reference (64 Mbp), 20,000 ONT reads, single-threaded.
WASM runs via wasmtime on the same x86_64 machine (AMD Ryzen 9 7900X).

| Stage | Native (x86 -t1) | WASM SIMD128 | Ratio |
|-------|------------------:|-------------:|------:|
| Index build | 0.9s | 1.3s | 1.4x |
| Sketching | 3.9s | 5.0s | 1.3x |
| Seeding | 7.3s | 9.5s | 1.3x |
| Chaining | 7.3s | 13.2s | 1.8x |
| Alignment (DP) | 97.5s | 164.7s | 1.7x |
| Post-chain | 7.0s | 10.7s | 1.5x |
| **Total mapping** | **123.4s** | **202.6s** | **1.64x** |

WASM SIMD128 is **1.64x slower** than native x86 single-threaded — competitive
for JIT-compiled code. The DP kernels (1.7x) have the largest gap due to
wasmtime's SIMD128 → x86 translation overhead. Non-SIMD stages (1.3-1.5x)
reflect general JIT overhead.

Without SIMD128 enabled (`-C target-feature=+simd128`), the DP kernels fall back
to scalar, and total slowdown increases to ~1.84x.

---

## Limitations

- **Single-threaded WASI**: wasmtime does not yet support WASI threads for Rust's
  `std::thread`, so WASI builds are always single-threaded.
- **Browser threading requires headers**: `SharedArrayBuffer` needs COOP/COEP HTTP
  headers. Embedding in iframes requires the parent page to also send these headers.
- **Memory**: WASM linear memory grows on demand but has a configured maximum (4 GB
  for threaded builds). Large reference genomes may approach this limit.
- **Gzip is decompressed on the JS side**: The WASM module itself takes raw bytes,
  but the demo's worker (`web/worker.js`) detects the gzip magic header and pipes
  the file through `DecompressionStream('gzip')` before chunking it into
  `AlignSession.append_ref` / `append_query`. So `.fa.gz` / `.fq.gz` files load
  end-to-end without an explicit decompress step. WASI consumers wanting the same
  behavior need to wire it up themselves.
- **Index serialization**: The WASM build does not support loading `.mmi` or `.rmmi`
  index files. Indices are built from FASTA at runtime.
- **Startup overhead**: wasmtime JIT compilation adds ~100ms startup time. Browser
  WASM loading depends on module size (~300-500 KB).
