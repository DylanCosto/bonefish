#!/bin/bash
# Build WASM without multithreading (single-threaded, no SharedArrayBuffer needed).
# Works on stable Rust.
#
# Usage: bash scripts/build-wasm.sh

set -euo pipefail

echo "[*] Building WASM (single-threaded + SIMD128)..."

# Build the bonefish-core lib crate. The lib's crate name is `bonefish`, so the
# wasm-pack output stays bonefish.js / bonefish_bg.wasm. wasm-pack writes to
# bonefish-core/pkg; relocate it to the repo-root pkg/ so the web/ demo's
# `../pkg/bonefish.js` import resolves. (`--out-dir` is avoided: wasm-pack forwards
# it to cargo as the nightly-only `--artifact-dir`.)
( cd bonefish-core && wasm-pack build --target web --release --no-default-features )
rm -rf pkg && mv bonefish-core/pkg pkg

echo "[*] Done. WASM package in pkg/"
echo "[*] Serve with: python3 -m http.server 8080 (or python3 scripts/serve.py for threads)"
