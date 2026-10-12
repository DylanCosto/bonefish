#!/bin/bash
# Build WASM with multithreading support (rayon via wasm-bindgen-rayon).
# Requires: nightly Rust, rust-src component, wasm-bindgen-cli
#
# Setup:
#   rustup component add rust-src --toolchain nightly
#   cargo install wasm-bindgen-cli
#
# Usage: bash scripts/build-wasm-threads.sh
# Serve: python3 scripts/serve.py

set -euo pipefail

echo "[*] Building WASM with threads + SIMD128..."

RUSTFLAGS="\
-C target-feature=+simd128,+atomics,+bulk-memory \
-C link-arg=--shared-memory \
-C link-arg=--max-memory=4294967296 \
-C link-arg=--import-memory \
-C link-arg=--export=__wasm_init_tls \
-C link-arg=--export=__tls_size \
-C link-arg=--export=__tls_align \
-C link-arg=--export=__tls_base" \
  rustup run nightly cargo build \
    -p bonefish-core \
    --target wasm32-unknown-unknown \
    --release \
    --features wasm-threads \
    --no-default-features \
    -Z build-std=panic_abort,std \
    --lib

echo "[*] Generating JS bindings..."
wasm-bindgen \
  --target web \
  --out-dir pkg \
  target/wasm32-unknown-unknown/release/bonefish.wasm

echo "[*] Patching workerHelpers.js for browser module loading..."
# wasm-bindgen-rayon's workerHelpers.js uses `import('../../..')` which resolves
# to a directory (pkg/), not a JS file. Browsers can't handle directory imports.
# Patch to use the explicit module path.
find pkg/snippets -name "workerHelpers.js" -exec \
  sed -i "s|import('../../..')|import('../../../bonefish.js')|g" {} \;

echo "[*] Done. WASM package in pkg/"
echo "[*] Serve with: python3 scripts/serve.py"
