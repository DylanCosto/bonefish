#!/usr/bin/env bash
# Build bonefish (release profile, jemalloc, stripped) and package it as
# dist/bonefish-v<version>-<platform>.tar.gz plus a .sha256 file.
#
# Usage: scripts/package-release.sh PLATFORM [RUST_TARGET]
#   e.g. scripts/package-release.sh linux-x86_64
#        scripts/package-release.sh macos-arm64 aarch64-apple-darwin
set -euo pipefail
cd "$(dirname "$0")/.."

platform=${1:?usage: package-release.sh PLATFORM [RUST_TARGET]}
target=${2:-}
version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
name=bonefish-v$version-$platform

CARGO_PROFILE_RELEASE_STRIP=true cargo build --release --locked --features jemalloc ${target:+--target "$target"}
bin=target/${target:+$target/}release/bonefish
"$bin" --version

rm -rf "dist/$name"
mkdir -p "dist/$name"
cp "$bin" README.md LICENSE "dist/$name/"
tar -C dist -czf "dist/$name.tar.gz" "$name"
rm -rf "dist/$name"
if command -v sha256sum >/dev/null; then sum=(sha256sum); else sum=(shasum -a 256); fi
(cd dist && "${sum[@]}" "$name.tar.gz" > "$name.tar.gz.sha256")
echo "dist/$name.tar.gz"
