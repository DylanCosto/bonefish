//! Cached environment-variable feature flags.
//!
//! `std::env::var` takes a process-wide lock and allocates a String on every
//! call. Reading these flags once at first access via `LazyLock` removes that
//! cost from per-read hot paths (DP/chain dispatch).

use std::sync::LazyLock;

#[inline(always)]
fn flag(name: &'static str) -> bool {
    std::env::var(name).is_ok()
}

pub static FORCE_SCALAR: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_FORCE_SCALAR"));
#[cfg(target_arch = "x86_64")]
pub static FORCE_SSE: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_FORCE_SSE"));
#[cfg(target_arch = "x86_64")]
pub static FORCE_AVX2: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_FORCE_AVX2"));
#[cfg(target_arch = "x86_64")]
pub static FORCE_AVX512: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_FORCE_AVX512"));
// wasm always uses the wasm chain kernel, so gate it out there to avoid a dead-code warning
#[cfg(not(target_arch = "wasm32"))]
pub static FORCE_SCALAR_CHAIN: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_FORCE_SCALAR_CHAIN"));
pub static COMPARE_SCALAR: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_COMPARE_SCALAR"));

// Phase 3 experiment knobs (extension work only; unset = minimap2 behaviour).
fn int_env(name: &'static str) -> Option<i32> {
    std::env::var(name).ok().and_then(|v| v.parse().ok())
}
/// RAMMAP_EXP_MAX_EXT: cap on the left/right extension window (default: max_gap).
pub static EXP_MAX_EXT: LazyLock<Option<i32>> = LazyLock::new(|| int_env("RAMMAP_EXP_MAX_EXT"));
/// RAMMAP_EXP_EXT_ZDROP: z-drop for left/right extensions only (default: opt.alignment.zdrop).
pub static EXP_EXT_ZDROP: LazyLock<Option<i32>> = LazyLock::new(|| int_env("RAMMAP_EXP_EXT_ZDROP"));
/// RAMMAP_EXP_EXT_BW: DP bandwidth for left/right extensions only (default: 1.5*w + 1).
pub static EXP_EXT_BW: LazyLock<Option<i32>> = LazyLock::new(|| int_env("RAMMAP_EXP_EXT_BW"));
/// RAMMAP_EXP_GAP_BW: cap on the long-gap-fill DP bandwidth (default: 1.5*w_long + 1).
pub static EXP_GAP_BW: LazyLock<Option<i32>> = LazyLock::new(|| int_env("RAMMAP_EXP_GAP_BW"));
/// RAMMAP_EXP_GAP_BW_ADAPT=margin: gap-fill bandwidth = min(default, |qlen - tlen| + margin),
/// so gaps holding a real large indel keep a band wide enough to cross it.
pub static EXP_GAP_BW_ADAPT: LazyLock<Option<i32>> = LazyLock::new(|| int_env("RAMMAP_EXP_GAP_BW_ADAPT"));
/// RAMMAP_EXP_RESCUE_COVER=frac: skip the RMQ rescue re-chain when the best chain already
/// spans at least frac of the read (minimap2 runs it whenever the span exceeds 10%).
pub static EXP_RESCUE_COVER: LazyLock<Option<f32>> =
    LazyLock::new(|| std::env::var("RAMMAP_EXP_RESCUE_COVER").ok().and_then(|v| v.parse().ok()));
/// RAMMAP_EXP_UNGAPPED: also use the short-read "ungapped first" gap fill for long reads
/// (equal-length gaps whose ungapped score beats any gapped alternative skip DP).
pub static EXP_UNGAPPED: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_EXP_UNGAPPED"));
/// RAMMAP_EXP_EXT_RESCUE: re-run a z-dropped narrow-band extension with the default band when
/// the unaligned tail still matches the reference at a diagonal shift wider than the band.
pub static EXP_EXT_RESCUE: LazyLock<bool> = LazyLock::new(|| flag("RAMMAP_EXP_EXT_RESCUE"));

/// RAMMAP_EXP_SEG_RMQ=1: RMQ chaining uses a segment tree over query-position ranks instead of
/// the AVL tree (same minimum; ties between equal priorities may resolve differently).
pub static EXP_SEG_RMQ: LazyLock<bool> =
    LazyLock::new(|| std::env::var("RAMMAP_EXP_SEG_RMQ").is_ok_and(|v| v != "0"));
/// RAMMAP_EXP_RMQ_CHECK=1: run both RMQ structures and count disagreements (debug).
pub static EXP_RMQ_CHECK: LazyLock<bool> =
    LazyLock::new(|| std::env::var("RAMMAP_EXP_RMQ_CHECK").is_ok_and(|v| v != "0"));
/// RAMMAP_EXP_EXT_ADAPT=W: adaptive extension band of W cells (see --ext-adaptive).
pub static EXP_EXT_ADAPT: LazyLock<Option<i32>> =
    LazyLock::new(|| std::env::var("RAMMAP_EXP_EXT_ADAPT").ok().and_then(|v| v.parse().ok()));
