pub mod align;
pub mod api;
pub mod fasta;
#[cfg(not(target_arch = "wasm32"))]
pub mod sorted_bam;

// Convenience re-exports for library users
pub use api::{Aligner, Mapping, MapResult, MapOpts, Preset, Strand, CigarOp};
pub use api::{dp_align, dp_global, dp_local, dp_extension, DpScoring, DpAlignment, encode_nt4};
// Runtime knob for the per-thread DP scratch-cache retention cap (peak-RSS control).
pub use align::dp::{set_dp_cache_cap_mb, DEFAULT_DP_CACHE_CAP_MB, DP_CACHE_UNBOUNDED};

// The CLI registers this target while drawing progress on stderr. Otherwise
// these longstanding diagnostic messages go directly to stderr.
fn terminal_message(level: log::Level, message: std::fmt::Arguments<'_>) {
    const TARGET: &str = "bonefish::terminal";
    if log::log_enabled!(target: TARGET, log::Level::Trace) {
        log::log!(target: TARGET, level, "{message}");
    } else {
        eprintln!("{message}");
    }
}
