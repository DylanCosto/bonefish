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

static TERMINAL_MESSAGES_TO_LOG: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Send the library's stderr diagnostics (index loaded, input warnings) to the `log` crate,
/// target `bonefish::terminal`, instead of printing them. The CLI does this while it draws
/// progress on stderr; other callers keep the plain stderr messages unless they opt in.
pub fn route_terminal_messages_to_log(enabled: bool) {
    TERMINAL_MESSAGES_TO_LOG.store(enabled, std::sync::atomic::Ordering::Relaxed);
}

fn terminal_message(level: log::Level, message: std::fmt::Arguments<'_>) {
    if TERMINAL_MESSAGES_TO_LOG.load(std::sync::atomic::Ordering::Relaxed) {
        log::log!(target: "bonefish::terminal", level, "{message}");
    } else {
        eprintln!("{message}");
    }
}
