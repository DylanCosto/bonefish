//! Mapping orchestration: sketch -> seed -> chain -> filter.
//!
//! `map_query` is the entry point: it sketches the query into minimizers,
//! collects seed hits from the index, chains them, and filters chains by
//! score and parent assignment. `MapOptions` holds all algorithm parameters
//! (seeding, chaining, scoring, flags). `MapContext` provides per-thread
//! reusable buffers to avoid allocation. Each surviving chain is returned as
//! a `Mapping` (reference coordinates, score, strand, and anchor list).

use bitflags::bitflags;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

// Sub-timers for post-chain breakdown (nanoseconds, summed across threads)
static PC_RESCUE_NS: AtomicU64 = AtomicU64::new(0);
static PC_RECHAIN_NS: AtomicU64 = AtomicU64::new(0);
static PC_RECHAIN_COUNT: AtomicU64 = AtomicU64::new(0);
static PC_RESCUE_SKIP: AtomicU64 = AtomicU64::new(0);
static PC_RESCUE_COUNT: AtomicU64 = AtomicU64::new(0);
static PC_RESCUE_ANCHORS: AtomicU64 = AtomicU64::new(0);
static PC_BUILD_REGS_NS: AtomicU64 = AtomicU64::new(0);
static PC_PARENT_NS: AtomicU64 = AtomicU64::new(0);
static PC_DIVEST_NS: AtomicU64 = AtomicU64::new(0);
static PC_SQUEEZED_NS: AtomicU64 = AtomicU64::new(0);

pub fn print_post_chain_breakdown() {
    eprintln!("--- Post-chain sub-timers ---");
    eprintln!("  RMQ rescue:      {:.3}s ({} calls, {} anchors, {} skipped: no joinable chains)", PC_RESCUE_NS.load(Relaxed) as f64 / 1e9, PC_RESCUE_COUNT.load(Relaxed), PC_RESCUE_ANCHORS.load(Relaxed), PC_RESCUE_SKIP.load(Relaxed));
    eprintln!("    of which high-occ re-chain: {:.3}s ({} calls)", PC_RECHAIN_NS.load(Relaxed) as f64 / 1e9, PC_RECHAIN_COUNT.load(Relaxed));
    eprintln!("  Build regs:      {:.3}s", PC_BUILD_REGS_NS.load(Relaxed) as f64 / 1e9);
    eprintln!("  Parent+filter:   {:.3}s", PC_PARENT_NS.load(Relaxed) as f64 / 1e9);
    eprintln!("  Div+strand:      {:.3}s", PC_DIVEST_NS.load(Relaxed) as f64 / 1e9);
    eprintln!("  Squeezed+bounds: {:.3}s", PC_SQUEEZED_NS.load(Relaxed) as f64 / 1e9);
    eprintln!("  Extension probe stops: {}", crate::align::extend::PROBE_STOPS.load(Relaxed));
    eprintln!("  Sub-anchored long gaps: {}", crate::align::extend::GAP_SUB_USED.load(Relaxed));
    eprintln!("  Narrow gap band: {} used ({} widened), {} too wide, {} edge fallbacks",
        crate::align::extend::GAP_NARROW_USED.load(Relaxed), crate::align::extend::GAP_NARROW_WIDENED.load(Relaxed),
        crate::align::extend::GAP_NARROW_KMER.load(Relaxed), crate::align::extend::GAP_NARROW_EDGE.load(Relaxed));
    if *crate::align::env_flags::EXP_RMQ_CHECK {
        use crate::align::chain_rmq::{RMQ_CHECK_N, RMQ_CHECK_DIFF, RMQ_CHECK_DIFF_PRI};
        eprintln!("  RMQ check: {} queries, {} different index ({} different priority)", RMQ_CHECK_N.load(Relaxed), RMQ_CHECK_DIFF.load(Relaxed), RMQ_CHECK_DIFF_PRI.load(Relaxed));
    }
    eprintln!("-----------------------------");
}
use crate::align::sketch::{sketch_sequence, Minimizer};
use crate::align::index::Index;
use crate::align::chain::{chain_anchors, chain_anchors_cdna};
#[cfg(feature = "parallel")]
use crate::align::chain::chain_anchors_partitioned;
use crate::align::chain_rmq::chain_anchors_rmq;
use crate::align::sort::radix_sort_128x;
use crate::align::filter::{FilterParams, ParentState, Filterable, check_secondary_filter, scale_alt_score};

use crate::align::seed::{hash64, compute_read_hash, collect_seed_hits, collect_seed_hits_with_occ, filter_minimizers_by_occ};

// Behavioral flag constants
bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct AlignFlags: u64 {
        const NO_DIAG          = 0x001;
        const NO_DUAL          = 0x002;
        const NO_QUAL          = 0x010;
        const OUT_CIGAR        = 0x020;
        /// Spliced alignment (RNA/cDNA presets): intron-aware chaining and DP.
        const SPLICE           = 0x080;
        /// Splice signals of the forward transcript strand (GT..AG on the read's strand).
        const SPLICE_FOR       = 0x100;
        /// Splice signals of the reverse transcript strand (CT..AC).
        const SPLICE_REV       = 0x200;
        const NO_LJOIN         = 0x400;
        const NO_PRINT_2ND     = 0x4000;
        const EQX              = 0x8000;
        const LONG_CIGAR       = 0x10000;
        const SOFTCLIP         = 0x20000;
        /// Score splice-site flanking bases (`--splice-flank`).
        const SPLICE_FLANK     = 0x40000;
        const FOR_ONLY         = 0x100000;
        const REV_ONLY         = 0x200000;
        const ALL_CHAINS       = 0x800000;
        const NO_END_FLT       = 0x1000000;
        const HARD_MASK_LEVEL  = 0x2000000;
        const SAM_HIT_ONLY     = 0x4000000;
        const PAF_NO_HIT       = 0x8000000;
        const NO_HASH_NAME     = 0x40000000;
        const RMQ_CHAIN        = 0x80000000;
        const QSTRAND          = 0x100000000;
        const NO_INV           = 0x200000000;
        const SECONDARY_SEQ    = 0x800000000;
        /// Write introns of spliced alignments as BED instead of SAM/PAF (`--write-junc`).
        const OUT_JUNC         = 0x10000000000;
        /// Partition anchors by ref_id and chain partitions in parallel.
        /// Beneficial for assembly-to-assembly (asm2asm) workloads where a
        /// single query has millions of anchors across many references.
        /// Causes minor non-determinism in chain selection (~0.01% of output).
        const PAR_CHAIN        = 0x20000000000;
    }
}

// Re-export seed encoding constants from sketch.rs (where Minimizer is defined)
pub use crate::align::sketch::{SEED_SEG_SHIFT, SEED_SEG_MASK, SEED_SELF};

#[derive(Debug, Clone)]
pub struct SeedingParams {
    pub mid_occ: usize,
    pub max_occ: usize,
    pub max_max_occ: usize,
    pub occ_dist: i32,
    pub q_occ_frac: f32,
    pub min_mid_occ: i32,
    pub max_mid_occ: i32,
    /// Keep 1/f of the query minimizer hash space for reads of at least min_len bp
    /// (`--query-subsample f[,min_len]`): fewer anchors, faster, loses some SV support.
    pub query_subsample: Option<(u32, u32)>,
}

#[derive(Debug, Clone)]
pub struct ChainingParams {
    pub min_cnt: i32,
    pub min_chain_score: i32,
    pub max_gap: i32,
    pub max_gap_ref: i32,
    pub max_dist_x: i32,
    pub max_dist_y: i32,
    pub bandwidth: i32,
    pub bandwidth_long: i32,
    pub max_chain_skip: i32,
    pub max_chain_iter: i32,
    pub chn_pen_gap: f32,
    pub chn_pen_skip: f32,
    pub chain_gap_scale: f32,
    pub rmq_rescue_size: i32,
    pub rmq_rescue_ratio: f32,
    pub rmq_inner_dist: i32,
    pub rmq_size_cap: i32,
    /// RMQ chaining with a segment tree instead of the AVL tree: same range minima, but a tie
    /// between exactly equal priorities may pick a different anchor (output not guaranteed
    /// byte-identical to minimap2).
    pub seg_rmq: bool,
    /// Run the RMQ rescue re-chain only when two first-pass chains could be joined
    /// (collinear, same strand, gaps within the long-join distance and bandwidth).
    pub rescue_join: bool,
    /// Skip the RMQ long-join rescue for reads with at least this many anchors per query base
    /// (tandem-repeat/satellite reads; 0 = never skip).
    pub rescue_max_density: u32,
}

#[derive(Debug, Clone)]
pub struct ScoringParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gap_open2: i32,
    pub gap_extend2: i32,
    pub transition: i32,
    pub ambig_penalty: i32,
    /// Spliced alignment: cost of a non-canonical splice site (`-C`).
    pub noncanon_penalty: i32,
    /// Spliced alignment: bonus for an annotated junction (`--junc-bonus`).
    pub junc_bonus: i32,
    /// Spliced alignment: penalty for an unannotated junction when annotations are loaded (`--junc-pen`).
    pub junc_pen: i32,
}

#[derive(Debug, Clone)]
pub struct AlignmentParams {
    pub zdrop: i32,
    pub zdrop_inv: i32,
    pub end_bonus: i32,
    pub max_sw_mat: i64,
    pub min_dp_max: i32,
    pub min_dp_len: i32,
    /// Sub-anchored fill of long gaps (see `extend::gap_subanchored`).
    pub gap_subanchor: bool,
    /// Sub-anchoring skips tandem repeats with a period up to this many bases (0 = none).
    pub gap_subanchor_tr: u32,
    pub anchor_ext_len: i32,
    pub anchor_ext_shift: i32,
    pub max_clip_ratio: f32,
    /// DP bandwidth for left/right extensions (None: 1.5*w + 1, as minimap2).
    pub ext_bw: Option<i32>,
    /// Gap-fill bandwidth = min(default, |qlen - tlen| + margin) (None: minimap2 bandwidth).
    pub gap_bw_margin: Option<i32>,
    /// Ungapped-first gap fill for long reads (as minimap2 does for short reads).
    pub ungapped_fill: bool,
    /// With a narrowed extension band: re-run z-dropped extensions with the default band when
    /// the unaligned tail shows an indel wider than the narrow band.
    pub ext_rescue: bool,
    /// Extension probe (len, min score per kb): for extension windows longer than 1.5 x len,
    /// first align only the first `len` bases; keep that result if it z-drops or scores below
    /// the minimum (a low-identity read end), otherwise align the full window.
    pub ext_probe: Option<(i32, i32, i32)>,
    /// Lower bound of the adaptive gap-fill band (None: the default bandwidth).
    pub gap_bw_floor: Option<i32>,
    /// Narrow gap-fill band |qlen - tlen| + margin, used unless read k-mers match the target
    /// only off that band, or the path nears its edge (then the regular band).
    pub gap_bw_narrow: Option<i32>,
    /// Extensions with an adaptive band of this many cells per anti-diagonal (AVX2).
    pub ext_adaptive: Option<i32>,
    /// Dual-affine DP with the scalar-semantics AVX2 kernel (`dp::dual_clean`) instead of the
    /// minimap2-identical one.
    pub clean_dp: bool,
}

#[derive(Debug, Clone)]
pub struct FilteringParams {
    pub best_n: i32,
    pub pri_ratio: f32,
    pub mask_level: f32,
    pub mask_len: i32,
    pub alt_drop: f32,
    pub seed: i32,
    pub chain_skip_scale: f32,
    pub max_qlen: i32,
}

#[derive(Debug, Clone)]
pub struct MapOptions {
    pub seeding: SeedingParams,
    pub chaining: ChainingParams,
    pub scoring: ScoringParams,
    pub alignment: AlignmentParams,
    pub filtering: FilteringParams,
    pub flags: AlignFlags,
    pub mini_batch_size: i64,
    /// Annotated splice junctions (`--junc-bed`) for the spliced presets, resolved against
    /// the reference names of the loaded index.
    pub junctions: Option<std::sync::Arc<crate::align::junc::JunctionDb>>,
}

impl Default for MapOptions {
    fn default() -> Self {
        MapOptions {
            seeding: SeedingParams {
                mid_occ: 0,
                max_occ: 0,
                max_max_occ: 4095,
                occ_dist: 500,
                q_occ_frac: 0.01,
                min_mid_occ: 10,
                max_mid_occ: 1000000,
                query_subsample: None,
            },
            chaining: ChainingParams {
                min_cnt: 3,
                min_chain_score: 40,
                max_gap: 5000,
                max_gap_ref: -1,
                max_dist_x: 5000,
                max_dist_y: 5000,
                bandwidth: 500,
                bandwidth_long: 20000,
                max_chain_skip: 25,
                max_chain_iter: 5000,
                chn_pen_gap: 0.12, // chain_gap_scale * 0.01 * k (0.8 * 0.01 * 15)
                chn_pen_skip: 0.0,
                chain_gap_scale: 0.8,
                rmq_rescue_size: 1000,
                rmq_rescue_ratio: 0.1,
                rmq_inner_dist: 1000,
                rmq_size_cap: 100000,
                seg_rmq: false,
                rescue_join: false,
                rescue_max_density: 0,
            },
            scoring: ScoringParams {
                match_score: 2,
                mismatch_penalty: 4,
                gap_open: 4,
                gap_extend: 2,
                gap_open2: 24,
                gap_extend2: 1,
                transition: 0,
                ambig_penalty: 1,
                noncanon_penalty: 0,
                junc_bonus: 0,
                junc_pen: 0,
            },
            alignment: AlignmentParams {
                zdrop: 400,
                zdrop_inv: 200,
                end_bonus: -1,
                max_sw_mat: 100_000_000,
                min_dp_max: 80, // min_chain_score * a = 40 * 2
                min_dp_len: 200,
                gap_subanchor: false,
                gap_subanchor_tr: 0,
                anchor_ext_len: 20,
                anchor_ext_shift: 6,
                max_clip_ratio: 1.0,
                ext_bw: None,
                gap_bw_margin: None,
                ungapped_fill: false,
                ext_rescue: false,
                ext_probe: None,
                gap_bw_floor: None,
                gap_bw_narrow: None,
                ext_adaptive: None,
                clean_dp: false,
            },
            filtering: FilteringParams {
                best_n: 5,
                pri_ratio: 0.8,
                mask_level: 0.5,
                mask_len: i32::MAX,
                alt_drop: 0.15,
                seed: 11,
                chain_skip_scale: 0.0,
                max_qlen: 0,
            },
            flags: AlignFlags::empty(),
            // 100 Mbp: as fast as 500 Mbp at 16 threads on a human genome, 2.7 GB less memory.
            mini_batch_size: 100_000_000,
            junctions: None,
        }
    }
}

/// Reusable buffers for the chaining stage. Separated from MapContext to
/// narrow the interface: chaining functions receive only what they need.
pub struct ChainingBuffers {
    // DP arrays (reused across reads to avoid per-read allocation)
    pub predecessors: Vec<i32>,
    pub scores: Vec<i32>,
    pub peak_scores: Vec<i32>,
    pub visited: Vec<i32>,
    // Backtrack candidate buffer (one (score, idx) per candidate chain)
    pub bt_candidates: Vec<Minimizer>,
    // SoA buffers for SIMD chaining (extracted from Minimizer AoS)
    pub soa_ref_pos: Vec<i32>,
    pub soa_query_pos: Vec<i32>,
    pub soa_query_span: Vec<i32>,
    pub soa_ref_id_strand: Vec<u32>,
    pub soa_scores_buf: Vec<i32>,
    /// ref_pos - query_pos per anchor (wrapping), for the one-load diagonal band test.
    pub soa_diag: Vec<i32>,
    /// Diagonal buckets for the DP chaining scan (see `chain_simd::DiagBuckets`).
    pub diag_buckets: crate::align::chain_simd::DiagBuckets,
}

impl Default for ChainingBuffers {
    fn default() -> Self {
        Self::new()
    }
}

impl ChainingBuffers {
    pub fn new() -> Self {
        ChainingBuffers {
            predecessors: Vec::with_capacity(4096),
            scores: Vec::with_capacity(4096),
            peak_scores: Vec::with_capacity(4096),
            visited: Vec::with_capacity(4096),
            bt_candidates: Vec::with_capacity(1024),
            soa_ref_pos: Vec::new(),
            soa_query_pos: Vec::new(),
            soa_query_span: Vec::new(),
            soa_ref_id_strand: Vec::new(),
            soa_scores_buf: Vec::new(),
            soa_diag: Vec::new(),
            diag_buckets: Default::default(),
        }
    }
}

// ---------------------------------------------------------------------------
// Trait definitions for swappable pipeline stages
// ---------------------------------------------------------------------------

/// Trait for anchor chaining strategies.
///
/// Given a reference-sorted array of anchors, produces chain descriptors and
/// reordered anchor sequences. Two implementations exist:
/// - [`crate::align::chain::chain_anchors`] — DP-based chaining (default)
/// - [`crate::align::chain_rmq::chain_anchors_rmq`] — RMQ-based chaining (for asm/hqae presets)
///
/// The output `Vec<u64>` encodes `(score << 32 | count)` per chain, and the
/// `Vec<Minimizer>` contains reordered anchors (chains packed contiguously).
pub trait Chainer {
    fn chain(
        &self,
        params: &ChainingParams,
        anchors: &mut [Minimizer],
        bufs: &mut ChainingBuffers,
    ) -> (Vec<u64>, Vec<Minimizer>);
}

/// Trait for the alignment/extension stage.
///
/// Takes a chain's anchors and nt4-encoded sequences, produces an alignment
/// result with CIGAR operations. The single implementation wraps the dp-based
/// [`crate::align::extend::align_anchors`] function.
///
/// All coordinates are region-relative: the caller extracts a target region from
/// the packed index and adjusts anchor ref_pos before calling, then restores
/// absolute coordinates from the result.
pub trait Aligner {
    fn align(
        &self,
        anchors: &mut [Minimizer],
        qseq: &[u8],
        tseq: &[u8],
        opt: &MapOptions,
        ctx: &mut crate::align::extend::AlignmentContext,
        call: &crate::align::extend::AlignAnchorContext,
    ) -> crate::align::extend::AlignResult;
}

/// Default DP-based chainer (wraps [`crate::align::chain::chain_anchors`]).
pub struct DpChainer {
    pub max_dist_x: i32,
    pub max_dist_y: i32,
}

impl Chainer for DpChainer {
    fn chain(
        &self,
        params: &ChainingParams,
        anchors: &mut [Minimizer],
        bufs: &mut ChainingBuffers,
    ) -> (Vec<u64>, Vec<Minimizer>) {
        chain_anchors(params, self.max_dist_x, self.max_dist_y, anchors, bufs)
    }
}

/// RMQ-based chainer (wraps [`crate::align::chain_rmq::chain_anchors_rmq`]).
pub struct RmqChainer;

impl Chainer for RmqChainer {
    fn chain(
        &self,
        params: &ChainingParams,
        anchors: &mut [Minimizer],
        bufs: &mut ChainingBuffers,
    ) -> (Vec<u64>, Vec<Minimizer>) {
        chain_anchors_rmq(params, anchors, bufs)
    }
}

/// Default aligner (wraps [`crate::align::extend::align_anchors`]).
pub struct RMAligner;

impl Aligner for RMAligner {
    fn align(
        &self,
        anchors: &mut [Minimizer],
        qseq: &[u8],
        tseq: &[u8],
        opt: &MapOptions,
        _ctx: &mut crate::align::extend::AlignmentContext,
        call: &crate::align::extend::AlignAnchorContext,
    ) -> crate::align::extend::AlignResult {
        // The trait gives `tseq: &[u8]` (immutable, pre-filled by the
        // caller). `align_anchors` now needs `&mut [u8]` for its lazy-fill
        // path, but with `lazy_extract = None` it treats tseq as pre-filled
        // and never writes. Copy to a local `Vec` here to satisfy the
        // signature; this is the API trait path, not the hot mapping path
        // (the hot path in `pipeline.rs` calls `align_anchors` directly with
        // a thread-local buffer and `lazy_extract = Some(...)`).
        let mut tseq_buf = tseq.to_vec();
        crate::align::extend::align_anchors(anchors, qseq, &mut tseq_buf, None, opt, call)
    }
}

/// Per-thread reusable buffers for the mapping pipeline (seeding + chaining).
pub struct MapContext {
    pub anchors: Vec<Minimizer>,
    pub minimizers: Vec<Minimizer>,
    /// Chaining buffers (passed to chain_anchors/chain_anchors_rmq).
    pub chain_bufs: ChainingBuffers,
    // Seeding buffer
    pub mini_pos: Vec<u64>,
    /// Scratch buffers reused across collect_seed_hits* calls.
    pub seed_scratch: crate::align::seed::SeedScratch,
}

impl Default for MapContext {
    fn default() -> Self {
        Self::new()
    }
}

impl MapContext {
    pub fn new() -> Self {
        MapContext {
            anchors: Vec::with_capacity(4096),
            minimizers: Vec::with_capacity(4096),
            chain_bufs: ChainingBuffers::new(),
            mini_pos: Vec::with_capacity(1024),
            seed_scratch: crate::align::seed::SeedScratch::new(),
        }
    }
}

#[derive(Clone)]
pub struct Mapping {
    pub ref_id: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub score: i32,
    pub initial_chain_score: i32,     // initial chain score before splitting
    pub anchor_count: i32,
    pub anchors: Vec<Minimizer>,
    pub is_reverse: bool,
    pub s2_score: Option<i32>,
    pub is_secondary: bool,
    pub hash: u32,
    /// Pre-alignment n_sub from parent assignment (persists for the root chain).
    pub pre_num_suboptimal: i32,
    /// Nearby-seed bounds for extension boundary tightening
    pub left_bound_rs1: i32,
    pub left_bound_qs1: i32,
    pub right_bound_re1: i32,
    pub right_bound_qe1: i32,
    /// Original offset into the sorted anchor array
    pub original_as: usize,
    /// Offset of this chain's anchors in the squeezed array (for z-drop split bound computation)
    pub sq_start: usize,
    /// Estimated divergence from estimate_divergence
    pub div: f32,
    /// Whether this chain was retained only due to opposite-strand retention
    /// in secondary filtering (would otherwise have been dropped).
    pub strand_retained: bool,
    /// Parent index in the compacted array (for filter_strand_retained)
    pub compact_parent: usize,
    /// True if this split part was created when zdrop detected an inversion.
    pub split_inv: bool,
    /// Whether this chain maps to an ALT contig
    pub is_alt: bool,
}

impl Filterable for Mapping {
    fn query_start(&self) -> usize { self.query_start }
    fn query_end(&self) -> usize { self.query_end }
    fn score(&self) -> i32 { self.score }
    fn is_reverse(&self) -> bool { self.is_reverse }
    fn is_alt(&self) -> bool { self.is_alt }
}

use crate::align::stats::AlignmentStats;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;
#[cfg(target_arch = "wasm32")]
use web_time::Instant;

/// Compute nearby-seed extension bounds from the squeezed anchor array.
/// Returns (left_rs1, left_qs1, right_re1, right_qe1).
/// Used for both initial chains and z-drop split sub-chains.
pub fn compute_bounds_from_squeezed(
    squeezed: &[Minimizer],
    sq_start: usize,
    sq_cnt: usize,
    tlen: i32,
    qlen: i32,
    min_cnt: i32,
) -> (i32, i32, i32, i32) {
    if sq_cnt == 0 { return (0, 0, tlen, qlen); }

    let first = &squeezed[sq_start];
    let last = &squeezed[sq_start + sq_cnt - 1];
    let rid_strand = first.ref_id_strand();
    let first_q_span = first.query_span();
    let rs0_pre = first.ref_pos() + 1 - first_q_span;
    let qs0_pre = first.query_pos() + 1 - first_q_span;


    // Left nearby-seed inspection.
    let mut left_rs1 = 0i32;
    let mut left_qs1 = 0i32;
    {
        let mut l = 0i32;
        let mut j = sq_start as i64 - 1;
        while j >= 0 {
            let aj = &squeezed[j as usize];
            if aj.ref_id_strand() != rid_strand {
                break;
            }
            let span_j = aj.query_span();
            let xj = aj.ref_pos() + 1 - span_j;
            let yj = aj.query_pos() + 1 - span_j;
            let qual = xj < rs0_pre && yj < qs0_pre;
            if qual {
                l += 1;
                if l > min_cnt {
                    let dist = std::cmp::max(rs0_pre - xj, qs0_pre - yj);
                    left_rs1 = rs0_pre - dist;
                    left_qs1 = qs0_pre - dist;
                    if left_rs1 < 0 { left_rs1 = 0; }
                    break;
                }
            }
            j -= 1;
        }
    }

    // Right nearby-seed inspection.
    let re0_pre = last.ref_pos() + 1;
    let qe0_pre = last.query_pos() + 1;
    let mut right_re1 = tlen;
    let mut right_qe1 = qlen;
    {
        let mut l = 0i32;
        let mut j = sq_start + sq_cnt;
        while j < squeezed.len() {
            let aj = &squeezed[j];
            if aj.ref_id_strand() != rid_strand { break; }
            let xj = aj.ref_pos() + 1;
            let yj = aj.query_pos() + 1;
            if xj > re0_pre && yj > qe0_pre {
                l += 1;
                if l > min_cnt {
                    let dist = std::cmp::max(xj - re0_pre, yj - qe0_pre);
                    right_re1 = re0_pre + dist;
                    right_qe1 = qe0_pre + dist;
                    break;
                }
            }
            j += 1;
        }
    }

    (left_rs1, left_qs1, right_re1, right_qe1)
}

/// Get forward-strand query position from an anchor
fn get_forward_query_pos(qlen: usize, a: &Minimizer) -> i32 {
    let x = a.query_pos(); // lower 32 bits = query position
    let q_span = a.query_span();
    if (a.x >> 63) != 0 {
        // Reverse strand: revert to forward
        (qlen as i32) - 1 - (x + 1 - q_span)
    } else {
        x
    }
}

/// Binary search for a forward query position in the sorted mini_pos array
fn get_minimizer_index(qlen: usize, a: &Minimizer, mini_pos: &[u64]) -> Option<usize> {
    let x = get_forward_query_pos(qlen, a);
    // Binary search on lower 32 bits of mini_pos
    let mut l: i64 = 0;
    let mut r: i64 = mini_pos.len() as i64 - 1;
    while l <= r {
        let m = ((l + r) as u64 >> 1) as i64;
        let y = mini_pos[m as usize] as i32; // lower 32 bits
        if y < x { l = m + 1; }
        else if y > x { r = m - 1; }
        else { return Some(m as usize); }
    }
    None
}

/// Estimate chain divergence from minimizer matching
fn estimate_divergence(regs: &mut [Mapping], _anchors: &[Minimizer], qlen: usize, mini_pos: &[u64], mi: &Index) {
    if mini_pos.is_empty() { return; }

    // Compute average k-mer span
    let sum_k: u64 = mini_pos.iter().map(|&p| (p >> 32) & 0xff).sum();
    let avg_k = sum_k as f32 / mini_pos.len() as f32;
    let n = mini_pos.len();

    for reg in regs.iter_mut() {
        // reg.div is already set to -1.0 in constructor
        if reg.anchor_count == 0 || reg.anchors.is_empty() { continue; }

        // Find the first anchor in forward query order
        let first_anchor = if reg.is_reverse {
            &reg.anchors[reg.anchors.len() - 1]
        } else {
            &reg.anchors[0]
        };

        let st = match get_minimizer_index(qlen, first_anchor, mini_pos) {
            Some(idx) => idx,
            None => continue,
        };

        // Walk through anchors and mini_pos to count matches
        let cnt = reg.anchors.len();
        let mut en = st;
        let mut n_match = 1i32;
        let mut k = 1usize;
        let mut j = st + 1;
        while j < n && k < cnt {
            let anchor = if reg.is_reverse {
                &reg.anchors[cnt - 1 - k]
            } else {
                &reg.anchors[k]
            };
            let x = get_forward_query_pos(qlen, anchor);
            if x == (mini_pos[j] as i32) {
                k += 1;
                en = j;
                n_match += 1;
            }
            j += 1;
        }

        let mut n_tot = (en - st + 1) as i32;
        if reg.query_start as f32 > avg_k && reg.ref_start as f32 > avg_k { n_tot += 1; }
            // Check divergence estimate validity
        let l_ref = mi.seqs[reg.ref_id].len;
        if (qlen - reg.query_start) as f32 > avg_k && (l_ref - reg.ref_end) as f32 > avg_k { n_tot += 1; }

        // Compute divergence (esterr.c:62)
        reg.div = if n_match >= n_tot {
            0.0f32
        } else {
            (1.0 - (n_match as f64 / n_tot as f64).powf(1.0 / avg_k as f64)) as f32
        };
    }
}

/// Filter strand-retained chains with high divergence
fn filter_strand_retained(regs: &mut Vec<Mapping>) {
    // Collect parent div values first (since we can't borrow regs mutably while iterating)
    let divs: Vec<f32> = regs.iter().map(|r| r.div).collect();
    let parents: Vec<usize> = regs.iter().map(|r| r.compact_parent).collect();
    let strand_retained: Vec<bool> = regs.iter().map(|r| r.strand_retained).collect();

    let mut keep = vec![true; regs.len()];
    for i in 0..regs.len() {
        if !strand_retained[i] { continue; }
        let p = parents[i];
        let parent_div = if p < divs.len() { divs[p] } else { 0.0 };
        // Keep if div < parent_div * 5.0 or div < 0.01.
        if divs[i] < parent_div * 5.0 || divs[i] < 0.01 {
            continue; // keep
        }
        keep[i] = false;
    }

    let mut j = 0;
    for i in 0..regs.len() {
        if keep[i] {
            if j != i {
                regs[j] = regs[i].clone();
            }
            j += 1;
        }
    }
    regs.truncate(j);
}

/// A chain's span on the target and query, plus the hash identifying it.
struct ChainCoords {
    rid: usize,
    rev: bool,
    rs: usize,
    re: usize,
    qs: usize,
    qe: usize,
    hash: u32,
}

/// Derive a chain's coordinates and identity hash from its anchors, which must
/// be exactly the chain's own (the hash mixes in their count). `qlen` flips
/// reverse-strand chains into query orientation; `is_qstrand` keeps them in
/// chain orientation instead.
fn chain_coords(
    anchors: &[Minimizer],
    qlen: usize,
    is_qstrand: bool,
    read_hash: u32,
) -> ChainCoords {
    let first = &anchors[0];
    let last = &anchors[anchors.len() - 1];
    let span = first.query_span() as usize;
    let rev = (first.x >> 63) == 1;
    // Guard against a seed span reaching before the target start.
    let rs = ((first.ref_pos() as usize) + 1).saturating_sub(span);
    let re = (last.ref_pos() as usize) + 1;
    let qs_fwd = (first.query_pos() as usize) + 1 - span;
    let qe_fwd = (last.query_pos() as usize) + 1;
    let (qs, qe) = if !rev || is_qstrand {
        (qs_fwd, qe_fwd)
    } else {
        (qlen.saturating_sub(qe_fwd), qlen.saturating_sub(qs_fwd))
    };
    let h = hash64(hash64(first.x).wrapping_add(hash64(first.y)) ^ (read_hash as u64)) as u32;
    ChainCoords {
        rid: first.ref_id(),
        rev, rs, re, qs, qe,
        hash: (anchors.len() as u32) ^ h,
    }
}

pub fn map_query(
    opt: &MapOptions,
    mi: &Index,
    qname: &str,
    qseq: &[u8],
    ctx: &mut MapContext,
) -> (Vec<Mapping>, usize, AlignmentStats, Vec<Minimizer>) {
    let mut stats = AlignmentStats { n_reads: 1, ..Default::default() };

    let qlen = qseq.len();
    let read_hash = compute_read_hash(qname, qlen, opt.filtering.seed as u32, opt.flags);

    let t0 = Instant::now();
    ctx.minimizers.clear();
    sketch_sequence(qseq, qlen, mi.window_size, mi.kmer_size, 0, &mut ctx.minimizers);
    let qsub = opt.seeding.query_subsample.or_else(|| crate::align::env_flags::EXP_QSUB.map(|f| (f.max(1) as u32,
        std::env::var("BONEFISH_EXP_QSUB_MIN").ok().and_then(|v| v.parse().ok()).unwrap_or(0))));
    if let Some((f, qsub_min)) = qsub && qlen >= qsub_min as usize {
        // Keep 1/f of the minimizer hash space (deterministic, as if the index were sparser).
        let f = f.max(1) as u64;
        ctx.minimizers.retain(|m| ((m.x >> 8).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 32) % f == 0);
    }

    // Drop query minimizers that appear too many times in the query itself.
    if opt.seeding.q_occ_frac > 0.0 && opt.seeding.mid_occ > 0 && ctx.minimizers.len() > opt.seeding.mid_occ {
        filter_minimizers_by_occ(&mut ctx.minimizers, opt.seeding.mid_occ, opt.seeding.q_occ_frac);
    }

    stats.t_sketch = t0.elapsed();
    stats.n_seeds = ctx.minimizers.len();

    let t1 = Instant::now();
    // collect_seed_hits clears ctx.anchors internally now
    ctx.mini_pos.clear();
    let qn_opt = Some(qname);
    let mut rep_len = collect_seed_hits(opt, mi, qlen, &ctx.minimizers, &mut ctx.anchors, &mut ctx.mini_pos, qn_opt, &mut ctx.seed_scratch);
    stats.t_seed = t1.elapsed();
    stats.n_anchors = ctx.anchors.len();

    let t2 = Instant::now();

    // Compute max chaining gaps
    let max_chain_gap_qry = opt.chaining.max_gap;
    let max_chain_gap_ref = if opt.chaining.max_gap_ref > 0 {
        opt.chaining.max_gap_ref
    } else {
        opt.chaining.max_gap
    };

    // Primary chaining: cDNA-mode DP for spliced presets, RMQ when RMQ_CHAIN is set (asm
    // presets), otherwise DP
    let is_splice = opt.flags.contains(AlignFlags::SPLICE);
    let mut anchors = std::mem::take(&mut ctx.anchors);
    let (mut u, mut chains) = if is_splice {
        chain_anchors_cdna(&opt.chaining, max_chain_gap_ref, max_chain_gap_qry, &mut anchors, &mut ctx.chain_bufs)
    } else if opt.flags.contains(AlignFlags::RMQ_CHAIN) {
        chain_anchors_rmq(&opt.chaining, &mut anchors, &mut ctx.chain_bufs)
    } else {
        // Use partitioned parallel chaining when there are enough anchors
        // to benefit from parallelism across reference sequences.
        // Partitioned parallel chaining: split anchors by ref_id and chain
        // each partition independently. Only when PAR_CHAIN flag is set
        // (asm2asm workloads with many anchors across many references).
        #[cfg(feature = "parallel")]
        if opt.flags.contains(AlignFlags::PAR_CHAIN) && anchors.len() >= 5000 {
            chain_anchors_partitioned(&opt.chaining, max_chain_gap_ref, max_chain_gap_qry, &mut anchors)
        } else {
            chain_anchors(&opt.chaining, max_chain_gap_ref, max_chain_gap_qry, &mut anchors, &mut ctx.chain_bufs)
        }
        #[cfg(not(feature = "parallel"))]
        chain_anchors(&opt.chaining, max_chain_gap_ref, max_chain_gap_qry, &mut anchors, &mut ctx.chain_bufs)
    };
    // A reference gap limit below -g (--max-gap-ref) can leave a read with sparse anchors
    // without any chain: chain it again with -g.
    if u.is_empty() && !anchors.is_empty() && opt.chaining.max_gap_ref > 0 && max_chain_gap_ref < opt.chaining.max_gap
        && !opt.flags.contains(AlignFlags::RMQ_CHAIN) && !is_splice {
        (u, chains) = chain_anchors(&opt.chaining, opt.chaining.max_gap, max_chain_gap_qry, &mut anchors, &mut ctx.chain_bufs);
    }
    ctx.anchors = anchors;

    stats.t_chain = t2.elapsed();
    stats.n_chains = chains.len();

    let t3 = Instant::now();
    let pc0 = Instant::now();

    // RMQ rescue re-chaining
    // Only when: bw_long > bw, not no_ljoin or splice, multiple chains
    // Does NOT exclude RMQ_CHAIN
    let skip_rescue = opt.flags.intersects(AlignFlags::NO_LJOIN | AlignFlags::SPLICE);
    if opt.chaining.bandwidth_long > opt.chaining.bandwidth && u.len() > 1 && !skip_rescue {
        let first_chain_cnt = (u[0] & 0xFFFFFFFF) as usize;
        if first_chain_cnt > 0 {
            let first_anchor = &chains[0];
            let last_anchor = &chains[first_chain_cnt - 1];
            let st = first_anchor.query_pos();
            let en = last_anchor.query_pos();
            let coverage = (en - st) as usize;

            // RMQ rescue: check if covered range is small or anchors span long query region
            let needs_rescue = ((qlen.saturating_sub(coverage) > opt.chaining.rmq_rescue_size as usize)
                || (coverage > (qlen as f32 * opt.chaining.rmq_rescue_ratio) as usize))
                && crate::align::env_flags::EXP_RESCUE_COVER.is_none_or(|f| (coverage as f32) < f * qlen as f32)
                && (opt.chaining.rescue_max_density == 0 || ctx.anchors.len() < opt.chaining.rescue_max_density as usize * qlen);
            let needs_rescue = needs_rescue && (!(opt.chaining.rescue_join || *crate::align::env_flags::EXP_RESCUE_JOIN) || {
                let max_dist = opt.chaining.max_gap.max(opt.chaining.bandwidth_long);
                let j = chains_joinable(&u, &chains, max_dist, opt.chaining.bandwidth_long);
                if !j { PC_RESCUE_SKIP.fetch_add(1, Relaxed); }
                j
            });

            if needs_rescue {
                PC_RESCUE_COUNT.fetch_add(1, Relaxed);
                PC_RESCUE_ANCHORS.fetch_add(chains.len() as u64, Relaxed);

                // chains is the output of the prior chain pass and is no longer
                // needed once we reassign below, so re-use it in place: truncate to
                // n_a, sort by (ref_id_strand, x), pass to the rescue call. Saves
                // an n_a × 16-byte clone per rescue trigger.
                let n_a: usize = u.iter().map(|&v| (v & 0xFFFFFFFF) as usize).sum();
                chains.truncate(n_a);
                radix_sort_128x(&mut chains);

                let rescue_params = ChainingParams { bandwidth: opt.chaining.bandwidth_long, ..opt.chaining.clone() };
                let (u2, chains2) = chain_anchors_rmq(
                    &rescue_params, &mut chains, &mut ctx.chain_bufs,
                );

                u = u2;
                chains = chains2;
            }
        }
    } else if opt.seeding.max_occ > opt.seeding.mid_occ && rep_len > 0 && !opt.flags.contains(AlignFlags::RMQ_CHAIN) {
        // High-occ re-chaining: only re-chain when no chains were found.
        let rechain = u.is_empty();

        if rechain {
            let rc0 = Instant::now();
            PC_RECHAIN_COUNT.fetch_add(1, Relaxed);
            rep_len = collect_seed_hits_with_occ(opt, mi, qlen, &ctx.minimizers, &mut ctx.anchors, &mut ctx.mini_pos, opt.seeding.max_occ, qn_opt, &mut ctx.seed_scratch);

            let mut anchors = std::mem::take(&mut ctx.anchors);
            let (u2, chains2) = if is_splice {
                chain_anchors_cdna(&opt.chaining, max_chain_gap_ref, max_chain_gap_qry, &mut anchors, &mut ctx.chain_bufs)
            } else {
                chain_anchors(&opt.chaining, max_chain_gap_ref, max_chain_gap_qry, &mut anchors, &mut ctx.chain_bufs)
            };
            ctx.anchors = anchors;
            u = u2;
            chains = chains2;
            PC_RECHAIN_NS.fetch_add(rc0.elapsed().as_nanos() as u64, Relaxed);
        }
    }

    let pc1 = Instant::now();
    PC_RESCUE_NS.fetch_add((pc1 - pc0).as_nanos() as u64, Relaxed);

    let mut regs = Vec::with_capacity(u.len());

    let mut k_offset = 0;
    for &u_val in u.iter() {
        let sc = (u_val >> 32) as i32;
        let cnt = (u_val & 0xFFFFFFFF) as i32;

        if cnt == 0 { continue; }

        // Collect chain anchors
        let chain_anchors_vec = chains[k_offset..k_offset + cnt as usize].to_vec();

        let is_qstrand = opt.flags.contains(AlignFlags::QSTRAND);
        let ChainCoords { rid, rev, rs, re, qs, qe, hash: chain_hash } =
            chain_coords(&chain_anchors_vec, qseq.len(), is_qstrand, read_hash);

        regs.push(Mapping {
            ref_id: rid,
            query_start: qs,
            query_end: qe,
            ref_start: rs,
            ref_end: re,
            score: sc,
            initial_chain_score: sc,  // initially same as score
            anchor_count: cnt,
            anchors: chain_anchors_vec,
            is_reverse: rev,
            s2_score: None,
            is_secondary: false,
            hash: chain_hash,
            pre_num_suboptimal: 0,
            left_bound_rs1: 0,
            left_bound_qs1: 0,
            right_bound_re1: mi.seqs[rid].len as i32,
            right_bound_qe1: qseq.len() as i32,
            original_as: k_offset,
            sq_start: 0, // set during squeezed array construction
            div: -1.0,
            strand_retained: false,
            compact_parent: 0,
            split_inv: false,
            is_alt: false,
        });

        k_offset += cnt as usize;
    }
    
    // Filter alignments
    if opt.chaining.min_chain_score > 0 {
        regs.retain(|r| r.score >= opt.chaining.min_chain_score);
    }

    // Strand filtering (FOR_ONLY / REV_ONLY)
    if opt.flags.contains(AlignFlags::FOR_ONLY) {
        regs.retain(|r| !r.is_reverse);
    }
    if opt.flags.contains(AlignFlags::REV_ONLY) {
        regs.retain(|r| r.is_reverse);
    }

    // Mark ALT-contig chains; then sort chains by descending score.
    let n_alt = mi.seqs.iter().filter(|s| s.is_alt).count();
    if n_alt > 0 {
        for r in regs.iter_mut() {
            if mi.seqs[r.ref_id].is_alt { r.is_alt = true; }
        }
        let alt_drop = opt.filtering.alt_drop;
        regs.sort_by(|a, b| {
            let sa = if a.is_alt { scale_alt_score(a.score, alt_drop) } else { a.score };
            let sb = if b.is_alt { scale_alt_score(b.score, alt_drop) } else { b.score };
            sb.cmp(&sa).then_with(|| b.hash.cmp(&a.hash))
        });
    } else {
        // Sort by (score, hash) descending, by (score, hash) descending
        regs.sort_by(|a, b| b.score.cmp(&a.score).then_with(|| b.hash.cmp(&a.hash)));
    }

    let pc2 = Instant::now();
    PC_BUILD_REGS_NS.fetch_add((pc2 - pc1).as_nanos() as u64, Relaxed);

    if regs.is_empty() { stats.t_post = t3.elapsed(); return (regs, rep_len, stats, Vec::new()); }

    // Skip parent-based filtering when ALL_CHAINS is set (e.g. ava-ont, ava-pb)
    if opt.flags.contains(AlignFlags::ALL_CHAINS) {
        // Nearby-seed computation on post-filter chains
        // Build squeezed anchor array from surviving chains, sorted by position
        let mut pos_order: Vec<usize> = (0..regs.len()).collect();
        pos_order.sort_by_key(|&i| regs[i].original_as);

        let mut squeezed: Vec<Minimizer> = Vec::new();
        let mut sq_offsets: Vec<(usize, usize)> = vec![(0, 0); regs.len()];
        for &pi in &pos_order {
            let start = squeezed.len();
            squeezed.extend_from_slice(&regs[pi].anchors);
            sq_offsets[pi] = (start, regs[pi].anchors.len());
            regs[pi].sq_start = start;
        }

        {
            for ri in 0..regs.len() {
                let (sq_start, sq_cnt) = sq_offsets[ri];
                if sq_cnt == 0 { continue; }
                let first = &squeezed[sq_start];
                let last = &squeezed[sq_start + sq_cnt - 1];
                let rid_strand = first.ref_id_strand();
                let first_q_span = first.query_span();
                let rs0_pre = first.ref_pos() + 1 - first_q_span;
                let qs0_pre = first.query_pos() + 1 - first_q_span;
                let last_re0 = last.ref_pos() + 1;
                let last_qe0 = last.query_pos() + 1;
                let tlen = mi.seqs[regs[ri].ref_id].len as i32;
                let qlen_i32 = qlen as i32;

                // Left nearby-seed
                let mut left_rs1 = 0i32;
                let mut left_qs1 = 0i32;
                {
                    let mut l = 0i32;
                    let mut j = sq_start as i64 - 1;
                    while j >= 0 {
                        let aj = &squeezed[j as usize];
                        if aj.ref_id_strand() != rid_strand { break; }
                        let span_j = aj.query_span();
                        let xj = aj.ref_pos() + 1 - span_j;
                        let yj = aj.query_pos() + 1 - span_j;
                        if xj < rs0_pre && yj < qs0_pre {
                            l += 1;
                            if l > opt.chaining.min_cnt {
                                let dist = std::cmp::max(rs0_pre - xj, qs0_pre - yj);
                                left_rs1 = rs0_pre - dist;
                                left_qs1 = qs0_pre - dist;
                                if left_rs1 < 0 { left_rs1 = 0; }
                                if left_qs1 < 0 { left_qs1 = 0; }
                                break;
                            }
                        }
                        j -= 1;
                    }
                }

                // Right nearby-seed
                let mut right_re1 = tlen;
                let mut right_qe1 = qlen_i32;
                {
                    let mut l = 0i32;
                    let sq_end = sq_start + sq_cnt;
                    let mut j = sq_end;
                    while j < squeezed.len() {
                        let aj = &squeezed[j];
                        if aj.ref_id_strand() != rid_strand { break; }
                        let xj = aj.ref_pos() + 1;
                        let yj = aj.query_pos() + 1;
                        if xj > last_re0 && yj > last_qe0 {
                            l += 1;
                            if l > opt.chaining.min_cnt {
                                let dist = std::cmp::max(xj - last_re0, yj - last_qe0);
                                right_re1 = std::cmp::min(last_re0 + dist, tlen);
                                right_qe1 = std::cmp::min(last_qe0 + dist, qlen_i32);
                                break;
                            }
                        }
                        j += 1;
                    }
                }

                regs[ri].left_bound_rs1 = left_rs1;
                regs[ri].left_bound_qs1 = left_qs1;
                regs[ri].right_bound_re1 = right_re1;
                regs[ri].right_bound_qe1 = right_qe1;
            }
        }

        // estimate_divergence runs unconditionally after chain_post, even for ALL_CHAINS.
        estimate_divergence(&mut regs, &ctx.anchors, qlen, &ctx.mini_pos, mi);

        stats.t_post = t3.elapsed();
        return (regs, rep_len, stats, squeezed);
    }

    // Parent-based filtering
    // Uses shared filter module to avoid code duplication with pipeline.rs
    let filter_params = FilterParams::new(opt, mi);
    let mut parent_state = ParentState::new(regs.len(), filter_params.mask_level, filter_params.mask_len, filter_params.hard_mask_level);
    parent_state.init_from_items(&regs);
    parent_state.assign_parents(&regs);

    // Compute subsc and n_sub per parent before filtering; these persist into
    // the post-alignment parent assignment call.
    let mut subsc_pre = vec![0i32; regs.len()];
    let mut n_sub_pre = vec![0i32; regs.len()];
    for i in 0..regs.len() {
        let pi = parent_state.parent[i];
        if pi != i {
            // ALT penalty for subsc (pre-alignment, chain score only).
            let sci = if !regs[pi].is_alt && regs[i].is_alt {
                scale_alt_score(regs[i].score, opt.filtering.alt_drop)
            } else {
                regs[i].score
            };
            if sci > subsc_pre[pi] {
                subsc_pre[pi] = sci;
            }
            // Pre-alignment cnt_sub uses anchor counts only (no dp_max yet).
            if regs[i].anchor_count >= regs[pi].anchor_count {
                n_sub_pre[pi] += 1;
            }
        }
    }
    // Write subsc and n_sub into regs (for parents)
    for i in 0..regs.len() {
        if parent_state.is_primary(i) {
            if subsc_pre[i] > 0 {
                regs[i].s2_score = Some(subsc_pre[i]);
            }
            regs[i].pre_num_suboptimal = n_sub_pre[i];
        }
    }

    // Two-pass secondary filtering. Parent score/strand lookups in the first pass
    // read from the original (unmodified) regs. Compaction happens in the second
    // pass so no parent read ever sees a slot overwritten by a later entry.
    let n = regs.len();
    let mut keep = vec![false; n];
    let mut n_second = 0i32;
    for i in 0..n {
        let p = parent_state.parent[i];
        if p == i {
            keep[i] = true;
            continue;
        }
        let p_score = regs[p].score;
        let p_rev = regs[p].is_reverse;
        let filter_result = check_secondary_filter(
            regs[i].score,
            regs[i].is_reverse,
            p_score,
            p_rev,
            &filter_params,
            true, // check_strand=true for pre-alignment filtering.
        );

        if filter_result.passes && n_second < opt.filtering.best_n {
            let p_qs = regs[p].query_start;
            let p_qe = regs[p].query_end;
            let p_rid = regs[p].ref_id;
            let p_rs = regs[p].ref_start;
            let p_re = regs[p].ref_end;
            let identical = regs[i].query_start == p_qs && regs[i].query_end == p_qe
                && regs[i].ref_id == p_rid && regs[i].ref_start == p_rs && regs[i].ref_end == p_re;
            if !identical {
                regs[i].is_secondary = true;
                if filter_result.passes_strand && !filter_result.passes_ratio && !filter_result.passes_min_diff {
                    regs[i].strand_retained = true;
                }
                keep[i] = true;
                n_second += 1;
            }
        }
    }

    // Compact regs based on `keep[]`. Build the orig→compacted index map first,
    // then populate compact_parent using parent_state (which indexes original slots).
    let mut orig_to_compact: Vec<usize> = vec![usize::MAX; n];
    let mut k = 0usize;
    for i in 0..n {
        if keep[i] { orig_to_compact[i] = k; k += 1; }
    }
    let mut compacted: Vec<Mapping> = Vec::with_capacity(k);
    for (orig, r) in regs.drain(..).enumerate() {
        if !keep[orig] { continue; }
        let mut r = r;
        let orig_parent = parent_state.parent[orig];
        r.compact_parent = orig_to_compact[orig_parent];
        compacted.push(r);
    }
    regs = compacted;
    let mut filtered_regs = regs;
    let pc3 = Instant::now();
    PC_PARENT_NS.fetch_add((pc3 - pc2).as_nanos() as u64, Relaxed);
    // estimate_divergence + filter_strand_retained
    estimate_divergence(&mut filtered_regs, &ctx.anchors, qlen, &ctx.mini_pos, mi);
    let has_strand_retained = filtered_regs.iter().any(|r| r.strand_retained);
    if has_strand_retained {
        // Need div values computed by estimate_divergence before filtering
        filter_strand_retained(&mut filtered_regs);
    }
    let pc4 = Instant::now();
    PC_DIVEST_NS.fetch_add((pc4 - pc3).as_nanos() as u64, Relaxed);
    // Nearby-seed computation on post-filter chains
    // Build squeezed anchor array from surviving chains, sorted by position
    // Sort surviving regs by original anchor array offset (matching mm_squeeze_a)
    // mm_squeeze_a sorts chains by regs[i].as, the position of their first anchor
    // in the original sorted anchor array. We stored this as original_as during chaining.
    let mut pos_order: Vec<usize> = (0..filtered_regs.len()).collect();
    pos_order.sort_by_key(|&i| filtered_regs[i].original_as);

    // Build squeezed anchor array and track offsets
    let mut squeezed: Vec<Minimizer> = Vec::new();
    let mut sq_offsets: Vec<(usize, usize)> = vec![(0, 0); filtered_regs.len()]; // (start, count) in squeezed
    for &pi in &pos_order {
        let start = squeezed.len();
        squeezed.extend_from_slice(&filtered_regs[pi].anchors);
        sq_offsets[pi] = (start, filtered_regs[pi].anchors.len());
        filtered_regs[pi].sq_start = start;
    }

    // Compute nearby-seed bounds for each surviving chain
    for ri in 0..filtered_regs.len() {
        let (sq_start, sq_cnt) = sq_offsets[ri];
        if sq_cnt == 0 { continue; }
        let rid = filtered_regs[ri].ref_id;
        let tlen_i32 = mi.seqs[rid].len as i32;
        let qlen_i32 = qseq.len() as i32;
        let (left_rs1, left_qs1, right_re1, right_qe1) =
            compute_bounds_from_squeezed(&squeezed, sq_start, sq_cnt, tlen_i32, qlen_i32, opt.chaining.min_cnt);
        filtered_regs[ri].left_bound_rs1 = left_rs1;
        filtered_regs[ri].left_bound_qs1 = left_qs1;
        filtered_regs[ri].right_bound_re1 = right_re1;
        filtered_regs[ri].right_bound_qe1 = right_qe1;
    }

    PC_SQUEEZED_NS.fetch_add((Instant::now() - pc4).as_nanos() as u64, Relaxed);
    stats.t_post = t3.elapsed();
    (filtered_regs, rep_len, stats, squeezed)
}

/// Whether two of the chains (`u`: score << 32 | count, anchors in `chains` chain by chain,
/// each in reference order) could be joined by long-gap chaining: same reference and strand,
/// the second starting after the first in both reference and query (overlap allowed) with
/// both gaps within `max_dist` and the diagonal shift within `bw`.
fn chains_joinable(u: &[u64], chains: &[Minimizer], max_dist: i32, bw: i32) -> bool {
    joinable_pairs(u, chains, max_dist, bw, true).is_some()
}

/// Shared pair search: chains sorted by (reference/strand, start, chain index) so each chain is
/// only compared with later-starting chains within `max_dist` of its end (the pair test needs
/// the second chain to start after the first and within `max_dist` of its end). With `first`,
/// returns as soon as one pair is found (the flags are then incomplete).
fn joinable_pairs(u: &[u64], chains: &[Minimizer], max_dist: i32, bw: i32, first: bool) -> Option<Vec<bool>> {
    let mut ext: Vec<(u64, i32, i32, i32, i32, usize)> = Vec::with_capacity(u.len()); // (rid_strand, rs, re, qs, qe, chain)
    let mut off = 0usize;
    for &v in u {
        let n = (v & 0xFFFFFFFF) as usize;
        if n == 0 { continue; }
        let (a0, a1) = (&chains[off], &chains[off + n - 1]);
        ext.push((a0.ref_id_strand(), a0.ref_pos(), a1.ref_pos(), a0.query_pos(), a1.query_pos(), ext.len()));
        off += n;
    }
    let n = ext.len();
    ext.sort_unstable_by_key(|e| (e.0, e.1, e.5));
    let mut member = vec![false; n];
    let mut any = false;
    for i in 0..n {
        let f = ext[i];
        for s in &ext[i + 1..] {
            // Sorted: same reference/strand and start >= f's start; later chains start no earlier.
            if s.0 != f.0 || s.1 as i64 - f.2 as i64 > max_dist as i64 { break; }
            // s must start after f starts and end after f ends, in both reference and query.
            if !(s.3 > f.3 && s.2 > f.2 && s.4 > f.4) { continue; }
            let (dr, dq) = (s.1 - f.2, s.3 - f.4);
            if dr <= max_dist && dq <= max_dist && (dr - dq).abs() <= bw {
                if first { return Some(member); }
                member[f.5] = true;
                member[s.5] = true;
                any = true;
            }
        }
    }
    any.then_some(member)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::index::Index;

    #[test]
    fn test_repetitive_filtering() {
        // Create an index where a k-mer appears many times
        let mut seqs = Vec::new();
        // 20 sequences, each containing "AAAAAAAAAAAAAAAA"
        for i in 0..20 {
             seqs.push((format!("seq{}", i), "AAAAAAAAAAAAAAAA".as_bytes().to_vec()));
        }
        // Build index with default parameters (w=10, k=15). 
        // 16 A's -> 2 k-mers of len 15 (overlapping). 
        // All seqs have identical content.
        let idx = Index::build(seqs, 10, 15, 50000);
        
        // A single A...A k-mer will appear 20+ times.
        
        let mut opt = MapOptions::default();
        opt.chaining.min_cnt = 1; // Allow single-anchor chains
        opt.chaining.min_chain_score = 10; // Allow low-scoring chains
        opt.seeding.mid_occ = 10; // Aggressive filter
        opt.seeding.max_occ = 10; // Same as mid_occ to prevent re-chaining path
        opt.seeding.occ_dist = 0; // Disable select_seeds, use simple mid_occ filtering

        // Need ctx
        let mut ctx = MapContext::new();

        let query = "AAAAAAAAAAAAAAAA";
        let (res, _, _, _) = map_query(&opt, &idx, "query", query.as_bytes(), &mut ctx);
        assert_eq!(res.len(), 0, "Should filter out repetitive k-mers with mid_occ=10");
        
        opt.seeding.mid_occ = 50; // Relaxed filter
        let (res, _, _, _) = map_query(&opt, &idx, "query", query.as_bytes(), &mut ctx);

        assert!(!res.is_empty(), "Should align repetitive k-mers with mid_occ=50 (found {})", res.len());
    }

    #[test]
    fn test_secondary_alignment() {
        // Create 2 identical sequences
        let seqs = vec![
            ("seq1".to_string(), "ACGTACGTACGTACGT".as_bytes().to_vec()),
            ("seq2".to_string(), "ACGTACGTACGTACGT".as_bytes().to_vec()),
        ];
        let idx = Index::build(seqs, 10, 15, 50000);
        
        let mut opt = MapOptions::default();
        opt.chaining.min_cnt = 1; 
        opt.chaining.min_chain_score = 10;
        opt.filtering.best_n = 5; // Allow secondaries
        
        let mut ctx = MapContext::new();
        let query = "ACGTACGTACGTACGT";
        let (res, _, _, _) = map_query(&opt, &idx, "query", query.as_bytes(), &mut ctx);
        
        // Should find 2 alignments
        assert_eq!(res.len(), 2, "Should find 2 alignments (1 primary, 1 secondary)");
        assert!(!res[0].is_secondary, "First should be primary");
        assert!(res[1].is_secondary, "Second should be secondary (overlapping)");
    }
}
