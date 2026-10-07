use clap::Parser;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::io::Write;
use std::time::Instant;

// Global allocator: Rust system allocator by default. Opt in to jemalloc or
// mimalloc via `--features jemalloc` / `--features mimalloc` (both need a C
// toolchain at build time). mimalloc wins precedence if both are set.
#[cfg(all(not(target_arch = "wasm32"), feature = "mimalloc"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[cfg(all(not(target_arch = "wasm32"), feature = "jemalloc", not(feature = "mimalloc")))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use rammap::align::index::Index;
use rammap::align::map::{MapOptions, AlignFlags};
use rammap::align::pipeline::{OutputConfig, ReadInfo};
use rammap::align::stats::AlignmentStats;
use rammap::align::junc::{self, JunctionDb};
use rammap::align::jump::JumpDb;
use rammap::align::split;

// ─── CLI definition ───

#[derive(Parser, Debug)]
#[command(name = "rammap")]
#[command(version, about = "A pure-Rust minimap2-compatible sequence aligner")]
struct Cli {
    #[command(flatten)]
    align: AlignArgs,

    /// Verbose output
    #[arg(short, long, global = true)]
    verbose: bool,
}

#[derive(Parser, Debug)]
struct AlignArgs {
    /// Target FASTA file (or .mmi/.idx index)
    pub target: String,

    /// Query FASTA/FASTQ file(s) (- for stdin). Multiple files supported.
    #[arg(trailing_var_arg = true)]
    pub queries: Vec<String>,

    /// Dump index to file
    #[arg(short = 'd', long)]
    pub dump_index: Option<String>,

    /// Preset applied before other options. Accepted values:
    ///   lr / map-ont     — Nanopore vs reference (default-like)
    ///   map-pb / map10k  — PacBio CLR vs reference (HPC seeds)
    ///   map-hifi / map-ccs — PacBio HiFi vs reference
    ///   map-iclr         — Illumina Complete Long Reads vs reference
    ///   map-iclr-prerender — ICLR pre-render variant
    ///   lr:hq            — accurate long reads (<1% error) vs reference
    ///   lr:hq-balanced   — lr:hq with near-output-preserving DP savings (--gap-bw-margin 500 --ungapped-fill)
    ///   lr:hq-fast       — lr:hq with cheaper DP and chaining (--ext-bw 200 --ext-adaptive 64 --ext-rescue --ext-probe 2000,1000 --gap-bw-margin 500 --gap-bw-narrow 32 --ungapped-fill --max-gap-ref 5000 --seg-rmq --clean-dp --rescue-join)
    ///   lr:hq-turbo      — lr:hq-fast plus --max-chain-skip 10 --max-gap-ref 2000 --rescue-subset --ext-probe 1000,1000 --gap-bw-margin 300 --gap-band-asym --stable-anchor-sort (faster; loses some SV support)
    ///   lr:hqae          — HQ long-read assembly eval (k=25 w=51, RMQ)
    ///   asm5 / asm10 / asm20 — asm-to-ref ~0.1/1/5% divergence
    ///   ava-ont / ava-pb — all-vs-all Nanopore / PacBio read overlap
    ///   short / sr       — short reads vs reference
    ///   splice / splice:hq — spliced long reads / accurate long reads
    ///   splice:sr        — spliced short RNA-seq reads
    ///   cdna             — cDNA / splice alias
    #[arg(short = 'x', long, verbatim_doc_comment)]
    pub preset: Option<String>,

    /// Minimizer k-mer length
    #[arg(short = 'k')]
    pub k: Option<usize>,

    /// Minimizer window size
    #[arg(short = 'w')]
    pub w: Option<usize>,

    /// Use homopolymer-compressed k-mer
    #[arg(short = 'H')]
    pub is_hpc: bool,

    /// Bandwidth used for chaining and long-join. Format: bw,bw_long
    #[arg(short = 'r', long)]
    pub bw: Option<String>,

    /// Minimal chaining score
    #[arg(short = 'm', long)]
    pub min_chain_score: Option<i32>,

    /// Minimal number of minimizers on a chain
    #[arg(short = 'n', long)]
    pub min_cnt: Option<i32>,

    /// Stop chain elongation if there are no minimizers in INT-bp
    #[arg(short = 'g', long)]
    pub max_gap: Option<i32>,

    /// Chain gap penalty scale
    #[arg(long)]
    pub chain_gap_scale: Option<f32>,

    /// Chain skip penalty scale
    #[arg(long)]
    pub chain_skip_scale: Option<f32>,

    /// Filter out top FLOAT fraction of repetitive minimizers (>=1 sets integer count)
    #[arg(short = 'f', long)]
    pub mid_occ_frac: Option<f32>,

    /// Ignore seed hits with more than INT occurrences
    #[arg(short = 'M', long)]
    pub mid_occ: Option<usize>,

    /// Output CIGAR in PAF
    #[arg(short = 'c', long)]
    pub output_cigar: bool,

    /// Z-drop score and inversion z-drop (format: zdrop[,zdrop_inv])
    #[arg(short = 'z', long = "z-drop")]
    pub z_drop: Option<String>,

    /// Output CS tag (short/long/none)
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "short")]
    pub cs: Option<Option<String>>,

    /// Output SAM format
    #[arg(short = 'a', long)]
    pub output_sam: bool,

    /// Matching score
    #[arg(short = 'A')]
    pub scor_a: Option<i32>,

    /// Mismatch penalty
    #[arg(short = 'B')]
    pub scor_b: Option<i32>,

    /// Gap open penalty (format: O1,O2)
    #[arg(short = 'O')]
    pub scor_o: Option<String>,

    /// Gap extension penalty (format: E1,E2)
    #[arg(short = 'E')]
    pub scor_e: Option<String>,

    /// Minimal peak DP alignment score
    #[arg(short = 's')]
    pub min_dp_max: Option<i32>,

    /// Max number of secondary alignments
    #[arg(short = 'N')]
    pub best_n: Option<i32>,

    /// Min secondary-to-primary score ratio
    #[arg(short = 'p')]
    pub pri_ratio: Option<f32>,

    /// SAM read group line in a format like '@RG\tID:foo\tSM:bar'
    #[arg(short = 'R')]
    pub sam_rg: Option<String>,

    /// Copy input FASTA/Q comments to output
    #[arg(short = 'y')]
    pub copy_comment: bool,

    /// Output =/X CIGAR operators
    #[arg(long)]
    pub eqx: bool,

    /// Output MD tag
    #[arg(long = "MD")]
    pub output_md: bool,

    /// Output ds tag (cs extension with INDEL position uncertainty)
    #[arg(long)]
    pub ds: bool,

    /// Number of threads for alignment
    #[arg(short = 't', long, default_value_t = 3)]
    pub align_threads: usize,

    /// Adaptive mid_occ floor,ceiling (format: min,max)
    #[arg(short = 'U')]
    pub mid_occ_range: Option<String>,

    /// Max fragment length for paired-end
    #[arg(short = 'F')]
    pub max_frag_len: Option<i32>,

    /// Max intron length (sets max_gap_ref + bw for splice)
    #[arg(short = 'G')]
    pub max_intron_len: Option<i32>,

    /// Only map to forward strand
    #[arg(long)]
    pub for_only: bool,

    /// Only map to reverse strand
    #[arg(long)]
    pub rev_only: bool,

    /// Disable long-gap joining
    #[arg(long)]
    pub no_long_join: bool,

    /// Output all chains (skip secondary filtering)
    #[arg(short = 'X')]
    pub all_chains: bool,

    /// Transition mismatch penalty
    #[arg(short = 'b')]
    pub transition: Option<i32>,

    /// Cost of non-canonical splicing
    #[arg(short = 'C')]
    pub noncan: Option<i32>,

    /// Splice model: 0=original, 1=miniprot [1]
    #[arg(short = 'J')]
    pub splice_model: Option<i32>,

    /// Splice site direction: f=forward, r=reverse, b=both, n=none
    #[arg(short = 'u')]
    pub splice_direction: Option<String>,

    /// Output to file instead of stdout
    #[arg(short = 'o')]
    pub output: Option<String>,

    /// Write CIGAR with >65535 ops at CG tag
    #[arg(short = 'L')]
    pub long_cigar: bool,

    /// Use soft clipping for supplementary alignments
    #[arg(short = 'Y')]
    pub softclip: bool,

    /// Don't output base quality in SAM
    #[arg(short = 'Q')]
    pub no_qual: bool,

    /// No exact diagonal for all-vs-all
    #[arg(short = 'D')]
    pub no_diag: bool,

    /// Seed occurrence distance
    #[arg(short = 'e')]
    pub occ_dist: Option<i32>,

    /// Minibatch size for mapping (query bases per batch; suffix K/M/G OK) [500M / preset dependent]
    #[arg(short = 'K', long = "mb-size")]
    pub mini_batch_size: Option<String>,

    /// Whether to output secondary alignments (yes/no)
    #[arg(long)]
    pub secondary: Option<String>,

    /// Max number of anchors to skip during chaining
    #[arg(long)]
    pub max_chain_skip: Option<i32>,

    /// Max chaining iterations
    #[arg(long)]
    pub max_chain_iter: Option<i32>,

    /// DP bandwidth for read-end extensions [default: 1.5*bandwidth+1]
    #[arg(long)]
    pub ext_bw: Option<i32>,

    /// Gap-fill bandwidth = min(default, |qlen-tlen| + INT), keeping large indels crossable
    #[arg(long)]
    pub gap_bw_margin: Option<i32>,

    /// Try an ungapped fill first for equal-length gaps (long reads; short reads always do)
    #[arg(long)]
    pub ungapped_fill: bool,

    /// With --ext-bw: re-run a z-dropped extension with the default band when its unaligned tail
    /// shows an indel wider than the narrow band
    #[arg(long)]
    pub ext_rescue: bool,

    /// Extension probe LEN,MIN_SCORE_PER_KB[,LOW_BW]: align the first LEN bases of a long
    /// extension window first and stop there if the alignment z-drops or scores below
    /// MIN_SCORE_PER_KB (with LOW_BW: align a low-scoring window in full with band LOW_BW)
    #[arg(long)]
    pub ext_probe: Option<String>,

    /// Max reference gap between chained anchors in the main (DP) chaining pass [-g]; the
    /// extension window and the long-join rescue keep -g
    #[arg(long)]
    pub max_gap_ref: Option<i32>,

    /// Gap fill: try band |qlen-tlen| + INT first unless read k-mers match the target only off
    /// that band; fall back to the regular band when the path nears its edge
    #[arg(long)]
    pub gap_bw_narrow: Option<i32>,

    /// Lower bound of the --gap-bw-margin band (default: the -r bandwidth x 1.5)
    #[arg(long)]
    pub gap_bw_floor: Option<i32>,

    /// Extensions: adaptive band of INT cells per anti-diagonal that follows the alignment
    /// (0 = static --ext-bw band)
    #[arg(long)]
    pub ext_adaptive: Option<i32>,

    /// Dual-affine DP with exact band boundaries (the scalar reference's semantics) instead of
    /// minimap2's SSE kernel behaviour; faster, output not minimap2-identical
    #[arg(long)]
    pub clean_dp: bool,

    /// RMQ chaining with a segment tree (faster; ties between equal-priority predecessors may
    /// resolve differently from minimap2)
    #[arg(long)]
    pub seg_rmq: bool,

    /// Run the long-gap (RMQ) rescue re-chain only when two first-pass chains could be joined [on in lr:hq-fast]
    #[arg(long)]
    pub rescue_join: bool,

    /// Asymmetric gap-fill bands: cover only the diagonals between the gap's ends (plus off-band
    /// k-mer evidence) [on in lr:hq-turbo; faster, loses some SV support]
    #[arg(long)]
    pub gap_band_asym: bool,

    /// Stable LSD radix sort for anchors (ties keep seed order; not minimap2's order) [on in lr:hq-turbo]
    #[arg(long)]
    pub stable_anchor_sort: bool,

    /// With --rescue-join, re-chain only the joinable chains and keep the others [on in lr:hq-turbo]
    #[arg(long)]
    pub rescue_subset: bool,

    /// Turn off --rescue-join (e.g. to undo it in lr:hq-fast)
    #[arg(long, conflicts_with = "rescue_join")]
    pub no_rescue_join: bool,

    /// Score for ambiguous bases
    #[arg(long = "score-N")]
    pub score_n: Option<i32>,

    /// Bonus for reaching end of query
    #[arg(long)]
    pub end_bonus: Option<i32>,

    /// End seed penalty (anchor_ext_shift)
    #[arg(long)]
    pub end_seed_pen: Option<i32>,

    /// Drop an alignment if BOTH ends are clipped above this ratio of qlen [1.0]
    #[arg(long)]
    pub max_clip_ratio: Option<f32>,

    /// Mask level for repeat filtering
    #[arg(long)]
    pub mask_level: Option<f32>,

    /// Mask length
    #[arg(long)]
    pub mask_len: Option<i32>,

    /// Hard masking level
    #[arg(long)]
    pub hard_mask_level: bool,

    /// Disable end filtering
    #[arg(long)]
    pub no_end_flt: bool,

    /// Cap Smith-Waterman matrix size
    #[arg(long)]
    pub cap_sw_mem: Option<i64>,

    /// Min DP alignment length
    #[arg(long)]
    pub min_dp_len: Option<i32>,

    /// Random seed
    #[arg(long)]
    pub seed: Option<i32>,

    /// Use RMQ chaining (yes/no)
    #[arg(long)]
    pub rmq: Option<String>,

    /// RMQ inner distance
    #[arg(long)]
    pub rmq_inner: Option<i32>,

    /// Dual mapping (yes/no)
    #[arg(long)]
    pub dual: Option<String>,

    /// Query occurrence fraction
    #[arg(long)]
    pub q_occ_frac: Option<f32>,

    /// Splice junction bonus
    #[arg(long)]
    pub junc_bonus: Option<i32>,

    /// Splice junction penalty
    #[arg(long)]
    pub junc_pen: Option<i32>,

    /// BED file with splice junction annotations
    #[arg(long)]
    pub junc_bed: Option<String>,

    /// Splice score file
    #[arg(long)]
    pub spsc: Option<String>,

    /// Scale factor for SPSC scores (default: 0.7)
    #[arg(long)]
    pub spsc_scale: Option<f32>,

    /// Output unmapped reads in PAF
    #[arg(long)]
    pub paf_no_hit: bool,

    /// Only output hits in SAM
    #[arg(long)]
    pub sam_hit_only: bool,

    /// Output SEQ for secondary alignments
    #[arg(long)]
    pub secondary_seq: bool,

    /// Number of bases loaded into memory for indexing [8G]
    #[arg(short = 'I')]
    pub batch_size: Option<String>,

    /// Write temporary files with this prefix for multi-part index merging
    #[arg(long)]
    pub split_prefix: Option<String>,

    /// File listing ALT contig names (one per line)
    #[arg(long)]
    pub alt: Option<String>,

    /// Score drop for ALT contigs [0.15]
    #[arg(long)]
    pub alt_drop: Option<f32>,

    /// Maximum query length (skip reads longer than this)
    #[arg(long)]
    pub max_qlen: Option<i32>,

    /// Don't store target sequences in index
    #[arg(long)]
    pub idx_no_seq: bool,

    /// Write the -d index in the raw format, which loads by memory-mapping (near-instant).
    /// With a pre-built index as the target, -d converts it.
    #[arg(long)]
    pub idx_raw: bool,

    /// Pairing mode: no, weak, or strong [determined by preset]
    #[arg(long)]
    pub pairing: Option<String>,

    /// Don't use query name in read hash
    #[arg(long)]
    pub no_hash_name: bool,

    /// Write junction BED instead of SAM/PAF
    #[arg(long)]
    pub write_junc: bool,

    /// Query-strand alignment mode (for SV detection)
    #[arg(long)]
    pub qstrand: bool,

    /// Parallel chaining across reference sequences (asm2asm workloads)
    #[arg(long)]
    pub par_chain: bool,

    /// Enable splice alignment mode
    #[arg(long)]
    pub splice: bool,

    /// Short-read mode: dna, rna, or no
    #[arg(long, num_args = 0..=1, default_missing_value = "dna")]
    pub sr: Option<String>,

    /// Fragment mode (yes/no)
    #[arg(long)]
    pub frag: Option<String>,

    /// Splice flank scoring (yes/no)
    #[arg(long)]
    pub splice_flank: Option<String>,

    /// Heap-based seed sorting (yes/no)
    #[arg(long)]
    pub heap_sort: Option<String>,

    /// BED12 junction file for jump extension
    #[arg(short = 'j', long = "junc-jump")]
    pub junc_jump: Option<String>,

    /// BED12 pass1 rescue junctions
    #[arg(long = "pass1")]
    pub pass1: Option<String>,

    /// Minimum matching bases for jump extension [3]
    #[arg(long)]
    pub jump_min_match: Option<i32>,
}

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();

    env_logger::Builder::new()
        .filter_level(if args.verbose {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Warn
        })
        .init();

    run(args.align)
}

// ─── Alignment logic ───

/// Read ALT contig names from file and mark matching sequences in index.
fn read_alt_list(mi: &mut Index, alt_path: &str) -> anyhow::Result<usize> {
    use std::io::BufRead;
    let file = std::fs::File::open(alt_path)
        .map_err(|e| anyhow::anyhow!("Cannot open alt file '{}': {}", alt_path, e))?;
    let reader = std::io::BufReader::new(file);
    // Collect names into owned set to avoid borrow conflict
    let mut alt_names: std::collections::HashSet<String> = std::collections::HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let name = line.split_whitespace().next().unwrap_or("").to_string();
        if !name.is_empty() { alt_names.insert(name); }
    }
    let mut n_alt = 0;
    for seq in mi.seqs.iter_mut() {
        if alt_names.contains(&seq.name) {
            seq.is_alt = true;
            n_alt += 1;
        }
    }
    if n_alt > 0 {
        eprintln!("[*] Found {} ALT contigs", n_alt);
    }
    Ok(n_alt)
}

/// Parse a number with optional K/M/G suffix.
/// Uses decimal multipliers (1e3/1e6/1e9) and supports floating point (e.g. "4.5G").
fn parse_num(s: &str) -> i64 {
    let s = s.trim();
    // Try to parse as float, then check suffix
    let (num_str, mult) = if s.ends_with('G') || s.ends_with('g') {
        (&s[..s.len() - 1], 1e9)
    } else if s.ends_with('M') || s.ends_with('m') {
        (&s[..s.len() - 1], 1e6)
    } else if s.ends_with('K') || s.ends_with('k') {
        (&s[..s.len() - 1], 1e3)
    } else {
        (s, 1.0)
    };
    let x: f64 = num_str.parse().unwrap_or(0.0);
    (x * mult + 0.499) as i64
}

/// Apply preset configurations.
/// The `is_hpc` output indicates whether homopolymer compression should be used.
fn apply_preset(opt: &mut MapOptions, k: &mut usize, w: &mut usize, is_hpc: &mut bool, preset: &str) -> anyhow::Result<()> {
    rammap::api::apply_preset_str(opt, k, w, is_hpc, preset).map_err(|e| anyhow::anyhow!("{}", e))
}

fn run(cli: AlignArgs) -> anyhow::Result<()> {
    let do_ds = cli.ds;
    let cs_mode: Option<&str> = cli.cs.as_ref().map(|v| v.as_deref().unwrap_or("short"));
    if let Some(m) = cs_mode && m != "short" && m != "long" && m != "none" {
        eprintln!("[WARNING] --cs only takes 'short', 'long', or 'none'. '{}' is treated as 'short'.", m);
    }
    let cs_is_none = cs_mode == Some("none");
    let cs_long = cs_mode == Some("long");
    let do_cigar = cli.output_cigar || cli.output_sam || cli.cs.is_some() || cli.output_md || do_ds || cli.write_junc;
    let do_cs = cli.cs.is_some() && !cs_is_none;
    let do_md = cli.output_md;
    let rg_id: Option<String> = if let Some(rg_line) = &cli.sam_rg {
        let expanded = rg_line.replace("\\t", "\t");
        expanded.split('\t').find(|token| token.starts_with("ID:")).map(|token| token[3..].to_string())
    } else {
        None
    };
    let out_cfg = OutputConfig {
        do_cigar,
        do_cs,
        cs_long,
        do_md,
        do_ds,
        eqx: cli.eqx,
        output_sam: cli.output_sam,
        rg_id: rg_id.clone(),
        split_mode: false,
    };
    let split_out_cfg = OutputConfig {
        split_mode: true,
        ..out_cfg.clone()
    };
    let copy_comment = cli.copy_comment;
    rammap::fasta::bam::set_copy_tags(copy_comment);

    #[cfg(not(feature = "parallel"))]
    if cli.align_threads > 1 {
        eprintln!("Warning: -t/--threads specified but executable was compiled without 'parallel' feature. Running in single-threaded mode.");
    }

    #[cfg(feature = "parallel")]
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.align_threads)
        .build_global()
        .map_err(|e| anyhow::anyhow!("Failed to initialize thread pool: {}", e))?;

    let mut k = 15;
    let mut w = 10;
    let mut is_hpc = cli.is_hpc;
    let mut opt = MapOptions::default();

    if let Some(preset) = &cli.preset {
        apply_preset(&mut opt, &mut k, &mut w, &mut is_hpc, preset)?;
    }

    // CLI overrides (applied after preset)
    if let Some(vk) = cli.k { k = vk; }
    if let Some(vw) = cli.w { w = vw; }

    if let Some(r_str) = &cli.bw {
        let parts: Vec<&str> = r_str.split(',').collect();
        if !parts.is_empty() && let Ok(val) = parts[0].parse() { opt.chaining.bandwidth = val; }
        if parts.len() > 1 && let Ok(val) = parts[1].parse() { opt.chaining.bandwidth_long = val; }
    }

    if let Some(v) = cli.min_chain_score { opt.chaining.min_chain_score = v; }
    if let Some(v) = cli.min_cnt { opt.chaining.min_cnt = v; }
    if let Some(v) = cli.max_gap { opt.chaining.max_gap = v; }
    if let Some(v) = cli.chain_gap_scale {
        opt.chaining.chain_gap_scale = v;
    }
    if let Some(v) = cli.mid_occ { opt.seeding.mid_occ = v; }
    if let Some(v) = cli.scor_a { opt.scoring.match_score = v; }
    if let Some(v) = cli.scor_b { opt.scoring.mismatch_penalty = v; }

    if let Some(s) = &cli.scor_o {
        let parts: Vec<&str> = s.split(',').collect();
        if let Ok(v) = parts[0].parse() { opt.scoring.gap_open = v; }
        if parts.len() > 1 && let Ok(v) = parts[1].parse() { opt.scoring.gap_open2 = v; }
    }
    if let Some(s) = &cli.scor_e {
        let parts: Vec<&str> = s.split(',').collect();
        if let Ok(v) = parts[0].parse() { opt.scoring.gap_extend = v; }
        if parts.len() > 1 && let Ok(v) = parts[1].parse() { opt.scoring.gap_extend2 = v; }
    }
    if let Some(v) = cli.best_n { opt.filtering.best_n = v; }
    if let Some(v) = cli.pri_ratio { opt.filtering.pri_ratio = v; }
    if let Some(s) = &cli.z_drop {
        let parts: Vec<&str> = s.split(',').collect();
        if let Ok(v) = parts[0].parse::<i32>() {
            opt.alignment.zdrop = v;
            opt.alignment.zdrop_inv = v; // zdrop_inv = zdrop by default
        }
        if parts.len() > 1 && let Ok(v) = parts[1].parse() { opt.alignment.zdrop_inv = v; }
    }
    if let Some(s) = &cli.mid_occ_range {
        let parts: Vec<&str> = s.split(',').collect();
        if let Ok(v) = parts[0].parse() { opt.seeding.min_mid_occ = v; }
        if parts.len() > 1 && let Ok(v) = parts[1].parse() { opt.seeding.max_mid_occ = v; }
    }
    if let Some(v) = cli.max_frag_len { opt.pairing.max_frag_len = v; }
    if let Some(v) = cli.max_intron_len {
        opt.chaining.max_gap_ref = v;
        opt.chaining.bandwidth = v; opt.chaining.bandwidth_long = v;
    }
    if cli.for_only { opt.flags.insert(AlignFlags::FOR_ONLY); }
    if cli.rev_only { opt.flags.insert(AlignFlags::REV_ONLY); }
    if cli.no_long_join { opt.flags.insert(AlignFlags::NO_LJOIN); }
    if cli.all_chains { opt.flags.insert(AlignFlags::ALL_CHAINS); }
    if let Some(v) = cli.min_dp_max { opt.alignment.min_dp_max = v; }
    if let Some(v) = cli.transition { opt.scoring.transition = v; }
    if let Some(v) = cli.noncan { opt.scoring.noncanon_penalty = v; }
    if let Some(j) = cli.splice_model {
        if j == 0 { opt.flags.insert(AlignFlags::SPLICE_OLD); }
        else { opt.flags.remove(AlignFlags::SPLICE_OLD); }
    }
    if let Some(v) = cli.score_n { opt.scoring.ambig_penalty = v; }
    if let Some(v) = cli.end_bonus { opt.alignment.end_bonus = v; }
    if let Some(v) = cli.end_seed_pen { opt.alignment.anchor_ext_shift = v; }
    if let Some(v) = cli.max_chain_skip { opt.chaining.max_chain_skip = v; }
    if let Some(v) = cli.max_chain_iter { opt.chaining.max_chain_iter = v; }
    if let Some(v) = cli.ext_bw { opt.alignment.ext_bw = Some(v); }
    if let Some(v) = cli.gap_bw_margin { opt.alignment.gap_bw_margin = Some(v); }
    if cli.ungapped_fill { opt.alignment.ungapped_fill = true; }
    if cli.ext_rescue { opt.alignment.ext_rescue = true; }
    if cli.seg_rmq { opt.chaining.seg_rmq = true; }
    if cli.rescue_join { opt.chaining.rescue_join = true; }
    if cli.no_rescue_join { opt.chaining.rescue_join = false; }
    if cli.gap_band_asym { opt.alignment.gap_band_asym = true; }
    if cli.stable_anchor_sort { opt.seeding.stable_anchor_sort = true; }
    if cli.rescue_subset { opt.chaining.rescue_subset = true; }
    if cli.clean_dp { opt.alignment.clean_dp = true; }
    if let Some(v) = cli.ext_adaptive { opt.alignment.ext_adaptive = if v > 0 { Some(v) } else { None }; }
    if let Some(v) = cli.gap_bw_floor { opt.alignment.gap_bw_floor = Some(v); }
    if let Some(v) = cli.gap_bw_narrow { opt.alignment.gap_bw_narrow = if v > 0 { Some(v) } else { None }; }
    if let Some(v) = cli.max_gap_ref { opt.chaining.max_gap_ref = v; }
    if let Some(s) = &cli.ext_probe {
        let v: Vec<i32> = s.split(',').filter_map(|x| x.trim().parse().ok()).collect();
        match v.as_slice() {
            [0] => opt.alignment.ext_probe = None,
            [len, min] if *len > 0 => opt.alignment.ext_probe = Some((*len, *min, 0)),
            [len, min, low_bw] if *len > 0 => opt.alignment.ext_probe = Some((*len, *min, *low_bw)),
            _ => eprintln!("Warning: --ext-probe expects LEN,MIN_SCORE_PER_KB[,LOW_BW]"),
        }
    }
    if let Some(v) = cli.max_clip_ratio { opt.alignment.max_clip_ratio = v; }
    if let Some(v) = cli.mask_level { opt.filtering.mask_level = v; }
    if let Some(v) = cli.mask_len { opt.filtering.mask_len = v; }
    if let Some(v) = cli.occ_dist { opt.seeding.occ_dist = v; }
    if let Some(v) = cli.q_occ_frac { opt.seeding.q_occ_frac = v; }
    if let Some(v) = cli.cap_sw_mem { opt.alignment.max_sw_mat = v; }
    if let Some(v) = cli.min_dp_len { opt.alignment.min_dp_len = v; }
    if let Some(v) = cli.seed { opt.filtering.seed = v; }
    if let Some(v) = cli.alt_drop { opt.filtering.alt_drop = v; }
    if let Some(v) = cli.max_qlen { opt.filtering.max_qlen = v; }
    if let Some(s) = &cli.pairing {
        match s.as_str() {
            "no" => { opt.flags.insert(AlignFlags::INDEPEND_SEG); }
            "weak" => { opt.flags.insert(AlignFlags::WEAK_PAIRING); opt.flags.remove(AlignFlags::INDEPEND_SEG); }
            "strong" => { opt.flags.remove(AlignFlags::INDEPEND_SEG | AlignFlags::WEAK_PAIRING); }
            _ => { eprintln!("[WARNING] unrecognized argument for --pairing; assuming 'strong'."); opt.flags.remove(AlignFlags::INDEPEND_SEG | AlignFlags::WEAK_PAIRING); }
        }
    }
    if cli.no_hash_name { opt.flags.insert(AlignFlags::NO_HASH_NAME); }
    if cli.write_junc { opt.flags.insert(AlignFlags::OUT_JUNC); }
    if cli.qstrand { opt.flags.insert(AlignFlags::QSTRAND | AlignFlags::NO_INV); }
    if cli.par_chain { opt.flags.insert(AlignFlags::PAR_CHAIN); }
    if cli.splice { opt.flags.insert(AlignFlags::SPLICE); }
    if let Some(ref val) = cli.sr {
        match val.as_str() {
            "dna" => { opt.flags.insert(AlignFlags::SHORT_READ); }
            "rna" => { opt.flags.insert(AlignFlags::SR_RNA); }
            "no" => { opt.flags.remove(AlignFlags::SHORT_READ | AlignFlags::SR_RNA); }
            v => { opt.flags.insert(AlignFlags::SHORT_READ); eprintln!("[WARNING] --sr only takes 'dna' or 'rna'. Invalid value '{}' assumed to be 'dna'.", v); }
        }
    }
    if let Some(ref s) = cli.frag {
        match s.as_str() {
            "yes" => { opt.flags.insert(AlignFlags::FRAG_MODE); }
            "no" => { opt.flags.remove(AlignFlags::FRAG_MODE); }
            _ => { eprintln!("[WARNING] --frag only takes 'yes' or 'no'."); }
        }
    }
    if let Some(ref s) = cli.splice_flank {
        match s.as_str() {
            "yes" => { opt.flags.insert(AlignFlags::SPLICE_FLANK); }
            "no" => { opt.flags.remove(AlignFlags::SPLICE_FLANK); }
            _ => { eprintln!("[WARNING] --splice-flank only takes 'yes' or 'no'."); }
        }
    }
    if let Some(ref s) = cli.heap_sort {
        match s.as_str() {
            "yes" => { opt.flags.insert(AlignFlags::HEAP_SORT); }
            "no" => { opt.flags.remove(AlignFlags::HEAP_SORT); }
            _ => { eprintln!("[WARNING] --heap-sort only takes 'yes' or 'no'."); }
        }
    }
    if let Some(v) = cli.jump_min_match { opt.filtering.jump_min_match = v; }
    if let Some(v) = cli.junc_bonus { opt.scoring.junc_bonus = v; }
    if let Some(v) = cli.junc_pen { opt.scoring.junc_pen = v; }
    if let Some(v) = cli.rmq_inner { opt.chaining.rmq_inner_dist = v; }
    if cli.no_diag { opt.flags.insert(AlignFlags::NO_DIAG); }
    if cli.hard_mask_level { opt.flags.insert(AlignFlags::HARD_MASK_LEVEL); }
    if cli.no_end_flt { opt.flags.insert(AlignFlags::NO_END_FLT); }
    if cli.softclip { opt.flags.insert(AlignFlags::SOFTCLIP); }
    if cli.no_qual { opt.flags.insert(AlignFlags::NO_QUAL); }
    if cli.long_cigar { opt.flags.insert(AlignFlags::LONG_CIGAR); }
    if cli.paf_no_hit { opt.flags.insert(AlignFlags::PAF_NO_HIT); }
    if cli.sam_hit_only { opt.flags.insert(AlignFlags::SAM_HIT_ONLY); }
    if cli.secondary_seq { opt.flags.insert(AlignFlags::SECONDARY_SEQ); }
    if cli.output_cigar || cli.cs.is_some() { opt.flags.insert(AlignFlags::OUT_CIGAR); }
    if let Some(s) = &cli.secondary {
        match s.as_str() {
            "yes" => { opt.flags.remove(AlignFlags::NO_PRINT_2ND); }
            "no" => { opt.flags.insert(AlignFlags::NO_PRINT_2ND); }
            _ => { eprintln!("Warning: --secondary must be 'yes' or 'no'"); }
        }
    }
    if let Some(s) = &cli.rmq {
        match s.as_str() {
            "yes" => { opt.flags.insert(AlignFlags::RMQ_CHAIN); }
            "no" => { opt.flags.remove(AlignFlags::RMQ_CHAIN); }
            _ => { eprintln!("Warning: --rmq must be 'yes' or 'no'"); }
        }
    }
    if let Some(s) = &cli.dual {
        match s.as_str() {
            "yes" => { opt.flags.remove(AlignFlags::NO_DUAL); }
            "no" => { opt.flags.insert(AlignFlags::NO_DUAL); }
            _ => { eprintln!("Warning: --dual must be 'yes' or 'no'"); }
        }
    }
    if let Some(s) = &cli.splice_direction {
        // Clear existing splice direction flags, then set based on user choice
        opt.flags.remove(AlignFlags::SPLICE_FOR | AlignFlags::SPLICE_REV);
        match s.as_str() {
            "f" => { opt.flags.insert(AlignFlags::SPLICE_FOR); }
            "r" => { opt.flags.insert(AlignFlags::SPLICE_REV); }
            "b" => { opt.flags.insert(AlignFlags::SPLICE_FOR | AlignFlags::SPLICE_REV); }
            "n" => { /* both cleared above */ }
            _ => { eprintln!("Warning: -u must be f, r, b, or n"); }
        }
    }

    if let Some(v) = cli.chain_skip_scale {
        opt.filtering.chain_skip_scale = v;
    }
    rammap::align::dp::set_clean_dual(opt.alignment.clean_dp);
    // Recompute chaining penalties after all CLI overrides
    rammap::api::finalize_options(&mut opt, k);

    let target_path = &cli.target;
    let batch_size: u64 = cli.batch_size.as_ref()
        .map(|s| parse_num(s) as u64)
        .unwrap_or(8_000_000_000); // 8G default

    // -K: overrides opt.mini_batch_size set by preset (defaults: 500M, 50M for sr/short, 100M for splice:sr)
    if let Some(s) = cli.mini_batch_size.as_ref() {
        opt.mini_batch_size = parse_num(s);
    }

    // Parse query files from positional args
    let frag_mode = opt.flags.contains(AlignFlags::FRAG_MODE);
    let no_query = cli.queries.is_empty() && cli.dump_index.is_some();
    let (query_files, query2_path, two_file_pe): (Vec<String>, Option<String>, bool) = {
        if cli.queries.is_empty() {
            if cli.dump_index.is_some() {
                // -d with no query files → index-only mode (build/dump and exit)
                (vec![], None, false)
            } else {
                // No query files, no -d → read from stdin.
                (vec!["-".to_string()], None, false)
            }
        } else if cli.queries.len() == 2 && frag_mode {
            // Exactly 2 files with frag_mode (sr preset) → two-file paired-end
            (vec![cli.queries[0].clone()], Some(cli.queries[1].clone()), true)
        } else {
            // Single or multiple independent query files
            (cli.queries.clone(), None, false)
        }
    };
    let pe_mode = two_file_pe || frag_mode;

    // Validate query file existence
    for qf in &query_files {
        if qf != "-" && !std::path::Path::new(qf).exists() {
            anyhow::bail!("Input file '{}' not found.", qf);
        }
    }
    if let Some(ref q2) = query2_path && !std::path::Path::new(q2).exists() {
        anyhow::bail!("Second query file '{}' not found.", q2);
    }
    let query_path = query_files.first().map(|s| s.as_str()).unwrap_or("-");

    // ─── Output setup ───
    let t_start = Instant::now();
    let output_file: Option<std::fs::File> = if let Some(path) = &cli.output {
        Some(std::fs::File::create(path)
            .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", path, e))?)
    } else {
        None
    };
    let stdout = std::io::stdout();
    let mut handle: Box<dyn Write + Send> = if let Some(f) = output_file {
        Box::new(std::io::BufWriter::new(f))
    } else {
        Box::new(std::io::BufWriter::new(stdout))
    };

    let mut total_stats = AlignmentStats::default();
    let mid_occ_frac = cli.mid_occ_frac.unwrap_or(2e-4);
    let run_cfg = RunConfig {
        mid_occ_frac,
        query2_path: query2_path.as_deref(),
        pe_mode,
        two_file_pe,
        copy_comment,
        cli: &cli,
        out_cfg: &out_cfg,
        split_out_cfg: &split_out_cfg,
    };
    // ─── Build index parts iterator ───
    // Each iteration yields one Index part. For small references or large batch_size,
    // this will be a single iteration (matching current behavior exactly).
    // Extension is just a naming convention and not always reliable (e.g. a real
    // minimap2/rammap index saved under a non-standard name) — fall back to sniffing
    // the file's magic bytes so a mis-named index isn't silently misread as sequence
    // input to build a new (bogus) index from.
    let is_idx_file = target_path.ends_with(".mmi") || target_path.ends_with(".idx") || target_path.ends_with(".rmmi")
        || Index::sniff_is_index_file(target_path);
    let mut n_parts = 0usize;
    let mut sam_header_written = false;
    let split_prefix = cli.split_prefix.clone();
    let mut saved_k = k;
    let mut saved_w = w;
    let mut saved_is_hpc = is_hpc;

    if is_idx_file {
        // ─── Load from .mmi/.idx: iterate over parts ───
        eprintln!("[*] Loading index from {}...", target_path);
        let t_load = Instant::now();
        // Raw indexes are memory-mapped (all parts at once; pages fault in on use);
        // other formats are deserialized one part at a time.
        #[cfg(not(target_arch = "wasm32"))]
        let mut mapped_parts: Option<std::collections::VecDeque<Index>> = if Index::sniff_is_raw_index(target_path) {
            Some(Index::load_raw_parts(target_path)
                .map_err(|e| anyhow::anyhow!("Error loading index: {}", e))?.into())
        } else {
            None
        };
        #[cfg(target_arch = "wasm32")]
        let mut mapped_parts: Option<std::collections::VecDeque<Index>> = None;
        let f = std::fs::File::open(target_path)
            .map_err(|e| anyhow::anyhow!("Error opening index '{}': {}", target_path, e))?;
        let mut idx_reader = std::io::BufReader::new(f);
        // -d with a pre-built index: convert it (e.g. to the raw format with --idx-raw).
        let mut dump_writer: Option<std::io::BufWriter<std::fs::File>> = if let Some(dump_path) = &cli.dump_index {
            Some(std::io::BufWriter::new(std::fs::File::create(dump_path)
                .map_err(|e| anyhow::anyhow!("Failed to create index '{}': {}", dump_path, e))?))
        } else {
            None
        };

        loop {
            let next = match mapped_parts.as_mut() {
                Some(parts) => parts.pop_front(),
                None => Index::load_part(&mut idx_reader)
                    .map_err(|e| anyhow::anyhow!("Error loading index part: {}", e))?,
            };
            let mut mi = match next {
                Some(idx) => idx,
                None => break,
            };
            n_parts += 1;
            if let Some(dw) = dump_writer.as_mut() {
                if cli.idx_no_seq {
                    let mut idx_copy = mi.clone();
                    idx_copy.strip_sequences();
                    idx_copy.save_part_as(dw, cli.idx_raw)
                } else {
                    mi.save_part_as(dw, cli.idx_raw)
                }.map_err(|e| anyhow::anyhow!("Failed to save index part: {}", e))?;
            }
            let total_len: u64 = mi.seqs.iter().map(|s| s.len as u64).sum();
            eprintln!("[*] Loaded part {} in {:.3}s — k: {}; w: {}; hpc: {}; #seq: {}; total len: {}",
                n_parts, t_load.elapsed().as_secs_f64(),
                mi.kmer_size, mi.window_size, if mi.homopolymer_compressed { 1 } else { 0 },
                mi.seqs.len(), total_len);
            if n_parts == 1 && mapped_parts.is_none() && cli.dump_index.is_none() {
                eprintln!("[*] Tip: `rammap -d INDEX.raw --idx-raw {}` converts this index to the raw format, which loads near-instantly", target_path);
            }
            if n_parts == 1 && do_cigar && !mi.has_sequences() {
                anyhow::bail!("The prebuilt index doesn't contain sequences. Cannot produce CIGAR/cs/MD output.");
            }
            saved_k = mi.kmer_size;
            saved_w = mi.window_size;
            saved_is_hpc = mi.homopolymer_compressed;
            if let Some(ref alt_path) = cli.alt {
                read_alt_list(&mut mi, alt_path)?;
            }

            if let Some(ref prefix) = split_prefix {
                // Split mode: write results to temp files
                if query_path == "-" && n_parts > 1 {
                    anyhow::bail!("Cannot use stdin for queries with multi-part index. Use a named file or increase -I.");
                }
                let mut split_writer = split::split_init(prefix, n_parts - 1, &mi)
                    .map_err(|e| anyhow::anyhow!("Failed to init split temp file: {}", e))?;
                for qf in &query_files {
                    let part_stats = map_one_part_split(
                        &mi, &mut opt, &run_cfg,
                        qf, n_parts, &mut split_writer,
                    )?;
                    total_stats = total_stats + part_stats;
                }
            } else {
                for qf in &query_files {
                    let part_stats = map_one_part(
                        &mi, &mut opt, &run_cfg,
                        qf, n_parts, &mut sam_header_written,
                        &mut handle,
                    )?;
                    total_stats = total_stats + part_stats;
                }
            }
        }

        if n_parts == 0 {
            anyhow::bail!("Empty index file: {}", target_path);
        }
    } else {
        // ─── Build from FASTA/FASTQ: batch reading ───
        eprintln!("[*] Reading target: {}", target_path);

        let is_fastq = target_path.ends_with(".fq") || target_path.ends_with(".fastq")
            || target_path.ends_with(".fq.gz") || target_path.ends_with(".fastq.gz");

        // For batch reading, we need a streaming reader (can't use mmap for multi-part)
        // But for single-part (which is the common case), mmap is faster.
        // Strategy: read first batch; if EOF, we have single-part. Otherwise, switch to streaming.
        #[cfg(not(target_arch = "wasm32"))]
        let use_streaming = is_fastq || target_path.ends_with(".gz") || std::fs::metadata(target_path).map(|m| m.len()).unwrap_or(0) > batch_size;
        #[cfg(target_arch = "wasm32")]
        let use_streaming = true;

        if !use_streaming {
            // Non-gz FASTA: try mmap for single-part fast path
            #[cfg(not(target_arch = "wasm32"))]
            {
                let all_seqs = rammap::fasta::read_fasta(target_path)
                    .map_err(|e| anyhow::anyhow!("Error reading target: {}", e))?;

                // Check if we need multiple parts
                let total_bases: u64 = all_seqs.iter().map(|(_, s)| s.len() as u64).sum();

                if total_bases <= batch_size {
                    // ─── Single-part fast path (identical to previous behavior) ───
                    n_parts = 1;
                    let t0 = Instant::now();
                    eprintln!("[*] Building index for {} sequences...", all_seqs.len());
                    let mut mi = Index::build(all_seqs, w, k, is_hpc, usize::MAX);
                    mi.index = 0;
                    eprintln!("[*] Index built in {:.3}s", t0.elapsed().as_secs_f64());
                    if let Some(ref alt_path) = cli.alt {
                        read_alt_list(&mut mi, alt_path)?;
                    }

                    // mid_occ per-part
                    if opt.seeding.mid_occ == 0 {
                        opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
                        eprintln!("[*] Calculated mid_occ: {} (frac={})", opt.seeding.mid_occ, mid_occ_frac);
                    } else {
                        eprintln!("[*] Using mid_occ: {}", opt.seeding.mid_occ);
                    }

                    if let Some(dump_path) = &cli.dump_index {
                        eprintln!("[*] Saving index to {}...", dump_path);
                        if cli.idx_no_seq {
                            let mut idx_copy = mi.clone();
                            idx_copy.strip_sequences();
                            if let Err(e) = idx_copy.save_as(dump_path, cli.idx_raw) {
                                eprintln!("Warning: Failed to save index: {}", e);
                            }
                        } else if let Err(e) = mi.save_as(dump_path, cli.idx_raw) {
                            eprintln!("Warning: Failed to save index: {}", e);
                        }
                    }

                    for qf in &query_files {
                        let part_stats = map_one_part(
                            &mi, &mut opt, &run_cfg,
                            qf, 1, &mut sam_header_written,
                            &mut handle,
                        )?;
                        total_stats = total_stats + part_stats;
                    }
                } else {
                    // ─── Multi-part from mmap: split all_seqs into batches ───
                    let mut batch_seqs: Vec<(String, Vec<u8>)> = Vec::new();
                    let mut batch_bases: u64 = 0;
                    let mut dump_writer: Option<std::io::BufWriter<std::fs::File>> = if let Some(dump_path) = &cli.dump_index {
                        Some(std::io::BufWriter::new(std::fs::File::create(dump_path)
                            .map_err(|e| anyhow::anyhow!("Failed to create index '{}': {}", dump_path, e))?))
                    } else {
                        None
                    };

                    // Helper closure for processing one batch
                    let mut process_batch = |batch_seqs: &mut Vec<(String, Vec<u8>)>, batch_bases: &mut u64,
                        n_parts: &mut usize, opt: &mut MapOptions, dump_writer: &mut Option<std::io::BufWriter<std::fs::File>>,
                        sam_header_written: &mut bool, total_stats: &mut AlignmentStats| -> anyhow::Result<()>
                    {
                        let seqs = std::mem::take(batch_seqs);
                        *batch_bases = 0;
                        *n_parts += 1;
                        let t0 = Instant::now();
                        eprintln!("[*] Building index part {} for {} sequences...", *n_parts, seqs.len());
                        let mut mi = Index::build(seqs, w, k, is_hpc, usize::MAX);
                        mi.index = *n_parts - 1;
                        eprintln!("[*] Index part {} built in {:.3}s", *n_parts, t0.elapsed().as_secs_f64());
                        if let Some(ref alt_path) = cli.alt {
                            read_alt_list(&mut mi, alt_path)?;
                        }

                        if opt.seeding.mid_occ == 0 {
                            opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
                            eprintln!("[*] Calculated mid_occ: {} (frac={}) [part {}]", opt.seeding.mid_occ, mid_occ_frac, *n_parts);
                        }

                        if let Some(dw) = dump_writer.as_mut() {
                            if cli.idx_no_seq {
                                let mut idx_copy = mi.clone();
                                idx_copy.strip_sequences();
                                idx_copy.save_part_as(dw, cli.idx_raw).map_err(|e| anyhow::anyhow!("Failed to save index part: {}", e))?;
                            } else {
                                mi.save_part_as(dw, cli.idx_raw).map_err(|e| anyhow::anyhow!("Failed to save index part: {}", e))?;
                            }
                        }

                        if let Some(ref prefix) = split_prefix {
                            let mut split_writer = split::split_init(prefix, *n_parts - 1, &mi)
                                .map_err(|e| anyhow::anyhow!("Failed to init split temp file: {}", e))?;
                            for qf in &query_files {
                                let part_stats = map_one_part_split(
                                    &mi, opt, &run_cfg,
                                    qf, *n_parts, &mut split_writer,
                                )?;
                                *total_stats = *total_stats + part_stats;
                            }
                        } else {
                            for qf in &query_files {
                                let part_stats = map_one_part(
                                    &mi, opt, &run_cfg,
                                    qf, *n_parts, sam_header_written,
                                    &mut handle,
                                )?;
                                *total_stats = *total_stats + part_stats;
                            }
                        }
                        Ok(())
                    };

                    for (name, seq) in all_seqs {
                        batch_bases += seq.len() as u64;
                        batch_seqs.push((name, seq));
                        if batch_bases > batch_size {
                            process_batch(&mut batch_seqs, &mut batch_bases,
                                &mut n_parts, &mut opt, &mut dump_writer,
                                &mut sam_header_written, &mut total_stats)?;
                        }
                    }
                    // Final batch
                    if !batch_seqs.is_empty() {
                        process_batch(&mut batch_seqs, &mut batch_bases,
                            &mut n_parts, &mut opt, &mut dump_writer,
                            &mut sam_header_written, &mut total_stats)?;
                    }
                }
            }
        } else {
            // ─── Streaming FASTA/FASTQ: batch reading ───
            let mut reader = rammap::fasta::open(target_path)
                .map_err(|e| anyhow::anyhow!("Error reading target: {}", e))?;
            let mut dump_writer: Option<std::io::BufWriter<std::fs::File>> = if let Some(dump_path) = &cli.dump_index {
                Some(std::io::BufWriter::new(std::fs::File::create(dump_path)
                    .map_err(|e| anyhow::anyhow!("Failed to create index '{}': {}", dump_path, e))?))
            } else {
                None
            };

            loop {
                let (batch_seqs, is_eof) = reader.read_batch(batch_size)
                    .map_err(|e| anyhow::anyhow!("Error reading target batch: {}", e))?;
                if batch_seqs.is_empty() { break; }

                n_parts += 1;
                let t0 = Instant::now();
                eprintln!("[*] Building index part {} for {} sequences...", n_parts, batch_seqs.len());
                let mut mi = Index::build(batch_seqs, w, k, is_hpc, usize::MAX);
                mi.index = n_parts - 1;
                eprintln!("[*] Index part {} built in {:.3}s", n_parts, t0.elapsed().as_secs_f64());
                if let Some(ref alt_path) = cli.alt {
                    read_alt_list(&mut mi, alt_path)?;
                }

                // mid_occ per-part: only calculate once
                if opt.seeding.mid_occ == 0 {
                    opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
                    eprintln!("[*] Calculated mid_occ: {} (frac={}) [part {}]", opt.seeding.mid_occ, mid_occ_frac, n_parts);
                }

                if let Some(dw) = dump_writer.as_mut() {
                    if cli.idx_no_seq {
                        let mut idx_copy = mi.clone();
                        idx_copy.strip_sequences();
                        idx_copy.save_part_as(dw, cli.idx_raw).map_err(|e| anyhow::anyhow!("Failed to save index part: {}", e))?;
                    } else {
                        mi.save_part_as(dw, cli.idx_raw).map_err(|e| anyhow::anyhow!("Failed to save index part: {}", e))?;
                    }
                }

                if let Some(ref prefix) = split_prefix {
                    let mut split_writer = split::split_init(prefix, n_parts - 1, &mi)
                        .map_err(|e| anyhow::anyhow!("Failed to init split temp file: {}", e))?;
                    for qf in &query_files {
                        let part_stats = map_one_part_split(
                            &mi, &mut opt, &run_cfg,
                            qf, n_parts, &mut split_writer,
                        )?;
                        total_stats = total_stats + part_stats;
                    }
                } else {
                    for qf in &query_files {
                        let part_stats = map_one_part(
                            &mi, &mut opt, &run_cfg,
                            qf, n_parts, &mut sam_header_written,
                            &mut handle,
                        )?;
                        total_stats = total_stats + part_stats;
                    }
                }
                drop(mi);

                if is_eof { break; }
            }
        }
    }

    if n_parts > 1 && split_prefix.is_none() {
        eprintln!("[WARNING] Index has {} parts. Without --split-prefix, MAPQ and secondary flags may be inaccurate.", n_parts);
    }

    // Merge split results if --split-prefix was used
    if let Some(ref prefix) = split_prefix && n_parts >= 1 {
        eprintln!("[*] Merging {} split parts...", n_parts);
        let pe_flip = if pe_mode {
            let pe_ori = opt.pairing.pe_ori;
            Some(((pe_ori >> 1) & 1 != 0, pe_ori & 1 != 0))
        } else {
            None
        };
        let merge_cfg = split::MergeConfig {
            prefix,
            n_parts,
            k: saved_k,
            w: saved_w,
            is_hpc: saved_is_hpc,
            pe_mode,
            pe_flip,
            query_path,
            query2_path: query2_path.as_deref(),
            two_file_pe,
            sam_rg: cli.sam_rg.as_deref(),
            copy_comment,
        };
        let merge_stats = split::merge_split_results(
            &merge_cfg, &opt, &out_cfg, &mut *handle,
        )?;
        total_stats = total_stats + merge_stats;
    }

    if no_query {
        let total_time = t_start.elapsed().as_secs_f64();
        eprintln!("[*] Done in {:.3}s ({} index part{}). No query file provided; skipping mapping.",
            total_time, n_parts, if n_parts != 1 { "s" } else { "" });
    } else {
        #[cfg(feature = "census")]
        rammap::align::dp::census::flush();
        #[cfg(feature = "gpu")]
        if let Some(line) = rammap::align::dp::gpu::report() { eprintln!("[*] {line}"); }
        let total_time = t_start.elapsed().as_secs_f64();
        eprintln!("[*] Mapping done in {:.3}s ({} index part{})", total_time, n_parts, if n_parts != 1 { "s" } else { "" });

        eprintln!("--- Alignment Breakdown ---");
        eprintln!("  Sequence Sketching: {:.3}s", total_stats.t_sketch.as_secs_f64());
        eprintln!("  Seeding (Lookup):   {:.3}s", total_stats.t_seed.as_secs_f64());
        eprintln!("  Chaining:           {:.3}s", total_stats.t_chain.as_secs_f64());
        eprintln!("  Alignment (Ext):    {:.3}s", total_stats.t_align.as_secs_f64());
        eprintln!("  Post-chain:         {:.3}s", total_stats.t_post.as_secs_f64());
        eprintln!("  Total Measured:     {:.3}s",
            (total_stats.t_sketch + total_stats.t_seed + total_stats.t_chain + total_stats.t_align + total_stats.t_post).as_secs_f64()
        );
        eprintln!("  Reads: {}, Seeds: {}, Anchors: {}, Chains: {}",
            total_stats.n_reads, total_stats.n_seeds, total_stats.n_anchors, total_stats.n_chains);
        eprintln!("---------------------------");
        rammap::align::map::print_post_chain_breakdown();
    }

    Ok(())
}

/// Run-level configuration fields that stay constant across index parts.
struct RunConfig<'a> {
    mid_occ_frac: f32,
    cli: &'a AlignArgs,
    query2_path: Option<&'a str>,
    pe_mode: bool,
    two_file_pe: bool,
    copy_comment: bool,
    out_cfg: &'a OutputConfig,
    split_out_cfg: &'a OutputConfig,
}

/// Map all queries against one index part. Re-opens query files for each part.
/// Returns (stats, n_reads_processed).
fn map_one_part(
    mi: &Index,
    opt: &mut MapOptions,
    run: &RunConfig,
    query_path: &str,
    part_number: usize,  // 1-based
    sam_header_written: &mut bool,
    handle: &mut Box<dyn Write + Send>,
) -> anyhow::Result<AlignmentStats> {
    let mid_occ_frac = run.mid_occ_frac;
    let cli = run.cli;
    let query2_path = run.query2_path;
    let pe_mode = run.pe_mode;
    let two_file_pe = run.two_file_pe;
    let copy_comment = run.copy_comment;
    let out_cfg = run.out_cfg;
    let mini_batch_size = opt.mini_batch_size as u64;
    // mid_occ per-part: only calculate if still <= 0.
    // Once calculated for the first part, it stays for subsequent parts (options.c:73).
    if opt.seeding.mid_occ == 0 {
        opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
    }

    // Load junction annotations per-part (rids are local to this part)
    let junc_db: Option<JunctionDb> = if cli.junc_bed.is_some() || cli.spsc.is_some() {
        let name_to_rid: std::collections::HashMap<String, usize> = mi.seqs.iter().enumerate()
            .map(|(i, s)| (s.name.clone(), i)).collect();
        if let Some(bed_path) = &cli.junc_bed {
            Some(junc::load_bed_junctions(bed_path, &name_to_rid, mi.seqs.len())?)
        } else if let Some(spsc_path) = &cli.spsc {
            let max_sc = junc::max_spsc_bonus(opt.scoring.gap_open2, opt.scoring.gap_open);
            let scale = cli.spsc_scale.unwrap_or(0.7);
            let seq_lens: Vec<usize> = mi.seqs.iter().map(|s| s.len).collect();
            Some(junc::load_spsc_scores(spsc_path, &name_to_rid, mi.seqs.len(), &seq_lens, max_sc, scale)?)
        } else {
            None
        }
    } else {
        None
    };

    // Load jump junction database (-j / --pass1) per-part
    let jump_db: Option<JumpDb> = {
        let mut db: Option<JumpDb> = None;
        if let Some(ref path) = cli.junc_jump {
            match JumpDb::load(mi, path, 0x1, -1) {
                Ok(j) => db = Some(j),
                Err(e) => eprintln!("[WARNING] failed to load jump BED file: {}", e),
            }
        }
        if let Some(ref path) = cli.pass1 {
            match JumpDb::load(mi, path, 0x2, 5) {
                Ok(j) => {
                    if let Some(existing) = &mut db {
                        existing.merge(&j);
                    } else {
                        db = Some(j);
                    }
                }
                Err(e) => eprintln!("[WARNING] failed to load pass1 BED file: {}", e),
            }
        }
        db
    };

    // SAM header: write once for single-part, or per-part for multi-part without split-prefix
    if cli.output_sam && !*sam_header_written {
        writeln!(handle, "@HD\tVN:1.6\tSO:unsorted\tGO:query")?;
        for s in &mi.seqs {
            writeln!(handle, "@SQ\tSN:{}\tLN:{}", s.name, s.len)?;
        }
        let version = env!("CARGO_PKG_VERSION");
        let cmd_line = std::env::args().collect::<Vec<String>>().join(" ");
        writeln!(handle, "@PG\tID:rammap\tPN:rammap\tVN:{}\tCL:{}", version, cmd_line)?;
        if let Some(rg_line) = &cli.sam_rg {
            let expanded = rg_line.replace("\\t", "\t");
            writeln!(handle, "{}", expanded)?;
        }
        *sam_header_written = true;
    }

    // Re-open query files for this part.
    // stdin can't be re-read for subsequent parts.
    if query_path == "-" && part_number > 1 {
        // If stdin is being used and this is not the first part, skip mapping silently.
        // This happens with -d (dump index only) or if stdin was already consumed by part 1.
        return Ok(AlignmentStats::default());
    }

    let mut reader = rammap::fasta::open(query_path)
        .map_err(|e| anyhow::anyhow!("Error reading query: {}", e))?;
    // Single-end reads only need the header comment for -y; skip copying it otherwise.
    if !pe_mode && !copy_comment { reader.set_keep_description(false); }
    let reader2 = if two_file_pe {
        let q2 = query2_path.unwrap();
        Some(rammap::fasta::open(q2)
            .map_err(|e| anyhow::anyhow!("Error reading query2: {}", e))?)
    } else {
        None
    };

    // Parallel branches below assign `total_stats` unconditionally from the
    // thread::scope result; only the non-parallel accumulation path needs the
    // default initial value.
    #[cfg(feature = "parallel")]
    let total_stats: AlignmentStats;
    #[cfg(not(feature = "parallel"))]
    let mut total_stats = AlignmentStats::default();
    let mut record_iter = reader.records();
    let mut record_iter2 = reader2.map(|r| r.records());

    if pe_mode {
        type PairData = (String, Vec<u8>, Option<String>, Option<String>, String, Vec<u8>, Option<String>, Option<String>);

        // Three-stage pipeline: reader thread → worker pool → writer thread.
        // sync_channel(1) gives reader 1 chunk lookahead. The writer thread
        // runs on its own and consumes results so the worker pool isn't idle
        // while output is serialized.
        #[cfg(feature = "parallel")]
        {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<PairData>>(1);
            let (tx_out, rx_out) = std::sync::mpsc::sync_channel::<Vec<(String, AlignmentStats)>>(2);

            total_stats = std::thread::scope(|s| -> anyhow::Result<AlignmentStats> {
                s.spawn(move || {
                    loop {
                        let mut chunk_data: Vec<PairData> = Vec::new();
                        let mut chunk_bases: u64 = 0;
                        loop {
                            let r1 = record_iter.next();
                            let rec1 = match r1 {
                                Some(Ok(r)) => r,
                                Some(Err(e)) => { eprintln!("Warning: Error reading record: {}", e); continue; },
                                None => break,
                            };

                            let r2 = if two_file_pe {
                                record_iter2.as_mut().unwrap().next()
                            } else {
                                record_iter.next()
                            };
                            let rec2 = match r2 {
                                Some(Ok(r)) => r,
                                Some(Err(e)) => { eprintln!("Warning: Error reading R2 record: {}", e); continue; },
                                None => {
                                    eprintln!("Warning: Odd number of reads in PE mode, last read ignored");
                                    break;
                                }
                            };

                            let q1 = rec1.quality().map(|qs| String::from_utf8_lossy(qs).to_string());
                            let c1 = if copy_comment { rec1.description().map(|s| s.to_string()) } else { None };
                            let q2 = rec2.quality().map(|qs| String::from_utf8_lossy(qs).to_string());
                            let c2 = if copy_comment { rec2.description().map(|s| s.to_string()) } else { None };
                            chunk_bases += (rec1.sequence().len() + rec2.sequence().len()) as u64;
                            chunk_data.push((
                                rec1.name().to_string(), rec1.sequence().to_vec(), q1, c1,
                                rec2.name().to_string(), rec2.sequence().to_vec(), q2, c2,
                            ));
                            if chunk_bases >= mini_batch_size { break; }
                        }
                        let done = chunk_data.is_empty();
                        if tx.send(chunk_data).is_err() { break; }
                        if done { break; }
                    }
                });

                let writer_handle = s.spawn(move || -> anyhow::Result<AlignmentStats> {
                    let mut acc = AlignmentStats::default();
                    while let Ok(results) = rx_out.recv() {
                        for (res, stats) in results {
                            if !res.is_empty() {
                                write!(handle, "{}", res)?;
                            }
                            acc = acc + stats;
                        }
                    }
                    Ok(acc)
                });

                while let Ok(chunk_data) = rx.recv() {
                    if chunk_data.is_empty() { break; }

                    let results: Vec<(String, AlignmentStats)> = chunk_data.par_iter().map_init(
                        || (rammap::align::extend::AlignmentContext::new(), rammap::align::map::MapContext::new()),
                        |(ctx, map_ctx), (qname1, qseq1, qual1, comment1, qname2, qseq2, qual2, comment2)| {
                            let r1 = ReadInfo { qname: qname1, qseq: qseq1, qual: qual1.as_deref(), comment: comment1.as_deref(), n_seg: 2, seg_idx: 0 };
                            let r2 = ReadInfo { qname: qname2, qseq: qseq2, qual: qual2.as_deref(), comment: comment2.as_deref(), n_seg: 2, seg_idx: 1 };
                            match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                rammap::align::pipeline::align_and_format_pair(
                                    opt, mi, &r1, &r2,
                                    ctx, map_ctx, junc_db.as_ref(), out_cfg,
                                )
                            })) {
                                Ok(result) => result,
                                Err(_) => {
                                    eprintln!("[WARNING] alignment panicked for pair: {}/{}", qname1, qname2);
                                    (String::new(), AlignmentStats::default())
                                }
                            }
                        }
                    ).collect();

                    if tx_out.send(results).is_err() { break; }
                }
                drop(tx_out);
                writer_handle.join().expect("writer thread panicked")
            })?;
        }

        #[cfg(not(feature = "parallel"))]
        {
            loop {
                let mut chunk_data: Vec<PairData> = Vec::new();
                let mut chunk_bases: u64 = 0;
                loop {
                    let r1 = if two_file_pe {
                        record_iter.next()
                    } else {
                        record_iter.next()
                    };
                    let rec1 = match r1 {
                        Some(Ok(r)) => r,
                        Some(Err(e)) => { eprintln!("Warning: Error reading record: {}", e); continue; },
                        None => break,
                    };

                    let r2 = if two_file_pe {
                        record_iter2.as_mut().unwrap().next()
                    } else {
                        record_iter.next()
                    };
                    let rec2 = match r2 {
                        Some(Ok(r)) => r,
                        Some(Err(e)) => { eprintln!("Warning: Error reading R2 record: {}", e); continue; },
                        None => {
                            eprintln!("Warning: Odd number of reads in PE mode, last read ignored");
                            break;
                        }
                    };

                    let q1 = rec1.quality().map(|qs| String::from_utf8_lossy(qs).to_string());
                    let c1 = if copy_comment { rec1.description().map(|s| s.to_string()) } else { None };
                    let q2 = rec2.quality().map(|qs| String::from_utf8_lossy(qs).to_string());
                    let c2 = if copy_comment { rec2.description().map(|s| s.to_string()) } else { None };
                    chunk_bases += (rec1.sequence().len() + rec2.sequence().len()) as u64;
                    chunk_data.push((
                        rec1.name().to_string(), rec1.sequence().to_vec(), q1, c1,
                        rec2.name().to_string(), rec2.sequence().to_vec(), q2, c2,
                    ));
                    if chunk_bases >= mini_batch_size { break; }
                }

                if chunk_data.is_empty() { break; }

                let results: Vec<(String, AlignmentStats)> = {
                    let mut ctx = rammap::align::extend::AlignmentContext::new();
                    let mut map_ctx = rammap::align::map::MapContext::new();
                    chunk_data.iter().map(|(qname1, qseq1, qual1, comment1, qname2, qseq2, qual2, comment2)| {
                        let r1 = ReadInfo { qname: qname1, qseq: qseq1, qual: qual1.as_deref(), comment: comment1.as_deref(), n_seg: 2, seg_idx: 0 };
                        let r2 = ReadInfo { qname: qname2, qseq: qseq2, qual: qual2.as_deref(), comment: comment2.as_deref(), n_seg: 2, seg_idx: 1 };
                        rammap::align::pipeline::align_and_format_pair(
                            opt, mi, &r1, &r2,
                            &mut ctx, &mut map_ctx, junc_db.as_ref(), out_cfg,
                        )
                    }).collect()
                };

                for (res, stats) in results {
                    write!(handle, "{}", res)?;
                    total_stats = total_stats + stats;
                }
            }
        }
    } else {
        // Three-stage pipeline: reader thread → worker pool → writer thread.
        // Output is decoupled from the worker pool so cores stay saturated
        // through chunk transitions.
        #[cfg(feature = "parallel")]
        {
            type SeData = (String, Vec<u8>, Option<String>, Option<String>);
            // Chunks the reader may run ahead: input decompression is steady but alignment speed
            // varies with the reads, so a deeper queue absorbs fast stretches. About 800 Mbp of
            // reads (~1.6 GB with qualities), 2 to 16 chunks; env RAMMAP_LOOKAHEAD overrides.
            let lookahead: usize = std::env::var("RAMMAP_LOOKAHEAD").ok().and_then(|v| v.parse().ok())
                .unwrap_or_else(|| (800_000_000 / mini_batch_size.max(1)).clamp(2, 16) as usize);
            let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<SeData>>(lookahead);
            let (tx_out, rx_out) = std::sync::mpsc::sync_channel::<Vec<(String, AlignmentStats)>>(2);

            total_stats = std::thread::scope(|s| -> anyhow::Result<AlignmentStats> {
                s.spawn(move || {
                    loop {
                        let mut chunk_data: Vec<SeData> = Vec::new();
                        let mut chunk_bases: u64 = 0;
                        loop {
                            match record_iter.next() {
                                Some(Ok(r)) => {
                                    let (name, desc, seq, qual) = r.into_parts();
                                    let q = qual.map(|v| String::from_utf8(v).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned()));
                                    let c = if copy_comment { desc } else { None };
                                    chunk_bases += seq.len() as u64;
                                    chunk_data.push((name, seq, q, c));
                                },
                                Some(Err(e)) => {
                                    eprintln!("Warning: Error reading record: {}", e);
                                    continue;
                                },
                                None => break,
                            }
                            if chunk_bases >= mini_batch_size { break; }
                        }
                        let done = chunk_data.is_empty();
                        if tx.send(chunk_data).is_err() { break; }
                        if done { break; }
                    }
                });

                let writer_handle = s.spawn(move || -> anyhow::Result<AlignmentStats> {
                    let mut acc = AlignmentStats::default();
                    while let Ok(results) = rx_out.recv() {
                        for (res, stats) in results {
                            if !res.is_empty() {
                                write!(handle, "{}", res)?;
                            }
                            acc = acc + stats;
                        }
                    }
                    Ok(acc)
                });

                // Each iteration aligns one new chunk *plus* the re-runs of earlier chunks' pending
                // reads, so with the GPU path on (`--features gpu`, RAMMAP_GPU=1) the GPU works on
                // one round's requests while the CPU works on the next round. A read that needs an
                // uncached GPU DP result finishes as a speculative dry run that queues its GPU
                // requests; its output is discarded and it is re-run in the next iteration, after
                // `wait_pending` has cached those results. Chunks are written strictly in order,
                // once all their reads are done. Without the GPU path every read completes in its
                // chunk's first iteration, which is the original one-chunk-at-a-time behaviour.
                #[cfg(feature = "gpu")]
                let gpu_on = rammap::align::dp::gpu::enabled();
                #[cfg(not(feature = "gpu"))]
                #[allow(unused_variables)]
                let gpu_on = false;
                let verbose = gpu_on && std::env::var_os("RAMMAP_GPU_VERBOSE").is_some();
                // A read still pending after this many speculative attempts (chains of z-drop
                // mispredictions, typically one ultra-long read) finishes on the CPU instead.
                const MAX_GPU_ATTEMPTS: u8 = 3;
                let align_one = |ctx: &mut rammap::align::extend::AlignmentContext,
                                 map_ctx: &mut rammap::align::map::MapContext,
                                 (qname, qseq, qual, comment): &SeData, use_gpu: bool| {
                    let ri = ReadInfo { qname, qseq, qual: qual.as_deref(), comment: comment.as_deref(), n_seg: 1, seg_idx: 0 };
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        #[allow(unused_mut)]
                        let mut run = || rammap::align::pipeline::align_and_format_query(
                            opt, mi, &ri, ctx, map_ctx, junc_db.as_ref(),
                            jump_db.as_ref(), out_cfg,
                        );
                        #[cfg(feature = "gpu")]
                        if use_gpu { return rammap::align::dp::gpu::retry_scope(run); }
                        let _ = use_gpu;
                        run()
                    }))
                };
                // CPU-only: stream reads into the pool chunk by chunk, so the stragglers at the
                // end of one chunk overlap with later chunks instead of leaving cores idle at a
                // per-chunk barrier. Two chunks are always admitted; beyond that a new chunk is
                // admitted whenever fewer than 2 x threads reads are left queued or running
                // (e.g. one long read still holding the oldest chunk), up to MAX_IN_FLIGHT.
                // Chunks are still written in order.
                if !gpu_on {
                    use std::sync::atomic::{AtomicUsize, Ordering};
                    struct Chunk {
                        data: Vec<SeData>,
                        done: Vec<std::sync::OnceLock<(String, AlignmentStats)>>,
                        remaining: AtomicUsize,
                    }
                    thread_local! {
                        static CTX: std::cell::RefCell<Option<(rammap::align::extend::AlignmentContext, rammap::align::map::MapContext)>> =
                            const { std::cell::RefCell::new(None) };
                    }
                    const MIN_IN_FLIGHT: usize = 2;
                    const MAX_IN_FLIGHT: usize = 8;
                    let low_water = 2 * rayon::current_num_threads();
                    // Reads admitted but not yet aligned, over all chunks in flight.
                    let outstanding = std::sync::Arc::new(AtomicUsize::new(0));
                    // usize::MAX = "outstanding fell below the low-water mark", else a finished chunk id.
                    let (tx_done, rx_done) = std::sync::mpsc::channel::<usize>();
                    let align_one = &align_one;
                    let mut send_failed = false;
                    rayon::in_place_scope(|sc| {
                        let mut in_flight: std::collections::VecDeque<(usize, std::sync::Arc<Chunk>)> = std::collections::VecDeque::new();
                        let mut finished: std::collections::HashSet<usize> = std::collections::HashSet::new();
                        let mut next_id = 0usize;
                        let mut input_done = false;
                        loop {
                            // Write out completed chunks in order.
                            while let Some((id, _)) = in_flight.front() {
                                if !finished.remove(id) { break; }
                                let (_, mut c) = in_flight.pop_front().unwrap();
                                // Workers drop their handle just after the final decrement.
                                let c = loop {
                                    match std::sync::Arc::try_unwrap(c) {
                                        Ok(c) => break c,
                                        Err(back) => { c = back; std::thread::yield_now(); }
                                    }
                                };
                                let results: Vec<(String, AlignmentStats)> = c.done.into_iter().map(|r| r.into_inner().expect("every read aligned")).collect();
                                if tx_out.send(results).is_err() { send_failed = true; break; }
                            }
                            if send_failed { break; }
                            let admit = in_flight.len() < MIN_IN_FLIGHT
                                || (in_flight.len() < MAX_IN_FLIGHT && outstanding.load(Ordering::Acquire) < low_water);
                            if !input_done && admit {
                                match rx.recv() {
                                    Ok(chunk_data) if !chunk_data.is_empty() => {
                                        let n = chunk_data.len();
                                        let id = next_id;
                                        next_id += 1;
                                        let c = std::sync::Arc::new(Chunk { data: chunk_data, done: (0..n).map(|_| std::sync::OnceLock::new()).collect(),
                                                                            remaining: AtomicUsize::new(n) });
                                        outstanding.fetch_add(n, Ordering::AcqRel);
                                        for ri in 0..n {
                                            let c = std::sync::Arc::clone(&c);
                                            let tx_done = tx_done.clone();
                                            let outstanding = std::sync::Arc::clone(&outstanding);
                                            sc.spawn(move |_| {
                                                let read = &c.data[ri];
                                                let r = CTX.with(|cell| {
                                                    let mut g = cell.borrow_mut();
                                                    let (ctx, map_ctx) = g.get_or_insert_with(|| (rammap::align::extend::AlignmentContext::new(), rammap::align::map::MapContext::new()));
                                                    align_one(ctx, map_ctx, read, false)
                                                });
                                                let r = r.unwrap_or_else(|_| {
                                                    eprintln!("[WARNING] alignment panicked for read: {}", read.0);
                                                    (String::new(), AlignmentStats::default())
                                                });
                                                let _ = c.done[ri].set(r);
                                                let last = c.remaining.fetch_sub(1, Ordering::AcqRel) == 1;
                                                drop(c);
                                                if last { let _ = tx_done.send(id); }
                                                if outstanding.fetch_sub(1, Ordering::AcqRel) == low_water { let _ = tx_done.send(usize::MAX); }
                                            });
                                        }
                                        in_flight.push_back((id, c));
                                        continue;
                                    }
                                    _ => input_done = true,
                                }
                            }
                            if in_flight.is_empty() { break; }
                            match rx_done.recv() {
                                Ok(usize::MAX) => {} // below the low-water mark: re-check admission
                                Ok(id) => { finished.insert(id); }
                                Err(_) => break,
                            }
                        }
                    });
                    drop(tx_out);
                    return writer_handle.join().expect("writer thread panicked");
                }
                struct OpenChunk {
                    data: Vec<SeData>,
                    done: Vec<Option<(String, AlignmentStats)>>,
                    pending: Vec<usize>,
                    attempts: Vec<u8>,
                }
                let mut open: std::collections::VecDeque<OpenChunk> = std::collections::VecDeque::new();
                let mut input_done = false;
                for iter in 0usize.. {
                    if !input_done {
                        match rx.recv() {
                            Ok(chunk_data) if !chunk_data.is_empty() => {
                                let n = chunk_data.len();
                                open.push_back(OpenChunk { data: chunk_data, done: (0..n).map(|_| None).collect(),
                                                           pending: (0..n).collect(), attempts: vec![0; n] });
                            }
                            _ => input_done = true,
                        }
                    }
                    if open.is_empty() { break; }
                    assert!(iter < 1_000_000, "GPU retry did not converge");

                    // Results for every GPU request queued in the previous iteration.
                    let t_wait = std::time::Instant::now();
                    #[cfg(feature = "gpu")]
                    if gpu_on { rammap::align::dp::gpu::wait_pending(); }
                    let wait_s = t_wait.elapsed().as_secs_f64();

                    let work: Vec<(usize, usize)> = open.iter().enumerate()
                        .flat_map(|(ci, o)| o.pending.iter().map(move |&ri| (ci, ri)))
                        .collect();
                    for o in open.iter_mut() { o.pending.clear(); }
                    let t_cpu = std::time::Instant::now();
                    let attempted: Vec<(usize, usize, Option<(String, AlignmentStats)>)> = work.par_iter().map_init(
                        || (rammap::align::extend::AlignmentContext::new(), rammap::align::map::MapContext::new()),
                        |(ctx, map_ctx), &(ci, ri)| {
                            let read = &open[ci].data[ri];
                            let use_gpu = gpu_on && open[ci].attempts[ri] < MAX_GPU_ATTEMPTS;
                            let attempt = align_one(ctx, map_ctx, read, use_gpu);
                            #[cfg(feature = "gpu")]
                            if use_gpu && rammap::align::dp::gpu::was_speculative() { return (ci, ri, None); }
                            match attempt {
                                Ok(result) => (ci, ri, Some(result)),
                                Err(_) => {
                                    eprintln!("[WARNING] alignment panicked for read: {}", read.0);
                                    (ci, ri, Some((String::new(), AlignmentStats::default())))
                                }
                            }
                        }
                    ).collect();
                    let cpu_s = t_cpu.elapsed().as_secs_f64();
                    for (ci, ri, r) in attempted {
                        match r {
                            Some(x) => open[ci].done[ri] = Some(x),
                            None => { open[ci].pending.push(ri); open[ci].attempts[ri] += 1; }
                        }
                    }
                    if verbose {
                        let n_pend: usize = open.iter().map(|o| o.pending.len()).sum();
                        eprintln!("[gpu] iter {iter}: gpu wait {wait_s:.2}s, cpu {cpu_s:.2}s on {} reads, {} still pending, {} open chunks",
                                  work.len(), n_pend, open.len());
                    }
                    let mut send_failed = false;
                    while open.front().is_some_and(|o| o.pending.is_empty()) {
                        let o = open.pop_front().unwrap();
                        let results: Vec<(String, AlignmentStats)> = o.done.into_iter().map(|r| r.expect("every read aligned")).collect();
                        if tx_out.send(results).is_err() { send_failed = true; break; }
                    }
                    if send_failed { break; }
                    #[cfg(feature = "gpu")]
                    if gpu_on && open.is_empty() { rammap::align::dp::gpu::clear_cache(); }
                }
                drop(tx_out);
                writer_handle.join().expect("writer thread panicked")
            })?;
        }

        #[cfg(not(feature = "parallel"))]
        {
            loop {
                let mut chunk_data: Vec<(String, Vec<u8>, Option<String>, Option<String>)> = Vec::new();
                let mut chunk_bases: u64 = 0;
                loop {
                    match record_iter.next() {
                        Some(Ok(r)) => {
                            let (name, desc, seq, qual) = r.into_parts();
                            let q = qual.map(|v| String::from_utf8(v).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned()));
                            let c = if copy_comment { desc } else { None };
                            chunk_bases += seq.len() as u64;
                            chunk_data.push((name, seq, q, c));
                        },
                        Some(Err(e)) => {
                            eprintln!("Warning: Error reading record: {}", e);
                            continue;
                        },
                        None => break,
                    }
                    if chunk_bases >= mini_batch_size { break; }
                }

                if chunk_data.is_empty() { break; }

                let results: Vec<(String, AlignmentStats)> = {
                    let mut ctx = rammap::align::extend::AlignmentContext::new();
                    let mut map_ctx = rammap::align::map::MapContext::new();
                    chunk_data.iter().map(|(qname, qseq, qual, comment)| {
                        let ri = ReadInfo { qname, qseq, qual: qual.as_deref(), comment: comment.as_deref(), n_seg: 1, seg_idx: 0 };
                         rammap::align::pipeline::align_and_format_query(
                            opt, mi, &ri, &mut ctx, &mut map_ctx, junc_db.as_ref(),
                            jump_db.as_ref(), out_cfg,
                        )
                    }).collect()
                };

                for (res, stats) in results {
                    write!(handle, "{}", res)?;
                    total_stats = total_stats + stats;
                }
            }
        }
    }

    Ok(total_stats)
}

/// Map all queries against one index part, writing raw results to a split temp file.
/// Used when --split-prefix is set for multi-part merge.
fn map_one_part_split(
    mi: &Index,
    opt: &mut MapOptions,
    run: &RunConfig,
    query_path: &str,
    part_number: usize,  // 1-based
    split_writer: &mut std::io::BufWriter<std::fs::File>,
) -> anyhow::Result<AlignmentStats> {
    use rammap::align::pipeline::{process_query, process_query_from_regs};
    use rammap::align::map::{MapContext, map_query_multi};
    use rammap::align::extend::rev_comp;
    let mid_occ_frac = run.mid_occ_frac;
    let cli = run.cli;
    let query2_path = run.query2_path;
    let pe_mode = run.pe_mode;
    let two_file_pe = run.two_file_pe;
    let out_cfg = run.split_out_cfg;

    // mid_occ per-part: only calculate if still <= 0
    if opt.seeding.mid_occ == 0 {
        opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
    }

    // Load junction annotations per-part
    let junc_db: Option<JunctionDb> = if cli.junc_bed.is_some() || cli.spsc.is_some() {
        let name_to_rid: std::collections::HashMap<String, usize> = mi.seqs.iter().enumerate()
            .map(|(i, s)| (s.name.clone(), i)).collect();
        if let Some(bed_path) = &cli.junc_bed {
            Some(junc::load_bed_junctions(bed_path, &name_to_rid, mi.seqs.len())?)
        } else if let Some(spsc_path) = &cli.spsc {
            let max_sc = junc::max_spsc_bonus(opt.scoring.gap_open2, opt.scoring.gap_open);
            let scale = cli.spsc_scale.unwrap_or(0.7);
            let seq_lens: Vec<usize> = mi.seqs.iter().map(|s| s.len).collect();
            Some(junc::load_spsc_scores(spsc_path, &name_to_rid, mi.seqs.len(), &seq_lens, max_sc, scale)?)
        } else {
            None
        }
    } else {
        None
    };

    // Load jump junction database (-j / --pass1) per-part
    let jump_db: Option<JumpDb> = {
        let mut db: Option<JumpDb> = None;
        if let Some(ref path) = cli.junc_jump {
            match JumpDb::load(mi, path, 0x1, -1) {
                Ok(j) => db = Some(j),
                Err(e) => eprintln!("[WARNING] failed to load jump BED file: {}", e),
            }
        }
        if let Some(ref path) = cli.pass1 {
            match JumpDb::load(mi, path, 0x2, 5) {
                Ok(j) => {
                    if let Some(existing) = &mut db {
                        existing.merge(&j);
                    } else {
                        db = Some(j);
                    }
                }
                Err(e) => eprintln!("[WARNING] failed to load pass1 BED file: {}", e),
            }
        }
        db
    };

    // stdin can't be re-read for subsequent parts
    if query_path == "-" && part_number > 1 {
        return Ok(AlignmentStats::default());
    }

    let reader = rammap::fasta::open(query_path)
        .map_err(|e| anyhow::anyhow!("Error reading query: {}", e))?;
    let reader2 = if two_file_pe {
        let q2 = query2_path.unwrap();
        Some(rammap::fasta::open(q2)
            .map_err(|e| anyhow::anyhow!("Error reading query2: {}", e))?)
    } else {
        None
    };

    let mut total_stats = AlignmentStats::default();
    let mut record_iter = reader.records();
    let mut record_iter2 = reader2.map(|r| r.records());
    let is_weak = opt.flags.contains(AlignFlags::WEAK_PAIRING);
    let is_independ = opt.flags.contains(AlignFlags::INDEPEND_SEG);
    let pe_ori = opt.pairing.pe_ori;
    let flip_r1 = (pe_ori >> 1) & 1 != 0;
    let flip_r2 = pe_ori & 1 != 0;

    if pe_mode {
        let mut ctx = rammap::align::extend::AlignmentContext::new();
        let mut map_ctx = MapContext::new();

        loop {
            let rec1 = match record_iter.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => { eprintln!("Warning: Error reading record: {}", e); continue; },
                None => break,
            };
            let rec2 = if two_file_pe {
                match record_iter2.as_mut().unwrap().next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => { eprintln!("Warning: Error reading R2 record: {}", e); continue; },
                    None => break,
                }
            } else {
                match record_iter.next() {
                    Some(Ok(r)) => r,
                    Some(Err(e)) => { eprintln!("Warning: Error reading R2 record: {}", e); continue; },
                    None => break,
                }
            };

            let qname1 = rec1.name().to_string();
            let qseq1 = rec1.sequence().to_vec();
            let qname2 = rec2.name().to_string();
            let qseq2 = rec2.sequence().to_vec();

            let qseq1_work = if flip_r1 { rev_comp(&qseq1) } else { qseq1.clone() };
            let qseq2_work = if flip_r2 { rev_comp(&qseq2) } else { qseq2.clone() };

            let (pq1, pq2, frag_gap) = if is_independ {
                let pq1 = process_query(opt, mi, &qname1, &qseq1_work, &mut ctx, &mut map_ctx,
                    junc_db.as_ref(), out_cfg);
                let pq2 = process_query(opt, mi, &qname2, &qseq2_work, &mut ctx, &mut map_ctx,
                    junc_db.as_ref(), out_cfg);
                (pq1, pq2, opt.chaining.max_gap_ref)
            } else if is_weak {
                let pq1 = process_query(opt, mi, &qname1, &qseq1_work, &mut ctx, &mut map_ctx,
                    junc_db.as_ref(), out_cfg);
                let pq2 = process_query(opt, mi, &qname1, &qseq2_work, &mut ctx, &mut map_ctx,
                    junc_db.as_ref(), out_cfg);
                (pq1, pq2, opt.chaining.max_gap_ref)
            } else {
                let seqs: Vec<&[u8]> = vec![&qseq1_work, &qseq2_work];
                let qlens = vec![qseq1_work.len(), qseq2_work.len()];
                let multi = map_query_multi(opt, mi, &qname1, &seqs, &qlens, &mut map_ctx);
                let shared_rep_len = multi.rep_len;
                let frag_gap = multi.frag_gap;
                let mut per_seg = multi.per_seg;
                let (regs2, _) = if per_seg.len() > 1 { per_seg.remove(1) } else { (Vec::new(), Vec::new()) };
                let (regs1, _) = if !per_seg.is_empty() { per_seg.remove(0) } else { (Vec::new(), Vec::new()) };
                let pq1 = process_query_from_regs(opt, mi, &qseq1_work, &mut ctx, &mut map_ctx,
                    junc_db.as_ref(), out_cfg, regs1, shared_rep_len);
                let pq2 = process_query_from_regs(opt, mi, &qseq2_work, &mut ctx, &mut map_ctx,
                    junc_db.as_ref(), out_cfg, regs2, shared_rep_len);
                (pq1, pq2, frag_gap)
            };

            // Write both segments to temp file
            split::split_write_query(split_writer, &pq1.results, pq1.rep_len, frag_gap)
                .map_err(|e| anyhow::anyhow!("Error writing split temp: {}", e))?;
            split::split_write_query(split_writer, &pq2.results, pq2.rep_len, frag_gap)
                .map_err(|e| anyhow::anyhow!("Error writing split temp: {}", e))?;

            total_stats = total_stats + pq1.stats + pq2.stats;
        }
    } else {
        // Single-end split write
        let mut ctx = rammap::align::extend::AlignmentContext::new();
        let mut map_ctx = MapContext::new();

        loop {
            let rec = match record_iter.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => { eprintln!("Warning: Error reading record: {}", e); continue; },
                None => break,
            };

            let qname = rec.name().to_string();
            let qseq = rec.sequence().to_vec();

            let mut pq = process_query(opt, mi, &qname, &qseq, &mut ctx, &mut map_ctx,
                junc_db.as_ref(), out_cfg);

            // Jump splice extension for single-end split.
            if let Some(ref jdb) = jump_db {
                let is_splice = opt.flags.contains(AlignFlags::SPLICE);
                if is_splice {
                    for r in pq.results.iter_mut() {
                        rammap::align::jump::jump_split(mi, opt, qseq.len(), &qseq, r, jdb);
                    }
                }
            }

            split::split_write_query(split_writer, &pq.results, pq.rep_len, 0)
                .map_err(|e| anyhow::anyhow!("Error writing split temp: {}", e))?;

            total_stats = total_stats + pq.stats;
        }
    }

    split_writer.flush()?;
    Ok(total_stats)
}
