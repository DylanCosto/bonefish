use clap::Parser;
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

use bonefish::align::index::Index;
use bonefish::align::map::{MapOptions, AlignFlags};
use bonefish::align::pipeline::{OutputConfig, ReadInfo};
use bonefish::align::stats::AlignmentStats;
use bonefish::align::split;

// ─── CLI definition ───

#[derive(Parser, Debug)]
#[command(name = "bonefish")]
#[command(version, about = "bonefish: a fast aligner for Oxford Nanopore DNA reads")]
struct Cli {
    #[command(flatten)]
    align: AlignArgs,

    /// Verbose output
    #[arg(short, long, global = true)]
    verbose: bool,
}

#[derive(Parser, Debug)]
struct AlignArgs {
    /// Reference FASTA (plain or gzipped) or an index built with -d (.idx; raw .raw with --idx-raw)
    pub target: String,

    /// Reads: FASTQ, FASTA or BAM (e.g. dorado's unaligned BAM), plain or gzipped; - for stdin.
    /// Several files are mapped one after another.
    pub queries: Vec<String>,

    /// Dump index to file
    #[arg(short = 'd', long)]
    pub dump_index: Option<String>,

    /// Preset applied before other options (names ignore case):
    ///   dna / ont        — ONT reads (R10.4.1 or R9.4.1) vs a reference genome (default)
    ///   hifi             — PacBio HiFi reads vs a reference genome; alias: map-hifi
    ///   lr:hq / lr:hq-fast / map-ont / lr — accepted as aliases of dna
    ///   rna              — direct RNA vs a reference genome (spliced; minimap2 -x splice -uf -k14)
    ///   cdna             — cDNA vs a reference genome (spliced; minimap2 -x splice); alias: splice
    ///   splice:hq        — accurate cDNA (minimap2 -x splice:hq)
    ///   tx               — RNA/cDNA vs a transcriptome, for quantification (minimap2 -x map-ont -N 100)
    /// dna = minimap2 lr:hq seeding and scoring plus --ext-bw 200 --ext-adaptive 64 --ext-rescue
    /// --ext-probe 2000,1000 --gap-bw-margin 500 --gap-bw-narrow 32 --ungapped-fill --max-gap-ref 5000
    /// --seg-rmq --clean-dp --rescue-join --rescue-max-density 2 --gap-subanchor-tr 150
    /// cdna = -k15 -w5 -g2k -G200k -A1 -B2 -O2,32 -E1,0 -C9 -z200,100 -ub --junc-bonus 9
    /// --junc-pen 5 --splice-flank yes; rna = cdna with -k14 -uf
    /// tx = dna with -k15 -w10 -g5000 -U10,1000000 -N100
    /// hifi = dna, with stretches of fragmented long gaps realigned under map-hifi scoring
    /// (-A1 -B4 -O6,26 -E2,1); uses dna indexes
    #[arg(short = 'x', long, verbatim_doc_comment)]
    pub preset: Option<String>,

    /// Minimizer k-mer length
    #[arg(short = 'k')]
    pub k: Option<usize>,

    /// Minimizer window size
    #[arg(short = 'w')]
    pub w: Option<usize>,


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

    /// Copy read tags to the output: SAM-style tags from FASTQ/FASTA comments, and the
    /// per-read tags of BAM input except alignment tags (dorado's MM/ML/MN, RG, qs, ch, ...)
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
    #[arg(short = 't', long = "threads", default_value_t = 3, value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..=4096))]
    pub align_threads: usize,

    /// Adaptive mid_occ floor,ceiling (format: min,max)
    #[arg(short = 'U')]
    pub mid_occ_range: Option<String>,

    /// Spliced presets: max intron length (sets the chaining reference gap and bandwidth) [200k]
    #[arg(short = 'G', value_name = "NUM")]
    pub max_intron_len: Option<String>,

    /// Only map to forward strand
    #[arg(long)]
    pub for_only: bool,

    /// Only map to reverse strand
    #[arg(long, conflicts_with = "for_only")]
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

    /// Spliced presets: cost of a non-canonical splice site
    #[arg(short = 'C')]
    pub noncan: Option<i32>,

    /// Spliced presets: transcript strands to score splice signals for: f (forward, direct
    /// RNA), r (reverse), b (both, the better one is kept; cDNA) or n (none)
    #[arg(short = 'u', value_name = "f|r|b|n")]
    pub splice_strand: Option<String>,

    /// Spliced presets: annotated junctions (BED12 transcripts, or one intron per BED6 line)
    #[arg(long, value_name = "FILE")]
    pub junc_bed: Option<String>,

    /// Spliced presets: score bonus for an annotated junction [9]
    #[arg(long, value_name = "INT")]
    pub junc_bonus: Option<i32>,

    /// Spliced presets: penalty for a junction absent from --junc-bed [5]
    #[arg(long, value_name = "INT")]
    pub junc_pen: Option<i32>,

    /// Spliced presets: score the bases flanking splice sites (yes/no) [yes]
    #[arg(long, value_name = "yes|no")]
    pub splice_flank: Option<String>,

    /// Spliced presets: write the introns of primary alignments (MAPQ >= 10) as BED instead of SAM/PAF
    #[arg(long)]
    pub write_junc: bool,

    /// Output to file instead of stdout (FILE.bam writes a coordinate-sorted BAM, as --sort)
    #[arg(short = 'o')]
    pub output: Option<String>,

    /// With -o FILE.bam: write a coordinate-sorted BAM and FILE.bam.bai (as `-a | samtools sort --write-index`)
    #[arg(long)]
    pub sort: bool,

    /// --sort: memory for record batches (two are alive at a time) [4G]
    #[arg(long, value_name = "NUM")]
    pub sort_mem: Option<String>,

    /// --sort: prefix of temporary run files [the output path]
    #[arg(long, value_name = "PREFIX")]
    pub sort_tmp: Option<String>,

    /// --sort: BGZF compression level of the BAM
    #[arg(long, value_name = "INT", default_value_t = 6)]
    pub sort_level: i32,

    /// Write CIGAR with >65535 ops at CG tag
    #[arg(short = 'L')]
    pub long_cigar: bool,

    /// Use soft clipping for supplementary alignments
    #[arg(short = 'Y')]
    pub softclip: bool,

    /// Don't output base quality in SAM
    #[arg(short = 'Q')]
    pub no_qual: bool,

    /// Skip seeds on the exact diagonal (self-mapping)
    #[arg(short = 'D', hide = true)]
    pub no_diag: bool,

    /// Seed occurrence distance
    #[arg(short = 'e')]
    pub occ_dist: Option<i32>,

    /// Minibatch size for mapping (query bases per batch; suffix K/M/G OK) [100M]
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

    /// DP bandwidth for read-end extensions [200]
    #[arg(long)]
    pub ext_bw: Option<i32>,

    /// Gap-fill bandwidth = min(default, |qlen-tlen| + INT), keeping large indels crossable
    #[arg(long)]
    pub gap_bw_margin: Option<i32>,

    /// Try an ungapped fill first for equal-length gaps [on]
    #[arg(long, hide = true)]
    pub ungapped_fill: bool,

    /// With --ext-bw: re-run a z-dropped extension with the default band when its unaligned tail
    /// shows an indel wider than the narrow band [on]
    #[arg(long, hide = true)]
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
    /// minimap2's SSE kernel behaviour [on]
    #[arg(long, hide = true)]
    pub clean_dp: bool,

    /// RMQ chaining with a segment tree (faster; ties between equal-priority predecessors may
    /// resolve differently from minimap2) [on]
    #[arg(long, hide = true)]
    pub seg_rmq: bool,

    /// Run the long-gap (RMQ) rescue re-chain only when two first-pass chains could be joined [on]
    #[arg(long, hide = true)]
    pub rescue_join: bool,

    /// Skip the long-join rescue for reads with >= INT anchors per query base (0 = never) [2 in dna]
    #[arg(long, value_name = "INT")]
    pub rescue_max_density: Option<u32>,

    /// Sub-anchored fill of long gaps (exact unique 15-mer segments, narrow DP between) [on]
    #[arg(long, hide = true)]
    pub gap_subanchor: bool,
    /// No sub-anchored gap fill (v0.1.4 behaviour; ~4% slower)
    #[arg(long, hide = true, conflicts_with = "gap_subanchor")]
    pub no_gap_subanchor: bool,
    /// Sub-anchoring skips tandem repeats of period up to INT bp; 0 = sub-anchor in them too (v0.1.3
    /// behaviour, where a segment on another repeat copy forces a deletion-insertion detour) [150]
    #[arg(long, hide = true, value_name = "INT")]
    pub gap_subanchor_tr: Option<u32>,

    /// Keep 1/F of the query minimizers (by hash) for reads of at least MIN_LEN bp [default 8000]:
    /// fewer anchors, much faster, loses some SV support (e.g. --query-subsample 2)
    #[arg(long, value_name = "F[,MIN_LEN]")]
    pub query_subsample: Option<String>,

    /// Turn off --rescue-join (e.g. to undo it in dna)
    #[arg(long, conflicts_with = "rescue_join")]
    pub no_rescue_join: bool,

    /// Score for ambiguous bases
    #[arg(long = "score-N")]
    pub score_n: Option<i32>,

    /// Bonus for reaching end of query
    #[arg(long)]
    pub end_bonus: Option<i32>,

    /// Penalty for seeds near read ends when extending
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

    /// Use the masking level as a hard filter for secondary chains
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
    #[arg(long, hide = true)]
    pub dual: Option<String>,

    /// Query occurrence fraction
    #[arg(long)]
    pub q_occ_frac: Option<f32>,






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

    /// Raw index: copy the seed hash table into transparent huge pages (Linux). About 1% less
    /// CPU from fewer TLB misses, for ~6.5 GB more memory (human) and a few seconds of loading.
    #[arg(long)]
    pub hugepages: bool,

    /// Write the -d index in the raw format, which loads by memory-mapping (near-instant).
    /// With a pre-built index as the target, -d converts it.
    #[arg(long)]
    pub idx_raw: bool,


    /// Don't use query name in read hash
    #[arg(long)]
    pub no_hash_name: bool,


    /// Query-strand alignment mode (for SV detection)
    #[arg(long)]
    pub qstrand: bool,

    /// Parallel chaining across reference sequences (reads with very many anchors)
    #[arg(long, hide = true)]
    pub par_chain: bool,








}

static VERBOSE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn main() -> anyhow::Result<()> {
    let args = Cli::parse();
    VERBOSE.store(args.verbose || std::env::var_os("BONEFISH_TIMERS").is_some(), std::sync::atomic::Ordering::Relaxed);

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
fn parse_num(opt: &str, s: &str) -> anyhow::Result<i64> {
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
    let x: f64 = num_str.parse().map_err(|_| anyhow::anyhow!("{opt}: '{s}' is not a number (suffix K, M or G allowed)"))?;
    let v = (x * mult + 0.499) as i64;
    if v <= 0 { anyhow::bail!("{opt}: '{s}' must be positive"); }
    Ok(v)
}

/// Parse a comma-separated list of 1..=max numbers for option `opt`.
fn parse_list<T: std::str::FromStr>(opt: &str, s: &str, max: usize) -> anyhow::Result<Vec<T>> {
    let v: Vec<T> = s.split(',').map(|x| x.trim().parse::<T>())
        .collect::<Result<_, _>>().map_err(|_| anyhow::anyhow!("{opt}: cannot parse '{s}'"))?;
    if v.is_empty() || v.len() > max {
        anyhow::bail!("{opt}: expected {} number{}, got '{s}'", if max == 1 { "1".to_string() } else { format!("1 to {max}") }, if max == 1 { "" } else { "s" });
    }
    Ok(v)
}

/// Apply preset configurations.
fn apply_preset(opt: &mut MapOptions, k: &mut usize, w: &mut usize, preset: &str) -> anyhow::Result<()> {
    bonefish::api::apply_preset_str(opt, k, w, preset).map_err(|e| anyhow::anyhow!("{}", e))
}

fn run(mut cli: AlignArgs) -> anyhow::Result<()> {
    // BAM output is always the coordinate-sorted BAM; it implies SAM records.
    if cli.output.as_deref().is_some_and(|p| p.ends_with(".bam")) { cli.sort = true; }
    if cli.sort { cli.output_sam = true; }
    if cli.write_junc && cli.output_sam { anyhow::bail!("--write-junc writes BED, not SAM/BAM (-a, --sort, -o FILE.bam)"); }
    if let Some(k) = cli.k && !(1..=28).contains(&k) { anyhow::bail!("-k must be between 1 and 28 (got {k})"); }
    if let Some(w) = cli.w && !(1..=255).contains(&w) { anyhow::bail!("-w must be between 1 and 255 (got {w})"); }
    if let Some(p) = cli.pri_ratio && !(0.0..=1.0).contains(&p) { anyhow::bail!("-p must be between 0 and 1 (got {p})"); }
    if let Some(rg) = &cli.sam_rg {
        let line = rg.replace("\\t", "\t");
        if !line.starts_with("@RG\t") || !line.split('\t').any(|t| t.len() > 3 && t.starts_with("ID:")) {
            anyhow::bail!("-R must be a read group line with an ID, like '@RG\\tID:sample1\\tSM:sample1' (got '{rg}')");
        }
    }
    if let Some(d) = &cli.dump_index && d.ends_with(".mmi") {
        eprintln!("Warning: -d writes bonefish's own index format, which minimap2 cannot read; consider a .idx name");
    }
    let do_ds = cli.ds;
    let cs_mode: Option<&str> = cli.cs.as_ref().map(|v| v.as_deref().unwrap_or("short"));
    if let Some(m) = cs_mode && m != "short" && m != "long" && m != "none" {
        anyhow::bail!("--cs takes 'short', 'long' or 'none' (got '{m}')");
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
    bonefish::fasta::bam::set_copy_tags(copy_comment);
    // -R replaces the read groups of BAM input.
    bonefish::fasta::bam::set_drop_rg(cli.sam_rg.is_some());

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
    let mut opt = MapOptions::default();

    apply_preset(&mut opt, &mut k, &mut w, cli.preset.as_deref().unwrap_or("dna"))?;

    // CLI overrides (applied after preset)
    if let Some(vk) = cli.k { k = vk; }
    if let Some(vw) = cli.w { w = vw; }

    let is_splice = opt.flags.contains(AlignFlags::SPLICE);
    if !is_splice && (cli.max_intron_len.is_some() || cli.noncan.is_some() || cli.splice_strand.is_some()
        || cli.junc_bed.is_some() || cli.junc_bonus.is_some() || cli.junc_pen.is_some()
        || cli.splice_flank.is_some() || cli.write_junc) {
        anyhow::bail!("-G, -C, -u, --junc-bed, --junc-bonus, --junc-pen, --splice-flank and --write-junc need a spliced preset (-x rna or -x cdna)");
    }
    if let Some(s) = &cli.max_intron_len {
        let v = parse_num("-G", s)?.min(i32::MAX as i64) as i32;
        opt.chaining.max_gap_ref = v;
        opt.chaining.bandwidth = v;
        opt.chaining.bandwidth_long = v;
    }
    if let Some(v) = cli.noncan { opt.scoring.noncanon_penalty = v; }
    if let Some(v) = cli.junc_bonus { opt.scoring.junc_bonus = v; }
    if let Some(v) = cli.junc_pen { opt.scoring.junc_pen = v; }
    if let Some(u) = &cli.splice_strand {
        opt.flags.remove(AlignFlags::SPLICE_FOR | AlignFlags::SPLICE_REV);
        match u.as_str() {
            "f" => opt.flags.insert(AlignFlags::SPLICE_FOR),
            "r" => opt.flags.insert(AlignFlags::SPLICE_REV),
            "b" => opt.flags.insert(AlignFlags::SPLICE_FOR | AlignFlags::SPLICE_REV),
            "n" => {}
            _ => anyhow::bail!("-u takes f, r, b or n (got '{u}')"),
        }
    }
    if let Some(f) = &cli.splice_flank {
        match f.as_str() {
            "yes" => opt.flags.insert(AlignFlags::SPLICE_FLANK),
            "no" => opt.flags.remove(AlignFlags::SPLICE_FLANK),
            _ => anyhow::bail!("--splice-flank takes yes or no (got '{f}')"),
        }
    }
    if cli.write_junc { opt.flags.insert(AlignFlags::OUT_JUNC); }

    if let Some(r_str) = &cli.bw {
        let v = parse_list("-r", r_str, 2)?;
        opt.chaining.bandwidth = v[0];
        if v.len() > 1 { opt.chaining.bandwidth_long = v[1]; }
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
        let v = parse_list("-O", s, 2)?;
        opt.scoring.gap_open = v[0];
        if v.len() > 1 { opt.scoring.gap_open2 = v[1]; }
    }
    if let Some(s) = &cli.scor_e {
        let v = parse_list("-E", s, 2)?;
        opt.scoring.gap_extend = v[0];
        if v.len() > 1 { opt.scoring.gap_extend2 = v[1]; }
    }
    if let Some(v) = cli.best_n { opt.filtering.best_n = v; }
    if let Some(v) = cli.pri_ratio { opt.filtering.pri_ratio = v; }
    if let Some(s) = &cli.z_drop {
        let v = parse_list::<i32>("-z", s, 2)?;
        opt.alignment.zdrop = v[0];
        opt.alignment.zdrop_inv = if v.len() > 1 { v[1] } else { v[0] }; // zdrop_inv = zdrop by default
    }
    if let Some(s) = &cli.mid_occ_range {
        let v = parse_list("-U", s, 2)?;
        opt.seeding.min_mid_occ = v[0];
        if v.len() > 1 { opt.seeding.max_mid_occ = v[1]; }
    }
    if cli.for_only { opt.flags.insert(AlignFlags::FOR_ONLY); }
    if cli.rev_only { opt.flags.insert(AlignFlags::REV_ONLY); }
    if cli.no_long_join { opt.flags.insert(AlignFlags::NO_LJOIN); }
    if cli.all_chains { opt.flags.insert(AlignFlags::ALL_CHAINS); }
    if let Some(v) = cli.min_dp_max { opt.alignment.min_dp_max = v; }
    if let Some(v) = cli.transition { opt.scoring.transition = v; }
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
    if let Some(d) = cli.rescue_max_density { opt.chaining.rescue_max_density = d; }
    if cli.gap_subanchor { opt.alignment.gap_subanchor = true; }
    if cli.no_gap_subanchor { opt.alignment.gap_subanchor = false; }
    if let Some(d) = cli.gap_subanchor_tr { opt.alignment.gap_subanchor_tr = d; }
    if let Some(sv) = &cli.query_subsample {
        let v: Vec<u32> = sv.split(',').filter_map(|x| x.trim().parse().ok()).collect();
        match v.as_slice() {
            [f] if *f >= 1 => opt.seeding.query_subsample = Some((*f, 8000)),
            [f, m] if *f >= 1 => opt.seeding.query_subsample = Some((*f, *m)),
            _ => anyhow::bail!("--query-subsample expects F[,MIN_LEN] (got '{sv}')"),
        }
    }
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
            _ => anyhow::bail!("--ext-probe expects LEN,MIN_SCORE_PER_KB[,LOW_BW] (got '{s}')"),
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
    if cli.no_hash_name { opt.flags.insert(AlignFlags::NO_HASH_NAME); }
    if cli.qstrand { opt.flags.insert(AlignFlags::QSTRAND | AlignFlags::NO_INV); }
    if cli.par_chain { opt.flags.insert(AlignFlags::PAR_CHAIN); }
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
            _ => anyhow::bail!("--secondary must be 'yes' or 'no' (got '{s}')"),
        }
    }
    if let Some(s) = &cli.rmq {
        match s.as_str() {
            "yes" => { opt.flags.insert(AlignFlags::RMQ_CHAIN); }
            "no" => { opt.flags.remove(AlignFlags::RMQ_CHAIN); }
            _ => anyhow::bail!("--rmq must be 'yes' or 'no' (got '{s}')"),
        }
    }
    if let Some(s) = &cli.dual {
        match s.as_str() {
            "yes" => { opt.flags.remove(AlignFlags::NO_DUAL); }
            "no" => { opt.flags.insert(AlignFlags::NO_DUAL); }
            _ => anyhow::bail!("--dual must be 'yes' or 'no' (got '{s}')"),
        }
    }
    if let Some(v) = cli.chain_skip_scale {
        opt.filtering.chain_skip_scale = v;
    }
    bonefish::align::dp::set_clean_dual(opt.alignment.clean_dp);
    // Recompute chaining penalties after all CLI overrides
    bonefish::api::finalize_options(&mut opt, k);

    let target_path = &cli.target;
    let batch_size: u64 = match cli.batch_size.as_ref() {
        Some(s) => parse_num("-I", s)? as u64,
        None => 8_000_000_000, // 8G default
    };

    // -K: overrides opt.mini_batch_size (default 100M)
    if let Some(s) = cli.mini_batch_size.as_ref() {
        opt.mini_batch_size = parse_num("-K", s)?;
    }

    // Parse query files from positional args
    let no_query = cli.queries.is_empty() && cli.dump_index.is_some();
    let query_files: Vec<String> = if cli.queries.is_empty() {
        // -d with no query files → index-only mode (build/dump and exit); otherwise stdin.
        if cli.dump_index.is_some() { vec![] } else { vec!["-".to_string()] }
    } else {
        // One or more query files, mapped one after another.
        cli.queries.clone()
    };

    // Validate query file existence
    for qf in &query_files {
        if qf != "-" && !std::path::Path::new(qf).exists() {
            anyhow::bail!("Input file '{}' not found.", qf);
        }
    }
    // Check the reference and the output combination before -o creates (truncates) its file.
    std::fs::File::open(target_path).map_err(|e| anyhow::anyhow!("cannot open reference '{}': {}", target_path, e))?;
    if cli.sort {
        if !cli.output.as_deref().is_some_and(|p| p.ends_with(".bam")) { anyhow::bail!("--sort needs -o FILE.bam"); }
        if let Some(m) = cli.sort_mem.as_deref() { parse_num("--sort-mem", m)?; }
    }

    // ─── Output setup ───
    let t_start = Instant::now();
    let output_file: Option<std::fs::File> = if let Some(path) = &cli.output {
        Some(std::fs::File::create(path)
            .map_err(|e| anyhow::anyhow!("Failed to create output file '{}': {}", path, e))?)
    } else {
        None
    };
    let stdout = std::io::stdout();
    let sorted_out = if cli.sort {
        let path = cli.output.as_deref().filter(|p| p.ends_with(".bam"))
            .ok_or_else(|| anyhow::anyhow!("--sort needs -o FILE.bam"))?;
        let cfg = bonefish::sorted_bam::SortConfig {
            mem: match cli.sort_mem.as_deref() { Some(m) => parse_num("--sort-mem", m)? as usize, None => 4 << 30 },
            tmp_prefix: cli.sort_tmp.clone().unwrap_or_else(|| format!("{}.tmp", path)),
            threads: cli.align_threads.max(1),
            level: cli.sort_level,
        };
        Some(std::sync::Arc::new(std::sync::Mutex::new(Some(bonefish::sorted_bam::SortedBamWriter::new(path, cfg)?))))
    } else {
        None
    };
    let mut handle: Box<dyn Write + Send> = if let Some(w) = &sorted_out {
        drop(output_file);
        Box::new(bonefish::sorted_bam::SharedSink(w.clone()))
    } else if let Some(f) = output_file {
        Box::new(std::io::BufWriter::new(f))
    } else {
        Box::new(std::io::BufWriter::new(stdout))
    };

    let mut total_stats = AlignmentStats::default();
    let mid_occ_frac = cli.mid_occ_frac.unwrap_or(2e-4);
    let mut input_rg: Vec<String> = Vec::new();
    if cli.sam_rg.is_none() {
        for qf in query_files.iter() {
            for l in bonefish::fasta::input_read_groups(qf).unwrap_or_default() {
                if !input_rg.contains(&l) { input_rg.push(l); }
            }
        }
    }
    let run_cfg = RunConfig {
        mid_occ_frac,
        input_rg,
        copy_comment,
        cli: &cli,
        out_cfg: &out_cfg,
        split_out_cfg: &split_out_cfg,
    };
    // ─── Build index parts iterator ───
    // Each iteration yields one Index part. For small references or large batch_size,
    // this will be a single iteration (matching current behavior exactly).
    // Extension is just a naming convention and not always reliable (e.g. a real
    // minimap2/bonefish index saved under a non-standard name) — fall back to sniffing
    // the file's magic bytes so a mis-named index isn't silently misread as sequence
    // input to build a new (bogus) index from.
    let is_idx_file = target_path.ends_with(".mmi") || target_path.ends_with(".idx") || target_path.ends_with(".rmmi")
        || Index::sniff_is_index_file(target_path);
    let mut n_parts = 0usize;
    let mut sam_header_written = false;
    let mut split_output = SplitOutput { prefix: cli.split_prefix.clone(), temp_dir: None };
    split_output.prepare(false, &query_files)?;
    let mut saved_k = k;
    let mut saved_w = w;

    if is_idx_file {
        // ─── Load from .mmi/.idx: iterate over parts ───
        eprintln!("[*] Loading index from {}...", target_path);
        let t_load = Instant::now();
        // Raw indexes are memory-mapped (all parts at once; pages fault in on use);
        // other formats are deserialized one part at a time.
        #[cfg(not(target_arch = "wasm32"))]
        let mut mapped_parts: Option<std::collections::VecDeque<Index>> = if Index::sniff_is_raw_index(target_path) {
            Some(Index::load_raw_parts_with(target_path, cli.hugepages)
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
            if n_parts == 1 {
                let has_more = match &mapped_parts {
                    Some(parts) => !parts.is_empty(),
                    None => !std::io::BufRead::fill_buf(&mut idx_reader)?.is_empty(),
                };
                split_output.prepare(cli.output_sam && has_more, &query_files)?;
            }
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
            eprintln!("[*] Loaded part {} in {:.3}s — k: {}; w: {}; #seq: {}; total len: {}",
                n_parts, t_load.elapsed().as_secs_f64(),
                mi.kmer_size, mi.window_size,
                mi.seqs.len(), total_len);
            if n_parts == 1 && mapped_parts.is_none() && cli.dump_index.is_none() {
                eprintln!("[*] Tip: `bonefish -d INDEX.raw --idx-raw {}` converts this index to the raw format, which loads near-instantly", target_path);
            }
            if n_parts == 1 && do_cigar && !mi.has_sequences() {
                anyhow::bail!("The prebuilt index doesn't contain sequences. Cannot produce CIGAR/cs/MD output.");
            }
            saved_k = mi.kmer_size;
            saved_w = mi.window_size;
            bonefish::api::finalize_options(&mut opt, mi.kmer_size);
            if n_parts == 1 && (cli.k.is_some_and(|v| v != mi.kmer_size) || cli.w.is_some_and(|v| v != mi.window_size)) {
                eprintln!("Warning: -k/-w are ignored with a prebuilt index (it has k={}, w={})", mi.kmer_size, mi.window_size);
            }
            if let Some(ref alt_path) = cli.alt {
                read_alt_list(&mut mi, alt_path)?;
            }

            if let Some(ref prefix) = split_output.prefix {
                // Split mode: write results to temp files
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
                let all_seqs = bonefish::fasta::read_fasta(target_path)
                    .map_err(|e| anyhow::anyhow!("cannot read reference '{}': {}", target_path, e))?;

                // Check if we need multiple parts
                let total_bases: u64 = all_seqs.iter().map(|(_, s)| s.len() as u64).sum();

                if total_bases <= batch_size {
                    // ─── Single-part fast path (identical to previous behavior) ───
                    n_parts = 1;
                    let t0 = Instant::now();
                    eprintln!("[*] Building index for {} sequences...", all_seqs.len());
                    let mut mi = Index::build(all_seqs, w, k, usize::MAX);
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

                    if let Some(ref prefix) = split_output.prefix {
                        let mut split_writer = split::split_init(prefix, 0, &mi)?;
                        for qf in &query_files {
                            total_stats = total_stats + map_one_part_split(
                                &mi, &mut opt, &run_cfg, qf, 1, &mut split_writer,
                            )?;
                        }
                    } else {
                        for qf in &query_files {
                            let part_stats = map_one_part(
                                &mi, &mut opt, &run_cfg,
                                qf, 1, &mut sam_header_written,
                                &mut handle,
                            )?;
                            total_stats = total_stats + part_stats;
                        }
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
                        split_output.prepare(cli.output_sam, &query_files)?;
                        let seqs = std::mem::take(batch_seqs);
                        *batch_bases = 0;
                        *n_parts += 1;
                        let t0 = Instant::now();
                        eprintln!("[*] Building index part {} for {} sequences...", *n_parts, seqs.len());
                        let mut mi = Index::build(seqs, w, k, usize::MAX);
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

                        if let Some(ref prefix) = split_output.prefix {
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
            let mut reader = bonefish::fasta::open(target_path)
                .map_err(|e| anyhow::anyhow!("cannot read reference '{}': {}", target_path, e))?;
            let mut dump_writer: Option<std::io::BufWriter<std::fs::File>> = if let Some(dump_path) = &cli.dump_index {
                Some(std::io::BufWriter::new(std::fs::File::create(dump_path)
                    .map_err(|e| anyhow::anyhow!("Failed to create index '{}': {}", dump_path, e))?))
            } else {
                None
            };

            loop {
                let (batch_seqs, is_eof) = reader.read_batch(batch_size)
                    .map_err(|e| anyhow::anyhow!("cannot read reference '{}': {}", target_path, e))?;
                if batch_seqs.is_empty() { break; }
                split_output.prepare(cli.output_sam && !is_eof, &query_files)?;

                n_parts += 1;
                let t0 = Instant::now();
                eprintln!("[*] Building index part {} for {} sequences...", n_parts, batch_seqs.len());
                let mut mi = Index::build(batch_seqs, w, k, usize::MAX);
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

                if let Some(ref prefix) = split_output.prefix {
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

    if n_parts > 1 && split_output.prefix.is_none() {
        eprintln!("[WARNING] Index has {} parts. Without --split-prefix, MAPQ and secondary flags may be inaccurate.", n_parts);
    }

    // Merge split results if --split-prefix was used
    if let Some(ref prefix) = split_output.prefix && n_parts >= 1 {
        eprintln!("[*] Merging {} split parts...", n_parts);
        let merge_cfg = split::MergeConfig {
            prefix,
            n_parts,
            k: saved_k,
            w: saved_w,
            query_paths: &query_files,
            sam_rg: cli.sam_rg.as_deref(),
            input_rg: &run_cfg.input_rg,
            copy_comment,
        };
        let merge_stats = split::merge_split_results(
            &merge_cfg, &opt, &out_cfg, &mut *handle,
        )?;
        total_stats = total_stats + merge_stats;
    }

    handle.flush()?;
    drop(handle);
    if let Some(w) = sorted_out {
        let t_sort = Instant::now();
        let w = std::sync::Arc::try_unwrap(w).map_err(|_| anyhow::anyhow!("sorted output still in use"))?
            .into_inner().unwrap().unwrap();
        let n = w.finish()?;
        eprintln!("[*] Sorted BAM: {} records, merged and written in {:.3}s after mapping", n, t_sort.elapsed().as_secs_f64());
    }

    if no_query {
        let total_time = t_start.elapsed().as_secs_f64();
        eprintln!("[*] Done in {:.3}s ({} index part{}). No query file provided; skipping mapping.",
            total_time, n_parts, if n_parts != 1 { "s" } else { "" });
    } else {
        #[cfg(feature = "census")]
        bonefish::align::dp::census::flush();
        let total_time = t_start.elapsed().as_secs_f64();
        eprintln!("[*] Mapped {} reads in {:.3}s ({} index part{})", total_stats.n_reads, total_time, n_parts, if n_parts != 1 { "s" } else { "" });

        // Per-stage timers (thread-seconds) for profiling: -v or BONEFISH_TIMERS=1.
        if VERBOSE.load(std::sync::atomic::Ordering::Relaxed) {
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
        bonefish::align::map::print_post_chain_breakdown();
        }
    }

    Ok(())
}

/// Automatically merge multi-part SAM/BAM output so its header contains every reference.
/// The user can still choose a temporary-file location with --split-prefix.
struct SplitOutput {
    prefix: Option<String>,
    temp_dir: Option<std::path::PathBuf>,
}

impl SplitOutput {
    fn prepare(&mut self, multi_part_sam: bool, queries: &[String]) -> anyhow::Result<()> {
        if self.prefix.is_none() && !multi_part_sam { return Ok(()); }
        if queries.iter().any(|q| q == "-" || std::fs::metadata(q).is_ok_and(|m| !m.is_file())) {
            anyhow::bail!("Split-index merging requires regular query files, not stdin or pipes. Increase -I to map in one part.");
        }
        if self.prefix.is_none() {
            let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?.as_nanos();
            let dir = std::env::temp_dir().join(format!("bonefish-split-{}-{}", std::process::id(), stamp));
            std::fs::create_dir(&dir)?;
            self.prefix = Some(dir.join("part").to_string_lossy().into_owned());
            self.temp_dir = Some(dir);
        }
        Ok(())
    }
}

impl Drop for SplitOutput {
    fn drop(&mut self) {
        if let Some(dir) = &self.temp_dir { let _ = std::fs::remove_dir_all(dir); }
    }
}

/// Run-level configuration fields that stay constant across index parts.
struct RunConfig<'a> {
    mid_occ_frac: f32,
    /// `@RG` header lines of the BAM inputs (kept unless -R sets the read group).
    input_rg: Vec<String>,
    cli: &'a AlignArgs,
    copy_comment: bool,
    out_cfg: &'a OutputConfig,
    split_out_cfg: &'a OutputConfig,
}

/// Load --junc-bed against this index part (reference ids are local to the part).
fn load_junctions(opt: &mut MapOptions, mi: &Index, cli: &AlignArgs) -> anyhow::Result<()> {
    if let Some(path) = &cli.junc_bed {
        let db = bonefish::align::junc::JunctionDb::load_bed_for_index(path, mi)
            .map_err(|e| anyhow::anyhow!("--junc-bed {}: {}", path, e))?;
        opt.junctions = Some(std::sync::Arc::new(db));
    }
    Ok(())
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
    let copy_comment = run.copy_comment;
    let out_cfg = run.out_cfg;
    let mini_batch_size = opt.mini_batch_size as u64;
    // mid_occ per-part: only calculate if still <= 0.
    // Once calculated for the first part, it stays for subsequent parts (options.c:73).
    if opt.seeding.mid_occ == 0 {
        opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
    }
    load_junctions(opt, mi, cli)?;

    // Single-part SAM header; multi-part SAM/BAM headers are written by the split merge.
    if cli.output_sam && !*sam_header_written {
        writeln!(handle, "@HD\tVN:1.6\tSO:unsorted\tGO:query")?;
        for s in &mi.seqs {
            writeln!(handle, "@SQ\tSN:{}\tLN:{}", s.name, s.len)?;
        }
        let version = env!("CARGO_PKG_VERSION");
        let cmd_line = std::env::args().collect::<Vec<String>>().join(" ");
        writeln!(handle, "@PG\tID:bonefish\tPN:bonefish\tVN:{}\tCL:{}", version, cmd_line)?;
        if let Some(rg_line) = &cli.sam_rg {
            let expanded = rg_line.replace("\\t", "\t");
            writeln!(handle, "{}", expanded)?;
        }
        for l in &run.input_rg {
            writeln!(handle, "{}", l)?;
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

    let mut reader = bonefish::fasta::open(query_path)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    // Single-end reads only need the header comment for -y; skip copying it otherwise.
    if !copy_comment { reader.set_keep_description(false); }

    // Parallel branches below assign `total_stats` unconditionally from the
    // thread::scope result; only the non-parallel accumulation path needs the
    // default initial value.
    #[cfg(feature = "parallel")]
    let total_stats: AlignmentStats;
    #[cfg(not(feature = "parallel"))]
    let mut total_stats = AlignmentStats::default();
    let mut record_iter = reader.records();

    {
        // Three-stage pipeline: reader thread → worker pool → writer thread.
        // Output is decoupled from the worker pool so cores stay saturated
        // through chunk transitions.
        #[cfg(feature = "parallel")]
        {
            type SeData = (String, Vec<u8>, Option<String>, Option<String>);
            // Chunks the reader may run ahead: input decompression is steady but alignment speed
            // varies with the reads, so a deeper queue absorbs fast stretches. About 800 Mbp of
            // reads (~1.6 GB with qualities), 2 to 16 chunks; env BONEFISH_LOOKAHEAD overrides.
            let lookahead: usize = std::env::var("BONEFISH_LOOKAHEAD").ok().and_then(|v| v.parse().ok())
                .unwrap_or_else(|| (800_000_000 / mini_batch_size.max(1)).clamp(2, 16) as usize);
            let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<SeData>>(lookahead);
            let (tx_out, rx_out) = std::sync::mpsc::sync_channel::<Vec<(String, AlignmentStats)>>(2);

            // A read error stops reading; it is reported (non-zero exit) once the reads before it
            // are written.
            let read_err: std::sync::Mutex<Option<String>> = std::sync::Mutex::new(None);
            let read_err_ref = &read_err;
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
                                    let c = if copy_comment { desc.and_then(|d| bonefish::fasta::bam::filter_comment_tags(&d)) } else { None };
                                    // Copied tags count towards the batch size: ONT's MM/ML tags take more
                                    // memory than the bases themselves.
                                    chunk_bases += (seq.len() + c.as_ref().map_or(0, |c| c.len())) as u64;
                                    chunk_data.push((name, seq, q, c));
                                },
                                Some(Err(e)) => {
                                    *read_err_ref.lock().unwrap() = Some(e.to_string());
                                    break;
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

                let align_one = |ctx: &mut bonefish::align::extend::AlignmentContext,
                                 map_ctx: &mut bonefish::align::map::MapContext,
                                 (qname, qseq, qual, comment): &SeData| {
                    let ri = ReadInfo { qname, qseq, qual: qual.as_deref(), comment: comment.as_deref() };
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        bonefish::align::pipeline::align_and_format_query(opt, mi, &ri, ctx, map_ctx, out_cfg)
                    }))
                };
                // Stream reads into the pool chunk by chunk, so the stragglers at the
                // end of one chunk overlap with later chunks instead of leaving cores idle at a
                // per-chunk barrier. Two chunks are always admitted; beyond that a new chunk is
                // admitted whenever fewer than 2 x threads reads are left queued or running
                // (e.g. one long read still holding the oldest chunk), up to MAX_IN_FLIGHT.
                // Chunks are still written in order.
                use std::sync::atomic::{AtomicUsize, Ordering};
                struct Chunk {
                    data: Vec<SeData>,
                    done: Vec<std::sync::OnceLock<(String, AlignmentStats)>>,
                    remaining: AtomicUsize,
                }
                thread_local! {
                    static CTX: std::cell::RefCell<Option<(bonefish::align::extend::AlignmentContext, bonefish::align::map::MapContext)>> =
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
                                                let (ctx, map_ctx) = g.get_or_insert_with(|| (bonefish::align::extend::AlignmentContext::new(), bonefish::align::map::MapContext::new()));
                                                align_one(ctx, map_ctx, read)
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
                // A failed writer stops the coordinator from draining the input queue.
                // Disconnect it before the scope waits for a reader blocked in send().
                drop(rx);
                drop(tx_out);
                writer_handle.join().expect("writer thread panicked")
            })?;
            if let Some(e) = read_err.into_inner().unwrap() {
                anyhow::bail!("{}: {}", query_path, e);
            }
        }

        #[cfg(not(feature = "parallel"))]
        {
            loop {
                #[allow(clippy::type_complexity)]
                let mut chunk_data: Vec<(String, Vec<u8>, Option<String>, Option<String>)> = Vec::new();
                let mut chunk_bases: u64 = 0;
                loop {
                    match record_iter.next() {
                        Some(Ok(r)) => {
                            let (name, desc, seq, qual) = r.into_parts();
                            let q = qual.map(|v| String::from_utf8(v).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned()));
                            let c = if copy_comment { desc.and_then(|d| bonefish::fasta::bam::filter_comment_tags(&d)) } else { None };
                            chunk_bases += (seq.len() + c.as_ref().map_or(0, |c| c.len())) as u64;
                            chunk_data.push((name, seq, q, c));
                        },
                        Some(Err(e)) => anyhow::bail!("{}: {}", query_path, e),
                        None => break,
                    }
                    if chunk_bases >= mini_batch_size { break; }
                }

                if chunk_data.is_empty() { break; }

                let results: Vec<(String, AlignmentStats)> = {
                    let mut ctx = bonefish::align::extend::AlignmentContext::new();
                    let mut map_ctx = bonefish::align::map::MapContext::new();
                    chunk_data.iter().map(|(qname, qseq, qual, comment)| {
                        let ri = ReadInfo { qname, qseq, qual: qual.as_deref(), comment: comment.as_deref() };
                         bonefish::align::pipeline::align_and_format_query(
                            opt, mi, &ri, &mut ctx, &mut map_ctx, out_cfg,
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
    use bonefish::align::pipeline::process_query;
    use bonefish::align::map::MapContext;
    let mid_occ_frac = run.mid_occ_frac;
    let out_cfg = run.split_out_cfg;

    // mid_occ per-part: only calculate if still <= 0
    if opt.seeding.mid_occ == 0 {
        opt.seeding.mid_occ = mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
    }
    load_junctions(opt, mi, run.cli)?;

    // stdin can't be re-read for subsequent parts
    if query_path == "-" && part_number > 1 {
        return Ok(AlignmentStats::default());
    }

    let reader = bonefish::fasta::open(query_path)
        .map_err(|e| anyhow::anyhow!("{}", e))?;

    let mut total_stats = AlignmentStats::default();
    let mut record_iter = reader.records();

    {
        // Single-end split write
        let mut ctx = bonefish::align::extend::AlignmentContext::new();
        let mut map_ctx = MapContext::new();

        loop {
            let rec = match record_iter.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => anyhow::bail!("{}: {}", query_path, e),
                None => break,
            };

            let qname = rec.name().to_string();
            let qseq = rec.sequence().to_vec();

            let pq = process_query(opt, mi, &qname, &qseq, &mut ctx, &mut map_ctx, out_cfg);

            split::split_write_query(split_writer, &pq.results, pq.rep_len)
                .map_err(|e| anyhow::anyhow!("Error writing split temp: {}", e))?;

            total_stats = total_stats + pq.stats;
        }
    }

    split_writer.flush()?;
    Ok(total_stats)
}
