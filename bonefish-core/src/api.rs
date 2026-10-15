//! High-level alignment API for library consumers.
//!
//! This module provides a simple interface for aligning sequences against a
//! reference, hiding the internal pipeline details (MapContext, AlignmentContext,
//! coordinate adjustments, etc.).
//!
//! # Quick Start
//!
//! ```no_run
//! use bonefish::api::{Aligner, Preset};
//! use bonefish::Strand;
//!
//! let aligner = Aligner::from_index("reference.mmi", Preset::Dna).unwrap();
//! let results = aligner.map_seq("read1", b"ACGTACGTACGT...");
//! for m in &results.mappings {
//!     println!("{}\t{}\t{}\t{}\tMAPQ={}",
//!         m.target_name, m.target_start, m.target_end,
//!         if m.strand == Strand::Forward { "+" } else { "-" },
//!         m.mapq);
//! }
//! ```
//!
//! # Thread Safety
//!
//! `Aligner` is `Send + Sync` — wrap it in `Arc<Aligner>` to share across threads.
//! Each call to `map_seq` allocates lightweight per-call buffers internally.

use std::io;
use std::sync::Arc;
use crate::align::index::Index;
use crate::align::map::{MapOptions, MapContext, AlignFlags};
use crate::align::extend::AlignmentContext;
use crate::align::pipeline::{self, OutputConfig};

/// Alignment preset (selected via the `-x` flag).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    /// ONT genomic DNA (the only preset).
    Dna,
}

impl Preset {
    /// Convert to the string form
    pub fn as_str(&self) -> &'static str {
        match self {
            Preset::Dna => "dna",
        }
    }
}

/// Strand orientation of an alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Structured CIGAR operation (length + op type).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// Operation length.
    pub len: u32,
    /// Operation type: 0=M, 1=I, 2=D, 3=N (intron skip).
    pub op: u8,
}

impl CigarOp {
    /// Operation type as a character (M, I, D, N, S, H, =, X).
    pub fn op_char(&self) -> char {
        match self.op { 0 => 'M', 1 => 'I', 2 => 'D', 3 => 'N', 4 => 'S', 5 => 'H', 7 => '=', 8 => 'X', _ => '?' }
    }
}

/// Per-call options for CS/MD tag generation, overriding the aligner default.
#[derive(Debug, Clone, Copy, Default)]
pub struct MapOpts {
    /// Generate CS tag for this call (overrides aligner setting if Some).
    pub cs: Option<bool>,
    /// Generate MD tag for this call (overrides aligner setting if Some).
    pub md: Option<bool>,
}

/// A single alignment result.
#[derive(Debug, Clone)]
pub struct Mapping {
    /// Target sequence name (shared via Arc to reduce allocation when many
    /// alignments reference the same target).
    pub target_name: Arc<str>,
    /// Target sequence numeric index in the reference.
    pub target_id: usize,
    /// Target sequence length.
    pub target_len: usize,
    /// Aligned region start on target (0-based).
    pub target_start: usize,
    /// Aligned region end on target (exclusive).
    pub target_end: usize,
    /// Aligned region start on query (0-based).
    pub query_start: usize,
    /// Aligned region end on query (exclusive).
    pub query_end: usize,
    /// Alignment strand.
    pub strand: Strand,
    /// Mapping quality (0-60).
    pub mapq: i32,
    /// Whether this is a primary alignment.
    pub is_primary: bool,
    /// Whether this is a supplementary alignment.
    pub is_supplementary: bool,
    /// Number of matching bases.
    pub matches: usize,
    /// Alignment block length.
    pub block_len: usize,
    /// Edit distance (NM tag).
    pub edit_distance: u32,
    /// CIGAR string (if CIGAR output enabled).
    pub cigar: Option<String>,
    /// Structured CIGAR (if CIGAR output enabled). Each element is (len, op).
    pub cigar_ops: Option<Vec<CigarOp>>,
    /// CS tag string (if requested).
    pub cs: Option<String>,
    /// MD tag string (if requested).
    pub md: Option<String>,
    /// Alignment score (AS tag).
    pub score: i32,
    /// Sequence divergence (0.0 = identical).
    pub divergence: f64,
}

/// Result of aligning one read (or read pair).
#[derive(Debug, Clone)]
pub struct MapResult {
    /// All alignments for this read, ordered by score (primary first).
    pub mappings: Vec<Mapping>,
}

/// High-level aligner wrapping an index and alignment parameters.
///
/// Construct via [`Aligner::from_index`], [`Aligner::from_fasta`], or
/// [`Aligner::from_seqs`], then call [`Aligner::map_seq`] per read.
///
/// # Builder-Style Preset Methods
///
/// For minimap2-rs compatibility, preset methods are also available:
///
/// ```no_run
/// use bonefish::api::Aligner;
/// let aligner = Aligner::from_fasta("ref.fa", Aligner::map_ont()).unwrap();
/// ```
pub struct Aligner {
    index: Index,
    options: MapOptions,
    out_cfg: OutputConfig,
}

/// Largest supported k-mer size.
pub const MAX_KMER: usize = 28;

/// Recompute the [`MapOptions`] fields derived from the k-mer size of the index
/// that will actually be used. Call after a preset and any overrides have been
/// applied.
pub fn finalize_options(opt: &mut MapOptions, k: usize) {
    opt.chaining.chn_pen_gap = (opt.chaining.chain_gap_scale as f64 * 0.01 * k as f64) as f32;
    opt.chaining.chn_pen_skip = (opt.filtering.chain_skip_scale as f64 * 0.01 * k as f64) as f32;
}

/// Calibrate `seeding.mid_occ` against the index, unless it has already been
/// set explicitly (by a preset or by the caller).
pub fn calibrate_mid_occ(opt: &mut MapOptions, mi: &Index, mid_occ_frac: f32) {
    if opt.seeding.mid_occ == 0 {
        opt.seeding.mid_occ =
            mi.cal_mid_occ(mid_occ_frac, opt.seeding.min_mid_occ, opt.seeding.max_mid_occ);
    }
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, msg)
}

/// Builder for [`Aligner`], overriding preset defaults at construction time.
///
/// Only parameters consumed while the index is built or while the seed
/// occurrence threshold is calibrated need to be set here. Anything else in
/// [`MapOptions`] can be adjusted after construction via
/// [`Aligner::options_mut`].
///
/// ```no_run
/// use bonefish::api::{Aligner, Preset};
/// // Dense seeding for short, degraded cDNA reads: w = 1 indexes every k-mer,
/// // at proportionally higher index memory cost.
/// let aligner = Aligner::builder(Preset::Dna).k(12).w(1).from_fasta("cdna.fa")?;
/// # Ok::<(), std::io::Error>(())
/// ```
#[derive(Debug, Clone)]
pub struct AlignerBuilder {
    preset: Preset,
    k: Option<usize>,
    w: Option<usize>,
    index_max_occ: usize,
    mid_occ: Option<usize>,
    mid_occ_frac: f32,
    min_mid_occ: Option<i32>,
    max_mid_occ: Option<i32>,
    cigar: bool,
}

impl AlignerBuilder {
    fn new(preset: Preset) -> Self {
        AlignerBuilder {
            preset,
            k: None,
            w: None,
            index_max_occ: usize::MAX,
            mid_occ: None,
            mid_occ_frac: 2e-4,
            min_mid_occ: None,
            max_mid_occ: None,
            cigar: true,
        }
    }

    /// k-mer size, in `1..=`[`MAX_KMER`]. Also scales the chaining gap penalty.
    pub fn k(mut self, k: usize) -> Self { self.k = Some(k); self }

    /// Minimizer window size. `1` indexes every k-mer.
    pub fn w(mut self, w: usize) -> Self { self.w = Some(w); self }


    /// Occurrence cap applied while building the index (default: uncapped).
    pub fn index_max_occ(mut self, max_occ: usize) -> Self { self.index_max_occ = max_occ; self }

    /// Seed occurrence threshold. Setting this skips calibration from the index.
    pub fn mid_occ(mut self, mid_occ: usize) -> Self { self.mid_occ = Some(mid_occ); self }

    /// Fraction of the most repetitive minimizers to filter when calibrating
    /// `mid_occ` (default: 2e-4).
    pub fn mid_occ_frac(mut self, frac: f32) -> Self { self.mid_occ_frac = frac; self }

    /// Floor for the calibrated `mid_occ`.
    pub fn min_mid_occ(mut self, min_mid_occ: i32) -> Self { self.min_mid_occ = Some(min_mid_occ); self }

    /// Ceiling for the calibrated `mid_occ`.
    pub fn max_mid_occ(mut self, max_mid_occ: i32) -> Self { self.max_mid_occ = Some(max_mid_occ); self }

    /// Produce CIGAR strings (default: true).
    pub fn cigar(mut self, cigar: bool) -> Self { self.cigar = cigar; self }

    /// Build the index from a FASTA reference file.
    pub fn from_fasta(self, path: &str) -> io::Result<Aligner> {
        let (opt, k, w) = self.resolve()?;
        let seqs = crate::fasta::read_fasta(path)?;
        Ok(self.finish(opt, Index::build(seqs, w, k, self.index_max_occ)))
    }

    /// Build the index from in-memory `(name, sequence)` pairs.
    pub fn from_seqs(self, seqs: Vec<(String, Vec<u8>)>) -> io::Result<Aligner> {
        let (opt, k, w) = self.resolve()?;
        Ok(self.finish(opt, Index::build(seqs, w, k, self.index_max_occ)))
    }

    /// Load a pre-built index file (bonefish `RMMI`, minimap2 `.mmi`, or legacy
    /// format; detected from the file's magic).
    pub fn from_index(self, path: &str) -> io::Result<Aligner> {
        let index = Index::load(path)?;
        self.from_loaded_index(index)
    }

    /// Load a pre-built index from any seekable reader positioned at its start.
    /// See [`Aligner::from_index_reader`].
    pub fn from_index_reader<R: std::io::Read + std::io::Seek>(
        self,
        reader: &mut R,
    ) -> io::Result<Aligner> {
        let index = Index::load_part(reader)?
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Empty index"))?;
        self.from_loaded_index(index)
    }

    /// Apply the preset to an already-loaded index. Its seeding geometry is
    /// fixed, so any `k` or `w` set here must agree with the stored values.
    #[allow(clippy::wrong_self_convention)] // named like the public from_* finishers above
    fn from_loaded_index(self, index: Index) -> io::Result<Aligner> {
        let (opt, _, _) = self.resolve()?;
        if let Some(k) = self.k && k != index.kmer_size {
            return Err(invalid(format!("k={} does not match index k={}", k, index.kmer_size)));
        }
        if let Some(w) = self.w && w != index.window_size {
            return Err(invalid(format!("w={} does not match index w={}", w, index.window_size)));
        }
        Ok(self.finish(opt, index))
    }

    /// Preset plus overrides, before any index-derived values are filled in.
    fn resolve(&self) -> io::Result<(MapOptions, usize, usize)> {
        let mut k = 15usize;
        let mut w = 10usize;
        let mut opt = MapOptions::default();
        apply_preset_str(&mut opt, &mut k, &mut w, self.preset.as_str())
            .expect("Preset enum variants are all recognized");
        let k = self.k.unwrap_or(k);
        let w = self.w.unwrap_or(w);

        if k == 0 || k > MAX_KMER {
            return Err(invalid(format!("k must be in 1..={}, got {}", MAX_KMER, k)));
        }
        if w == 0 {
            return Err(invalid("w must be at least 1".to_string()));
        }
        if !(0.0..1.0).contains(&self.mid_occ_frac) {
            return Err(invalid(format!("mid_occ_frac must be in [0, 1), got {}", self.mid_occ_frac)));
        }

        if let Some(v) = self.mid_occ { opt.seeding.mid_occ = v; }
        if let Some(v) = self.min_mid_occ { opt.seeding.min_mid_occ = v; }
        if let Some(v) = self.max_mid_occ { opt.seeding.max_mid_occ = v; }
        if opt.seeding.min_mid_occ > opt.seeding.max_mid_occ {
            return Err(invalid(format!(
                "min_mid_occ ({}) exceeds max_mid_occ ({})",
                opt.seeding.min_mid_occ, opt.seeding.max_mid_occ)));
        }
        Ok((opt, k, w))
    }

    fn finish(&self, mut opt: MapOptions, index: Index) -> Aligner {
        finalize_options(&mut opt, index.kmer_size);
        calibrate_mid_occ(&mut opt, &index, self.mid_occ_frac);
        if self.cigar {
            opt.flags.insert(AlignFlags::OUT_CIGAR);
        } else {
            opt.flags.remove(AlignFlags::OUT_CIGAR);
        }
        let out_cfg = OutputConfig {
            do_cigar: self.cigar, do_cs: false, cs_long: false, do_md: false, do_ds: false,
            eqx: false, output_sam: false, rg_id: None, split_mode: false,
        };
        Aligner { index, options: opt, out_cfg }
    }
}

impl Aligner {
    // --- Preset convenience methods ---

    /// ONT genomic DNA preset (the only preset).
    pub fn dna() -> Preset { Preset::Dna }
    /// Alias of [`Aligner::dna`] for minimap2-rs compatibility.
    pub fn map_ont() -> Preset { Preset::Dna }

    // --- Constructors ---

    /// Start a builder for overriding the preset's seeding parameters.
    ///
    /// ```no_run
    /// use bonefish::api::{Aligner, Preset};
    /// let aligner = Aligner::builder(Preset::Dna).k(12).w(1).from_fasta("ref.fa")?;
    /// # Ok::<(), std::io::Error>(())
    /// ```
    pub fn builder(preset: Preset) -> AlignerBuilder {
        AlignerBuilder::new(preset)
    }

    /// Load an aligner from a pre-built index file (bonefish `RMMI`, minimap2
    /// `.mmi`, or legacy format; the format is detected from the file's magic).
    pub fn from_index(path: &str, preset: Preset) -> io::Result<Self> {
        Aligner::builder(preset).from_index(path)
    }

    /// Load an aligner from any seekable reader positioned at the start of a
    /// bonefish (`RMMI`), minimap2 (`.mmi`), or legacy index. Unlike
    /// [`Aligner::from_index`], the caller owns the reader and may bound it to
    /// exactly the index bytes — e.g. a [`std::io::Cursor`] over a slice whose
    /// trailing footer/metadata has been removed, or an index embedded in a
    /// larger container — so loading never depends on the index format ignoring
    /// trailing bytes.
    pub fn from_index_reader<R: std::io::Read + std::io::Seek>(
        reader: &mut R,
        preset: Preset,
    ) -> io::Result<Self> {
        Aligner::builder(preset).from_index_reader(reader)
    }

    /// Build an aligner from a FASTA reference file.
    pub fn from_fasta(path: &str, preset: Preset) -> io::Result<Self> {
        Aligner::builder(preset).from_fasta(path)
    }

    /// Build an aligner from in-memory sequences (name, sequence pairs).
    ///
    /// ```no_run
    /// use bonefish::api::{Aligner, Preset};
    /// let seqs = vec![
    ///     ("chr1".to_string(), b"ACGTACGT...".to_vec()),
    ///     ("chr2".to_string(), b"TGCATGCA...".to_vec()),
    /// ];
    /// let aligner = Aligner::from_seqs(seqs, Preset::Dna);
    /// ```
    pub fn from_seqs(seqs: Vec<(String, Vec<u8>)>, preset: Preset) -> Self {
        Aligner::builder(preset)
            .from_seqs(seqs)
            .expect("Preset defaults are always valid")
    }

    // --- Mapping ---

    /// Align a single-end read against the reference.
    pub fn map_seq(&self, name: &str, seq: &[u8]) -> MapResult {
        self.map_seq_with(name, seq, MapOpts::default())
    }

    /// Align a single-end read with per-call options (CS/MD toggles).
    pub fn map_seq_with(&self, name: &str, seq: &[u8], opts: MapOpts) -> MapResult {
        let mut ctx = AlignmentContext::new();
        let mut map_ctx = MapContext::new();
        let out = self.resolve_out_cfg(&opts);
        let pq = pipeline::process_query(
            &self.options, &self.index, name, seq,
            &mut ctx, &mut map_ctx, &out,
        );
        to_map_result(&pq, &self.index, &out)
    }

    // --- Output ---

    /// Format a pre-computed MapResult as a PAF string.
    pub fn format_paf(&self, name: &str, query_len: usize, result: &MapResult) -> String {
        let mut lines = Vec::new();
        for aln in &result.mappings {
            let strand = if aln.strand == Strand::Forward { '+' } else { '-' };
            let mut line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                name, query_len, aln.query_start, aln.query_end,
                strand, aln.target_name, aln.target_len,
                aln.target_start, aln.target_end,
                aln.matches, aln.block_len, aln.mapq,
            );
            if let Some(ref cigar) = aln.cigar {
                line.push_str(&format!("\tcg:Z:{}", cigar));
            }
            lines.push(line);
        }
        lines.join("\n")
    }

    // --- Index management ---

    /// Save the index to a file (RMMI format).
    pub fn save_index(&self, path: &str) -> io::Result<()> {
        self.index.save(path).map_err(io::Error::other)
    }

    /// Access the underlying index.
    pub fn index(&self) -> &Index { &self.index }

    /// Access the current MapOptions.
    pub fn options(&self) -> &MapOptions { &self.options }

    /// Mutably access MapOptions for fine-tuning after construction.
    pub fn options_mut(&mut self) -> &mut MapOptions { &mut self.options }

    /// Mutably access OutputConfig for toggling CIGAR/CS/MD/SAM defaults.
    pub fn output_config_mut(&mut self) -> &mut OutputConfig { &mut self.out_cfg }

    // --- Internal ---

    fn resolve_out_cfg(&self, opts: &MapOpts) -> OutputConfig {
        let mut out = self.out_cfg.clone();
        if let Some(cs) = opts.cs { out.do_cs = cs; }
        if let Some(md) = opts.md { out.do_md = md; }
        out
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn to_map_result(
    pq: &pipeline::ProcessedQuery,
    mi: &Index,
    out: &OutputConfig,
) -> MapResult {
    // Cache Arc<str> per target to avoid cloning the name for each alignment
    let mut name_cache: Vec<Option<Arc<str>>> = vec![None; mi.seqs.len()];

    let mappings = pq.results.iter().zip(pq.mapqs.iter()).map(|(r, &mapq)| {
        let cigar_str = if out.do_cigar && !r.cigar_str.is_empty() { Some(r.cigar_str.clone()) } else { None };
        let cigar_ops = cigar_str.as_ref().map(|s| parse_cigar_string(s));
        let target_name = name_cache[r.ref_id].get_or_insert_with(|| {
            Arc::from(mi.seqs[r.ref_id].name.as_str())
        }).clone();
        Mapping {
            target_name,
            target_id: r.ref_id,
            target_len: mi.seqs[r.ref_id].len,
            target_start: r.ref_start,
            target_end: r.ref_end,
            query_start: r.query_start,
            query_end: r.query_end,
            strand: if r.is_reverse { Strand::Reverse } else { Strand::Forward },
            mapq,
            is_primary: !r.is_secondary,
            is_supplementary: r.split != 0,
            matches: r.matches,
            block_len: r.block_len,
            edit_distance: r.edit_distance,
            cigar: cigar_str,
            cigar_ops,
            cs: if out.do_cs && !r.cs_str.is_empty() { Some(r.cs_str.clone()) } else { None },
            md: if out.do_md && !r.md_str.is_empty() { Some(r.md_str.clone()) } else { None },
            score: r.align_score,
            divergence: r.divergence,
        }
    }).collect();
    MapResult { mappings }
}

/// Parse a CIGAR string like "10M2I5M3D8M" into structured CigarOps.
fn parse_cigar_string(s: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num = 0u32;
    for c in s.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as u32 - '0' as u32);
        } else {
            let op = match c {
                'M' => 0, 'I' => 1, 'D' => 2, 'N' => 3, 'S' => 4, 'H' => 5, '=' => 7, 'X' => 8,
                _ => continue,
            };
            if num > 0 { ops.push(CigarOp { len: num, op }); }
            num = 0;
        }
    }
    ops
}

// ---------------------------------------------------------------------------
// Preset application (moved from main.rs)
// ---------------------------------------------------------------------------

/// Apply a preset by string name ("dna", or one of its minimap2-name aliases).
///
/// This is the same function used by the CLI's `-x` flag. Modifies `opt`,
/// `k` and `w` in place. Returns an error if the preset name
/// is not recognized so the caller can abort rather than silently using
/// the default configuration.
pub fn apply_preset_str(opt: &mut MapOptions, k: &mut usize, w: &mut usize, preset: &str) -> Result<(), String> {
    // Reset SIMD cap each time so switching presets has predictable behavior.
    // Specific presets below opt back into a cap if their workload is hurt by
    // AVX-512's per-core license throttling (chain/seed-heavy presets where DP
    // is a small fraction of CPU time).
    crate::align::dp::set_simd_cap(crate::align::dp::SimdCap::Auto);
    match preset {
        // The one preset: ONT genomic DNA (formerly lr:hq-fast). minimap2's lr:hq seeding (k=19,
        // w=19) and scoring, with cheaper DP and chaining (ont-aligner Phase 3 onwards): narrower
        // extension band (with a wide-band retry when an extension's tail shows a wider indel), a
        // 2 kb probe that stops extending into low-identity read ends, an adaptive 64-cell
        // extension band that follows the alignment, adaptive gap-fill band with a narrow first
        // pass (widened where read k-mers match only off it), ungapped-first gap fill, 5 kb
        // reference gap limit in the main chaining pass (the long-join rescue keeps -g, and skips
        // repeat-dense reads), segment-tree RMQ, the clean DP kernel, and long gaps sub-anchored
        // on unique 15-mers outside tandem repeats of period up to 150 bp (inside them a segment can
        // land on another repeat copy; v0.1.3 sub-anchored there too, v0.1.4 not at all).
        // Validated against GIAB truth sets on R10.4.1 (sup, HAC) and R9.4.1:
        // SNP/indel accuracy at minimap2 level, SV recall at or above it. The minimap2 names
        // people habitually pass are accepted for the same preset.
        "dna" | "lr:hq" | "lr:hq-fast" | "map-ont" | "lr" => {
            *k = 19; *w = 19;
            opt.chaining.max_gap = 10000;
            opt.chaining.max_gap_ref = 5000;
            opt.chaining.seg_rmq = true;
            opt.chaining.rescue_join = true;
            opt.chaining.rescue_max_density = 2;
            opt.seeding.min_mid_occ = 50; opt.seeding.max_mid_occ = 500;
            opt.alignment.ext_bw = Some(200);
            opt.alignment.ext_rescue = true;
            opt.alignment.ext_probe = Some((2000, 1000, 0));
            opt.alignment.gap_bw_margin = Some(500);
            opt.alignment.gap_bw_narrow = Some(32);
            opt.alignment.ext_adaptive = Some(64);
            opt.alignment.clean_dp = true;
            opt.alignment.ungapped_fill = true;
            opt.alignment.gap_subanchor = true;
            opt.alignment.gap_subanchor_tr = 150;
        },
        _ => {
            return Err(format!("unrecognized preset '{}' (ONT genomic DNA only: dna, the default)", preset));
        }
    }
    finalize_options(opt, *k);
    Ok(())
}

// ---------------------------------------------------------------------------
// Low-level pairwise DP alignment API
// ---------------------------------------------------------------------------

/// Scoring parameters for pairwise DP alignment.
#[derive(Debug, Clone)]
pub struct DpScoring {
    /// Match score (positive, e.g., 2).
    pub match_score: i32,
    /// Mismatch penalty (positive, applied as negative, e.g., 4).
    pub mismatch: i32,
    /// Gap open penalty (positive, e.g., 4).
    pub gap_open: i32,
    /// Gap extend penalty (positive, e.g., 2).
    pub gap_extend: i32,
    /// Second gap open penalty for dual-affine model (0 to disable).
    pub gap_open2: i32,
    /// Second gap extend penalty for dual-affine model (0 to disable).
    pub gap_extend2: i32,
}

impl Default for DpScoring {
    fn default() -> Self {
        DpScoring { match_score: 2, mismatch: 4, gap_open: 4, gap_extend: 2, gap_open2: 0, gap_extend2: 0 }
    }
}

/// Result of a pairwise DP alignment.
#[derive(Debug, Clone)]
pub struct DpAlignment {
    /// Alignment score.
    pub score: i32,
    /// CIGAR string (M/I/D operations).
    pub cigar: String,
    /// Query start position (0-based). For extension mode, may be > 0.
    pub query_start: usize,
    /// Query end position (exclusive).
    pub query_end: usize,
    /// Target start position (0-based).
    pub target_start: usize,
    /// Target end position (exclusive).
    pub target_end: usize,
}

/// Align two nt4-encoded sequences using the SIMD-optimized DP engine.
///
/// Performs semi-global extension alignment used to fill the gap between two
/// adjacent chained anchors: aligns `query` against `target` with affine gap
/// penalties, returning the best-scoring alignment with CIGAR traceback.
///
/// Sequences must be nt4-encoded (0=A, 1=C, 2=G, 3=T, 4=N). Use
/// [`encode_nt4`] to convert from ASCII.
///
/// # Parameters
/// - `query`: nt4-encoded query sequence
/// - `target`: nt4-encoded target sequence
/// - `scoring`: gap penalties and match/mismatch scores
/// - `bandwidth`: band width for banded alignment (-1 for full matrix)
///
/// # Example
///
/// ```
/// use bonefish::api::{dp_align, DpScoring, encode_nt4};
///
/// let query = encode_nt4(b"ACGTACGTACGT");
/// let target = encode_nt4(b"ACGTACCTACGT");
/// let result = dp_align(&query, &target, &DpScoring::default(), -1);
/// println!("score={} cigar={}", result.score, result.cigar);
/// ```
pub fn dp_align(query: &[u8], target: &[u8], scoring: &DpScoring, bandwidth: i32) -> DpAlignment {
    use crate::align::dp;
    use crate::align::extend::build_scoring_matrix;

    let mat = build_scoring_matrix(scoring.match_score, scoring.mismatch);
    let mut ez = dp::DpResult::default();

    if scoring.gap_open2 > 0 || scoring.gap_extend2 > 0 {
        dp::extend_dual_affine(
            query, target, 5, &mat,
            scoring.gap_open as i8, scoring.gap_extend as i8,
            scoring.gap_open2 as i8, scoring.gap_extend2 as i8,
            bandwidth, -1, 0, dp::APPROX_MAX, &mut ez,
        );
    } else {
        dp::extend_single_affine(
            query, target, 5, &mat,
            scoring.gap_open as i8, scoring.gap_extend as i8,
            bandwidth, -1, 0, dp::APPROX_MAX, &mut ez,
        );
    }

    let cigar_raw = &ez.cigar;
    let cigar = raw_cigar_to_string(cigar_raw);
    let (qs, qe, ts, te) = cigar_bounds(cigar_raw);
    // Recompute score from CIGAR for consistency across SIMD variants
    let score = compute_score_from_cigar(
        cigar_raw, query, target, &mat,
        scoring.gap_open, scoring.gap_extend, scoring.gap_open2, scoring.gap_extend2,
    );
    DpAlignment {
        score,
        cigar,
        query_start: qs,
        query_end: qe,
        target_start: ts,
        target_end: te,
    }
}

/// Global alignment (Needleman-Wunsch): align full query to full target.
///
/// True banded NW using the Suzuki-Kasahara formulation
/// Produces a CIGAR spanning both sequences end-to-end with
/// proper NW boundary conditions (gap penalties on first row/column).
///
/// Uses single-affine gap model (gap_open2/gap_extend2 are ignored).
///
/// # Parameters
/// - `bandwidth`: band width for banded alignment (-1 for full matrix)
pub fn dp_global(query: &[u8], target: &[u8], scoring: &DpScoring, bandwidth: i32) -> DpAlignment {
    use crate::align::dp;
    use crate::align::extend::build_scoring_matrix;

    let mat = build_scoring_matrix(scoring.match_score, scoring.mismatch);
    let mut ez = dp::DpResult::default();
    dp::global_align(
        query, target, 5, &mat,
        scoring.gap_open, scoring.gap_extend,
        bandwidth, &mut ez,
    );

    let cigar = raw_cigar_to_string(&ez.cigar);
    let (qs, qe, ts, te) = cigar_bounds(&ez.cigar);
    DpAlignment { score: ez.score, cigar, query_start: qs, query_end: qe, target_start: ts, target_end: te }
}

/// Local alignment (Smith-Waterman): find the best-scoring local region.
///
/// Two-phase approach:
/// 1. Forward Smith-Waterman to find the optimal endpoint (score, qe, te)
/// 2. Reverse extension from endpoint to find start position + CIGAR
///
/// The reverse extension uses `EXTENSION_ONLY | REV_CIGAR` on reversed
/// prefixes — it extends backward from the endpoint, stopping at the
/// optimal start, and outputs the CIGAR in forward order. Both phases
/// use SIMD-accelerated kernels.
pub fn dp_local(query: &[u8], target: &[u8], scoring: &DpScoring) -> DpAlignment {
    use crate::align::dp;
    use crate::align::extend::build_scoring_matrix;

    let mat = build_scoring_matrix(scoring.match_score, scoring.mismatch);

    // Phase 1: Forward Smith-Waterman to find the best endpoint
    let mut qp = dp::lightweight_profile_init(query.len() as i32, query, 5, &mat);
    let (score, q_end, t_end) = dp::lightweight_align_i16(
        &mut qp, target.len() as i32, target, scoring.gap_open, scoring.gap_extend,
    );
    if score <= 0 || q_end < 0 || t_end < 0 {
        return DpAlignment { score: 0, cigar: String::new(), query_start: 0, query_end: 0, target_start: 0, target_end: 0 };
    }

    let qe = q_end as usize + 1;
    let te = t_end as usize + 1;

    // Phase 2: Reverse extension → start position + CIGAR in one pass
    let q_rev: Vec<u8> = query[..qe].iter().rev().copied().collect();
    let t_rev: Vec<u8> = target[..te].iter().rev().copied().collect();
    let mut ez = dp::DpResult::default();
    dp::extend_single_affine(
        &q_rev, &t_rev, 5, &mat,
        scoring.gap_open as i8, scoring.gap_extend as i8,
        -1, score, 0,
        dp::EXTENSION_ONLY | dp::REV_CIGAR,
        &mut ez,
    );

    let cigar = raw_cigar_to_string(&ez.cigar);
    let (_, cig_qlen, _, cig_tlen) = cigar_bounds(&ez.cigar);
    let qs = qe - cig_qlen;
    let ts = te - cig_tlen;

    DpAlignment {
        score,
        cigar,
        query_start: qs,
        query_end: qe,
        target_start: ts,
        target_end: te,
    }
}

/// Extension alignment: align from position 0 outward, stopping at the best score.
///
/// Useful for extending a seed match left or right. Returns the optimal-scoring
/// prefix alignment (may not reach the end of either sequence).
pub fn dp_extension(query: &[u8], target: &[u8], scoring: &DpScoring, bandwidth: i32) -> DpAlignment {
    use crate::align::dp::{self, EXTENSION_ONLY};
    use crate::align::extend::build_scoring_matrix;

    let mat = build_scoring_matrix(scoring.match_score, scoring.mismatch);
    let mut ez = dp::DpResult::default();

    let zdrop = 400; // default z-drop threshold for extension
    if scoring.gap_open2 > 0 || scoring.gap_extend2 > 0 {
        dp::extend_dual_affine(
            query, target, 5, &mat,
            scoring.gap_open as i8, scoring.gap_extend as i8,
            scoring.gap_open2 as i8, scoring.gap_extend2 as i8,
            bandwidth, zdrop, 0, EXTENSION_ONLY, &mut ez,
        );
    } else {
        dp::extend_single_affine(
            query, target, 5, &mat,
            scoring.gap_open as i8, scoring.gap_extend as i8,
            bandwidth, zdrop, 0, EXTENSION_ONLY, &mut ez,
        );
    }

    let cigar = raw_cigar_to_string(&ez.cigar);
    let (qs, qe, ts, te) = cigar_bounds(&ez.cigar);
    let ext_qe = if ez.max_score_query_pos >= 0 { std::cmp::min(qe, ez.max_score_query_pos as usize + 1) } else { qe };
    let ext_te = if ez.max_score_target_pos >= 0 { std::cmp::min(te, ez.max_score_target_pos as usize + 1) } else { te };
    DpAlignment {
        score: ez.score,
        cigar,
        query_start: qs,
        query_end: ext_qe,
        target_start: ts,
        target_end: ext_te,
    }
}

/// Encode an ASCII DNA sequence to nt4 (A=0, C=1, G=2, T=3, N=4).
///
/// Required for DP alignment functions which operate on nt4-encoded sequences.
pub fn encode_nt4(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| crate::align::extend::encode_nt4_byte(b)).collect()
}

/// Compute true alignment score from raw CIGAR + sequences.
/// Walks the CIGAR, scoring matches/mismatches from the matrix and gaps with
/// affine penalties. If gap_open2/gap_extend2 are nonzero, uses dual-affine
/// model (min of both gap costs).
#[allow(clippy::too_many_arguments)]
fn compute_score_from_cigar(
    cigar: &[u32], query: &[u8], target: &[u8],
    mat: &[i8], gapo: i32, gape: i32, gapo2: i32, gape2: i32,
) -> i32 {
    let dual = gapo2 > 0 || gape2 > 0;
    let mut score = 0i32;
    let mut qi = 0usize;
    let mut ti = 0usize;
    for &c in cigar {
        let len = (c >> 4) as usize;
        match c & 0xf {
            0 => { // M: match/mismatch
                for _ in 0..len {
                    if qi < query.len() && ti < target.len() {
                        score += mat[target[ti].min(4) as usize * 5 + query[qi].min(4) as usize] as i32;
                    }
                    qi += 1;
                    ti += 1;
                }
            }
            1 => { // I: insertion (query gap)
                let cost = gapo + gape * len as i32;
                if dual {
                    let cost2 = gapo2 + gape2 * len as i32;
                    score -= cost.min(cost2);
                } else {
                    score -= cost;
                }
                qi += len;
            }
            2 => { // D: deletion (target gap)
                let cost = gapo + gape * len as i32;
                if dual {
                    let cost2 = gapo2 + gape2 * len as i32;
                    score -= cost.min(cost2);
                } else {
                    score -= cost;
                }
                ti += len;
            }
            _ => {}
        }
    }
    score
}

/// Convert raw u32 CIGAR (len<<4 | op) to a human-readable string.
fn raw_cigar_to_string(cigar: &[u32]) -> String {
    let mut s = String::new();
    for &c in cigar {
        let len = c >> 4;
        let op = match c & 0xf {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            _ => '?',
        };
        s.push_str(&format!("{}{}", len, op));
    }
    s
}

/// Compute aligned region bounds from raw CIGAR.
fn cigar_bounds(cigar: &[u32]) -> (usize, usize, usize, usize) {
    let mut qlen = 0usize;
    let mut tlen = 0usize;
    for &c in cigar {
        let len = (c >> 4) as usize;
        match c & 0xf {
            0 => { qlen += len; tlen += len; } // M
            1 => { qlen += len; }               // I
            2 => { tlen += len; }               // D
            _ => {}
        }
    }
    (0, qlen, 0, tlen)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Small reference with enough distinct k-mers to build a usable index.
    fn test_seqs() -> Vec<(String, Vec<u8>)> {
        let bases = b"ACGTTGCAAGGCTTAACCGGATCGATCGGATTACAGGCCTTAACGTATCG";
        let seq: Vec<u8> = (0..4000).map(|i| bases[i % bases.len()]).collect();
        vec![("chr1".to_string(), seq)]
    }

    /// chain_gap_scale * 0.01 * k, with chain_gap_scale at its 0.8 default.
    /// Widened from f32 exactly as `finalize_options` does.
    fn expected_chn_pen_gap(k: usize) -> f32 {
        (0.8f32 as f64 * 0.01 * k as f64) as f32
    }

    #[test]
    fn test_builder_k_override_rescales_chn_pen_gap() {
        let a = Aligner::builder(Preset::Dna).k(12).w(1).from_seqs(test_seqs()).unwrap();
        assert_eq!(a.index().kmer_size, 12);
        assert_eq!(a.index().window_size, 1);
        assert_eq!(a.options().chaining.chn_pen_gap, expected_chn_pen_gap(12));

        let b = Aligner::builder(Preset::Dna).k(15).from_seqs(test_seqs()).unwrap();
        assert_eq!(b.options().chaining.chn_pen_gap, expected_chn_pen_gap(15));
    }

    #[test]
    fn test_from_index_derives_chn_pen_gap_from_index_k() {
        let path = std::env::temp_dir().join("bonefish_api_from_index_k25.rmmi");
        let p = path.to_str().unwrap();
        Aligner::builder(Preset::Dna).k(25).w(20).from_seqs(test_seqs()).unwrap()
            .save_index(p).unwrap();

        let a = Aligner::from_index(p, Preset::Dna).unwrap();
        assert_eq!(a.index().kmer_size, 25);
        assert_eq!(a.options().chaining.chn_pen_gap, expected_chn_pen_gap(25));

        // Geometry of a pre-built index is fixed; a conflicting override is an error.
        assert!(Aligner::builder(Preset::Dna).k(15).from_index(p).is_err());
        assert!(Aligner::builder(Preset::Dna).k(25).w(20).from_index(p).is_ok());

        let _ = std::fs::remove_file(p);
    }

    #[test]
    fn test_mid_occ_only_calibrated_when_unset() {
        // dna leaves it at 0, so it is calibrated from the index.
        let ont = Aligner::from_seqs(test_seqs(), Preset::Dna);
        assert!(ont.options().seeding.mid_occ > 0);

        let explicit = Aligner::builder(Preset::Dna).mid_occ(77).from_seqs(test_seqs()).unwrap();
        assert_eq!(explicit.options().seeding.mid_occ, 77);
    }

    #[test]
    fn test_out_cigar_flag_matches_output_config() {
        for a in [
            Aligner::from_seqs(test_seqs(), Preset::Dna),
            Aligner::builder(Preset::Dna).cigar(true).from_seqs(test_seqs()).unwrap(),
        ] {
            assert!(a.options().flags.contains(AlignFlags::OUT_CIGAR));
        }
        let off = Aligner::builder(Preset::Dna).cigar(false).from_seqs(test_seqs()).unwrap();
        assert!(!off.options().flags.contains(AlignFlags::OUT_CIGAR));
    }

    #[test]
    fn test_builder_rejects_invalid_params() {
        assert!(Aligner::builder(Preset::Dna).k(0).from_seqs(test_seqs()).is_err());
        assert!(Aligner::builder(Preset::Dna).k(MAX_KMER + 1).from_seqs(test_seqs()).is_err());
        assert!(Aligner::builder(Preset::Dna).w(0).from_seqs(test_seqs()).is_err());
        assert!(Aligner::builder(Preset::Dna).mid_occ_frac(1.0).from_seqs(test_seqs()).is_err());
        assert!(Aligner::builder(Preset::Dna).min_mid_occ(100).max_mid_occ(10)
            .from_seqs(test_seqs()).is_err());
        assert!(Aligner::builder(Preset::Dna).k(MAX_KMER).from_seqs(test_seqs()).is_ok());
    }

    #[test]
    fn test_from_index_reader_bounded_and_trailing() {
        // Build a tiny index and serialize it (RMMI magic + bincode) into a
        // buffer. Then append junk "footer" bytes after the index.
        let seqs = vec![
            ("s1".to_string(), b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec()),
            ("s2".to_string(), b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA".to_vec()),
            ("s3".to_string(), b"ACGTNNNNACGTACGTACGTACGTNNNNACGT".to_vec()),
        ];
        let index = Index::build(seqs, 10, 15, usize::MAX);
        let mut buf: Vec<u8> = Vec::new();
        index.save_part(&mut buf).expect("serialize index");
        let index_len = buf.len();
        buf.extend_from_slice(b"OARFISHSIG-LIKE-TRAILING-FOOTER-BYTES");

        // (a) A reader bounded to exactly the index bytes loads correctly.
        let mut bounded = std::io::Cursor::new(&buf[..index_len]);
        let a_bounded = Aligner::from_index_reader(&mut bounded, Preset::Dna)
            .expect("from_index_reader (bounded) should succeed");
        assert_eq!(a_bounded.index.seqs.len(), 3);

        // (b) A reader that also contains trailing bytes loads the same index
        // (the loader stops at the end of the serialized structure).
        let mut with_trailing = std::io::Cursor::new(&buf[..]);
        let a_trailing = Aligner::from_index_reader(&mut with_trailing, Preset::Dna)
            .expect("from_index_reader (with trailing bytes) should succeed");
        assert_eq!(a_trailing.index.seqs.len(), 3);

        // Both reader paths agree with each other on the sequence names.
        let names_b: Vec<&str> = a_bounded.index.seqs.iter().map(|t| t.name.as_str()).collect();
        let names_t: Vec<&str> = a_trailing.index.seqs.iter().map(|t| t.name.as_str()).collect();
        assert_eq!(names_b, names_t);
        assert_eq!(names_b, vec!["s1", "s2", "s3"]);
    }

    #[test]
    fn test_preset_as_str_roundtrip() {
        let presets = [Preset::Dna];
        for p in &presets {
            let s = p.as_str();
            assert!(!s.is_empty(), "preset {:?} should have a string form", p);
        }
    }

    #[test]
    fn test_apply_preset_aliases_and_unknown() {
        let apply = |name: &str| {
            let mut opt = MapOptions::default();
            let (mut k, mut w) = (15, 10);
            apply_preset_str(&mut opt, &mut k, &mut w, name).map(|_| (format!("{:?}", opt), k, w))
        };
        let dna = apply("dna").unwrap();
        assert_eq!((dna.1, dna.2), (19, 19));
        // Sub-anchoring must skip tandem repeats: inside them it misaligns reads (fixed in v0.1.4).
        assert!(dna.0.contains("gap_subanchor: true, gap_subanchor_tr: 150"), "dna preset: sub-anchor outside tandem repeats");
        for alias in ["lr:hq", "lr:hq-fast", "map-ont", "lr"] {
            assert_eq!(apply(alias).unwrap(), dna, "{alias} should be the dna preset");
        }
        for gone in ["lr:hq-turbo", "sr", "splice", "map-pb", "asm5"] {
            assert!(apply(gone).is_err(), "{gone} should be rejected");
        }
    }

    #[test]
    fn test_map_result_empty() {
        let result = MapResult { mappings: Vec::new() };
        assert!(result.mappings.is_empty());
    }

    #[test]
    fn test_dp_align_identical() {
        let seq = encode_nt4(b"ACGTACGTACGT");
        let result = dp_align(&seq, &seq, &DpScoring::default(), -1);
        assert!(result.score > 0, "identical seqs should have positive score");
        assert_eq!(result.cigar, "12M");
    }

    #[test]
    fn test_dp_align_with_mismatch() {
        let query = encode_nt4(b"ACGTACGTACGT");
        let target = encode_nt4(b"ACGTACCTACGT");
        //                              ^ mismatch
        let result = dp_align(&query, &target, &DpScoring::default(), -1);
        assert!(result.score > 0);
        assert!(result.cigar.contains('M'), "should contain M ops");
    }

    #[test]
    fn test_dp_align_with_insertion() {
        let query  = encode_nt4(b"ACGTAAAACGT");
        let target = encode_nt4(b"ACGTACGT");
        let result = dp_align(&query, &target, &DpScoring::default(), -1);
        assert!(result.cigar.contains('I'), "should contain insertion: {}", result.cigar);
    }

    #[test]
    fn test_dp_global_identical() {
        let query = encode_nt4(b"ACGTACGT");
        let target = encode_nt4(b"ACGTACGT");
        let result = dp_global(&query, &target, &DpScoring::default(), -1);
        assert_eq!(result.query_end, 8);
        assert_eq!(result.target_end, 8);
        assert_eq!(result.cigar, "8M");
        assert_eq!(result.score, 16); // 8 * 2
    }

    #[test]
    fn test_dp_global_with_indel() {
        let query  = encode_nt4(b"ACGTAAACGT");
        let target = encode_nt4(b"ACGTACGT");
        let result = dp_global(&query, &target, &DpScoring::default(), -1);
        // Must cover both sequences end-to-end
        let total_q: usize = result.cigar.chars().filter(|c| !c.is_ascii_digit())
            .zip(result.cigar.split(|c: char| !c.is_ascii_digit()).filter(|s| !s.is_empty()))
            .map(|(op, len)| { let l: usize = len.parse().unwrap(); if op == 'D' { 0 } else { l } }).sum();
        let total_t: usize = result.cigar.chars().filter(|c| !c.is_ascii_digit())
            .zip(result.cigar.split(|c: char| !c.is_ascii_digit()).filter(|s| !s.is_empty()))
            .map(|(op, len)| { let l: usize = len.parse().unwrap(); if op == 'I' { 0 } else { l } }).sum();
        assert_eq!(total_q, 10, "CIGAR should consume full query (10bp): cigar={}", result.cigar);
        assert_eq!(total_t, 8, "CIGAR should consume full target (8bp): cigar={}", result.cigar);
    }

    #[test]
    fn test_dp_local_finds_best_region() {
        // Embed a matching region in noise
        let query  = encode_nt4(b"TTTTACGTACGTACGTTTTT");
        let target = encode_nt4(b"GGGGACGTACGTACGTGGGG");
        let result = dp_local(&query, &target, &DpScoring::default());
        assert!(result.score > 0);
        // Local alignment should find the ACGTACGTACGT match, not the flanking noise
        assert!(result.query_start >= 3, "should skip noise prefix: qs={}", result.query_start);
        assert!(result.query_end <= 16, "should skip noise suffix: qe={}", result.query_end);
    }

    #[test]
    fn test_dp_extension() {
        // Extension aligns from position 0, stopping at the best score
        let query  = encode_nt4(b"ACGTACGTACGTACGTACGTACGTACGTACGTCGATCGATCGATCGATCG");
        let target = encode_nt4(b"ACGTACGTACGTACGTACGTACGTACGTACGTTTTTTTTTTTTTTTTTTTT");
        let result = dp_extension(&query, &target, &DpScoring::default(), -1);
        // Should produce a valid alignment with CIGAR
        assert!(!result.cigar.is_empty(), "extension should produce CIGAR");
    }

    #[test]
    fn test_dp_dual_affine() {
        let query = encode_nt4(b"ACGTACGT");
        let target = encode_nt4(b"ACGTACGT");
        let scoring = DpScoring {
            match_score: 2, mismatch: 4,
            gap_open: 4, gap_extend: 2,
            gap_open2: 24, gap_extend2: 1,
        };
        let result = dp_align(&query, &target, &scoring, -1);
        assert_eq!(result.score, 16); // 8 * 2
        assert_eq!(result.cigar, "8M");
    }

    #[test]
    fn test_encode_nt4() {
        let encoded = encode_nt4(b"ACGTNacgtn");
        assert_eq!(encoded, vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4]);
    }
}
