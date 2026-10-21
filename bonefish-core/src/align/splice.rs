//! Spliced alignment of a chain (RNA/cDNA presets).
//!
//! [`align_anchors_splice`] is the spliced counterpart of [`extend::align_anchors`]: the same
//! left extension / gap fill / right extension skeleton, with every DP segment aligned by the
//! intron-aware kernel ([`dp::extend_splice`], unbanded, so introns up to the chaining
//! bandwidth `-G` fit), junction annotations from `--junc-bed`, and minimap2's splice-specific
//! end trimming ([`fix_bad_ends_splice`]). The genomic speed heuristics of the `dna` preset
//! (narrow and adaptive extension bands, the low-identity probe, ungapped-first and
//! sub-anchored gap filling, the narrow gap-fill band) are not used here.
//!
//! [`extend::align_anchors`]: crate::align::extend::align_anchors

use std::cmp;

use crate::align::dp::{self, DpResult, APPROX_MAX, EXTENSION_ONLY, GENERIC_SCORING, REV_CIGAR, RIGHT_ALIGN,
    SPLICE_COMPLEX, SPLICE_FLANK, SPLICE_FORWARD, SPLICE_REVERSE};
use crate::align::extend::{append_cigar, build_scoring_matrix_full, compute_left_boundary, compute_right_boundary,
    filter_bad_seeds, filter_bad_seeds_alt, finalize_cigar, test_zdrop, AlignAnchorContext, AlignResult,
    SEED_IGNORE, SEED_LONG_JOIN, SEED_TANDEM};
use crate::align::index::Index;
use crate::align::junc::JunctionDb;
use crate::align::map::{AlignFlags, MapOptions};
use crate::align::sketch::Minimizer;

/// Score of the extension of `anchor` by `ext_len` bases on both sides (lightweight local DP),
/// as minimap2's `mm_seed_ext_score`. With `extract` (index, rid, region start), the reference
/// window is read from the index (anchor coordinates are region-relative); otherwise it is
/// sliced from `tseq`.
#[allow(clippy::too_many_arguments)]
fn seed_extension_score(
    qseq: &[u8],
    tseq: &[u8],
    anchor: &Minimizer,
    ext_len: i32,
    mat: &[i8],
    gapo: i32,
    gape: i32,
    extract: Option<(&Index, usize, usize)>,
) -> i32 {
    let q_span = anchor.query_span();
    let re_a = anchor.ref_pos() + 1;
    let rs_a = re_a - q_span;
    let qe_a = anchor.query_pos() + 1;
    let qs_a = qe_a - q_span;
    let qs = cmp::max(0, qs_a - ext_len) as usize;
    let qe = cmp::min(qseq.len() as i32, qe_a + ext_len) as usize;
    if qe <= qs { return 0; }

    let mut owned_t: Vec<u8> = Vec::new();
    let t_slice: &[u8] = if let Some((mi, rid, rgn_start)) = extract {
        let chrom_len = mi.seqs[rid].len as i32;
        let abs_rs = cmp::max(0, rs_a + rgn_start as i32 - ext_len) as usize;
        let abs_re = cmp::min(chrom_len, re_a + rgn_start as i32 + ext_len) as usize;
        if abs_re <= abs_rs { return 0; }
        owned_t.resize(abs_re - abs_rs, 0);
        mi.extract_nt4_into(rid, abs_rs, abs_re, &mut owned_t);
        &owned_t
    } else {
        let rs = cmp::max(0, rs_a - ext_len) as usize;
        let re = cmp::min(tseq.len() as i32, re_a + ext_len) as usize;
        if re <= rs { return 0; }
        &tseq[rs..re]
    };

    let mut qp = dp::lightweight_profile_init((qe - qs) as i32, &qseq[qs..qe], 5, mat);
    dp::lightweight_align_i16(&mut qp, t_slice.len() as i32, t_slice, gapo, gape).0
}

/// minimap2's `mm_fix_bad_ends_splice`: drop the first / last anchor when it is a short,
/// weakly supported exon across a long gap. Returns (offset, count) of the kept anchors.
#[allow(clippy::too_many_arguments)]
fn fix_bad_ends_splice(
    anchors: &[Minimizer],
    qseq: &[u8],
    tseq: &[u8],
    mat: &[i8],
    opt: &MapOptions,
    extract: Option<(&Index, usize, usize)>,
) -> (usize, usize) {
    let n = anchors.len();
    if n < 3 { return (0, n); }
    let ext_len = opt.alignment.anchor_ext_len;
    let ext_shift = opt.alignment.anchor_ext_shift as f64;
    let (gapo, gape) = (opt.scoring.gap_open, opt.scoring.gap_extend);
    let weak = |a: &Minimizer, gap: i32| -> bool {
        if gap <= 0 { return false; }
        let min_len = (gap as f64).ln() + ext_shift;
        if a.query_span() as f64 >= min_len { return false; }
        let score = seed_extension_score(qseq, tseq, a, ext_len, mat, gapo, gape, extract);
        (score as f64) / (mat[0] as f64) < min_len
    };
    let (mut off, mut cnt) = (0, n);
    if weak(&anchors[0], anchors[1].ref_pos() - anchors[0].ref_pos()) {
        off += 1;
        cnt -= 1;
    }
    if weak(&anchors[n - 1], anchors[n - 1].ref_pos() - anchors[n - 2].ref_pos()) {
        cnt -= 1;
    }
    (off, cnt)
}

/// Splice-site DP flags for one alignment pass. `strand` holds `SPLICE_FOR` and/or
/// `SPLICE_REV` (transcript strand relative to the read); on a reverse-strand chain the
/// reference-side orientation flips.
fn splice_dp_flags(opt: &MapOptions, strand: AlignFlags, rev: bool) -> i32 {
    let mut f = SPLICE_COMPLEX;
    if strand.contains(AlignFlags::SPLICE_FOR) {
        f |= if rev { SPLICE_REVERSE } else { SPLICE_FORWARD };
    }
    if strand.contains(AlignFlags::SPLICE_REV) {
        f |= if rev { SPLICE_FORWARD } else { SPLICE_REVERSE };
    }
    if opt.flags.contains(AlignFlags::SPLICE_FLANK) {
        f |= SPLICE_FLANK;
    }
    f
}

/// One spliced DP call; `max_sw_mat` (0 = no limit) skips oversized matrices as a z-drop.
#[allow(clippy::too_many_arguments)]
fn splice_dp(qs: &[u8], ts: &[u8], mat: &[i8; 25], opt: &MapOptions, zdrop: i32, end_bonus: i32, flags: i32,
             junc: Option<&[u8]>, ez: &mut DpResult) {
    let max_sw_mat = opt.alignment.max_sw_mat;
    if max_sw_mat > 0 && qs.len() as i64 * ts.len() as i64 > max_sw_mat {
        ez.max_score_query_pos = -1;
        ez.max_score_target_pos = -1;
        ez.max_query_end_target_pos = -1;
        ez.zdropped = 1;
        return;
    }
    let sc = &opt.scoring;
    dp::extend_splice(qs, ts, 5, mat, sc.gap_open as i8, sc.gap_extend as i8, sc.gap_open2 as i8,
        sc.noncanon_penalty as i8, zdrop, end_bonus, sc.junc_bonus as i8, sc.junc_pen as i8, flags, junc, ez);
}

/// Annotated splice sites of reference window [st, en) (region-relative; `ref_offset` makes
/// them absolute), or None without annotations.
fn junctions(db: Option<&JunctionDb>, rid: usize, ref_offset: i32, st: i32, en: i32) -> Option<Vec<u8>> {
    db.map(|db| {
        let mut buf = vec![0u8; (en - st) as usize];
        db.fill(rid, st + ref_offset, en + ref_offset, &mut buf);
        buf
    })
}

/// Spliced alignment of a chain. Arguments and result as [`crate::align::extend::align_anchors`],
/// plus `strand` (`SPLICE_FOR` / `SPLICE_REV`: which transcript strand's splice signals to score)
/// and the optional junction annotations.
#[allow(clippy::too_many_arguments)]
pub fn align_anchors_splice(
    anchors: &mut [Minimizer],
    qseq: &[u8],
    tseq: &mut [u8],
    lazy_extract: Option<(&Index, usize, usize)>,
    opt: &MapOptions,
    call: &AlignAnchorContext,
    strand: AlignFlags,
    junc_db: Option<&JunctionDb>,
) -> AlignResult {
    let empty = AlignResult { cigar_ops: Vec::new(), query_start: 0, query_end: 0, ref_start: 0, ref_end: 0, dp_score: 0,
        split_right_anchors: None, split_offset_in_orig: None, split_inv: false, probe_claim_qs: None, probe_claim_qe: None };
    if anchors.is_empty() { return empty; }

    let sc = &opt.scoring;
    let (zdrop, zdrop_inv, end_bonus) = (opt.alignment.zdrop, opt.alignment.zdrop_inv, opt.alignment.end_bonus);
    let min_dp_len = opt.alignment.min_dp_len;
    let seed_bounds = call.seed_bounds;
    let rid = call.rid;
    let ref_offset = call.ref_offset as i32;
    let mat = build_scoring_matrix_full(sc.match_score, sc.mismatch_penalty, sc.transition, sc.ambig_penalty);
    let generic_sc_flag = if sc.transition != 0 && sc.transition != sc.mismatch_penalty { GENERIC_SCORING } else { 0 };
    let splice_flag = splice_dp_flags(opt, strand, call.rev);
    let qlen = qseq.len() as i32;
    let tlen = tseq.len() as i32;
    let k_half = (call.k >> 1) as i32;

    // Fill a segment of `tseq` from the index (lazy mode: `tseq` arrives uninitialized; with
    // `lazy_extract = None` the caller has filled it).
    let fill = |tseq: &mut [u8], a: i32, b: i32| {
        if let Some((mi, rid, rgn_start)) = lazy_extract && b > a {
            let (a, b) = (a as usize, b as usize);
            mi.extract_nt4_into(rid, rgn_start + a, rgn_start + b, &mut tseq[a..b]);
        }
    };

    // --- Anchor trimming (minimap2's splice end fix instead of the genomic one) ---
    let (as1_offset, cnt1) = if opt.flags.contains(AlignFlags::NO_END_FLT) {
        (0, anchors.len())
    } else {
        fix_bad_ends_splice(anchors, qseq, tseq, &mat, opt, lazy_extract)
    };
    if cnt1 == 0 { return empty; }
    let max_gap = opt.chaining.max_gap;
    let mut work_anchors = anchors[as1_offset..as1_offset + cnt1].to_vec();
    filter_bad_seeds(&mut work_anchors, 10, 40, max_gap >> 1, 10);
    filter_bad_seeds_alt(&mut work_anchors, 30, max_gap >> 1);
    for i in 0..cnt1 {
        let fl = work_anchors[i].y & (SEED_IGNORE | SEED_LONG_JOIN);
        if fl != 0 { anchors[as1_offset + i].y |= fl; }
    }
    let (rs, qs) = (work_anchors[0].ref_pos() - k_half, work_anchors[0].query_pos() - k_half);
    let (re, qe) = (work_anchors[cnt1 - 1].ref_pos() - k_half, work_anchors[cnt1 - 1].query_pos() - k_half);

    let mut raw_cigar: Vec<u32> = Vec::new();
    let mut dp_score: i32 = 0;

    // --- Left extension ---
    let orig_first = &anchors[0];
    let (rs0, qs0) = compute_left_boundary(rs, qs, orig_first.ref_pos(), orig_first.query_pos(), orig_first.query_span(),
        seed_bounds.0, seed_bounds.1, opt);
    let (mut rs1, mut qs1) = (rs, qs);
    if qs > 0 && rs > 0 && qs > qs0 && rs > rs0 {
        fill(tseq, rs0, rs);
        let q_rev: Vec<u8> = qseq[qs0 as usize..qs as usize].iter().rev().copied().collect();
        let t_rev: Vec<u8> = tseq[rs0 as usize..rs as usize].iter().rev().copied().collect();
        let junc = junctions(junc_db, rid, ref_offset, rs0, rs).map(|mut j| { j.reverse(); j });
        let left_zdrop = if call.split_inv { zdrop_inv } else { zdrop };
        let mut ez = DpResult::default();
        splice_dp(&q_rev, &t_rev, &mat, opt, left_zdrop, end_bonus,
            splice_flag | EXTENSION_ONLY | RIGHT_ALIGN | REV_CIGAR | generic_sc_flag, junc.as_deref(), &mut ez);
        if !ez.cigar.is_empty() {
            append_cigar(&mut raw_cigar, &ez.cigar);
            dp_score += ez.max;
        }
        rs1 = rs - if ez.reach_end != 0 { ez.max_query_end_target_pos + 1 } else { ez.max_score_target_pos + 1 };
        qs1 = qs - if ez.reach_end != 0 { qs - qs0 } else { ez.max_score_query_pos + 1 };
    }

    // --- Gap filling between anchors ---
    let (mut re1, mut qe1) = (rs, qs);
    let (mut gap_rs, mut gap_qs) = (rs, qs);
    let mut dropped = false;
    let mut split_right_anchors: Option<Vec<Minimizer>> = None;
    let mut split_offset_in_orig: Option<usize> = None;
    let mut split_inv_flag = false;
    for i in 1..cnt1 {
        if (work_anchors[i].y & (SEED_IGNORE | SEED_TANDEM)) != 0 && i != cnt1 - 1 {
            continue;
        }
        let (anchor_re, anchor_qe) = (work_anchors[i].ref_pos() - k_half, work_anchors[i].query_pos() - k_half);
        re1 = anchor_re;
        qe1 = anchor_qe;
        let should_align = i == cnt1 - 1 || (work_anchors[i].y & SEED_LONG_JOIN) != 0
            || (anchor_qe - gap_qs >= min_dp_len && anchor_re - gap_rs >= min_dp_len);
        if !should_align { continue; }

        let (q_start, q_end) = (cmp::max(0, gap_qs), cmp::min(qlen, anchor_qe));
        let (t_start, t_end) = (cmp::max(0, gap_rs), cmp::min(tlen, anchor_re));
        fill(tseq, t_start, t_end);
        if q_end > q_start && t_end > t_start {
            let q_sub = &qseq[q_start as usize..q_end as usize];
            let t_sub = &tseq[t_start as usize..t_end as usize];
            let junc = junctions(junc_db, rid, ref_offset, gap_rs, anchor_re);
            // Gap fill: end_bonus -1; first pass with approximate z-drop, then a full pass if the
            // path shows a (possibly inversion) z-drop.
            let mut ez = DpResult::default();
            splice_dp(q_sub, t_sub, &mat, opt, zdrop, -1, splice_flag | APPROX_MAX | generic_sc_flag, junc.as_deref(), &mut ez);
            let zdrop_code = if !ez.cigar.is_empty() { test_zdrop(q_sub, t_sub, &ez.cigar, &mat, opt) } else { 0 };
            if zdrop_code != 0 {
                ez = DpResult::default();
                let zdrop2 = if zdrop_code == 2 { zdrop_inv } else { zdrop };
                splice_dp(q_sub, t_sub, &mat, opt, zdrop2, -1, splice_flag | generic_sc_flag, junc.as_deref(), &mut ez);
            }
            if !ez.cigar.is_empty() {
                append_cigar(&mut raw_cigar, &ez.cigar);
            }
            if ez.zdropped != 0 {
                dp_score += ez.max;
                re1 = gap_rs + ez.max_score_target_pos + 1;
                qe1 = gap_qs + ez.max_score_query_pos + 1;
                // Split the chain after the last anchor before the z-drop point.
                let max_re_split = gap_rs + ez.max_score_target_pos;
                let j_split = (0..i).rev().find(|&j| work_anchors[j].ref_pos() <= max_re_split).unwrap_or(0);
                if cnt1 as i32 - (j_split as i32 + 1) >= opt.chaining.min_cnt {
                    let split_off = as1_offset + j_split + 1;
                    if split_off < anchors.len() {
                        split_right_anchors = Some(anchors[split_off..].to_vec());
                        split_offset_in_orig = Some(split_off);
                        split_inv_flag = zdrop_code == 2;
                    }
                }
                dropped = true;
                break;
            }
            dp_score += ez.score;
        }
        gap_rs = anchor_re;
        gap_qs = anchor_qe;
    }

    // --- Right extension ---
    let (mut final_re, mut final_qe) = (re1, qe1);
    if !dropped && qe < qlen && re < tlen {
        let orig_last = &anchors[anchors.len() - 1];
        let (re0, qe0) = compute_right_boundary(re, qe, orig_last.ref_pos() + 1, orig_last.query_pos() + 1,
            seed_bounds.2, seed_bounds.3, opt, qlen, tlen, seed_bounds);
        if qe0 > qe && re0 > re {
            fill(tseq, re, re0);
            let q_sub = &qseq[qe as usize..qe0 as usize];
            let t_sub = &tseq[re as usize..re0 as usize];
            let junc = junctions(junc_db, rid, ref_offset, re, re0);
            let mut ez = DpResult::default();
            splice_dp(q_sub, t_sub, &mat, opt, zdrop, end_bonus, splice_flag | EXTENSION_ONLY | generic_sc_flag,
                junc.as_deref(), &mut ez);
            if !ez.cigar.is_empty() {
                append_cigar(&mut raw_cigar, &ez.cigar);
                dp_score += ez.max;
            }
            final_re = re + if ez.reach_end != 0 { ez.max_query_end_target_pos + 1 } else { ez.max_score_target_pos + 1 };
            final_qe = qe + if ez.reach_end != 0 { qe0 - qe } else { ez.max_score_query_pos + 1 };
        }
    }
    if final_re <= 0 { final_re = re1; }
    if final_qe <= 0 { final_qe = qe1; }

    // The alignment range must be filled for CIGAR finalization and tag generation.
    let a = cmp::max(0, rs1);
    fill(tseq, a, cmp::max(a, final_re));
    let (cigar_ops, final_qs, final_rs) = finalize_cigar(raw_cigar, qseq, tseq, qs1, rs1);

    AlignResult {
        cigar_ops,
        query_start: cmp::max(0, final_qs) as usize,
        query_end: cmp::max(0, final_qe) as usize,
        ref_start: cmp::max(0, final_rs) as usize,
        ref_end: cmp::max(0, final_re) as usize,
        dp_score,
        split_right_anchors,
        split_offset_in_orig,
        split_inv: split_inv_flag,
        probe_claim_qs: None,
        probe_claim_qe: None,
    }
}
