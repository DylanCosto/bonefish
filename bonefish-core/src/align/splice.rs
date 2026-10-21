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

/// Reference window margin (beyond the query gap) kept on each side of a compressed intron gap.
const GAP_MARGIN: usize = 64;
/// N bases separating the two windows of a compressed gap: crossing them other than by the
/// intron (long-gap) state costs far more than an intron.
const SEAM: usize = 256;
/// A junction may not lie this close to the seam (the splice-site scores read 3 bases around it).
const SEAM_CLEAR: usize = 8;
/// Query bases inspected on each side of a widened intron, and the most errors (inserted,
/// deleted or mismatched bases) allowed there. More suggests an exon inside the removed
/// reference (a missing anchor), which only the full DP can place.
const FLANK_LEN: i32 = 32;
/// Query bases on each side of a rough intron searched for in the removed reference.
const RESCUE_Q: i32 = 48;
/// Most exon candidates added per rescue round.
const MAX_CAND: usize = 6;
/// Longest exon assumed on either side of a rescue candidate (an unanchored exon this long is
/// very unlikely; the window would otherwise reach the neighbouring anchors).
const RESCUE_EXON_MAX: i32 = 300;
/// Rounds of exon rescue (each adds windows and aligns again).
const MAX_RESCUE_ROUNDS: usize = 2;
const FLANK_MAX: i32 = 6;

/// Gap fill through long reference gaps (introns) without the full query x reference matrix.
/// An intron costs the same at any length (second gap extension 0), so reference stretches far
/// from every chain anchor only matter as intron interiors. Between consecutive anchor points
/// (`pts`: (query, reference) positions in the gap, gap ends included) whose reference distance
/// exceeds the query distance by more than two windows, the middle is replaced by an N seam;
/// the query is aligned to the compressed reference and every intron crossing a seam is widened
/// back to its real length. When the query aligns poorly next to such an intron, an exon may lie
/// in the removed reference without an anchor: exact 10-mer matches of that query stretch are
/// searched there, and a matching stretch is kept as a new point and the gap aligned again.
/// Returns false (the caller runs the full DP) when nothing is cut, the alignment z-drops, or
/// the path crosses a seam other than inside one intron clear of it.
#[allow(clippy::too_many_arguments)]
fn splice_gap_compressed(q: &[u8], t: &[u8], pts: &[(i32, i32)], mat: &[i8; 25], opt: &MapOptions, zdrop: i32,
                         flags: i32, junc: Option<&[u8]>, ez: &mut DpResult) -> bool {
    let m = GAP_MARGIN as i32;
    let tl = t.len() as i32;
    // Reference intervals kept: around each anchor point, as far as the query distance to its
    // neighbours (plus a margin) could carry the alignment.
    let mut keeps: Vec<(i32, i32)> = Vec::with_capacity(2 * pts.len());
    for w in pts.windows(2) {
        let ((qa, ta), (qb, tb)) = (w[0], w[1]);
        let dq = (qb - qa).max(0);
        keeps.push((ta, ta + dq + m));
        keeps.push((tb - dq - m, tb));
    }
    for round in 0..=MAX_RESCUE_ROUNDS {
        keeps.sort_unstable();
        let mut cuts: Vec<(usize, usize)> = Vec::new();
        let mut covered = 0i32;
        for &(a, b) in &keeps {
            let a = a.clamp(0, tl);
            if a - covered > 2 * SEAM as i32 {
                cuts.push((covered as usize, a as usize));
            }
            covered = covered.max(b.clamp(0, tl));
        }
        if cuts.is_empty() {
            return false;
        }
        // Compressed reference; `seams` holds each seam's start in `tc` and the bases it replaces.
        let mut tc = Vec::with_capacity(t.len().min(q.len() * 4 + cuts.len() * (SEAM + 4 * GAP_MARGIN)));
        let mut jc = junc.map(|_| Vec::with_capacity(tc.capacity()));
        let mut seams: Vec<(usize, usize)> = Vec::with_capacity(cuts.len());
        let mut prev = 0;
        for &(a, b) in &cuts {
            tc.extend_from_slice(&t[prev..a]);
            if let (Some(v), Some(j)) = (jc.as_mut(), junc) { v.extend_from_slice(&j[prev..a]); }
            seams.push((tc.len(), b - a));
            tc.resize(tc.len() + SEAM, 4);
            if let Some(v) = jc.as_mut() { v.resize(v.len() + SEAM, 0); }
            prev = b;
        }
        tc.extend_from_slice(&t[prev..]);
        if let (Some(v), Some(j)) = (jc.as_mut(), junc) { v.extend_from_slice(&j[prev..]); }
        *ez = DpResult::default();
        splice_dp(q, &tc, mat, opt, zdrop, -1, flags, jc.as_deref(), ez);
        if ez.zdropped != 0 || ez.cigar.is_empty() {
            return false;
        }
        // Each seam must lie inside one N op, SEAM_CLEAR bases from both of its ends.
        let mut tpos = 0usize;
        let mut si = 0usize;
        for c in ez.cigar.iter_mut() {
            let (len, op) = ((*c >> 4) as usize, *c & 0xf);
            if !matches!(op, 0 | 2 | 3 | 7 | 8) {
                continue;
            }
            let end = tpos + len;
            let mut extra = 0usize;
            while si < seams.len() && seams[si].0 < end {
                let (s0, cut) = seams[si];
                if op != 3 || tpos + SEAM_CLEAR > s0 || end < s0 + SEAM + SEAM_CLEAR {
                    return false;
                }
                extra += cut - SEAM;
                si += 1;
            }
            *c += (extra as u32) << 4;
            tpos = end;
        }
        if si != seams.len() {
            return false;
        }
        if round == MAX_RESCUE_ROUNDS {
            return true;
        }
        let Some((seam_lo, seam_hi, at)) = rough_flank(q, &tc, &ez.cigar, &seams) else { return true };
        // An exon missing from that intron lies in one of the reference stretches it skips, and
        // its query bases lie between the anchor points on either side of those stretches.
        let (ra, rb) = (cuts[seam_lo].0 as i32, cuts[seam_hi - 1].1 as i32);
        let qa = pts.iter().filter(|p| p.1 <= ra).map(|p| p.0).max().unwrap_or(0).max(0);
        let qb = pts.iter().filter(|p| p.1 >= rb).map(|p| p.0).min().unwrap_or(q.len() as i32).min(q.len() as i32);
        let (q0, q1) = ((at - RESCUE_Q).max(qa) as usize, (at + RESCUE_Q).min(qb).max(0) as usize);
        if q1 <= q0 {
            return true;
        }
        // Candidates from every skipped stretch, best supported first, at most MAX_CAND in all.
        let mut cand: Vec<(usize, i32, i32)> = Vec::new();
        for &(a, b) in &cuts[seam_lo..seam_hi] {
            cand.extend(find_exons(&q[q0..q1], &t[a..b]).into_iter()
                .map(|(qo, to, n)| (n, (q0 + qo) as i32, (a + to) as i32)));
        }
        if cand.is_empty() {
            return true;
        }
        cand.sort_by_key(|c| std::cmp::Reverse(c.0));
        for &(_, xq, xt) in cand.iter().take(MAX_CAND) {
            keeps.push((xt - (xq - qa).min(RESCUE_EXON_MAX) - m, xt + (qb - xq).min(RESCUE_EXON_MAX) + m));
        }
    }
    true
}

/// The first intron crossing seams that has more than FLANK_MAX errors (inserted, deleted or
/// mismatched bases) within FLANK_LEN query bases on either side, as (first seam, one past its
/// last seam, query position of the intron), or None. `cigar` is already widened; positions in
/// `tc` are recovered by removing the widening.
fn rough_flank(q: &[u8], tc: &[u8], cigar: &[u32], seams: &[(usize, usize)]) -> Option<(usize, usize, i32)> {
    let n = cigar.len();
    let mut qp = Vec::with_capacity(n + 1);
    let mut tp = Vec::with_capacity(n + 1);
    let (mut qi, mut ti) = (0i32, 0i32);
    let mut si = 0usize;
    let mut seam_ops: Vec<(usize, usize, usize)> = Vec::new();
    for (k, &c) in cigar.iter().enumerate() {
        qp.push(qi);
        tp.push(ti);
        let (len, op) = ((c >> 4) as i32, c & 0xf);
        let mut tl = len;
        if op == 3 {
            let first = si;
            let mut extra = 0i32;
            while si < seams.len() && (seams[si].0 as i32) < ti + tl - extra {
                extra += (seams[si].1 - SEAM) as i32;
                si += 1;
            }
            if si > first { seam_ops.push((k, first, si)); }
            tl -= extra;
        }
        if matches!(op, 0 | 1 | 4 | 7 | 8) { qi += len; }
        if matches!(op, 0 | 2 | 3 | 7 | 8) { ti += tl; }
    }
    let errors = |lo: i32, hi: i32| -> i32 {
        let mut e = 0;
        for (k, &c) in cigar.iter().enumerate() {
            let (len, op) = ((c >> 4) as i32, c & 0xf);
            let (q0, t0) = (qp[k], tp[k]);
            if q0 >= hi { break; }
            match op {
                0 | 7 | 8 => {
                    for x in q0.max(lo)..(q0 + len).min(hi) {
                        let tq = (t0 + x - q0) as usize;
                        if q[x as usize] != tc[tq] || q[x as usize] > 3 { e += 1; }
                    }
                }
                1 => { e += ((q0 + len).min(hi) - q0.max(lo)).max(0); }
                2 if q0 > lo && q0 < hi => { e += len; }
                _ => {}
            }
        }
        e
    };
    let fmax = FLANK_MAX;
    seam_ops.iter().find(|&&(k, _, _)| {
        let at = qp[k];
        errors(at - FLANK_LEN, at) > fmax || errors(at, at + FLANK_LEN) > fmax
    }).map(|&(k, lo, hi)| (lo, hi, qp[k]))
}

/// Candidate exon locations for query stretch `qs` inside reference stretch `ts`: exact 10-mer
/// matches grouped by diagonal (within 8), largest groups first, at most MAX_CAND, as (query,
/// reference) centre and the group's number of matches. The aligner decides which, if any, is used, so a spurious candidate only
/// costs a slightly larger reference window.
fn find_exons(qs: &[u8], ts: &[u8]) -> Vec<(usize, usize, usize)> {
    const K: usize = 10;
    const MASK: u32 = (1 << (2 * K)) - 1;
    let mut out = Vec::new();
    if qs.len() < 16 || ts.len() < K {
        return out;
    }
    let mut qk: Vec<(u32, u32)> = Vec::with_capacity(qs.len());
    let (mut h, mut l) = (0u32, 0usize);
    for (i, &c) in qs.iter().enumerate() {
        if c > 3 { l = 0; continue; }
        h = ((h << 2) | c as u32) & MASK;
        l += 1;
        if l >= K { qk.push((h, (i + 1 - K) as u32)); }
    }
    qk.sort_unstable();
    let mut hits: Vec<(i32, u32)> = Vec::new(); // (diagonal, query position)
    let (mut h, mut l) = (0u32, 0usize);
    for (j, &c) in ts.iter().enumerate() {
        if c > 3 { l = 0; continue; }
        h = ((h << 2) | c as u32) & MASK;
        l += 1;
        if l < K { continue; }
        let tpos = (j + 1 - K) as i32;
        let mut i = qk.partition_point(|x| x.0 < h);
        while i < qk.len() && qk[i].0 == h {
            hits.push((tpos - qk[i].1 as i32, qk[i].1));
            i += 1;
        }
        if hits.len() > 4096 { return out; }
    }
    hits.sort_unstable();
    // Groups of hits with diagonals within 8 of the group's first.
    let mut groups: Vec<(usize, usize, usize)> = Vec::new(); // (count, start, end)
    let mut lo = 0;
    while lo < hits.len() {
        let mut hi = lo + 1;
        while hi < hits.len() && hits[hi].0 - hits[lo].0 <= 8 { hi += 1; }
        groups.push((hi - lo, lo, hi));
        lo = hi;
    }
    groups.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    for &(n, a, b) in groups.iter().take(MAX_CAND) {
        let qm = hits[a..b].iter().map(|x| x.1 as usize).sum::<usize>() / n + K / 2;
        let tm = qm as i32 + hits[a + n / 2].0;
        if qm < qs.len() && tm >= 0 && (tm as usize) < ts.len() {
            out.push((qm, tm as usize, n));
        }
    }
    out
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
    let mut gap_i = 0usize;
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
            // Anchor points inside this gap, in gap coordinates (gap ends included).
            let mut pts = Vec::with_capacity(i - gap_i + 1);
            pts.push((0, 0));
            for a in &work_anchors[gap_i + 1..i] {
                let (aq, ar) = (a.query_pos() - k_half - q_start, a.ref_pos() - k_half - t_start);
                if aq > 0 && ar > 0 && aq < q_end - q_start && ar < t_end - t_start { pts.push((aq, ar)); }
            }
            pts.push((q_end - q_start, t_end - t_start));
            if !splice_gap_compressed(q_sub, t_sub, &pts, &mat, opt, zdrop, splice_flag | APPROX_MAX | generic_sc_flag, junc.as_deref(), &mut ez) {
                ez = DpResult::default();
                splice_dp(q_sub, t_sub, &mat, opt, zdrop, -1, splice_flag | APPROX_MAX | generic_sc_flag, junc.as_deref(), &mut ez);
            }
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
        gap_i = i;
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

#[cfg(test)]
mod tests {
    use super::*;

    /// xorshift64 generator for reproducible synthetic sequences (nt4).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
        fn range(&mut self, lo: usize, hi: usize) -> usize { lo + self.below(hi - lo + 1) }
        fn base(&mut self) -> u8 { (self.next() >> 62) as u8 }
        fn seq(&mut self, n: usize) -> Vec<u8> { (0..n).map(|_| self.base()).collect() }
    }

    /// The cdna preset's options, scoring matrix and gap-fill DP flags (as `align_anchors_splice`
    /// builds them for a forward chain; `strand` as `splice_dp_flags`).
    fn cdna(strand: AlignFlags) -> (MapOptions, [i8; 25], i32) {
        let mut opt = MapOptions::default();
        let (mut k, mut w) = (15, 10);
        crate::api::apply_preset_str(&mut opt, &mut k, &mut w, "cdna").unwrap();
        let sc = &opt.scoring;
        let mat = build_scoring_matrix_full(sc.match_score, sc.mismatch_penalty, sc.transition, sc.ambig_penalty);
        assert_eq!(sc.transition, 0, "no GENERIC_SCORING flag needed");
        let flags = splice_dp_flags(&opt, strand, false) | APPROX_MAX;
        (opt, mat, flags)
    }

    fn both() -> AlignFlags { AlignFlags::SPLICE_FOR | AlignFlags::SPLICE_REV }

    /// Reference (first exon start .. last exon end, nt4) of a gene with exons of the given
    /// lengths separated by GT..AG introns of the given lengths with random interiors, and the
    /// exon intervals.
    fn gene(rng: &mut Rng, exons: &[usize], introns: &[usize]) -> (Vec<u8>, Vec<(usize, usize)>) {
        assert_eq!(exons.len(), introns.len() + 1);
        let mut t = Vec::new();
        let mut iv = Vec::new();
        for (i, &e) in exons.iter().enumerate() {
            iv.push((t.len(), t.len() + e));
            t.extend(rng.seq(e));
            if let Some(&l) = introns.get(i) {
                assert!(l >= 4);
                let mut x = rng.seq(l);
                x[..2].copy_from_slice(&[2, 3]);
                x[l - 2..].copy_from_slice(&[0, 2]);
                t.extend(x);
            }
        }
        (t, iv)
    }

    /// The exons of `t` concatenated, with each base substituted / followed by an insertion /
    /// deleted with probability `err` (per mille, each), except within 16 bases of the
    /// reference positions in `marks` and of both ends (anchors are exact k-mer matches there).
    /// Returns the query and, for each mark, its (query, reference) point.
    fn transcribe(rng: &mut Rng, t: &[u8], exons: &[(usize, usize)], err: [usize; 3], marks: &[usize])
        -> (Vec<u8>, Vec<(i32, i32)>) {
        let (first, last) = (exons[0].0, exons[exons.len() - 1].1);
        let guarded = |p: usize| p < first + 16 || p + 16 >= last || marks.iter().any(|&m| p + 16 > m && p < m + 16);
        let mut q = Vec::new();
        let mut pts = Vec::new();
        for &(a, b) in exons {
            for (p, &tp) in t.iter().enumerate().take(b).skip(a) {
                if marks.contains(&p) { pts.push((q.len() as i32, p as i32)); }
                let free = !guarded(p);
                if free && rng.below(1000) < err[2] { continue; }
                let c = if free && rng.below(1000) < err[0] { (tp + 1 + rng.below(3) as u8) & 3 } else { tp };
                q.push(c);
                if free && rng.below(1000) < err[1] { q.push(rng.base()); }
            }
        }
        assert_eq!(pts.len(), marks.len());
        (q, pts)
    }

    fn cigar_str(cigar: &[u32]) -> String {
        cigar.iter().map(|&c| format!("{}{}", c >> 4, b"MIDNSHP=X"[(c & 0xf) as usize] as char)).collect()
    }

    fn introns(cigar: &[u32]) -> usize { cigar.iter().filter(|&&c| c & 0xf == 3).count() }

    /// Full spliced DP of the gap, and the compressed attempt with anchor points `marks`
    /// (gap ends added): (compressed succeeded, full result, compressed result).
    fn run(q: &[u8], t: &[u8], marks: &[(i32, i32)], strand: AlignFlags, junc: Option<&[u8]>)
        -> (bool, DpResult, DpResult) {
        let (opt, mat, flags) = cdna(strand);
        let zdrop = opt.alignment.zdrop;
        let mut full = DpResult::default();
        splice_dp(q, t, &mat, &opt, zdrop, -1, flags, junc, &mut full);
        let mut pts = vec![(0, 0)];
        pts.extend_from_slice(marks);
        pts.push((q.len() as i32, t.len() as i32));
        let mut comp = DpResult::default();
        let ok = splice_gap_compressed(q, t, &pts, &mat, &opt, zdrop, flags, junc, &mut comp);
        (ok, full, comp)
    }

    fn assert_same(ok: bool, full: &DpResult, comp: &DpResult, what: &str) {
        assert!(ok, "{what}: compressed gap fill fell back (full {})", cigar_str(&full.cigar));
        assert_eq!(cigar_str(&comp.cigar), cigar_str(&full.cigar), "{what}: CIGAR");
        assert_eq!(comp.score, full.score, "{what}: score");
        assert_eq!(comp.zdropped, 0, "{what}");
    }

    #[test]
    fn test_compressed_gap_matches_full_dp() {
        for seed in 0..12u64 {
            let mut rng = Rng(0x9e37_79b9_7f4a_7c15 ^ (seed + 1).wrapping_mul(0x2545_f491_4f6c_dd1d));
            let ex = [rng.range(80, 150), rng.range(80, 150), rng.range(80, 150)];
            let ix = [rng.range(2_000, 50_000), rng.range(2_000, 50_000)];
            let (t, iv) = gene(&mut rng, &ex, &ix);
            // One anchor point inside every exon (the gap ends lie in the first and last).
            let marks: Vec<usize> = iv.iter().map(|&(a, b)| rng.range(a + 20, b - 20)).collect();
            // cDNA (both strands' splice signals) and a single-strand pass.
            let strand = if seed % 3 == 2 { AlignFlags::SPLICE_FOR } else { both() };
            for (err, label) in [([0, 0, 0], "exact"), ([20, 5, 5], "3% errors")] {
                let (q, pts) = transcribe(&mut rng, &t, &iv, err, &marks);
                let (ok, full, comp) = run(&q, &t, &pts, strand, None);
                let what = format!("seed {seed} {label} exons {ex:?} introns {ix:?}");
                if err[0] == 0 {
                    assert_eq!(introns(&full.cigar), 2, "{what}: full DP {}", cigar_str(&full.cigar));
                }
                assert_same(ok, &full, &comp, &what);
            }
        }
    }

    #[test]
    fn test_compressed_gap_matches_full_dp_with_junctions() {
        // Annotated junctions are compressed along with the reference.
        let mut rng = Rng(0x51ed_270b_a3c1_9e55);
        let (t, iv) = gene(&mut rng, &[120, 90, 140], &[7_000, 23_000]);
        let mut junc = vec![0u8; t.len()];
        for w in iv.windows(2) {
            junc[w[0].1] |= 1; // donor (first intron base)
            junc[w[1].0 - 1] |= 2; // acceptor (last intron base)
        }
        let marks: Vec<usize> = iv.iter().map(|&(a, b)| (a + b) / 2).collect();
        let (q, pts) = transcribe(&mut rng, &t, &iv, [20, 5, 5], &marks);
        for strand in [both(), AlignFlags::SPLICE_FOR] {
            let (ok, full, comp) = run(&q, &t, &pts, strand, Some(&junc));
            assert_same(ok, &full, &comp, "annotated");
        }
    }

    #[test]
    fn test_compressed_gap_rescues_unanchored_exon() {
        // A 60 bp exon with no anchor point lies in the removed reference: the compressed
        // alignment is rough next to the widened intron, and the exon is found by its 10-mers.
        for (exons, introns_len, anchored, n_introns) in [
            (vec![120, 60, 130], vec![8_000, 12_000], vec![true, false, true], 2),
            (vec![110, 60, 90, 120], vec![6_000, 9_000, 4_000], vec![true, false, true, true], 3),
            (vec![140, 60, 100], vec![20_000, 3_000], vec![true, false, true], 2),
        ] {
            for seed in 0..4u64 {
                let mut rng = Rng(0xd1b5_4a32_d192_ed03 ^ (seed + 7).wrapping_mul(0x9e37_79b9_7f4a_7c15));
                let (t, iv) = gene(&mut rng, &exons, &introns_len);
                let marks: Vec<usize> = iv.iter().zip(&anchored).filter(|x| *x.1).map(|(&(a, b), _)| (a + b) / 2).collect();
                let (q, pts) = transcribe(&mut rng, &t, &iv, [0, 0, 0], &marks);
                let (ok, full, comp) = run(&q, &t, &pts, both(), None);
                let what = format!("exons {exons:?} seed {seed}");
                assert_eq!(introns(&full.cigar), n_introns, "{what}: full DP {}", cigar_str(&full.cigar));
                assert_same(ok, &full, &comp, &what);
                // With one substitution in the middle of the unanchored exon its 10-mers still match.
                let mut q1 = q.clone();
                let mid = (iv[0].1 - iv[0].0) + 30; // exon 2 starts at query exon1_len
                q1[mid] = (q1[mid] + 1) & 3;
                let (ok, full, comp) = run(&q1, &t, &pts, both(), None);
                assert_eq!(introns(&full.cigar), n_introns, "{what} +1 mismatch: full DP {}", cigar_str(&full.cigar));
                assert_same(ok, &full, &comp, &format!("{what} +1 mismatch"));
            }
        }
    }

    #[test]
    fn test_compressed_gap_nothing_cut() {
        let mut rng = Rng(0x0123_4567_89ab_cdef);
        // Contiguous alignment: the anchor windows cover the whole reference.
        let t = rng.seq(500);
        let (ok, _, _) = run(&t, &t, &[(250, 250)], both(), None);
        assert!(!ok, "contiguous gap: nothing to cut");
        // Two 100 bp exons and points only at the gap ends: windows of 200 + 64 bases on each
        // side, so a cut needs a reference longer than 2 * 200 + 2 * 64 + 2 * SEAM = 1040.
        for (intron, cut) in [(500, false), (840, false), (841, true), (1_500, true)] {
            let (t, iv) = gene(&mut rng, &[100, 100], &[intron]);
            assert_eq!(t.len(), 200 + intron);
            let (q, _) = transcribe(&mut rng, &t, &iv, [0, 0, 0], &[]);
            let (ok, full, comp) = run(&q, &t, &[], both(), None);
            assert_eq!(introns(&full.cigar), 1, "intron {intron}: {}", cigar_str(&full.cigar));
            if cut {
                assert_same(ok, &full, &comp, &format!("intron {intron}"));
            } else {
                assert!(!ok, "intron {intron}: reference {} too short to cut", t.len());
            }
        }
    }

    #[test]
    fn test_compressed_gap_rejects_path_through_seam() {
        // Exon 2 (no anchor point) straddles the edge of the window kept right of exon 1's
        // point: exons [0,150), [234,294), [10294,10444); points at exon 1's last and exon 3's
        // first base keep [149, 274) and [10169, 10294), so the seam starts 40 bases into exon 2.
        // The best compressed path aligns those 40 bases and then needs an intron starting at
        // (or inside) the seam, which the seam rule rejects: the caller must run the full DP.
        let mut rng = Rng(0x7f4a_7c15_9e37_79b9);
        let (t, iv) = gene(&mut rng, &[150, 60, 150], &[84, 10_000]);
        assert_eq!(iv, vec![(0, 150), (234, 294), (10_294, 10_444)]);
        let (q, pts) = transcribe(&mut rng, &t, &iv, [0, 0, 0], &[149, 10_294]);
        assert_eq!(pts, vec![(149, 149), (210, 10_294)]);
        let (ok, full, _) = run(&q, &t, &pts, both(), None);
        assert_eq!(cigar_str(&full.cigar), "150M84N60M10000N150M");
        assert!(!ok, "a path through the seam must fall back to the full DP");
    }

    #[test]
    fn test_rough_flank() {
        // Compressed reference: exon 1 (100), 20 kept bases, the seam, 20 kept bases, exon 2
        // (100); the seam replaced 5000 bases, so the widened intron is 40 + 5000 long.
        let mut rng = Rng(0xa076_1d64_78bd_642f);
        let (e1, e2) = (rng.seq(100), rng.seq(100));
        let tc: Vec<u8> = [e1.clone(), rng.seq(20), vec![4; SEAM], rng.seq(20), e2.clone()].concat();
        let seams = [(120usize, 5_000usize)];
        let n = (40 + 5_000) << 4 | 3;
        let clean = [e1.clone(), e2.clone()].concat();
        assert_eq!(rough_flank(&clean, &tc, &[100 << 4, n, 100 << 4], &seams), None);
        // FLANK_MAX mismatches within FLANK_LEN bases of the intron are tolerated, one more is not.
        for side in [0usize, 1] {
            for k in [FLANK_MAX, FLANK_MAX + 1] {
                let mut q = clean.clone();
                for i in 0..k as usize {
                    let x = if side == 0 { 99 - 4 * i } else { 100 + 4 * i };
                    q[x] = (q[x] + 1) & 3;
                }
                let want = if k > FLANK_MAX { Some((0, 1, 100)) } else { None };
                assert_eq!(rough_flank(&q, &tc, &[100 << 4, n, 100 << 4], &seams), want, "side {side}, {k} errors");
            }
        }
        // Mismatches further than FLANK_LEN from the intron do not count.
        let mut q = clean.clone();
        for b in &mut q[..10] { *b = (*b + 1) & 3; }
        assert_eq!(rough_flank(&q, &tc, &[100 << 4, n, 100 << 4], &seams), None);
        // Inserted query bases next to the intron (an exon the compressed reference lacks).
        let q: Vec<u8> = [e1.clone(), rng.seq(10), e2.clone()].concat();
        assert_eq!(rough_flank(&q, &tc, &[100 << 4, 10 << 4 | 1, n, 100 << 4], &seams), Some((0, 1, 110)));
        // An intron that does not cross a seam is not inspected.
        let tc2: Vec<u8> = [e1.clone(), rng.seq(40), e2.clone()].concat();
        assert_eq!(rough_flank(&q, &tc2, &[100 << 4, 10 << 4 | 1, 40 << 4 | 3, 100 << 4], &[]), None);
    }

    #[test]
    fn test_find_exons() {
        let mut rng = Rng(0xbf58_476d_1ce4_e5b9);
        let mut ts = rng.seq(100_000);
        let qs = rng.seq(96);
        // Exact 30 bp stretch q[30..60) planted at reference 54321.
        ts[54_321..54_351].copy_from_slice(&qs[30..60]);
        let cand = find_exons(&qs, &ts);
        assert!(!cand.is_empty() && cand.len() <= 6, "{cand:?}");
        // Largest group first: the 21 10-mers of the stretch, centred at query 40 + 5.
        assert_eq!(cand[0], (45, 54_321 - 30 + 45, 21), "{cand:?}");

        // Ten planted copies (ten diagonals): at most six candidates, all planted copies.
        let mut ts = rng.seq(100_000);
        let starts: Vec<usize> = (0..10).map(|i| 3_000 + i * 9_000).collect();
        for &s in &starts { ts[s..s + 30].copy_from_slice(&qs[30..60]); }
        let cand = find_exons(&qs, &ts);
        assert_eq!(cand.len(), 6, "{cand:?}");
        for &(qm, tm, _) in &cand {
            assert!(starts.iter().any(|&s| tm == s - 30 + qm), "({qm}, {tm}) is not a planted copy");
        }

        // Too short a query (under 16 bases), or a reference shorter than a 10-mer: nothing.
        let mut ts = rng.seq(10_000);
        ts[5_000..5_015].copy_from_slice(&qs[0..15]);
        assert!(find_exons(&qs[0..15], &ts).is_empty());
        assert!(find_exons(&qs, &qs[0..9]).is_empty());
        // N bases break 10-mers.
        let mut qn = qs[30..60].to_vec();
        for i in (5..30).step_by(9) { qn[i] = 4; }
        let mut ts = rng.seq(10_000);
        ts[4_000..4_030].copy_from_slice(&qs[30..60]);
        assert!(find_exons(&qn, &ts).is_empty());
        // Repeats: more than 4096 hits give up (the aligner then keeps the rough alignment).
        assert!(find_exons(&[0u8; 48], &[0u8; 1_000]).is_empty());
    }

    /// An anchor-free exact 114 bp exon between two anchored exons: the window kept around the
    /// rescued exon spans the query between the neighbouring anchors, so the whole exon fits
    /// (a window of RESCUE_Q query bases once cut it short and shifted the next junction).
    #[test]
    fn test_compressed_gap_rescued_exon_truncated_by_window() {
        let seed = 2u64;
        let mut rng = Rng(0x3c6e_f372_fe94_f82b ^ (seed + 1).wrapping_mul(0x2545_f491_4f6c_dd1d));
        let (t, iv) = gene(&mut rng, &[150, 114, 150], &[30_000, 20_000]);
        let (q, pts) = transcribe(&mut rng, &t, &iv, [0, 0, 0], &[75, iv[2].0 + 75]);
        let (ok, full, comp) = run(&q, &t, &pts, both(), None);
        assert_eq!(cigar_str(&full.cigar), "150M30000N114M20000N150M");
        if ok { assert_same(ok, &full, &comp, "anchor-free 114 bp exon"); }
    }

    /// Random genes (2-5 exons of 30-200 bp, introns 100 bp - 60 kb, 0 / 3% / 8% errors), every
    /// other one with a middle exon lacking an anchor point. Whenever the compressed fill
    /// succeeds it must equal the full DP. With all exons anchored it always does. With an
    /// unanchored exon the rescue needs an exact 10-mer of it within RESCUE_Q query bases of
    /// the rough intron; at ~8% errors one case of 225 (seed 77) has none and keeps a worse
    /// alignment, the accepted cost of the heuristic.
    #[test]
    fn test_compressed_gap_stress_unanchored_exons() {
        let mut bad = Vec::new();
        for seed in 0..600u64 {
            let mut rng = Rng(0x1234_5678_9abc_def1 ^ (seed + 1).wrapping_mul(0x2545_f491_4f6c_dd1d));
            let ne = rng.range(2, 5);
            let ex: Vec<usize> = (0..ne).map(|_| rng.range(30, 200)).collect();
            let ix: Vec<usize> = (0..ne - 1)
                .map(|_| if rng.below(3) == 0 { rng.range(100, 2_000) } else { rng.range(2_000, 60_000) }).collect();
            let (t, iv) = gene(&mut rng, &ex, &ix);
            let unanchored = seed % 2 == 1 && ne >= 3;
            let skip = if unanchored { 1 + rng.below(ne - 2) } else { usize::MAX };
            let marks: Vec<usize> = iv.iter().enumerate().filter(|(i, _)| *i != skip)
                .map(|(_, &(a, b))| if b - a > 40 { rng.range(a + 20, b - 20) } else { (a + b) / 2 }).collect();
            let err = [[0, 0, 0], [20, 5, 5], [50, 15, 15]][rng.below(3)];
            let (q, pts) = transcribe(&mut rng, &t, &iv, err, &marks);
            let strand = if seed % 5 == 0 { AlignFlags::SPLICE_FOR } else { both() };
            let (ok, full, comp) = run(&q, &t, &pts, strand, None);
            if ok && (comp.cigar != full.cigar || comp.score != full.score) {
                bad.push(format!("seed {seed} (unanchored exon {unanchored}, errors {err:?}): full {} {} / compressed {} {}",
                    full.score, cigar_str(&full.cigar), comp.score, cigar_str(&comp.cigar)));
            }
        }
        let anchored_bad = bad.iter().filter(|b| b.contains("unanchored exon false")).count();
        assert!(anchored_bad == 0 && bad.len() <= 1, "compressed gap fill differs from the full DP:\n{}", bad.join("\n"));
    }
}
