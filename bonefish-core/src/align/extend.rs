//! Alignment extension stage: converts chained anchors into base-level alignments.
//!
//! Entry point is `align_anchors`, which takes region-relative nt4-encoded query/target
//! sequences and a set of anchors, and produces an `AlignResult` containing CIGAR,
//! score, and coordinate mappings. The core flow for each anchor chain is:
//! left-extension, inter-anchor gap-fill, right-extension — all delegated to the
//! DP kernels in `dp.rs`. After extension, z-drop splitting (`split_zdrop`) may
//! break a single alignment into multiple sub-alignments at score dips.
//!
//! Key types: `AlignAnchorContext` bundles per-region inputs (sequences, anchors,
//! options); `AlignResult` carries the final CIGAR, scores, and query/target ranges.
//! CS, MD, and DS tag formatting functions convert CIGAR + sequences into SAM tags.

use crate::align::sketch::Minimizer;
use crate::align::dp::{self, DpResult, APPROX_MAX, RIGHT_ALIGN, EXTENSION_ONLY, REV_CIGAR,
    GENERIC_SCORING};
use crate::align::index::Index;
use crate::align::map::{AlignFlags, MapOptions};
use std::cmp;

/// Reverse complement a DNA sequence (ASCII encoding).
// Table lookup (the per-base match compiled to mispredicting branches). Every byte maps to an
// ASCII base, so the output is always valid UTF-8.
const RC: [u8; 256] = {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'T'; t[b'a' as usize] = b'T';
    t[b'C' as usize] = b'G'; t[b'c' as usize] = b'G';
    t[b'G' as usize] = b'C'; t[b'g' as usize] = b'C';
    t[b'T' as usize] = b'A'; t[b't' as usize] = b'A';
    t[b'U' as usize] = b'A'; t[b'u' as usize] = b'A';
    t
};

pub fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| RC[b as usize]).collect()
}

/// Append the reverse complement of `seq` (as `rev_comp`) to `out`.
pub fn push_rev_comp(out: &mut String, seq: &[u8]) {
    // SAFETY: RC maps every byte to an ASCII letter, so `out` stays valid UTF-8.
    let v = unsafe { out.as_mut_vec() };
    v.reserve(seq.len());
    v.extend(seq.iter().rev().map(|&b| RC[b as usize]));
}

/// Reverse complement a nt4-encoded sequence (0=A,1=C,2=G,3=T,4=N).
pub fn rev_comp_nt4(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| if b < 4 { 3 - b } else { 4 }).collect()
}

/// Convert ASCII sequence to nt4 + reverse complement in one pass.
pub fn encode_nt4_rc(ascii_seq: &[u8]) -> Vec<u8> {
    ascii_seq.iter().rev().map(|&b| {
        let nt4 = match b { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'=>3, _=>4 };
        if nt4 < 4 { 3 - nt4 } else { 4 }
    }).collect()
}

/// nt4 value to lowercase ASCII base.
const NT4_TO_LOWER: [u8; 5] = [b'a', b'c', b'g', b't', b'n'];

/// nt4 value to uppercase ASCII base.
const NT4_TO_UPPER: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

// Re-export anchor flags from sketch.rs (canonical home for Minimizer bit-fields)
pub use crate::align::sketch::{SEED_IGNORE, SEED_LONG_JOIN, SEED_TANDEM};


#[cfg(any(target_arch = "aarch64", target_arch = "x86_64", target_arch = "wasm32"))]
/// Encode a single ASCII base to nt4.
#[inline]
pub fn encode_nt4_byte(b: u8) -> u8 {
    match b { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'=>3, _=>4 }
}

/// Find the longest stretch of consecutively matching anchors.
/// For SR mode, replaces fix_bad_ends + filter_bad_seeds + adjust_minier.
/// Returns (as1_offset, cnt1) into the anchors slice.
pub fn max_stretch(anchors: &[Minimizer]) -> (usize, usize) {
    let n = anchors.len();
    if n < 2 {
        return (0, n);
    }

    let mut max_score: i32 = -1;
    let mut max_i: usize = 0;
    let mut max_len: usize = 0;

    let mut score = anchors[0].query_span();
    let mut len: usize = 1;

    for i in 1..n {
        let q_span = anchors[i].query_span();
        let lr = anchors[i].ref_pos() - anchors[i - 1].ref_pos();
        let lq = anchors[i].query_pos() - anchors[i - 1].query_pos();
        if lq == lr {
            score += std::cmp::min(lq, q_span);
            len += 1;
        } else {
            if score > max_score {
                max_score = score;
                max_len = len;
                max_i = i - len;
            }
            score = q_span;
            len = 1;
        }
    }
    if score > max_score {
        max_len = len;
        max_i = n - len;
    }
    (max_i, max_len)
}

#[derive(Debug, Clone)]
pub struct CigarOp {
    pub op: char,
    pub len: u32,
}

impl std::fmt::Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.len, self.op)
    }
}

pub struct AlignmentContext {
    /// Reusable buffer for the extracted target region.
    /// Allocated once per thread, grown as needed, never shrunk —
    /// avoids per-mapping `malloc`/`free` churn for the multi-megabyte
    /// regions produced by chromosome-scale chains (asm presets).
    /// Filled via `set_len + extract_nt4_into` to skip the zero-fill
    /// that `vec![0u8; N]` would otherwise perform.
    pub target_buf: Vec<u8>,
}

impl Default for AlignmentContext {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignmentContext {
    pub fn new() -> Self {
        AlignmentContext { target_buf: Vec::new() }
    }
}

/// Per-call context for align_anchors (values that vary per mapping).
/// Context for a single [`align_anchors`] call. Carries per-chain metadata
/// that doesn't belong in [`MapOptions`] (which is shared across all chains).
pub struct AlignAnchorContext {
    /// (rs1, qs1, re1, qe1): nearby-seed bounds from the squeezed array.
    /// Used to compute left/right extension boundaries. Region-relative.
    pub seed_bounds: (i32, i32, i32, i32),
    /// Whether this chain is on the reverse strand.
    pub rev: bool,
    /// Reference sequence ID.
    pub rid: usize,
    /// Whether this chain was produced by a z-drop split (affects zdrop threshold).
    pub split_inv: bool,
    /// K-mer size from the index.
    pub k: usize,
    /// Offset from region-relative to absolute chromosome coordinates.
    /// Set to `rgn_start` when tseq is a region extraction; 0 for full chromosome.
    pub ref_offset: usize,
}


/// Build 5x5 scoring matrix.
/// mat[i*5+j] = a if i==j (match), -b if i!=j (mismatch), -sc_ambi for N
pub fn build_scoring_matrix(a: i32, b: i32) -> [i8; 25] {
    build_scoring_matrix_full(a, b, 0, 1)
}

/// Build 5x5 scoring matrix with optional transition scoring.
/// Encoding: A=0, C=1, G=2, T=3, N=4
/// Transitions: A↔G (0↔2), C↔T (1↔3) — get -transition instead of -b when transition > 0
/// sc_ambi: ambiguity (N) penalty (default: 1, stored as positive, negated here)
pub fn build_scoring_matrix_full(a: i32, b: i32, transition: i32, sc_ambi: i32) -> [i8; 25] {
    let mut mat = [0i8; 25];
    let ambi = -(sc_ambi as i8);
    for i in 0..4 {
        for j in 0..4 {
            if i == j {
                mat[i * 5 + j] = a as i8;
            } else if transition > 0 && ((i == 0 && j == 2) || (i == 2 && j == 0) || (i == 1 && j == 3) || (i == 3 && j == 1)) {
                mat[i * 5 + j] = -(transition as i8);
            } else {
                mat[i * 5 + j] = -(b as i8);
            }
        }
        mat[i * 5 + 4] = ambi;
    }
    for j in 0..5 {
        mat[4 * 5 + j] = ambi;
    }
    mat
}

/// Collect positions of gaps longer than min_gap between consecutive anchors.
fn collect_long_gaps(anchors: &[Minimizer], min_gap: i32) -> Vec<usize> {
    let mut gaps = Vec::new();
    for i in 1..anchors.len() {
        let gap = (anchors[i].query_pos() - anchors[i-1].query_pos())
                - (anchors[i].ref_pos() - anchors[i-1].ref_pos());
        if gap < -min_gap || gap > min_gap {
            gaps.push(i);
        }
    }
    gaps
}

/// Filter bad seeds between anchors.
/// Marks seeds with SEED_IGNORE flag.
fn filter_bad_seeds(anchors: &mut [Minimizer], min_gap: i32, diff_thres: i32, max_ext_len: i32, max_ext_cnt: usize) {
    let k_vec = collect_long_gaps(anchors, min_gap);
    let n = k_vec.len();
    if n <= 1 { return; }

    let mut max_val = 0i32;
    let mut max_st: isize = -1;
    let mut max_en: isize = -1;
    let mut k = 0;
    loop {
        if k == n || k as isize >= max_en {
            if max_en > 0 {
                for anchor in &mut anchors[k_vec[max_st as usize]..k_vec[max_en as usize]] {
                    anchor.y |= SEED_IGNORE;
                }
            }
            max_val = 0;
            max_st = -1;
            max_en = -1;
            if k == n { break; }
        }
        let i = k_vec[k];
        let gap = (anchors[i].query_pos() - anchors[i-1].query_pos())
                - (anchors[i].ref_pos() - anchors[i-1].ref_pos());
        let mut n_ins = if gap > 0 { gap } else { 0 };
        let mut n_del = if gap < 0 { -gap } else { 0 };
        let qs_start = anchors[i-1].query_pos();
        let rs_start = anchors[i-1].ref_pos();
        let mut max_diff = 0i32;
        let mut max_diff_l: isize = -1;
        let end_l = cmp::min(n, k + 1 + max_ext_cnt);
        for (l, &j) in k_vec.iter().enumerate().take(end_l).skip(k+1) {
            if anchors[j].query_pos() - qs_start > max_ext_len || anchors[j].ref_pos() - rs_start > max_ext_len { break; }
            let gap2 = (anchors[j].query_pos() - anchors[j-1].query_pos())
                     - (anchors[j].ref_pos() - anchors[j-1].ref_pos());
            if gap2 > 0 { n_ins += gap2; } else { n_del += -gap2; }
            let diff = n_ins + n_del - (n_ins - n_del).abs();
            if max_diff < diff {
                max_diff = diff;
                max_diff_l = l as isize;
            }
        }
        if max_diff > diff_thres && max_diff > max_val {
            max_val = max_diff;
            max_st = k as isize;
            max_en = max_diff_l;
        }
        k += 1;
    }
}

/// Alternative bad seed filtering.
/// Marks seeds with SEED_IGNORE and SEED_LONG_JOIN flags.
fn filter_bad_seeds_alt(anchors: &mut [Minimizer], min_gap: i32, max_ext: i32) {
    let k_vec = collect_long_gaps(anchors, min_gap);
    let n = k_vec.len();
    if n == 0 { return; }

    let mut k = 0;
    while k < n {
        let i = k_vec[k];
        let gap1_raw = (anchors[i].query_pos() - anchors[i-1].query_pos())
                     - (anchors[i].ref_pos() - anchors[i-1].ref_pos());
        let mut gap1 = gap1_raw.abs();
        let mut re1 = anchors[i].ref_pos();
        let mut qe1 = anchors[i].query_pos();
        let mut l = k + 1;
        while l < n {
            let j = k_vec[l];
            if anchors[j].query_pos() - qe1 > max_ext || anchors[j].ref_pos() - re1 > max_ext { break; }
            let gap2_raw = (anchors[j].query_pos() - anchors[j-1].query_pos())
                         - (anchors[j].ref_pos() - anchors[j-1].ref_pos());
            let q_span_pre = anchors[j-1].query_span();
            let rs2 = anchors[j-1].ref_pos() + q_span_pre;
            let qs2 = anchors[j-1].query_pos() + q_span_pre;
            let m = cmp::min(rs2 - re1, qs2 - qe1);
            let gap2 = gap2_raw.abs();
            if m > gap1 + gap2 { break; }
            re1 = anchors[j].ref_pos();
            qe1 = anchors[j].query_pos();
            gap1 = gap2;
            l += 1;
        }
        if l > k + 1 {
            let end = k_vec[l - 1];
            for anchor in &mut anchors[k_vec[k]..end] {
                anchor.y |= SEED_IGNORE;
            }
            anchors[end].y |= SEED_LONG_JOIN;
        }
        k = l;
    }
}

/// Fix bad chain endpoints by trimming weak boundary alignments.
pub fn fix_bad_ends(anchors: &[Minimizer], bw: i32, min_match: i32) -> (usize, usize) {
    let cnt = anchors.len();
    let mut as_idx = 0usize;
    let mut cnt_out = cnt;

    if cnt < 3 { return (as_idx, cnt_out); }

    // Compute mlen for the chain (fuzzy matched length)
    let mlen = {
        let mut ml = anchors[0].query_span();
        for i in 1..cnt {
            let span = anchors[i].query_span();
            let tl = anchors[i].ref_pos() - anchors[i-1].ref_pos();
            let ql = anchors[i].query_pos() - anchors[i-1].query_pos();
            // Fuzzy length: tl>span && ql>span ? span : tl<ql ? tl : ql
            ml += if tl > span && ql > span { span } else if tl < ql { tl } else { ql };
        }
        ml
    };

    // Trim from left
    let mut m = anchors[0].query_span();
    let mut l = m;
    for i in 1..(cnt - 1) {
        let q_span = anchors[i].query_span();
        if anchors[i].y & SEED_LONG_JOIN != 0 { break; }
        let lr = anchors[i].ref_pos() - anchors[i-1].ref_pos();
        let lq = anchors[i].query_pos() - anchors[i-1].query_pos();
        let min_lrlq = cmp::min(lr, lq);
        let max_lrlq = cmp::max(lr, lq);
        if max_lrlq - min_lrlq > l >> 1 { as_idx = i; }
        l += min_lrlq;
        m += cmp::min(min_lrlq, q_span);
        if l >= bw << 1 || (m >= min_match && m >= bw) || m >= mlen >> 1 { break; }
    }
    cnt_out = cnt - as_idx;

    // Trim from right
    m = anchors[cnt - 1].query_span();
    l = m;
    for i in (as_idx+1..(cnt - 1)).rev() {
        let q_span = anchors[i+1].query_span();
        if anchors[i+1].y & SEED_LONG_JOIN != 0 { break; }
        let lr = anchors[i+1].ref_pos() - anchors[i].ref_pos();
        let lq = anchors[i+1].query_pos() - anchors[i].query_pos();
        let min_lrlq = cmp::min(lr, lq);
        let max_lrlq = cmp::max(lr, lq);
        if max_lrlq - min_lrlq > l >> 1 { cnt_out = i + 1 - as_idx; }
        l += min_lrlq;
        m += cmp::min(min_lrlq, q_span);
        if l >= bw << 1 || (m >= min_match && m >= bw) || m >= mlen >> 1 { break; }
    }

    (as_idx, cnt_out)
}

/// CIGAR post-processing: merge ops, fix boundaries, apply right-alignment.
/// Returns (qshift, tshift) for coordinate adjustment
fn fix_cigar(cigar: &mut Vec<u32>, qseq: &[u8], tseq: &[u8], qs: i32, rs: i32) -> (i32, i32) {
    let (mut qshift, mut tshift) = (0i32, 0i32);
    if cigar.len() <= 1 { return (qshift, tshift); }

    // Phase 1: Indel left-alignment
    let mut toff = 0usize;
    let mut qoff = 0usize;
    let mut to_shrink = false;
    let n = cigar.len();
    for k in 0..n {
        let op = cigar[k] & 0xf;
        let len = (cigar[k] >> 4) as usize;
        if len == 0 { to_shrink = true; }
        if op == 0 { // M
            toff += len;
            qoff += len;
        } else if op == 1 || op == 2 { // I or D
            if k > 0 && k < n - 1 && (cigar[k-1] & 0xf) == 0 && (cigar[k+1] & 0xf) == 0 {
                let prev_len = (cigar[k-1] >> 4) as usize;
                let mut l = 0usize;
                if op == 1 { // INS: compare query bases
                    let qi = qs as usize + qoff;
                    while l < prev_len {
                        if qi < l + 1 || qi + len < l + 1 { break; }
                        if qseq[qi - 1 - l] != qseq[qi + len - 1 - l] { break; }
                        l += 1;
                    }
                } else { // DEL: compare target bases
                    let ti = rs as usize + toff;
                    while l < prev_len {
                        if ti < l + 1 || ti + len < l + 1 { break; }
                        if tseq[ti - 1 - l] != tseq[ti + len - 1 - l] { break; }
                        l += 1;
                    }
                }
                if l > 0 {
                    cigar[k-1] -= (l as u32) << 4;
                    cigar[k+1] += (l as u32) << 4;
                    qoff -= l;
                    toff -= l;
                }
                if l == prev_len { to_shrink = true; }
            }
            if op == 1 { qoff += len; }
            else { toff += len; }
        } else if op == 3 { // N_SKIP
            toff += len;
        }
    }

    // Phase 2: Fix adjacent I/D sequences like 5I6D7I
    // Uses unconditional k += 1 to match C for-loop semantics (++k always runs)
    if cigar.len() >= 3 {
        let mut k = 0;
        while k + 2 < cigar.len() {
            let op_k = cigar[k] & 0xf;
            let op_k1 = cigar[k+1] & 0xf;
            if op_k > 0 && op_k + op_k1 == 3 { // I+D or D+I
                let mut s = [0u32; 3]; // s[0]=M, s[1]=I, s[2]=D
                let mut l = k;
                while l < cigar.len() {
                    let op = cigar[l] & 0xf;
                    if op == 1 || op == 2 || (cigar[l] >> 4) == 0 {
                        s[op as usize] += cigar[l] >> 4;
                        l += 1;
                    } else {
                        break;
                    }
                }
                if s[1] > 0 && s[2] > 0 && l - k > 2 {
                    cigar[k] = (s[1] << 4) | 1; // combined I
                    cigar[k+1] = (s[2] << 4) | 2; // combined D
                    for c in &mut cigar[k+2..l] {
                        *c &= 0xf; // zero out length
                    }
                    to_shrink = true;
                }
                k = l;
            }
            k += 1;
        }
    }

    // Phase 3: Squeeze zero-length ops and merge adjacent same ops
    if to_shrink {
        // Remove zero-length ops
        cigar.retain(|&c| (c >> 4) != 0);
        // Merge adjacent same-op entries
        let mut merged: Vec<u32> = Vec::with_capacity(cigar.len());
        for &c in cigar.iter() {
            if let Some(last) = merged.last_mut() {
                if (*last & 0xf) == (c & 0xf) {
                    *last += (c >> 4) << 4;
                } else {
                    merged.push(c);
                }
            } else {
                merged.push(c);
            }
        }
        *cigar = merged;
    }

    // Phase 4: Remove leading I or D
    if !cigar.is_empty() {
        let first_op = cigar[0] & 0xf;
        if first_op == 1 || first_op == 2 { // leading I or D
            let l = (cigar[0] >> 4) as i32;
            if first_op == 1 { // INS
                qshift = l;
            } else { // DEL
                tshift = l;
            }
            cigar.remove(0);
        }
    }

    (qshift, tshift)
}

/// Test z-drop on a CIGAR string.
/// Returns 0 if no zdrop, 1 if zdrop detected, 2 if inversion zdrop detected
fn test_zdrop(
    qseq: &[u8], tseq: &[u8], cigar: &[u32], mat: &[i8; 25],
    opt: &MapOptions,
) -> i32 {
    let q_gap = opt.scoring.gap_open;
    let e_gap = opt.scoring.gap_extend;
    let zdrop = opt.alignment.zdrop;
    let zdrop_inv = opt.alignment.zdrop_inv;
    let max_gap = opt.chaining.max_gap;
    let min_chain_score = opt.chaining.min_chain_score;
    let min_dp_max = opt.alignment.min_dp_max;
    let a = opt.scoring.match_score;
    let flags = opt.flags;

    // Exact early out: every drop below the running maximum is paid for by penalties on the
    // path after it, so max_zdrop <= the path's total penalty. If even an upper bound of that
    // (mismatched or ambiguous bases x the largest mismatch penalty, plus the gap costs) is at
    // most both thresholds, neither the z-drop nor the inversion test can fire.
    let lim = zdrop.min(zdrop_inv);
    if lim >= 0 {
        let max_pen = mat.iter().map(|&v| -(v as i32)).max().unwrap_or(0).max(0) as i64;
        let (mut ti, mut qi) = (0usize, 0usize);
        let mut neg: i64 = 0;
        for &c in cigar {
            let op = c & 0xf;
            let len = (c >> 4) as usize;
            if op == 0 {
                let n = len.min(tseq.len().saturating_sub(ti)).min(qseq.len().saturating_sub(qi));
                let mm = tseq[ti..ti + n].iter().zip(&qseq[qi..qi + n]).filter(|&(&t, &q)| (t != q) | (t > 3) | (q > 3)).count();
                neg += mm as i64 * max_pen;
                ti += len;
                qi += len;
            } else if op == 1 || op == 2 || op == 3 {
                neg += (q_gap + e_gap * len as i32) as i64;
                if op == 1 { qi += len; } else { ti += len; }
            }
            if neg > lim as i64 { break; }
        }
        if neg <= lim as i64 { return 0; }
    }

    let mut score: i32 = 0;
    let mut max: i32 = i32::MIN;
    let mut max_i: i32 = -1;
    let mut max_j: i32 = -1;
    let mut i: i32 = 0;
    let mut j: i32 = 0;
    let mut max_zdrop: i32 = 0;
    // pos[0] = (t_start, t_end), pos[1] = (q_start, q_end) of worst z-drop region
    let mut pos: [[i32; 2]; 2] = [[-1, -1], [-1, -1]];

    for &c in cigar {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        if op == 0 { // M
            // Bases inside both sequences: no bounds checks (the general loop below handles any
            // part of the run past a sequence end, where the score does not change).
            let n_in = (len.min(tseq.len() as i32 - i).min(qseq.len() as i32 - j)).max(0);
            let (tp, qp) = (&tseq[i as usize..(i + n_in) as usize], &qseq[j as usize..(j + n_in) as usize]);
            for (&tc, &qc) in tp.iter().zip(qp) {
                // SAFETY: nt4 codes are < 5, so the index is < 25.
                score += unsafe { *mat.get_unchecked(tc as usize * 5 + qc as usize) } as i32;
                if score < max {
                    let diff = ((i - max_i) - (j - max_j)).abs();
                    let z = max - score - diff * e_gap;
                    if z > max_zdrop {
                        max_zdrop = z;
                        pos[0][0] = max_i; pos[0][1] = i;
                        pos[1][0] = max_j; pos[1][1] = j;
                    }
                } else {
                    max = score;
                    max_i = i;
                    max_j = j;
                }
                i += 1;
                j += 1;
            }
            for _l in n_in..len {
                let ti = i as usize;
                let qi = j as usize;
                if ti < tseq.len() && qi < qseq.len() {
                    score += mat[tseq[ti] as usize * 5 + qseq[qi] as usize] as i32;
                }
                // update_max_zdrop inline
                if score < max {
                    let li = i - max_i;
                    let lj = j - max_j;
                    let diff = (li - lj).abs();
                    let z = max - score - diff * e_gap;
                    if z > max_zdrop {
                        max_zdrop = z;
                        pos[0][0] = max_i; pos[0][1] = i;
                        pos[1][0] = max_j; pos[1][1] = j;
                    }
                } else {
                    max = score;
                    max_i = i;
                    max_j = j;
                }
                i += 1;
                j += 1;
            }
        } else if op == 1 || op == 2 || op == 3 { // I, D, N
            score -= q_gap + e_gap * len;
            if op == 1 { j += len; }
            else { i += len; }
            // update_max_zdrop inline
            if score < max {
                let li = i - max_i;
                let lj = j - max_j;
                let diff = (li - lj).abs();
                let z = max - score - diff * e_gap;
                if z > max_zdrop {
                    max_zdrop = z;
                    pos[0][0] = max_i; pos[0][1] = i;
                    pos[1][0] = max_j; pos[1][1] = j;
                }
            } else {
                max = score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Test if there is an inversion in the most dropped region
    let q_len = pos[1][1] - pos[1][0];
    let t_len = pos[0][1] - pos[0][0];
    if !flags.intersects(AlignFlags::FOR_ONLY | AlignFlags::REV_ONLY)
        && max_zdrop > zdrop_inv
        && q_len < max_gap
        && t_len < max_gap
    {
        // Reverse complement the query region
        let mut qseq2 = vec![0u8; q_len as usize];
        for k in 0..q_len as usize {
            let c = qseq[pos[1][1] as usize - k - 1];
            qseq2[k] = if c >= 4 { 4 } else { 3 - c };
        }
        let t0 = pos[0][0] as usize;
        let inv_score = dp::lightweight_score(&qseq2, &tseq[t0..t0 + t_len as usize], mat, q_gap, e_gap);
        if inv_score >= min_chain_score * a && inv_score >= min_dp_max {
            return 2; // potential inversion
        }
    }

    if max_zdrop > zdrop { 1 } else { 0 }
}

/// Append CIGAR operations with merging.
fn append_cigar(cigar: &mut Vec<u32>, new_ops: &[u32]) {
    if new_ops.is_empty() { return; }
    if !cigar.is_empty() && (cigar.last().unwrap() & 0xf) == (new_ops[0] & 0xf) {
        // Same op at boundary — merge
        *cigar.last_mut().unwrap() += (new_ops[0] >> 4) << 4;
        cigar.extend_from_slice(&new_ops[1..]);
    } else {
        cigar.extend_from_slice(new_ops);
    }
}

pub fn fmt_cigar(ops: &[CigarOp], eqx: bool) -> String {
    // One string, digits written directly (as `CigarOp`'s Display: length then op letter).
    fn push(out: &mut String, len: u32, op: char) {
        let mut d = [0u8; 10];
        let (mut n, mut k) = (len, d.len());
        loop { k -= 1; d[k] = b'0' + (n % 10) as u8; n /= 10; if n == 0 { break; } }
        // SAFETY: ASCII digits.
        out.push_str(unsafe { std::str::from_utf8_unchecked(&d[k..]) });
        out.push(op);
    }
    let mut out = String::with_capacity(ops.len() * 5);
    if eqx {
        for op in ops { push(&mut out, op.len, op.op); }
    } else {
        // Merge = and X into M
        let mut cur: Option<(char, u32)> = None;
        for op in ops {
            let ch = if op.op == '=' || op.op == 'X' { 'M' } else { op.op };
            match cur {
                Some((c, l)) if c == ch => cur = Some((c, l + op.len)),
                Some((c, l)) => { push(&mut out, l, c); cur = Some((ch, op.len)); }
                None => cur = Some((ch, op.len)),
            }
        }
        if let Some((c, l)) = cur { push(&mut out, l, c); }
    }
    out
}

pub fn calculate_alignment_score(ops: &[CigarOp], a: i32, b: i32, q: i32, e: i32, q2: i32, e2: i32) -> i32 {
    let mut score = 0;
    for op in ops {
        match op.op {
            '=' => score += a * op.len as i32,
            'X' => score -= b * op.len as i32,
            'M' => {
                 score += a * op.len as i32; 
            },
            'I' | 'D' => {
                let cost1 = q + e * op.len as i32;
                let cost2 = q2 + e2 * op.len as i32;
                score -= std::cmp::min(cost1, cost2);
            },
            _ => {} 
        }
    }
    score
}

/// Visitor pattern for CIGAR-based tag formatting.
/// Each visitor accumulates formatted output as `walk_cigar_ops` dispatches CIGAR operations.
trait CigarTagVisitor {
    /// Called for consecutive match/mismatch bases. `qi` and `ti` are 0-based offsets
    /// into the aligned query and target regions respectively.
    fn on_aligned(&mut self, op: char, len: usize, qi: usize, ti: usize);
    /// Called for insertion (query gap) operations.
    fn on_insertion(&mut self, len: usize, qi: usize);
    /// Called for deletion operations.
    fn on_deletion(&mut self, len: usize, ti: usize);
    /// Called for intron skip (N) operations.
    fn on_intron(&mut self, len: usize, ti: usize);
    /// Flush any accumulated state and return the formatted string.
    fn finish(&mut self) -> String;
}

/// Walk CIGAR ops, dispatching to the visitor for each operation type.
fn walk_cigar_ops(ops: &[CigarOp], visitor: &mut impl CigarTagVisitor) -> String {
    let mut qi = 0usize;
    let mut ti = 0usize;
    for op in ops {
        let len = op.len as usize;
        match op.op {
            '=' | 'X' | 'M' => {
                visitor.on_aligned(op.op, len, qi, ti);
                qi += len;
                ti += len;
            }
            'I' => {
                visitor.on_insertion(len, qi);
                qi += len;
            }
            'D' => {
                visitor.on_deletion(len, ti);
                ti += len;
            }
            'N' => {
                visitor.on_intron(len, ti);
                ti += len;
            }
            _ => {}
        }
    }
    visitor.finish()
}

// ---------------------------------------------------------------------------
// MdVisitor — produces MD:Z: tag
// ---------------------------------------------------------------------------

struct MdVisitor<'a> {
    tseq: &'a [u8],
    rs: usize,
    md: String,
    match_len: u32,
}

impl<'a> MdVisitor<'a> {
    fn new(tseq: &'a [u8], rs: usize) -> Self {
        Self { tseq, rs, md: String::new(), match_len: 0 }
    }
}

impl CigarTagVisitor for MdVisitor<'_> {
    #[inline]
    fn on_aligned(&mut self, op: char, len: usize, _qi: usize, ti: usize) {
        match op {
            '=' => {
                self.match_len += len as u32;
            }
            'X' => {
                let r_idx = self.rs + ti;
                for i in 0..len {
                    self.md.push_str(&self.match_len.to_string());
                    self.match_len = 0;
                    if r_idx + i < self.tseq.len() {
                        self.md.push(Index::NT4_TO_ASCII[self.tseq[r_idx + i].min(4) as usize] as char);
                    }
                }
            }
            'M' => {
                // Should not happen if coming from align_anchors producing =/X, but handle safely
            }
            _ => {}
        }
    }

    #[inline]
    fn on_insertion(&mut self, _len: usize, _qi: usize) {
        // Ignored in MD. Ref index does not move.
    }

    #[inline]
    fn on_deletion(&mut self, len: usize, ti: usize) {
        let r_idx = self.rs + ti;
        // Deletion from Reference
        // (match_len even when 0, e.g., after a mismatch: "0^ACGT")
        self.md.push_str(&self.match_len.to_string());
        self.match_len = 0;
        self.md.push('^');
        for i in 0..len {
            if r_idx + i < self.tseq.len() {
                self.md.push(Index::NT4_TO_ASCII[self.tseq[r_idx + i].min(4) as usize] as char);
            }
        }
    }

    #[inline]
    fn on_intron(&mut self, _len: usize, _ti: usize) {
        // Intron skip: not included in MD (ref index advanced by walk_cigar_ops)
    }

    fn finish(&mut self) -> String {
        if self.match_len > 0 {
            self.md.push_str(&self.match_len.to_string());
        }
        std::mem::take(&mut self.md)
    }
}

// ---------------------------------------------------------------------------
// CsVisitor — produces cs:Z: tag
// ---------------------------------------------------------------------------

struct CsVisitor<'a> {
    qseq: &'a [u8],
    tseq: &'a [u8],
    qs: usize,
    rs: usize,
    long: bool,
    cs: String,
}

impl<'a> CsVisitor<'a> {
    fn new(qseq: &'a [u8], tseq: &'a [u8], qs: usize, rs: usize, long: bool) -> Self {
        Self { qseq, tseq, qs, rs, long, cs: String::new() }
    }

    #[inline]
    fn flush_match(&mut self, run: &mut Vec<u8>) {
        if run.is_empty() { return; }
        if self.long {
            self.cs.push('=');
            self.cs.extend(run.iter().map(|&b| b as char));
        } else {
            self.cs.push_str(&format!(":{}", run.len()));
        }
        run.clear();
    }
}

impl CigarTagVisitor for CsVisitor<'_> {
    #[inline]
    fn on_aligned(&mut self, op: char, len: usize, qi: usize, ti: usize) {
        let q_idx = self.qs + qi;
        let r_idx = self.rs + ti;
        match op {
            '=' => {
                if self.long {
                    self.cs.push('=');
                    for j in 0..len {
                        if q_idx + j < self.qseq.len() {
                            self.cs.push(NT4_TO_UPPER[self.qseq[q_idx + j].min(4) as usize] as char);
                        }
                    }
                } else {
                    self.cs.push_str(&format!(":{}", len));
                }
            }
            'X' => {
                for j in 0..len {
                    let qb = NT4_TO_LOWER[self.qseq[q_idx + j].min(4) as usize];
                    let rb = NT4_TO_LOWER[self.tseq[r_idx + j].min(4) as usize];
                    self.cs.push_str(&format!("*{}{}", rb as char, qb as char));
                }
            }
            'M' => {
                let mut run: Vec<u8> = Vec::new();
                for j in 0..len {
                    if (q_idx + j) < self.qseq.len() && (r_idx + j) < self.tseq.len() {
                        let qb = self.qseq[q_idx + j];
                        let rb = self.tseq[r_idx + j];
                        if qb == rb {
                            run.push(NT4_TO_UPPER[qb.min(4) as usize]);
                        } else {
                            self.flush_match(&mut run);
                            self.cs.push_str(&format!(
                                "*{}{}",
                                NT4_TO_LOWER[rb.min(4) as usize] as char,
                                NT4_TO_LOWER[qb.min(4) as usize] as char
                            ));
                        }
                    }
                }
                self.flush_match(&mut run);
            }
            _ => {}
        }
    }

    #[inline]
    fn on_insertion(&mut self, len: usize, qi: usize) {
        let q_idx = self.qs + qi;
        self.cs.push('+');
        for i in 0..len {
            if q_idx + i < self.qseq.len() {
                self.cs.push(NT4_TO_LOWER[self.qseq[q_idx + i].min(4) as usize] as char);
            }
        }
    }

    #[inline]
    fn on_deletion(&mut self, len: usize, ti: usize) {
        let r_idx = self.rs + ti;
        self.cs.push('-');
        for i in 0..len {
            if r_idx + i < self.tseq.len() {
                self.cs.push(NT4_TO_LOWER[self.tseq[r_idx + i].min(4) as usize] as char);
            }
        }
    }

    #[inline]
    fn on_intron(&mut self, len: usize, ti: usize) {
        let r_idx = self.rs + ti;
        if len >= 2 && r_idx + len <= self.tseq.len() {
            let d1 = NT4_TO_LOWER[self.tseq[r_idx].min(4) as usize] as char;
            let d2 = NT4_TO_LOWER[self.tseq[r_idx + 1].min(4) as usize] as char;
            let a1 = NT4_TO_LOWER[self.tseq[r_idx + len - 2].min(4) as usize] as char;
            let a2 = NT4_TO_LOWER[self.tseq[r_idx + len - 1].min(4) as usize] as char;
            self.cs.push_str(&format!("~{}{}{}{}{}", d1, d2, len, a1, a2));
        }
    }

    fn finish(&mut self) -> String {
        std::mem::take(&mut self.cs)
    }
}

// ---------------------------------------------------------------------------
// DsVisitor — produces ds:Z: tag (homology-aware indel brackets)
// ---------------------------------------------------------------------------

struct DsVisitor<'a> {
    qseq: &'a [u8],
    tseq: &'a [u8],
    qs: usize,
    rs: usize,
    q_len: usize,
    t_len: usize,
    ds: String,
}

impl<'a> DsVisitor<'a> {
    fn new(qseq: &'a [u8], tseq: &'a [u8], qs: usize, rs: usize, q_len: usize, t_len: usize) -> Self {
        Self { qseq, tseq, qs, rs, q_len, t_len, ds: String::new() }
    }
}

impl CigarTagVisitor for DsVisitor<'_> {
    #[inline]
    fn on_aligned(&mut self, op: char, len: usize, qi: usize, ti: usize) {
        let q_idx = self.qs + qi;
        let r_idx = self.rs + ti;
        match op {
            '=' => {
                self.ds.push_str(&format!(":{}", len));
            }
            'X' => {
                for j in 0..len {
                    let qb = NT4_TO_LOWER[self.qseq[q_idx + j].min(4) as usize];
                    let rb = NT4_TO_LOWER[self.tseq[r_idx + j].min(4) as usize];
                    self.ds.push_str(&format!("*{}{}", rb as char, qb as char));
                }
            }
            'M' => {
                let mut match_len = 0u32;
                for j in 0..len {
                    if (q_idx + j) < self.qseq.len() && (r_idx + j) < self.tseq.len() {
                        let qb = self.qseq[q_idx + j];
                        let rb = self.tseq[r_idx + j];
                        if qb == rb {
                            match_len += 1;
                        } else {
                            if match_len > 0 {
                                self.ds.push_str(&format!(":{}", match_len));
                                match_len = 0;
                            }
                            self.ds.push_str(&format!(
                                "*{}{}",
                                NT4_TO_LOWER[rb.min(4) as usize] as char,
                                NT4_TO_LOWER[qb.min(4) as usize] as char
                            ));
                        }
                    }
                }
                if match_len > 0 {
                    self.ds.push_str(&format!(":{}", match_len));
                }
            }
            _ => {}
        }
    }

    #[inline]
    fn on_insertion(&mut self, len: usize, qi: usize) {
        // Homology-aware insertion formatting
        let y = qi; // current query offset (relative to aligned region start)
        // Right homology: count z from 1..=len where qseq[y+len-z] == qseq[y-z]
        let mut lr: usize = 0;
        for z in 1..=len {
            if y < z { break; } // y - z < 0
            if self.qseq[self.qs + y + len - z] != self.qseq[self.qs + y - z] { break; }
            lr += 1;
        }
        // Left homology: count z from 0..len where qseq[y+len+z] == qseq[y+z]
        let mut ll: usize = 0;
        for z in 0..len {
            if y + len + z >= self.q_len { break; }
            if self.qseq[self.qs + y + len + z] != self.qseq[self.qs + y + z] { break; }
            ll += 1;
        }
        self.ds.push('+');
        write_indel_ds(&mut self.ds, len, &self.qseq[self.qs + y..], ll, lr);
    }

    #[inline]
    fn on_deletion(&mut self, len: usize, ti: usize) {
        // Homology-aware deletion formatting
        let x = ti; // current target offset (relative to aligned region start)
        // Right homology: count z from 1..=len where tseq[x+len-z] == tseq[x-z]
        let mut lr: usize = 0;
        for z in 1..=len {
            if x < z { break; } // x - z < 0
            if self.tseq[self.rs + x + len - z] != self.tseq[self.rs + x - z] { break; }
            lr += 1;
        }
        // Left homology: count z from 0..len where tseq[x+len+z] == tseq[x+z]
        let mut ll: usize = 0;
        for z in 0..len {
            if x + len + z >= self.t_len { break; }
            if self.tseq[self.rs + x + z] != self.tseq[self.rs + x + len + z] { break; }
            ll += 1;
        }
        self.ds.push('-');
        write_indel_ds(&mut self.ds, len, &self.tseq[self.rs + x..], ll, lr);
    }

    #[inline]
    fn on_intron(&mut self, len: usize, ti: usize) {
        let r_idx = self.rs + ti;
        if len >= 2 && r_idx + len <= self.tseq.len() {
            let d1 = NT4_TO_LOWER[self.tseq[r_idx].min(4) as usize] as char;
            let d2 = NT4_TO_LOWER[self.tseq[r_idx + 1].min(4) as usize] as char;
            let a1 = NT4_TO_LOWER[self.tseq[r_idx + len - 2].min(4) as usize] as char;
            let a2 = NT4_TO_LOWER[self.tseq[r_idx + len - 1].min(4) as usize] as char;
            self.ds.push_str(&format!("~{}{}{}{}{}", d1, d2, len, a1, a2));
        }
    }

    fn finish(&mut self) -> String {
        std::mem::take(&mut self.ds)
    }
}

// ---------------------------------------------------------------------------
// Public formatting functions — thin wrappers over walk_cigar_ops
// ---------------------------------------------------------------------------

pub fn fmt_md(
    ops: &[CigarOp],
    tseq: &[u8],
    rs: usize,
) -> String {
    let mut visitor = MdVisitor::new(tseq, rs);
    walk_cigar_ops(ops, &mut visitor)
}

/// Alignment result with additional statistics.
/// Output of [`align_anchors`]. All coordinates are region-relative (the caller
/// adds `ref_offset` to convert ref coords back to absolute chromosome positions).
///
/// If z-drop splitting occurred, `split_right_anchors` contains the right-side
/// anchors (also region-relative) for a separate alignment pass.
pub struct AlignResult {
    /// CIGAR operations in =/X/I/D/N format (M ops are split into =/X).
    pub cigar_ops: Vec<CigarOp>,
    /// Aligned query start (region-relative for query).
    pub query_start: usize,
    /// Aligned query end (exclusive).
    pub query_end: usize,
    /// Aligned reference start (region-relative).
    pub ref_start: usize,
    /// Aligned reference end (exclusive, region-relative).
    pub ref_end: usize,
    /// Sum of DP scores from all alignment segments (left ext + gap fills + right ext).
    pub dp_score: i32,
    /// Right-side anchors from z-drop split (region-relative ref_pos). Caller must
    /// restore to absolute coords before using in a subsequent alignment.
    pub split_right_anchors: Option<Vec<Minimizer>>,
    /// Start / end of an extension window the low-identity probe left unaligned (query
    /// coordinates as `query_start`), so callers can treat it as covered by this alignment.
    pub probe_claim_qs: Option<usize>,
    pub probe_claim_qe: Option<usize>,
    /// Offset of the split point within the chain's anchor array.
    pub split_offset_in_orig: Option<usize>,
    /// True if z-drop split detected inversion (zdrop_code == 2).
    pub split_inv: bool,
}

/// Result of anchor trimming and coordinate setup.
struct TrimResult {
    as1_offset: usize,
    cnt1: usize,
    work_anchors: Vec<Minimizer>,
    rs: i32,
    qs: i32,
    re: i32,
    qe: i32,
}

/// Trim anchors and compute alignment coordinates (fix_bad_ends + filter_bad_seeds).
/// Returns None if no valid anchors remain after trimming.
fn trim_and_prepare_anchors(
    anchors: &mut [Minimizer],
    opt: &MapOptions,
    k_half: i32,
) -> Option<TrimResult> {
    let w = opt.chaining.bandwidth;
    let min_chain_score = opt.chaining.min_chain_score;
    let max_gap = opt.chaining.max_gap;
    let no_end_flt = opt.flags.contains(AlignFlags::NO_END_FLT);
    let (off, cnt) = if no_end_flt {
        (0, anchors.len())
    } else {
        fix_bad_ends(anchors, w, min_chain_score * 2)
    };
    if cnt == 0 { return None; }

    let mut wa = anchors[off..off + cnt].to_vec();

    // Filter bad seeds
    filter_bad_seeds(&mut wa, 10, 40, max_gap >> 1, 10);
    filter_bad_seeds_alt(&mut wa, 30, max_gap >> 1);

    // Propagate SEED_IGNORE and SEED_LONG_JOIN flags back to original anchors
    for i in 0..cnt {
        let fl = wa[i].y & (SEED_IGNORE | SEED_LONG_JOIN);
        if fl != 0 {
            anchors[off + i].y |= fl;
        }
    }

    // Adjust anchor coordinates (k-mer center offset)
    let first = &wa[0];
    let last = &wa[cnt - 1];
    let (rs, qs, re, qe) = (
        first.ref_pos() - k_half,
        first.query_pos() - k_half,
        last.ref_pos() - k_half,
        last.query_pos() - k_half,
    );

    Some(TrimResult {
        as1_offset: off,
        cnt1: cnt,
        work_anchors: wa,
        rs, qs, re, qe,
    })
}

/// Dual-affine extension with the optional narrow band (`ext_bw`), wide-band rescue and
/// low-identity probe. With `probe = Some((len, min_per_kb))` and a window longer than
/// 1.5 x len, the first `len` x `len` block is aligned first: up to anti-diagonal `len` its
/// cells equal the full window's, so an early z-drop is the full window's result. If the probe
/// z-drops, or its best score is below `min_per_kb` per 1,000 probe bases (a low-identity read
/// end that the chain could not seed), the probe result is kept and the rest of the window is
/// left unaligned, or, with a third value `low_bw` > 0, a low-scoring window is aligned in full
/// with band `low_bw`; otherwise the full window is aligned.
#[allow(clippy::too_many_arguments)]
fn extend_dual_ext(qs: &[u8], ts: &[u8], mat: &[i8], q: i8, e: i8, q2: i8, e2: i8, ext_bw: i32, bw: i32,
                   zdrop: i32, end_bonus: i32, flag: i32, ext_rescue: bool, probe: Option<(i32, i32, i32)>,
                   adaptive: Option<i32>, ez: &mut DpResult) -> bool {
    // With `adaptive`, the probe and the full window use a band of that many cells per
    // anti-diagonal that follows the alignment (indel drift over a long read end) instead of
    // the static band; an indel wider than about half the window shows up as a tail shift and
    // triggers the wide-band rescue below.
    let adaptive = adaptive.filter(|&w| w < ext_bw);
    let run = |qs: &[u8], ts: &[u8], band: i32, ez: &mut DpResult| match adaptive {
        Some(w) => dp::extend_dual_affine_adaptive(qs, ts, 5, mat, q, e, q2, e2, band, w, zdrop, end_bonus, flag, ez),
        None => dp::extend_dual_affine_dispatch(qs, ts, 5, mat, q, e, q2, e2, band, zdrop, end_bonus, flag, ez),
    };
    let mut done = false;
    let mut low_stop = false;
    let mut full_bw = ext_bw;
    if let Some((plen, min_per_kb, low_bw)) = probe {
        let p = plen as usize;
        if qs.len() > p + p / 2 && ts.len() > p {
            run(&qs[..p], &ts[..p], ext_bw, ez);
            let low = (ez.max as i64) * 1000 < min_per_kb as i64 * p as i64;
            if ez.reach_end == 0 && (ez.zdropped != 0 || (low && low_bw <= 0)) {
                PROBE_STOPS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                done = true;
                low_stop = ez.zdropped == 0;
            } else {
                if low { full_bw = low_bw.min(ext_bw); PROBE_STOPS.fetch_add(1, std::sync::atomic::Ordering::Relaxed); }
                *ez = DpResult::default();
            }
        }
    }
    if !done {
        run(qs, ts, full_bw, ez);
    }
    let min_shift = adaptive.map_or(ext_bw * 3 / 4, |w| w / 2);
    if ext_rescue && ext_bw < bw && ez.zdropped != 0
        && tail_indel_shift(qs, ts, (ez.max_score_query_pos + 1).max(0) as usize,
                            (ez.max_score_target_pos + 1).max(0) as usize, min_shift, bw) {
        *ez = DpResult::default();
        dp::extend_dual_affine_dispatch(qs, ts, 5, mat, q, e, q2, e2, bw, zdrop, end_bonus, flag, ez);
        low_stop = false;
    }
    low_stop
}

pub static PROBE_STOPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);



/// Gap fill with an optional narrow first pass (band |qlen - tlen| + `margin`). A banded
/// optimum can miss a path that leaves the band and comes back (a deletion and an insertion in
/// one VNTR) without touching the band edge, so the band is widened to cover every sampled read
/// 15-mer that matches the target only off it; the result is kept unless it z-drops or its path
/// comes within `EDGE` cells of the band edge (then the regular band `bw1`).
#[allow(clippy::too_many_arguments)]
/// Gap fill: `gap_dual`, or for long gaps with `subanchor` the sub-anchored fill.
#[allow(clippy::too_many_arguments)]
fn gap_fill(qs: &[u8], ts: &[u8], mat: &[i8], q: i8, e: i8, q2: i8, e2: i8, bw1: i32, narrow: Option<i32>,
            subanchor: bool, zdrop: i32, end_bonus: i32, flag: i32, ez: &mut DpResult) {
    let sub_min_len = *crate::align::env_flags::SUB_MIN;
    if subanchor && qs.len().max(ts.len()) >= sub_min_len
        && let Some(r) = gap_subanchored(qs, ts, mat, q, e, q2, e2, bw1, narrow.unwrap_or(32), zdrop, end_bonus, flag) {
        GAP_SUB_USED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        *ez = r;
        return;
    }
    gap_dual(qs, ts, mat, q, e, q2, e2, bw1, narrow, zdrop, end_bonus, flag, ez);
}

pub static GAP_SUB_USED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Sub-anchored fill of a long gap: exact 15-mer matches unique in both sides
/// are chained collinearly (longest increasing subsequence), merged into exact segments, and
/// only the pieces between segments are aligned (narrow bands). Turns one wide-band DP over
/// kilobases into a few small ones; None if there is too little to anchor on.
#[allow(clippy::too_many_arguments)]
fn gap_subanchored(qs: &[u8], ts: &[u8], mat: &[i8], q: i8, e: i8, q2: i8, e2: i8, bw1: i32, m: i32,
                   zdrop: i32, end_bonus: i32, flag: i32) -> Option<DpResult> {
    const K: usize = 15;
    const MIN_SEG: usize = 20;
    if qs.len() < K || ts.len() < K { return None; }
    let mask: u64 = (1u64 << (2 * K)) - 1;
    let kmers = |s: &[u8]| -> Vec<(u32, u32)> { // (k-mer, start), N-free k-mers only
        let mut v = Vec::with_capacity(s.len());
        let (mut h, mut l) = (0u64, 0usize);
        for (i, &c) in s.iter().enumerate() {
            if c > 3 { l = 0; h = 0; continue; }
            h = ((h << 2) | c as u64) & mask;
            l += 1;
            if l >= K { v.push((h as u32, (i + 1 - K) as u32)); }
        }
        v
    };
    // k-mers occurring exactly once on each side
    let unique = |mut v: Vec<(u32, u32)>| -> Vec<(u32, u32)> {
        v.sort_unstable();
        let mut out = Vec::with_capacity(v.len());
        let mut i = 0;
        while i < v.len() {
            let mut j = i + 1;
            while j < v.len() && v[j].0 == v[i].0 { j += 1; }
            if j == i + 1 { out.push(v[i]); }
            i = j;
        }
        out
    };
    let (tu, qu) = (unique(kmers(ts)), unique(kmers(qs)));
    // matches (qpos, tpos), by merge of the two sorted lists
    let mut hits: Vec<(u32, u32)> = Vec::new();
    let (mut a, mut b) = (0, 0);
    while a < qu.len() && b < tu.len() {
        match qu[a].0.cmp(&tu[b].0) {
            std::cmp::Ordering::Less => a += 1,
            std::cmp::Ordering::Greater => b += 1,
            std::cmp::Ordering::Equal => { hits.push((qu[a].1, tu[b].1)); a += 1; b += 1; }
        }
    }
    if hits.len() < 2 { return None; }
    hits.sort_unstable();
    // Longest chain increasing in both coordinates (patience sorting on tpos).
    let mut tails: Vec<usize> = Vec::new(); // index of the hit ending a chain of each length
    let mut prev = vec![usize::MAX; hits.len()];
    for (i, h) in hits.iter().enumerate() {
        let pos = tails.partition_point(|&j| hits[j].1 < h.1);
        if pos > 0 { prev[i] = tails[pos - 1]; }
        if pos == tails.len() { tails.push(i); } else { tails[pos] = i; }
    }
    let mut chain = Vec::with_capacity(tails.len());
    let mut c = *tails.last()?;
    while c != usize::MAX { chain.push(hits[c]); c = prev[c]; }
    chain.reverse();
    // Merge into exact segments [q0, q1) x [t0, t1) on one diagonal, non-overlapping.
    let mut segs: Vec<(usize, usize, usize)> = Vec::new(); // (q0, t0, len)
    for &(hq, ht) in &chain {
        let (hq, ht) = (hq as usize, ht as usize);
        if let Some(last) = segs.last_mut() && hq - last.0 == ht - last.1 && hq <= last.0 + last.2 {
            last.2 = last.2.max(hq + K - last.0);
            continue;
        }
        let (mut q0, mut t0, mut len) = (hq, ht, K);
        if let Some(&(pq, pt, pl)) = segs.last() {
            let d = (pq + pl).saturating_sub(q0).max((pt + pl).saturating_sub(t0));
            if d >= len { continue; }
            q0 += d; t0 += d; len -= d;
        }
        segs.push((q0, t0, len));
    }
    segs.retain(|s| s.2 >= MIN_SEG);
    if segs.is_empty() { return None; }
    // Align the pieces between segments; concatenate.
    let mut out = DpResult::default();
    let mut cigar: Vec<u32> = Vec::new();
    let mut score = 0i32;
    let (mut qp, mut tp) = (0usize, 0usize);
    let piece = |q0: usize, q1: usize, t0: usize, t1: usize, cigar: &mut Vec<u32>, score: &mut i32| -> bool {
        let (ql, tl) = (q1 - q0, t1 - t0);
        let gap_cost = |l: i32| (q as i32 + e as i32 * l).min(q2 as i32 + e2 as i32 * l);
        if ql == 0 && tl == 0 { return true; }
        if ql == 0 { append_cigar(cigar, &[(tl as u32) << 4 | 2]); *score -= gap_cost(tl as i32); return true; }
        if tl == 0 { append_cigar(cigar, &[(ql as u32) << 4 | 1]); *score -= gap_cost(ql as i32); return true; }
        let d = (tl as i32 - ql as i32).abs();
        let bw = bw1.min(d + 300).max(64);
        let narrow = if 2 * (d + m) < (ql.min(tl) as i32) { Some(m) } else { None };
        let mut ez = DpResult::default();
        gap_dual(&qs[q0..q1], &ts[t0..t1], mat, q, e, q2, e2, bw, narrow, zdrop, end_bonus, flag, &mut ez);
        if ez.zdropped != 0 || ez.cigar.is_empty() { return false; }
        append_cigar(cigar, &ez.cigar);
        *score += ez.score;
        true
    };
    for &(q0, t0, len) in &segs {
        if !piece(qp, q0, tp, t0, &mut cigar, &mut score) { return None; }
        append_cigar(&mut cigar, &[(len as u32) << 4]);
        score += mat[0] as i32 * len as i32;
        qp = q0 + len;
        tp = t0 + len;
    }
    if !piece(qp, qs.len(), tp, ts.len(), &mut cigar, &mut score) { return None; }
    out.cigar = cigar;
    out.score = score;
    Some(out)
}

fn gap_dual(qs: &[u8], ts: &[u8], mat: &[i8], q: i8, e: i8, q2: i8, e2: i8, bw1: i32, narrow: Option<i32>,
            zdrop: i32, end_bonus: i32, flag: i32, ez: &mut DpResult) {
    use std::sync::atomic::Ordering::Relaxed;
    const EDGE: i32 = 8;
    if let Some(m) = narrow {
        let mut nb = (ts.len() as i32 - qs.len() as i32).abs() + m;
        let evidence = offband_only_kmers(qs, ts, nb - EDGE);
        // Largest |diagonal| of the evidence bins, rounded out to the 16-wide bin.
        let need = evidence.map_or(0, |(lo, hi)| (hi + 1).max(16 - lo));
        if need > 0 {
            // Widen the band to cover those k-mers' diagonals (plus the margin).
            nb = need + m;
        }
        if nb < bw1 {
            // An anti-diagonal of a band nb holds up to nb + 1 cells, processed 32 at a time: trim
            // the band to the last whole vector (31, 63, ...) when that still leaves a margin of
            // MIN_MARGIN diagonals around the gap's own shift.
            const MIN_MARGIN: i32 = 20;
            let d = (ts.len() as i32 - qs.len() as i32).abs();
            let trimmed = (nb + 1) / 32 * 32 - 1;
            if need == 0 && trimmed < nb && trimmed >= d + MIN_MARGIN { nb = trimmed; }
            dp::extend_dual_affine_dispatch(qs, ts, 5, mat, q, e, q2, e2, nb, zdrop, end_bonus, flag, ez);
            if ez.zdropped == 0 && !ez.cigar.is_empty() && path_max_offdiag(&ez.cigar) < nb - EDGE {
                GAP_NARROW_USED.fetch_add(1, Relaxed);
                return;
            }
            GAP_NARROW_EDGE.fetch_add(1, Relaxed);
            *ez = DpResult::default();
        } else {
            GAP_NARROW_KMER.fetch_add(1, Relaxed);
        }
        if need > 0 { GAP_NARROW_WIDENED.fetch_add(1, Relaxed); }
    }
    dp::extend_dual_affine_dispatch(qs, ts, 5, mat, q, e, q2, e2, bw1, zdrop, end_bonus, flag, ez);
}

/// Largest |target_pos - query_pos| along a CIGAR path (ops len<<4|op; 0 M, 1 I, 2 D).
fn path_max_offdiag(cigar: &[u32]) -> i32 {
    let (mut d, mut m) = (0i32, 0i32);
    for &c in cigar {
        let len = (c >> 4) as i32;
        match c & 0xf {
            1 => { d -= len; }
            2 | 3 => { d += len; }
            _ => {}
        }
        m = m.max(d.abs());
    }
    m
}

thread_local! {
    /// Scratch for `offband_only_kmers`: target k-mers, off-band-only hits, in-band hit positions.
    static GAP_KMERS: std::cell::RefCell<(Vec<u64>, Vec<u64>, Vec<i32>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new())) };
}

/// Widest diagonal the band must cover: sampled read 15-mers (every 4th position) that occur in
/// `t` only off the band |t_pos - q_pos| <= `w` are binned by the diagonal of their nearest
/// occurrence (16-wide bins), and a bin with at least 3 of them is a path the band would miss
/// (a run of read sequence matching the target at that shift, as in a VNTR deletion plus
/// insertion), unless in-band hits outnumber it 2:1 over the same query stretch. Isolated
/// off-band hits, typically another copy of a repeat matched where the in-band copy has a
/// sequencing error, do not count. Returns the lowest and highest diagonal (t - q) covered by
/// such bins, if any. Sequences are nt4 (0-3; 4 = N).
fn offband_only_kmers(q: &[u8], t: &[u8], w: i32) -> Option<(i32, i32)> {
    const K: usize = 15;
    const MIN_HITS: usize = 3;
    if q.len() < K || t.len() < K {
        return None;
    }
    let mask: u64 = (1u64 << (2 * K)) - 1;
    GAP_KMERS.with(|cell| {
        let mut buf = cell.borrow_mut();
        let (tk, off, inb) = &mut *buf;
        off.clear();
        inb.clear();
        // Target k-mers in an open-addressing table of (k-mer << 32 | position) entries.
        const EMPTY: u64 = u64::MAX;
        let n_k = t.len() + 1 - K;
        let bits = (2 * n_k).next_power_of_two().trailing_zeros().max(4);
        let tmask = (1usize << bits) - 1;
        let slot = |kmer: u64| (kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - bits)) as usize;
        tk.clear();
        tk.resize(1 << bits, EMPTY);
        let (mut h, mut l) = (0u64, 0usize);
        for (i, &c) in t.iter().enumerate() {
            if c > 3 { l = 0; h = 0; continue; }
            h = ((h << 2) | c as u64) & mask;
            l += 1;
            if l >= K {
                let mut j = slot(h);
                while tk[j] != EMPTY { j = (j + 1) & tmask; }
                tk[j] = h << 32 | (i + 1 - K) as u64;
            }
        }
        let (mut h, mut l) = (0u64, 0usize);
        for (i, &c) in q.iter().enumerate() {
            if c > 3 { l = 0; h = 0; continue; }
            h = ((h << 2) | c as u64) & mask;
            l += 1;
            if l < K || (i + 1 - K) % 4 != 0 { continue; }
            let qp = (i + 1 - K) as i32;
            // Nearest occurrence; at equal distance the lower target position (as when the
            // occurrences were scanned in position order).
            let (mut nearest, mut nd) = (i32::MAX, 0i32);
            let mut j = slot(h);
            while tk[j] != EMPTY {
                let x = tk[j];
                if (x >> 32) == h {
                    let d = (x & 0xffff_ffff) as i32 - qp;
                    if d.abs() < nearest || (d.abs() == nearest && d < nd) { nearest = d.abs(); nd = d; }
                }
                j = (j + 1) & tmask;
            }
            if nearest == i32::MAX { continue; }
            if nearest <= w {
                inb.push(qp);
            } else {
                // (diagonal bin, query position) of an off-band-only hit
                off.push((((nd >> 4) + (1 << 30)) as u64) << 32 | qp as u64);
            }
        }
        if off.len() < MIN_HITS {
            return None;
        }
        off.sort_unstable();
        let mut range: Option<(i32, i32)> = None;
        let mut i = 0;
        while i < off.len() {
            let bin = off[i] >> 32;
            let mut j = i;
            while j < off.len() && off[j] >> 32 == bin { j += 1; }
            if j - i >= MIN_HITS {
                // Only a shifted path that the read supports about as well as the in-band one
                // over the same stretch can be optimal: ignore bins outnumbered 2:1 by in-band hits.
                let (q0, q1) = ((off[i] & 0xffff_ffff) as i32, (off[j - 1] & 0xffff_ffff) as i32);
                let n_in = inb.partition_point(|&x| x <= q1) - inb.partition_point(|&x| x < q0);
                if 2 * (j - i) >= n_in {
                    let d = (bin as i64 - (1 << 30)) as i32; // diagonals 16d .. 16d + 15
                    let (lo, hi) = range.unwrap_or((i32::MAX, i32::MIN));
                    range = Some((lo.min(16 * d), hi.max(16 * d + 15)));
                }
            }
            i = j;
        }
        range
    })
}

pub static GAP_NARROW_USED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static GAP_NARROW_KMER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static GAP_NARROW_EDGE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static GAP_NARROW_WIDENED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Does the unaligned tail of a z-dropped extension still match the target at one consistent
/// diagonal shift wider than `min_shift` (up to `max_shift`)? That is the signature of an indel
/// too wide for a narrowed extension band, as opposed to divergent or junk sequence. Compares
/// 15-mers of the first `TAIL` query bases past (qs, ts) against the target window and looks
/// for >= 4 hits in one 8-wide diagonal bucket. Sequences are nt4 (0-3; 4 = N).
fn tail_indel_shift(q: &[u8], t: &[u8], qs: usize, ts: usize, min_shift: i32, max_shift: i32) -> bool {
    const K: usize = 15;
    const TAIL: usize = 400;
    const MAX_OCC: usize = 8;
    if qs >= q.len() || ts >= t.len() {
        return false;
    }
    let span = TAIL + max_shift as usize;
    let qw = &q[qs..q.len().min(qs + span)];
    let tw = &t[ts..t.len().min(ts + span)];
    if qw.len() < K + 20 || tw.len() < K + 20 {
        return false;
    }
    let mask: u32 = (1u32 << (2 * K)) - 1;
    let kmers = |s: &[u8]| -> Vec<(usize, u32)> {
        let mut out = Vec::with_capacity(s.len());
        let (mut h, mut l) = (0u32, 0usize);
        for (i, &c) in s.iter().enumerate() {
            if c > 3 { l = 0; h = 0; continue; }
            h = ((h << 2) | c as u32) & mask;
            l += 1;
            if l >= K { out.push((i + 1 - K, h)); }
        }
        out
    };
    let mut index: std::collections::HashMap<u32, Vec<i32>> = std::collections::HashMap::new();
    for (j, h) in kmers(tw) {
        index.entry(h).or_default().push(j as i32);
    }
    let nb = (2 * max_shift as usize) / 8 + 2;
    let mut buckets = vec![0u16; nb];
    for (i, h) in kmers(qw).into_iter().take_while(|&(i, _)| i < TAIL) {
        if let Some(js) = index.get(&h) {
            if js.len() > MAX_OCC { continue; }   // repeat k-mer: no positional information
            for &j in js {
                let d = j - i as i32;
                if d.abs() <= max_shift {
                    buckets[((d + max_shift) / 8) as usize] += 1;
                }
            }
        }
    }
    let (best, cnt) = buckets.iter().enumerate().max_by_key(|&(_, &c)| c).map(|(b, &c)| (b, c)).unwrap();
    let d = (best as i32) * 8 + 4 - max_shift;
    cnt >= 4 && d.abs() > min_shift
}

/// Compute left extension boundary (rs0, qs0) for non-SR mode.
/// `orig_anchor_rpos`/`orig_anchor_qpos`: raw coordinates of original chain start anchor.
/// `seed_r_bound`/`seed_q_bound`: rs1/qs1 from nearby-seed inspection.
/// Returns (rs0, qs0).
fn compute_left_boundary(
    rs: i32, qs: i32,
    orig_anchor_rpos: i32, orig_anchor_qpos: i32, orig_q_span: i32,
    seed_r_bound: i32, seed_q_bound: i32,
    opt: &MapOptions,
) -> (i32, i32) {
    let max_gap = crate::align::env_flags::EXP_MAX_EXT.map_or(opt.chaining.max_gap, |m| m.min(opt.chaining.max_gap));
    let a = opt.scoring.match_score;
    let q = opt.scoring.gap_open;
    let e = opt.scoring.gap_extend;
    let mut rs0 = cmp::max(0, orig_anchor_rpos + 1 - orig_q_span);
    let mut qs0 = orig_anchor_qpos + 1 - orig_q_span;
    let (mut rs1_bound, mut qs1_bound) = if seed_r_bound > 0 || seed_q_bound > 0 {
        (seed_r_bound, seed_q_bound)
    } else {
        (0i32, 0i32)
    };

    if qs > 0 && rs > 0 {
        // Extend boundary by gap allowance
        let mut l = cmp::min(qs, max_gap);
        qs1_bound = cmp::max(qs1_bound, qs - l);
        qs0 = cmp::min(qs0, qs1_bound);
        l += if l * a > q { (l * a - q) / e } else { 0 };
        l = cmp::min(l, max_gap);
        l = cmp::min(l, rs);
        rs1_bound = cmp::max(rs1_bound, rs - l);
        rs0 = cmp::min(rs0, rs1_bound);
        rs0 = cmp::min(rs0, rs);
    } else {
        rs0 = rs;
        qs0 = qs;
    }
    (cmp::max(0, rs0), cmp::max(0, qs0))
}

/// Compute right extension boundary (re0, qe0) for non-SR mode.
/// `orig_anchor_rpos`/`orig_anchor_qpos`: raw coordinates of original chain end anchor.
/// `seed_r_bound`/`seed_q_bound`: re1/qe1 from nearby-seed inspection.
/// Returns (re0, qe0).
fn compute_right_boundary(
    re: i32, qe: i32,
    orig_anchor_rend: i32, orig_anchor_qend: i32,
    seed_r_bound: i32, seed_q_bound: i32,
    opt: &MapOptions,
    qlen: i32, tlen: i32,
    _seed_bounds: (i32, i32, i32, i32),
) -> (i32, i32) {
    let max_gap = crate::align::env_flags::EXP_MAX_EXT.map_or(opt.chaining.max_gap, |m| m.min(opt.chaining.max_gap));
    let a = opt.scoring.match_score;
    let q = opt.scoring.gap_open;
    let e = opt.scoring.gap_extend;
    let mut re0 = orig_anchor_rend;
    let mut qe0 = orig_anchor_qend;
    let (mut re1_bound, mut qe1_bound) = if seed_r_bound < tlen || seed_q_bound < qlen {
        (seed_r_bound, seed_q_bound)
    } else {
        (tlen, qlen)
    };

    let mut l = cmp::min(qlen - qe, max_gap);
    qe1_bound = cmp::min(qe1_bound, qe + l);
    qe0 = cmp::max(qe0, qe1_bound);
    l += if l * a > q { (l * a - q) / e } else { 0 };
    l = cmp::min(l, max_gap);
    l = cmp::min(l, tlen - re);
    re1_bound = cmp::min(re1_bound, re + l);
    re0 = cmp::max(re0, re1_bound);

    (cmp::min(re0, tlen), cmp::min(qe0, qlen))
}

/// CIGAR post-processing: fix_cigar + convert M→=/X + condense adjacent ops.
/// Returns (final_cigar_ops, adjusted_qs, adjusted_rs).
fn finalize_cigar(
    mut raw_cigar: Vec<u32>,
    qseq: &[u8],
    tseq: &[u8],
    qs: i32,
    rs: i32,
) -> (Vec<CigarOp>, i32, i32) {
    let mut final_qs = qs;
    let mut final_rs = rs;
    let (qshift, tshift) = fix_cigar(&mut raw_cigar, qseq, tseq, final_qs, final_rs);
    if qshift > 0 {
        final_qs += qshift;
    }
    if tshift > 0 {
        final_rs += tshift;
    }

    let cigar_ops = convert_cigar_to_eqx(&raw_cigar, qseq, tseq, final_qs as usize, final_rs as usize);

    // Condense adjacent same-type ops
    let mut condensed: Vec<CigarOp> = Vec::new();
    for op in cigar_ops {
        if let Some(last) = condensed.last_mut() {
            if last.op == op.op {
                last.len += op.len;
            } else {
                condensed.push(op);
            }
        } else {
            condensed.push(op);
        }
    }

    (condensed, final_qs, final_rs)
}

/// Align a chain of anchors to produce CIGAR operations.
///
/// This is the monolithic alignment entry point.
/// It performs left extension, per-gap fill, and right extension using DP
/// kernels, with z-drop splitting and anchor trimming.
///
/// # Inputs
/// - `anchors`: sorted by ref_pos, **region-relative** coordinates (caller subtracts
///   `rgn_start` before calling). May be mutated in place (SEED_IGNORE/SEED_LONG_JOIN
///   flags set during anchor trimming).
/// - `qseq`: nt4-encoded query sequence (0=A,1=C,2=G,3=T,4=N), on the correct strand.
/// - `tseq`: nt4-encoded target region (region-relative, typically 10-200 KB).
/// - `opt`: scoring/chaining/alignment parameters.
/// - `call`: per-chain context (seed_bounds, junction db, ref_offset for absolute coords).
///
/// # Outputs
/// [`AlignResult`] with region-relative coordinates. The caller adds `ref_offset` to
/// convert `ref_start`/`ref_end` back to absolute chromosome positions. Split anchors
/// (if any) also have region-relative ref_pos.
/// `lazy_extract` parameter (`Option<(mi, rid, rgn_start)>`): when `Some`,
/// segments of `tseq` are filled on demand from `mi` (each DP call fills its
/// own segment; post-DP processing fills the final alignment range). When
/// `None`, `tseq` is assumed pre-filled (qstrand+reverse path, where the
/// buffer holds a rev-comp'd full chromosome).
pub fn align_anchors(
    anchors: &mut [Minimizer],
    qseq: &[u8],
    tseq: &mut [u8],
    lazy_extract: Option<(&crate::align::index::Index, usize, usize)>,
    opt: &MapOptions,
    call: &AlignAnchorContext,
) -> AlignResult {
    // Local bindings from opt/call to avoid modifying the function body
    let a = opt.scoring.match_score;
    let b = opt.scoring.mismatch_penalty;
    let q = opt.scoring.gap_open;
    let e = opt.scoring.gap_extend;
    let q2 = opt.scoring.gap_open2;
    let e2 = opt.scoring.gap_extend2;
    let w = opt.chaining.bandwidth;
    let w_long = opt.chaining.bandwidth_long;
    let zdrop = opt.alignment.zdrop;
    let zdrop_inv = opt.alignment.zdrop_inv;
    let end_bonus = opt.alignment.end_bonus;
    let _k = call.k;
    let seed_bounds = call.seed_bounds;
    let min_cnt = opt.chaining.min_cnt;
    let max_sw_mat = opt.alignment.max_sw_mat;
    let transition = opt.scoring.transition;
    let sc_ambi = opt.scoring.ambig_penalty;
    let min_dp_len = opt.alignment.min_dp_len;
    let ref_offset = call.ref_offset as i32;
    let split_inv = call.split_inv;

    // --- Function body (formerly align_anchors_full) ---
    let empty = AlignResult { cigar_ops: Vec::new(), query_start: 0, query_end: 0, ref_start: 0, ref_end: 0, dp_score: 0, split_right_anchors: None, split_offset_in_orig: None, split_inv: false, probe_claim_qs: None, probe_claim_qe: None };
    // Query ends of extension windows the probe left unaligned as low-identity (if any).
    let mut probe_claim_qs: Option<i32> = None;
    let mut probe_claim_qe: Option<i32> = None;
    if anchors.is_empty() { return empty; }


    let mat = build_scoring_matrix_full(a, b, transition, sc_ambi);

    // Bandwidth scaling
    let bw = (w as f64 * 1.5 + 1.0) as i32;
    let mut bw_long = (w_long as f64 * 1.5 + 1.0) as i32;
    if let Some(cap) = *crate::align::env_flags::EXP_GAP_BW { bw_long = bw_long.min(cap); }
    if bw_long < bw { bw_long = bw; }
    let ext_bw = opt.alignment.ext_bw.or(*crate::align::env_flags::EXP_EXT_BW).unwrap_or(bw);
    let ext_rescue = opt.alignment.ext_rescue || *crate::align::env_flags::EXP_EXT_RESCUE;
    let ext_probe = opt.alignment.ext_probe;
    let ext_adapt = opt.alignment.ext_adaptive.or(*crate::align::env_flags::EXP_EXT_ADAPT);

    let qlen = qseq.len() as i32;
    let tlen = tseq.len() as i32;
    let k_half = (_k >> 1) as i32;

    // Generic scoring flag: needed when transition != b (non-uniform mismatch penalties)
    let generic_sc_flag = if transition != 0 && transition != b { GENERIC_SCORING } else { 0 };

    // --- Steps 1-3: Anchor trimming + coordinate computation ---

    // Lazy fill: `tseq` arrives uninitialized and we fill segments on demand
    // before each DP call + once for the final alignment range before CIGAR
    // finalization. When `lazy_extract` is `None` (qstrand+reverse path), the
    // caller has already filled `tseq`.
    let use_lazy = lazy_extract.is_some();

    let trim = match trim_and_prepare_anchors(
        anchors, opt, k_half,
    ) {
        Some(t) => t,
        None => return empty,
    };
    let TrimResult { as1_offset, cnt1, work_anchors, rs, qs, re, qe } = trim;



    // Raw CIGAR accumulator (u32 format: len<<4 | op)
    let mut raw_cigar: Vec<u32> = Vec::new();
    let mut dp_score: i32 = 0;

    // --- Step 4: Compute extension boundaries ---
    let (rs0, qs0);
    {
        let orig_first = &anchors[0];
        let (rs0_tmp, qs0_tmp) = compute_left_boundary(
            rs, qs,
            orig_first.ref_pos(), orig_first.query_pos(), orig_first.query_span(),
            seed_bounds.0, seed_bounds.1,
            opt,
        );
        rs0 = rs0_tmp;
        qs0 = qs0_tmp;
    }

    // SEED_SELF: clamp leftward extension distance to |q-r| of the chain's
    // FIRST original anchor (pre-trim) in absolute coords. `pre_qs` is in
    // absolute alignment-space query coords. `pre_rs` (after the earlier
    // rebase) is in region-relative ref coords; shift by `ref_offset` so
    // max_ext is computed in matching coord systems.
    let (mut rs0, mut qs0) = (rs0, qs0);
    if !anchors.is_empty() && (anchors[0].y & crate::align::map::SEED_SELF) != 0 {
        let pre_q_span = anchors[0].query_span();
        let pre_qs = anchors[0].query_pos() + 1 - pre_q_span;
        let pre_rs = anchors[0].ref_pos() + 1 - pre_q_span;
        let pre_rs_abs = pre_rs + ref_offset;
        let max_ext = (pre_qs - pre_rs_abs).unsigned_abs() as i32;
        if pre_rs - rs0 > max_ext { rs0 = pre_rs - max_ext; }
        if pre_qs - qs0 > max_ext { qs0 = pre_qs - max_ext; }
    }


    // --- Step 5: Left extension ---
    let (mut rs1, mut qs1) = (rs, qs);
    if qs > 0 && rs > 0 && qs > qs0 && rs > rs0 {
        let q_ext_len = (qs - qs0) as usize;

        // Lazy fill the left-extension segment [rs0, rs) in tseq.
        if use_lazy && let Some((mi, rid, rgn_start)) = lazy_extract {
            let a = rs0 as usize;
            let b = rs as usize;
            mi.extract_nt4_into(rid, rgn_start + a, rgn_start + b, &mut tseq[a..b]);
        }

        // Reverse sequences for left extension
        let mut q_rev: Vec<u8> = qseq[qs0 as usize..qs as usize].to_vec();
        let mut t_rev: Vec<u8> = tseq[rs0 as usize..rs as usize].to_vec();
        q_rev.reverse();
        t_rev.reverse();

        // Left extension flags: EXTZ_ONLY | RIGHT | REV_CIGAR
        let flag = EXTENSION_ONLY | RIGHT_ALIGN | REV_CIGAR | generic_sc_flag;

        // Left extension zdrop: use zdrop_inv for split_inv chains
        let ext_zdrop = crate::align::env_flags::EXP_EXT_ZDROP.unwrap_or(zdrop);
        let left_zdrop = if split_inv { zdrop_inv } else { ext_zdrop };

        let mut ez = DpResult::default();
        #[cfg(feature = "census")]
        let census_t0 = dp::census::start();
        // max_sw_mat check: skip alignment if matrix too large
        let sw_mat_size = q_rev.len() as i64 * t_rev.len() as i64;
        if max_sw_mat > 0 && sw_mat_size > max_sw_mat {
            ez.max_score_query_pos = -1;
            ez.max_score_target_pos = -1;
            ez.max_query_end_target_pos = -1;
            ez.zdropped = 1;
        } else if q == q2 && e == e2 {
            dp::extend_single_affine(&q_rev, &t_rev, 5, &mat, q as i8, e as i8, bw, left_zdrop, end_bonus, flag, &mut ez);
        } else {
            if extend_dual_ext(&q_rev, &t_rev, &mat, q as i8, e as i8, q2 as i8, e2 as i8, ext_bw, bw, left_zdrop, end_bonus, flag, ext_rescue, ext_probe, ext_adapt, &mut ez) {
                probe_claim_qs = Some(qs0);
            }
        }
        #[cfg(feature = "census")]
        dp::census::finish(census_t0, dp::census::Kind::LeftExt, &q_rev, &t_rev,
            &dp::census::Params { mat: &mat, q: q as i8, e: e as i8, q2: q2 as i8, e2: e2 as i8, bandwidth: bw, zdrop: left_zdrop, end_bonus, flags: flag }, &ez);

        if !ez.cigar.is_empty() {
            append_cigar(&mut raw_cigar, &ez.cigar);
            dp_score += ez.max;
        }
        // Compute consumed bases
        rs1 = rs - if ez.reach_end != 0 { ez.max_query_end_target_pos + 1 } else { ez.max_score_target_pos + 1 };
        qs1 = qs - if ez.reach_end != 0 { q_ext_len as i32 } else { ez.max_score_query_pos + 1 };
    }
    let (mut re1, mut qe1) = (rs, qs);

    // Gap filling: iterate seeds, accumulate, align when appropriate
    let mut gap_rs = rs; // current gap-fill start (reference)
    let mut gap_qs = qs; // current gap-fill start (query)
    let mut dropped = false;
    let mut split_right_anchors: Option<Vec<Minimizer>> = None;
    let mut split_offset_in_orig: Option<usize> = None;
    let mut split_inv_flag = false;

    for i in 1..cnt1 {
        if (work_anchors[i].y & (SEED_IGNORE | SEED_TANDEM)) != 0 && i != cnt1 - 1 {
            continue;
        }

        // Anchor coordinate computation
        let (anchor_re, anchor_qe) = (work_anchors[i].ref_pos() - k_half,
                                      work_anchors[i].query_pos() - k_half);
        re1 = anchor_re;
        qe1 = anchor_qe;

        let is_last = i == cnt1 - 1;
        let is_long_join = (work_anchors[i].y & SEED_LONG_JOIN) != 0;
        let should_align = is_last || is_long_join
            || (anchor_qe - gap_qs >= min_dp_len && anchor_re - gap_rs >= min_dp_len);
        // gap-fill always uses end_bonus = -1
        let gap_end_bonus: i32 = -1;

        if should_align {
            let mut bw1 = bw_long;
            if is_long_join {
                bw1 = cmp::max(anchor_qe - gap_qs, anchor_re - gap_rs);
            }
            if let Some(margin) = opt.alignment.gap_bw_margin.or(*crate::align::env_flags::EXP_GAP_BW_ADAPT) {
                let diff = ((anchor_qe - gap_qs) - (anchor_re - gap_rs)).abs();
                bw1 = cmp::max(opt.alignment.gap_bw_floor.unwrap_or(bw), cmp::min(bw1, diff + margin));
            }
            // Narrow first pass only where it saves real work: band well inside the matrix.
            let gap_narrow = opt.alignment.gap_bw_narrow.filter(|&m| {
                let (gq, gt) = (anchor_qe - gap_qs, anchor_re - gap_rs);
                let nb = (gq - gt).abs() + m;
                nb < bw1 && 2 * nb < gq.min(gt)
            });

            let q_start = cmp::max(0, gap_qs) as usize;
            let q_end = cmp::min(qlen, anchor_qe) as usize;
            let t_start = cmp::max(0, gap_rs) as usize;
            let t_end = cmp::min(tlen, anchor_re) as usize;

            // Lazy fill the gap segment [t_start, t_end) in tseq.
            if use_lazy && t_end > t_start && let Some((mi, rid, rgn_start)) = lazy_extract {
                mi.extract_nt4_into(rid, rgn_start + t_start, rgn_start + t_end, &mut tseq[t_start..t_end]);
            }

            if q_end > q_start && t_end > t_start {
                let q_sub = &qseq[q_start..q_end];
                let t_sub = &tseq[t_start..t_end];

                let mut ez = DpResult::default();
                let mut zdrop_code = 0i32;
                // max_sw_mat check: skip alignment if matrix too large
                // Reset: max_q=max_t=mqe_t=mte_q=-1, max=0, zdropped=1
                let sw_mat_size = q_sub.len() as i64 * t_sub.len() as i64;
                if max_sw_mat > 0 && sw_mat_size > max_sw_mat {
                    ez.max_score_query_pos = -1;
                    ez.max_score_target_pos = -1;
                    ez.max_query_end_target_pos = -1;
                    ez.zdropped = 1;
                } else if (opt.alignment.ungapped_fill || *crate::align::env_flags::EXP_UNGAPPED)
                    && q_sub.len() == t_sub.len() {
                    // Ungapped alignment, only when query and reference gap lengths match (no indels).
                    // Try ungapped score first; only use gapped alignment if ungapped is poor
                    let gap_len = q_sub.len() as i32;
                    let max_gapped_score = (gap_len - 2) * a - 2 * (q + e);
                    let mut ungapped_score = 0i32;
                    for j in 0..q_sub.len() {
                        if q_sub[j] >= 4 || t_sub[j] >= 4 {
                            ungapped_score += if sc_ambi > 0 { -sc_ambi } else { sc_ambi };
                        } else if q_sub[j] == t_sub[j] {
                            ungapped_score += a;
                        } else {
                            ungapped_score -= b;
                        }
                    }
                    if ungapped_score > max_gapped_score {
                        // Emit single M op (will be converted to =/X later)
                        ez.score = ungapped_score;
                        ez.cigar = vec![((gap_len as u32) << 4)]; // M op
                    } else {
                        // Fall through to gapped alignment
                        if q == q2 && e == e2 {
                            dp::extend_single_affine(q_sub, t_sub, 5, &mat, q as i8, e as i8, bw1, zdrop, gap_end_bonus, APPROX_MAX | generic_sc_flag, &mut ez);
                        } else {
                            gap_fill(q_sub, t_sub, &mat, q as i8, e as i8, q2 as i8, e2 as i8, bw1, gap_narrow, opt.alignment.gap_subanchor, zdrop, gap_end_bonus, APPROX_MAX | generic_sc_flag, &mut ez);
                        }
                    }
                    // Z-drop test applies to both paths
                    zdrop_code = if !ez.cigar.is_empty() && !dp::speculating() {
                        test_zdrop(q_sub, t_sub, &ez.cigar, &mat, opt)
                    } else { 0 };
                    if zdrop_code != 0 {
                        ez = DpResult::default();
                        let zdrop2 = if zdrop_code == 2 { zdrop_inv } else { zdrop };
                        if q == q2 && e == e2 {
                            dp::extend_single_affine(q_sub, t_sub, 5, &mat, q as i8, e as i8, bw1, zdrop2, gap_end_bonus, generic_sc_flag, &mut ez);
                        } else {
                            dp::extend_dual_affine_dispatch(q_sub, t_sub, 5, &mat, q as i8, e as i8, q2 as i8, e2 as i8, bw1, zdrop2, gap_end_bonus, generic_sc_flag, &mut ez);
                        }
                    }
                } else {
                    // First pass: APPROX_MAX
                    {
                        #[cfg(feature = "census")]
                        let census_t0 = dp::census::start();
                        if q == q2 && e == e2 {
                            dp::extend_single_affine(q_sub, t_sub, 5, &mat, q as i8, e as i8, bw1, zdrop, gap_end_bonus, APPROX_MAX | generic_sc_flag, &mut ez);
                        } else {
                            gap_fill(q_sub, t_sub, &mat, q as i8, e as i8, q2 as i8, e2 as i8, bw1, gap_narrow, opt.alignment.gap_subanchor, zdrop, gap_end_bonus, APPROX_MAX | generic_sc_flag, &mut ez);
                        }
                        #[cfg(feature = "census")]
                        dp::census::finish(census_t0, dp::census::Kind::Gap, q_sub, t_sub,
                            &dp::census::Params { mat: &mat, q: q as i8, e: e as i8, q2: q2 as i8, e2: e2 as i8, bandwidth: bw1, zdrop, end_bonus: gap_end_bonus, flags: APPROX_MAX | generic_sc_flag }, &ez);
                    }

                    // Z-drop test
                    zdrop_code = if !ez.cigar.is_empty() && !dp::speculating() {
                        test_zdrop(q_sub, t_sub, &ez.cigar, &mat, opt)
                    } else { 0 };
                    if zdrop_code != 0 {
                        ez = DpResult::default();
                        let zdrop2 = if zdrop_code == 2 { zdrop_inv } else { zdrop };
                        #[cfg(feature = "census")]
                        let census_t0 = dp::census::start();
                        if q == q2 && e == e2 {
                            dp::extend_single_affine(q_sub, t_sub, 5, &mat, q as i8, e as i8, bw1, zdrop2, gap_end_bonus, generic_sc_flag, &mut ez);
                        } else {
                            dp::extend_dual_affine_dispatch(q_sub, t_sub, 5, &mat, q as i8, e as i8, q2 as i8, e2 as i8, bw1, zdrop2, gap_end_bonus, generic_sc_flag, &mut ez);
                        }
                        #[cfg(feature = "census")]
                        dp::census::finish(census_t0, dp::census::Kind::GapRerun, q_sub, t_sub,
                            &dp::census::Params { mat: &mat, q: q as i8, e: e as i8, q2: q2 as i8, e2: e2 as i8, bandwidth: bw1, zdrop: zdrop2, end_bonus: gap_end_bonus, flags: generic_sc_flag }, &ez);
                    }
                }

                if !ez.cigar.is_empty() {
                    append_cigar(&mut raw_cigar, &ez.cigar);
                }

                if ez.zdropped != 0 {
                    dp_score += ez.max;
                    re1 = gap_rs + ez.max_score_target_pos + 1;
                    qe1 = gap_qs + ez.max_score_query_pos + 1;
                    // Find last anchor before z-drop truncation point
                    let max_re_split = gap_rs + ez.max_score_target_pos;
                    let mut j_split = i as i32 - 1;
                    while j_split >= 0 {
                        if work_anchors[j_split as usize].ref_pos() <= max_re_split {
                            break;
                        }
                        j_split -= 1;
                    }
                    if j_split < 0 { j_split = 0; }
                    // Check if remaining anchors >= min_cnt
                    let remaining = cnt1 as i32 - (j_split + 1);
                    if remaining >= min_cnt {
                        // Split: save right part anchors from original chain
                        // In our terms: split_offset = as1_offset + (j_split + 1) as usize
                        let split_off = as1_offset + (j_split + 1) as usize;
                        if split_off < anchors.len() {
                            split_right_anchors = Some(anchors[split_off..].to_vec());
                            split_offset_in_orig = Some(split_off);
                            // if zdrop_code == 2, mark split as inversion
                            if zdrop_code == 2 {
                                split_inv_flag = true;
                            }
                        }
                    }
                    dropped = true;
                    break;
                } else {
                    dp_score += ez.score;
                }
            }

            // Update gap-fill start for next segment
            gap_rs = anchor_re;
            gap_qs = anchor_qe;
        }
    }

    // --- Step 7: Right extension ---
    let mut final_re = re1;
    let mut final_qe = qe1;

    if !dropped && qe < qlen && re < tlen {
        // Compute right extension boundary
        let (mut qe0, mut re0);
        {
            let orig_last = &anchors[anchors.len() - 1];
            let (re0_tmp, qe0_tmp) = compute_right_boundary(
                re, qe,
                orig_last.ref_pos() + 1, orig_last.query_pos() + 1,
                seed_bounds.2, seed_bounds.3,
                opt,
                qlen, tlen,
                seed_bounds,
            );
            re0 = re0_tmp;
            qe0 = qe0_tmp;
        }

        // SEED_SELF: clamp rightward extension to |q-r| of LAST pre-trim
        // anchor. Same coord-system issue as the left clamp: `pre_qe` is
        // absolute, `pre_re` is region-relative after the anchor rebase;
        // compute max_ext in absolute coords.
        if !anchors.is_empty() && (anchors[0].y & crate::align::map::SEED_SELF) != 0 {
            let last_idx = anchors.len() - 1;
            let pre_qe = anchors[last_idx].query_pos() + 1;
            let pre_re = anchors[last_idx].ref_pos() + 1;
            let pre_re_abs = pre_re + ref_offset;
            let max_ext = (pre_qe - pre_re_abs).unsigned_abs() as i32;
            if re0 - pre_re > max_ext { re0 = pre_re + max_ext; }
            if qe0 - pre_qe > max_ext { qe0 = pre_qe + max_ext; }
        }

        if qe0 > qe && re0 > re {
            // Lazy fill the right-extension segment [re, re0) in tseq.
            if use_lazy && let Some((mi, rid, rgn_start)) = lazy_extract {
                let a = re as usize;
                let b = re0 as usize;
                mi.extract_nt4_into(rid, rgn_start + a, rgn_start + b, &mut tseq[a..b]);
            }
            let q_sub = &qseq[qe as usize..qe0 as usize];
            let t_sub = &tseq[re as usize..re0 as usize];

            // Right extension flags: EXTZ_ONLY
            let flag = EXTENSION_ONLY | generic_sc_flag;
            let mut ez = DpResult::default();
            #[cfg(feature = "census")]
            let census_t0 = dp::census::start();
            // max_sw_mat check: skip alignment if matrix too large
            let sw_mat_size = q_sub.len() as i64 * t_sub.len() as i64;
            if max_sw_mat > 0 && sw_mat_size > max_sw_mat {
                ez.max_score_query_pos = -1;
                ez.max_score_target_pos = -1;
                ez.max_query_end_target_pos = -1;
                ez.zdropped = 1;
            } else if q == q2 && e == e2 {
                dp::extend_single_affine(q_sub, t_sub, 5, &mat, q as i8, e as i8, bw, zdrop, end_bonus, flag, &mut ez);
            } else {
                let ext_zdrop = crate::align::env_flags::EXP_EXT_ZDROP.unwrap_or(zdrop);
                if extend_dual_ext(q_sub, t_sub, &mat, q as i8, e as i8, q2 as i8, e2 as i8, ext_bw, bw, ext_zdrop, end_bonus, flag, ext_rescue, ext_probe, ext_adapt, &mut ez) {
                    probe_claim_qe = Some(qe0);
                }
            }
            #[cfg(feature = "census")]
            dp::census::finish(census_t0, dp::census::Kind::RightExt, q_sub, t_sub,
                &dp::census::Params { mat: &mat, q: q as i8, e: e as i8, q2: q2 as i8, e2: e2 as i8, bandwidth: bw, zdrop, end_bonus, flags: flag }, &ez);

            if !ez.cigar.is_empty() {
                append_cigar(&mut raw_cigar, &ez.cigar);
                dp_score += ez.max;
            }

            // Consumed bases
            final_re = re + if ez.reach_end != 0 { ez.max_query_end_target_pos + 1 } else { ez.max_score_target_pos + 1 };
            final_qe = qe + if ez.reach_end != 0 { qe0 - qe } else { ez.max_score_query_pos + 1 };
        }
    }

    if final_re <= 0 { final_re = re1; }
    if final_qe <= 0 { final_qe = qe1; }

    // Lazy mode: fill the final alignment range [rs1, final_re) so that
    // `finalize_cigar` (indel right-alignment, op merging, =/X conversion)
    // and any downstream MD/cs tag generators can read tseq there.
    if use_lazy && let Some((mi, rid, rgn_start)) = lazy_extract {
        let a = cmp::max(0, rs1) as usize;
        let b = cmp::max(a as i32, final_re) as usize;
        if b > a {
            mi.extract_nt4_into(rid, rgn_start + a, rgn_start + b, &mut tseq[a..b]);
        }
    }

    // --- Steps 8-9: fix_cigar + convert M→=/X + condense ---
    let (condensed, final_qs, final_rs) = finalize_cigar(raw_cigar, qseq, tseq, qs1, rs1);


    AlignResult {
        cigar_ops: condensed,
        query_start: cmp::max(0, final_qs) as usize,
        query_end: cmp::max(0, final_qe) as usize,
        ref_start: cmp::max(0, final_rs) as usize,
        ref_end: cmp::max(0, final_re) as usize,
        dp_score,
        split_right_anchors,
        split_offset_in_orig,
        split_inv: split_inv_flag,
        probe_claim_qs: probe_claim_qs.map(|v| cmp::max(0, v) as usize),
        probe_claim_qe: probe_claim_qe.map(|v| cmp::max(0, v) as usize),
    }
}

/// Public wrapper for convert_cigar_to_eqx (used by pipeline.rs for inversion alignment)
pub fn convert_cigar_to_eqx_pub(raw_cigar: &[u32], qseq: &[u8], tseq: &[u8], qs: usize, rs: usize) -> Vec<CigarOp> {
    convert_cigar_to_eqx(raw_cigar, qseq, tseq, qs, rs)
}

/// Public wrapper for fix_cigar (used by pipeline.rs for inversion alignment)
pub fn fix_cigar_pub(cigar: &mut Vec<u32>, qseq: &[u8], tseq: &[u8], qs: i32, rs: i32) -> (i32, i32) {
    fix_cigar(cigar, qseq, tseq, qs, rs)
}

/// Convert M ops to =/X ops (EQX mode).
fn convert_cigar_to_eqx(raw_cigar: &[u32], qseq: &[u8], tseq: &[u8], qs: usize, rs: usize) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut qi = qs;
    let mut ti = rs;

    for &c in raw_cigar {
        let len = (c >> 4) as usize;
        let op = c & 0xf;

        match op {
            0 => { // M -> split into = and X
                let mut pos = 0;
                while pos < len {
                    // Run of matches
                    let start = pos;
                    while pos < len && qi + pos < qseq.len() && ti + pos < tseq.len()
                        && qseq[qi + pos] == tseq[ti + pos] {
                        pos += 1;
                    }
                    if pos > start {
                        ops.push(CigarOp { op: '=', len: (pos - start) as u32 });
                    }
                    // Run of mismatches
                    let start = pos;
                    while pos < len && qi + pos < qseq.len() && ti + pos < tseq.len()
                        && qseq[qi + pos] != tseq[ti + pos] {
                        pos += 1;
                    }
                    if pos > start {
                        ops.push(CigarOp { op: 'X', len: (pos - start) as u32 });
                    }
                }
                qi += len;
                ti += len;
            },
            1 => { // I
                ops.push(CigarOp { op: 'I', len: len as u32 });
                qi += len;
            },
            2 => { // D
                ops.push(CigarOp { op: 'D', len: len as u32 });
                ti += len;
            },
            3 => { // N
                ops.push(CigarOp { op: 'N', len: len as u32 });
                ti += len;
            },
            7 => { // = (already EQX)
                ops.push(CigarOp { op: '=', len: len as u32 });
                qi += len;
                ti += len;
            },
            8 => { // X (already EQX)
                ops.push(CigarOp { op: 'X', len: len as u32 });
                qi += len;
                ti += len;
            },
            _ => {
                ops.push(CigarOp { op: 'M', len: len as u32 });
                qi += len;
                ti += len;
            }
        }
    }
    ops
}

pub fn fmt_cs(
    ops: &[CigarOp],
    qseq: &[u8],
    tseq: &[u8],
    qs: usize,
    rs: usize,
    long: bool,
) -> String {
    let mut visitor = CsVisitor::new(qseq, tseq, qs, rs, long);
    walk_cigar_ops(ops, &mut visitor)
}

/// Generate ds:Z difference string from CIGAR and sequences.
pub fn fmt_ds(
    ops: &[CigarOp],
    qseq: &[u8],
    tseq: &[u8],
    qs: usize,
    rs: usize,
) -> String {
    // First pass: compute q_len and t_len (total aligned lengths for bounds checking)
    let mut q_len: usize = 0;
    let mut t_len: usize = 0;
    for op in ops {
        match op.op {
            '=' | 'X' | 'M' => { q_len += op.len as usize; t_len += op.len as usize; }
            'I' => { q_len += op.len as usize; }
            'D' | 'N' => { t_len += op.len as usize; }
            _ => {}
        }
    }
    let mut visitor = DsVisitor::new(qseq, tseq, qs, rs, q_len, t_len);
    walk_cigar_ops(ops, &mut visitor)
}

/// Format indel portion of ds:Z string.
fn write_indel_ds(out: &mut String, len: usize, seq: &[u8], ll: usize, lr: usize) {
    if ll + lr >= len {
        // Entire indel is ambiguous
        out.push('[');
        for &b in &seq[..len] {
            out.push(NT4_TO_LOWER[b.min(4) as usize] as char);
        }
        out.push(']');
    } else {
        let mut k: usize = 0;
        if ll > 0 {
            out.push('[');
            for i in 0..ll {
                out.push(NT4_TO_LOWER[seq[k + i].min(4) as usize] as char);
            }
            out.push(']');
            k += ll;
        }
        for i in 0..(len - lr - ll) {
            out.push(NT4_TO_LOWER[seq[k + i].min(4) as usize] as char);
        }
        k += len - lr - ll;
        if lr > 0 {
            out.push('[');
            for i in 0..lr {
                out.push(NT4_TO_LOWER[seq[k + i].min(4) as usize] as char);
            }
            out.push(']');
        }
    }
}
