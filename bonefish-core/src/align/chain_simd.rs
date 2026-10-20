//! SIMD-optimized chaining for x86_64 (AVX2), aarch64 (NEON), and wasm32 (SIMD128).
//!
//! Two-phase approach:
//! - Phase 1 (SIMD): Batch-compute chain scores for all predecessors using SIMD
//! - Phase 2 (Scalar): Apply max_chain_skip / visited logic over precomputed scores
//!
//! This guarantees bit-exact output with the scalar implementation while
//! getting 3-5x speedup on the score computation hot path.

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

#[cfg(target_arch = "wasm32")]
use std::arch::wasm32::*;

use crate::align::chain::{chain_backtrack, compute_chain_score, fast_log2};
use crate::align::map::{ChainingParams, ChainingBuffers};
use crate::align::sketch::Minimizer;
use crate::align::sort::radix_sort_128x;

/// Vectorized fast_log2 for 8 floats (AVX2).
/// Matches the scalar `fast_log2` exactly (same polynomial, same rounding).
/// MUST NOT use FMA — use explicit mul+add to match scalar rounding.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn simd_fast_log2_avx2(x: __m256) -> __m256 {
    let z = _mm256_castps_si256(x);

    // exponent = ((z >> 23) & 255) - 128
    let exp_raw = _mm256_and_si256(_mm256_srli_epi32(z, 23), _mm256_set1_epi32(255));
    let exp = _mm256_sub_epi32(exp_raw, _mm256_set1_epi32(128));
    let log_2 = _mm256_cvtepi32_ps(exp);

    // mantissa: clear exponent bits, set exponent to 127 (1.0 bias)
    let mant_bits = _mm256_or_si256(
        _mm256_and_si256(z, _mm256_set1_epi32(!(255i32 << 23))),
        _mm256_set1_epi32(127i32 << 23),
    );
    let z_f = _mm256_castsi256_ps(mant_bits);

    // Horner's: ((-0.34484843 * z_f) + 2.0246658) * z_f + (-0.6748776)
    let c0 = _mm256_set1_ps(-0.34484843f32);
    let c1 = _mm256_set1_ps(2.024_665_8f32);
    let c2 = _mm256_set1_ps(-0.674_877_6f32);

    let t1 = _mm256_mul_ps(c0, z_f);        // c0 * z_f
    let t2 = _mm256_add_ps(t1, c1);          // c0 * z_f + c1
    let t3 = _mm256_mul_ps(t2, z_f);         // (c0 * z_f + c1) * z_f
    let poly = _mm256_add_ps(t3, c2);        // ... + c2

    _mm256_add_ps(log_2, poly)
}

/// Batch-compute chain scores for predecessors j in [start, end) against anchor i.
/// Writes scores to sc_buf[start..end]. Invalid pairs get i32::MIN.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn compute_chain_scores_batch_avx2(
    qi: i32,
    ri: i32,
    rid_strand_i: u32,
    soa_ref_pos: &[i32],
    soa_query_pos: &[i32],
    soa_query_span: &[i32],
    soa_rid_strand: &[u32],
    max_dist_x: i32,
    max_dist_y: i32,
    bw: i32,
    chn_pen_gap: f32,
    chn_pen_skip: f32,
    sc_buf: &mut [i32],
    start: usize,
    end: usize,
) { unsafe {
    if start >= end {
        return;
    }

    let qi_v = _mm256_set1_epi32(qi);
    let ri_v = _mm256_set1_epi32(ri);
    let rid_i_v = _mm256_set1_epi32(rid_strand_i as i32);
    let max_dx_v = _mm256_set1_epi32(max_dist_x);
    let max_dy_v = _mm256_set1_epi32(max_dist_y);
    let bw_v = _mm256_set1_epi32(bw);
    let zero_v = _mm256_setzero_si256();
    let one_v = _mm256_set1_epi32(1);
    let min_val_v = _mm256_set1_epi32(i32::MIN);
    let pen_gap_v = _mm256_set1_ps(chn_pen_gap);
    let pen_skip_v = _mm256_set1_ps(chn_pen_skip);
    let half_v = _mm256_set1_ps(0.5f32);

    // Process 8 predecessors at a time
    let aligned_end = start + ((end - start) / 8) * 8;
    let mut j = start;

    while j < aligned_end {
        // Load 8 predecessors' SoA data
        let qj_v = _mm256_loadu_si256(soa_query_pos.as_ptr().add(j) as *const __m256i);
        let rj_v = _mm256_loadu_si256(soa_ref_pos.as_ptr().add(j) as *const __m256i);
        let span_v = _mm256_loadu_si256(soa_query_span.as_ptr().add(j) as *const __m256i);
        let rid_j_v = _mm256_loadu_si256(soa_rid_strand.as_ptr().add(j) as *const __m256i);

        // query_diff = qi - qj (wrapping subtraction, matches Rust wrapping_sub)
        let qdiff = _mm256_sub_epi32(qi_v, qj_v);
        // ref_diff = ri - rj
        let rdiff = _mm256_sub_epi32(ri_v, rj_v);

        // --- Build validity mask ---
        // 1. query_diff > 0
        let mut valid = _mm256_cmpgt_epi32(qdiff, zero_v);

        // 2. query_diff <= max_dist_x  (i.e., NOT (query_diff > max_dist_x))
        let qdiff_gt_max_dx = _mm256_cmpgt_epi32(qdiff, max_dx_v);
        valid = _mm256_andnot_si256(qdiff_gt_max_dx, valid);

        // 3. Same ref_id_strand
        let same_rid = _mm256_cmpeq_epi32(rid_i_v, rid_j_v);
        valid = _mm256_and_si256(valid, same_rid);

        // 4. ref_diff != 0 (same_rid already ensures same target)
        let rdiff_zero = _mm256_cmpeq_epi32(rdiff, zero_v);
        valid = _mm256_andnot_si256(rdiff_zero, valid);

        // 5. query_diff <= max_dist_y  (NOT (qdiff > max_dy))
        let qdiff_gt_max_dy = _mm256_cmpgt_epi32(qdiff, max_dy_v);
        valid = _mm256_andnot_si256(qdiff_gt_max_dy, valid);

        // 6. gap_width = abs(ref_diff - query_diff); gap_width <= bw
        let diff = _mm256_sub_epi32(rdiff, qdiff);
        let gap_w = _mm256_abs_epi32(diff);
        let gw_gt_bw = _mm256_cmpgt_epi32(gap_w, bw_v);
        valid = _mm256_andnot_si256(gw_gt_bw, valid);


        // No valid lane (usually: off the band): the scores are all i32::MIN.
        if _mm256_testz_si256(valid, valid) != 0 {
            _mm256_storeu_si256(sc_buf.as_mut_ptr().add(j) as *mut __m256i, min_val_v);
            j += 8;
            continue;
        }

        // --- Score computation ---
        // min_diff = min(ref_diff, query_diff)
        let min_diff = _mm256_min_epi32(rdiff, qdiff);

        // score = min(q_span_j, min_diff)
        let mut score = _mm256_min_epi32(span_v, min_diff);

        // --- Penalty computation ---
        // needs_penalty = (gap_width > 0) || (min_diff > q_span)
        let gw_gt_zero = _mm256_cmpgt_epi32(gap_w, zero_v);
        let md_gt_span = _mm256_cmpgt_epi32(min_diff, span_v);
        let needs_pen = _mm256_or_si256(gw_gt_zero, md_gt_span);

        // lin_pen = chn_pen_gap * gap_width + chn_pen_skip * min_diff
        let gap_w_f = _mm256_cvtepi32_ps(gap_w);
        let min_diff_f = _mm256_cvtepi32_ps(min_diff);
        let mul1 = _mm256_mul_ps(pen_gap_v, gap_w_f);
        let mul2 = _mm256_mul_ps(pen_skip_v, min_diff_f);
        let lin_pen = _mm256_add_ps(mul1, mul2);

        // log_pen = fast_log2(gap_width + 1) if gap_width >= 1, else 0
        let gw_plus1 = _mm256_add_epi32(gap_w, one_v);
        let gw_plus1_f = _mm256_cvtepi32_ps(gw_plus1);
        let log2_val = simd_fast_log2_avx2(gw_plus1_f);
        // Mask: zero out log_pen where gap_width == 0
        let log_pen = _mm256_and_ps(_mm256_castsi256_ps(gw_gt_zero), log2_val);

        // total_pen = lin_pen + 0.5 * log_pen
        let half_log = _mm256_mul_ps(half_v, log_pen);
        let total_pen = _mm256_add_ps(lin_pen, half_log);

        // Convert to i32 (truncation toward zero, matches Rust `as i32`)
        let pen_i32 = _mm256_cvttps_epi32(total_pen);

        // Apply penalty only where needs_pen is set
        let pen_masked = _mm256_and_si256(needs_pen, pen_i32);
        score = _mm256_sub_epi32(score, pen_masked);

        // Mask invalid scores to i32::MIN
        score = _mm256_blendv_epi8(min_val_v, score, valid);

        // Store 8 scores
        _mm256_storeu_si256(sc_buf.as_mut_ptr().add(j) as *mut __m256i, score);

        j += 8;
    }

    // Scalar remainder
    while j < end {
        // We can't call compute_chain_score directly because it takes Minimizer refs,
        // but we have SoA data. Recompute from SoA fields to match exactly.
        let qj = soa_query_pos[j];
        let rj = soa_ref_pos[j];
        let span_j = soa_query_span[j];
        let rid_j = soa_rid_strand[j];

        let query_diff = qi.wrapping_sub(qj);
        if query_diff <= 0 || query_diff > max_dist_x {
            sc_buf[j] = i32::MIN;
            j += 1;
            continue;
        }

        let ref_diff = ri.wrapping_sub(rj);
        if rid_strand_i == rid_j && (ref_diff == 0 || query_diff > max_dist_y) {
            sc_buf[j] = i32::MIN;
            j += 1;
            continue;
        }

        let gap_width = (ref_diff - query_diff).abs();
        if rid_strand_i == rid_j && gap_width > bw {
            sc_buf[j] = i32::MIN;
            j += 1;
            continue;
        }

        let min_diff = ref_diff.min(query_diff);
        let mut sc = span_j.min(min_diff);

        if gap_width > 0 || min_diff > span_j {
            let lin_pen = chn_pen_gap * (gap_width as f32) + chn_pen_skip * (min_diff as f32);
            let log_pen = if gap_width >= 1 {
                fast_log2((gap_width + 1) as f32)
            } else {
                0.0f32
            };
            sc -= (lin_pen + 0.5f32 * log_pen) as i32;
        }

        sc_buf[j] = sc;
        j += 1;
    }
}}

/// Predecessors of anchor i only score when their diagonal (ref - query position) is within
/// `bw` of i's; the rest of the DP chaining window is a no-op scan (most of it in repeats:
/// ~86% of candidates genome-wide, ~97% in tandem-repeat reads). Anchors are grouped into
/// diagonal buckets of width 2^shift >= bw, so i's in-band predecessors all lie in buckets
/// key-1..=key+1. Merging those buckets' anchor lists from the back visits the same in-band
/// predecessors in the same (descending) order as the full scan, and the side effects of the
/// scan (best predecessor, skip count, `visited` marks) only happen for in-band ones, so the
/// result is identical.
#[derive(Default)]
pub struct DiagBuckets {
    /// Bucket id and diagonal of each anchor.
    id: Vec<u32>,
    diag: Vec<i64>,
    /// Per bucket id: ids of buckets key-1, key, key+1 (u32::MAX if absent).
    nbr: Vec<[u32; 3]>,
    /// Per bucket: indices of the anchors processed so far, ascending.
    lists: Vec<Vec<u32>>,
    map: std::collections::HashMap<i64, u32, std::hash::BuildHasherDefault<DiagIntegerHasher>>,
    keys: Vec<i64>,
    g_idx: Vec<u32>,
    g_ref: Vec<i32>,
    g_qry: Vec<i32>,
    g_span: Vec<i32>,
    g_rid: Vec<u32>,
    g_sc: Vec<i32>,
}

/// Hasher for the bucket map's diagonal keys: a splitmix64 finalizer, much cheaper than the
/// default SipHash. Bucket ids follow insertion order and the map is never iterated, so the
/// hasher cannot change any result.
#[derive(Default)]
struct DiagIntegerHasher(u64);

impl std::hash::Hasher for DiagIntegerHasher {
    fn finish(&self) -> u64 {
        let mut x = self.0;
        x ^= x >> 30;
        x = x.wrapping_mul(0xbf58476d1ce4e5b9);
        x ^= x >> 27;
        x = x.wrapping_mul(0x94d049bb133111eb);
        x ^ (x >> 31)
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 ^ b as u64).wrapping_mul(0x100000001b3);
        }
    }
    fn write_i64(&mut self, value: i64) {
        self.0 = value as u64;
    }
    fn write_u64(&mut self, value: u64) {
        self.0 = value;
    }
}

/// Build the buckets for reads with at least this many anchors and at least
/// DIAG_BUCKET_MIN_DENSITY anchors per query base (repeat-dense reads; typical reads have
/// ~0.4 and a mostly in-band window, where the contiguous SIMD scan is cheaper).
const DIAG_BUCKET_MIN_N: usize = 1024;
const DIAG_BUCKET_MIN_DENSITY: usize = 2;
/// Use the buckets for anchors whose window holds more than this many predecessors.
const DIAG_BUCKET_MIN_WINDOW: usize = 128;

impl DiagBuckets {
    fn build(&mut self, ref_pos: &[i32], query_pos: &[i32], bw: i32) {
        let n = ref_pos.len();
        let shift = if bw > 1 { 32 - ((bw - 1) as u32).leading_zeros() } else { 0 };
        self.map.clear();
        self.keys.clear();
        self.id.clear();
        self.diag.clear();
        for j in 0..n {
            let d = ref_pos[j] as i64 - query_pos[j] as i64;
            self.diag.push(d);
            let key = d >> shift;
            let next = self.keys.len() as u32;
            let id = *self.map.entry(key).or_insert(next);
            if id == next { self.keys.push(key); }
            self.id.push(id);
        }
        let nb = self.keys.len();
        if self.lists.len() < nb { self.lists.resize_with(nb, Vec::new); }
        for l in &mut self.lists[..nb] { l.clear(); }
        self.nbr.clear();
        for &key in &self.keys {
            let get = |k: i64| self.map.get(&k).copied().unwrap_or(u32::MAX);
            self.nbr.push([get(key - 1), get(key), get(key + 1)]);
        }
    }
}

/// `compute_chain_scores_batch_avx2` for one scan chunk of at most 64 candidates `[start, end)`,
/// returning a bitmask of the valid ones (bit b: candidate start + b; exactly the candidates
/// whose score is not i32::MIN). Most groups are off the diagonal band, so the band test comes
/// first, from the precomputed diagonals (ref - query): |diag_i - diag_j| (wrapping) is the batch
/// kernel's |ref_diff - query_diff|, and a group where every lane exceeds `bw` costs one load and
/// writes nothing. Scores are written only for groups with a candidate lane; others must not be
/// read.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
#[allow(clippy::needless_range_loop)] // tail loop over sc_buf[j..end] builds a lane mask from k
unsafe fn chain_scores_masked_avx2(
    qi: i32, ri: i32, rid_strand_i: u32, diag_i: i32,
    soa_ref_pos: &[i32], soa_query_pos: &[i32], soa_query_span: &[i32], soa_rid_strand: &[u32], soa_diag: &[i32],
    max_dist_x: i32, max_dist_y: i32, bw: i32, chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], start: usize, end: usize,
) -> u64 { unsafe {
    debug_assert!(end - start <= 64);
    let diag_i_v = _mm256_set1_epi32(diag_i);
    let qi_v = _mm256_set1_epi32(qi);
    let ri_v = _mm256_set1_epi32(ri);
    let rid_i_v = _mm256_set1_epi32(rid_strand_i as i32);
    let max_dx_v = _mm256_set1_epi32(max_dist_x);
    let max_dy_v = _mm256_set1_epi32(max_dist_y);
    let bw_v = _mm256_set1_epi32(bw);
    let zero_v = _mm256_setzero_si256();
    let one_v = _mm256_set1_epi32(1);
    let min_val_v = _mm256_set1_epi32(i32::MIN);
    let pen_gap_v = _mm256_set1_ps(chn_pen_gap);
    let pen_skip_v = _mm256_set1_ps(chn_pen_skip);
    let half_v = _mm256_set1_ps(0.5f32);
    let mut mask = 0u64;
    let aligned_end = start + ((end - start) / 8) * 8;
    let mut j = start;
    while j < aligned_end {
        let dj = _mm256_loadu_si256(soa_diag.as_ptr().add(j) as *const __m256i);
        // Same signed test as the batch kernel (abs(i32::MIN) stays negative: not "too wide").
        let gap_w = _mm256_abs_epi32(_mm256_sub_epi32(diag_i_v, dj));
        let gw_gt_bw = _mm256_cmpgt_epi32(gap_w, bw_v);
        if _mm256_movemask_epi8(gw_gt_bw) == -1 {
            j += 8;
            continue;
        }
        let qj_v = _mm256_loadu_si256(soa_query_pos.as_ptr().add(j) as *const __m256i);
        let rj_v = _mm256_loadu_si256(soa_ref_pos.as_ptr().add(j) as *const __m256i);
        let span_v = _mm256_loadu_si256(soa_query_span.as_ptr().add(j) as *const __m256i);
        let rid_j_v = _mm256_loadu_si256(soa_rid_strand.as_ptr().add(j) as *const __m256i);
        let qdiff = _mm256_sub_epi32(qi_v, qj_v);
        let rdiff = _mm256_sub_epi32(ri_v, rj_v);
        let mut valid = _mm256_cmpgt_epi32(qdiff, zero_v);
        valid = _mm256_andnot_si256(_mm256_cmpgt_epi32(qdiff, max_dx_v), valid);
        valid = _mm256_and_si256(valid, _mm256_cmpeq_epi32(rid_i_v, rid_j_v));
        valid = _mm256_andnot_si256(_mm256_cmpeq_epi32(rdiff, zero_v), valid);
        valid = _mm256_andnot_si256(_mm256_cmpgt_epi32(qdiff, max_dy_v), valid);
        valid = _mm256_andnot_si256(gw_gt_bw, valid);
        let vm = _mm256_movemask_ps(_mm256_castsi256_ps(valid)) as u32;
        if vm == 0 {
            j += 8;
            continue;
        }
        let min_diff = _mm256_min_epi32(rdiff, qdiff);
        let mut score = _mm256_min_epi32(span_v, min_diff);
        let gw_gt_zero = _mm256_cmpgt_epi32(gap_w, zero_v);
        let md_gt_span = _mm256_cmpgt_epi32(min_diff, span_v);
        let needs_pen = _mm256_or_si256(gw_gt_zero, md_gt_span);
        let lin_pen = _mm256_add_ps(_mm256_mul_ps(pen_gap_v, _mm256_cvtepi32_ps(gap_w)),
                                    _mm256_mul_ps(pen_skip_v, _mm256_cvtepi32_ps(min_diff)));
        let log2_val = simd_fast_log2_avx2(_mm256_cvtepi32_ps(_mm256_add_epi32(gap_w, one_v)));
        let log_pen = _mm256_and_ps(_mm256_castsi256_ps(gw_gt_zero), log2_val);
        let total_pen = _mm256_add_ps(lin_pen, _mm256_mul_ps(half_v, log_pen));
        let pen_i32 = _mm256_cvttps_epi32(total_pen);
        score = _mm256_sub_epi32(score, _mm256_and_si256(needs_pen, pen_i32));
        score = _mm256_blendv_epi8(min_val_v, score, valid);
        _mm256_storeu_si256(sc_buf.as_mut_ptr().add(j) as *mut __m256i, score);
        mask |= (vm as u64) << (j - start);
        j += 8;
    }
    if j < end {
        compute_chain_scores_batch_avx2(qi, ri, rid_strand_i, soa_ref_pos, soa_query_pos, soa_query_span, soa_rid_strand,
            max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, sc_buf, j, end);
        for k in j..end {
            if sc_buf[k] != i32::MIN { mask |= 1u64 << (k - start); }
        }
    }
    mask
}}

/// AVX2 SIMD-optimized chaining.
/// Phase 1: Batch-compute scores with SIMD. Phase 2: Scalar skip/visited logic.
/// Produces bit-exact output matching `chain_anchors_scalar`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn chain_anchors_avx2(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) { unsafe {
    let n = a.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    let bw = opt.bandwidth;
    let mut max_dist_x = max_dist_x;
    let mut max_dist_y = max_dist_y;
    if max_dist_x < bw {
        max_dist_x = bw;
    }
    if max_dist_y < bw {
        max_dist_y = bw;
    }

    let real_max_drop = bw; // non-cdna path
    let no_chain_skip = opt.max_chain_skip <= 0;

    // --- Extract SoA fields from AoS Minimizer array ---
    let mut soa_ref_pos = std::mem::take(&mut ctx.soa_ref_pos);
    let mut soa_query_pos = std::mem::take(&mut ctx.soa_query_pos);
    let mut soa_query_span = std::mem::take(&mut ctx.soa_query_span);
    let mut soa_rid_strand = std::mem::take(&mut ctx.soa_ref_id_strand);
    let mut soa_sc_buf = std::mem::take(&mut ctx.soa_scores_buf);
    let mut soa_diag = std::mem::take(&mut ctx.soa_diag);
    soa_diag.resize(n, 0);

    soa_ref_pos.resize(n, 0);
    soa_query_pos.resize(n, 0);
    soa_query_span.resize(n, 0);
    soa_rid_strand.resize(n, 0);
    soa_sc_buf.resize(n, 0);

    for j in 0..n {
        soa_ref_pos[j] = a[j].ref_pos();
        soa_query_pos[j] = a[j].query_pos();
        soa_query_span[j] = a[j].query_span();
        soa_rid_strand[j] = (a[j].x >> 32) as u32;
        soa_diag[j] = soa_ref_pos[j].wrapping_sub(soa_query_pos[j]);
    }

    // --- DP buffers ---
    let mut predecessors = std::mem::take(&mut ctx.predecessors); let mut bt_candidates = std::mem::take(&mut ctx.bt_candidates);
    let mut scores = std::mem::take(&mut ctx.scores);
    let mut peak_scores = std::mem::take(&mut ctx.peak_scores);
    let mut visited = std::mem::take(&mut ctx.visited);
    predecessors.resize(n, 0i32);
    scores.resize(n, 0i32);
    peak_scores.resize(n, 0i32);
    visited.clear();
    visited.resize(n, 0i32);

    let mut global_max_score = 0;
    let mut window_start: usize = 0;
    let mut best_anchor_idx: i64 = -1;

    let mut bk = std::mem::take(&mut ctx.diag_buckets);
    let use_buckets = !no_chain_skip && n >= DIAG_BUCKET_MIN_N && {
        let (lo, hi) = soa_query_pos.iter().fold((i32::MAX, i32::MIN), |(l, h), &q| (l.min(q), h.max(q)));
        n >= DIAG_BUCKET_MIN_DENSITY * (hi - lo + 1).max(1) as usize
    };
    if use_buckets { bk.build(&soa_ref_pos, &soa_query_pos, bw); }

    for i in 0..n {
        let mut best_predecessor: i64 = -1;
        let mut best_score = a[i].query_span();
        let mut skip_count: i32 = 0;

        while window_start < i
            && (a[i].ref_id_strand() != a[window_start].ref_id_strand()
                || (a[i].x > a[window_start].x + max_dist_x as u64))
        {
            window_start += 1;
        }

        if !no_chain_skip && (i - window_start) > opt.max_chain_iter as usize {
            window_start = i - opt.max_chain_iter as usize;
        }

        let mut end_j: i64 = if window_start > 0 {
            window_start as i64 - 1
        } else {
            -1
        };

        if no_chain_skip {
            // No early termination — compute all scores, then scan
            if i > window_start {
                compute_chain_scores_batch_avx2(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, window_start, i,
                );
            }
            for j in window_start..i {
                let sc = soa_sc_buf[j];
                if sc == i32::MIN { continue; }
                let total_sc = sc + scores[j];
                if total_sc > best_score { best_score = total_sc; best_predecessor = j as i64; }
            }
        } else if use_buckets && i - window_start > DIAG_BUCKET_MIN_WINDOW {
            // In-band predecessors only, merged from the three diagonal buckets, descending.
            const CHUNK_SIZE: usize = 64;
            let nb = bk.nbr[bk.id[i] as usize];
            // Merge cursors: list base pointers and remaining lengths (empty for absent buckets).
            let mut lp: [*const u32; 3] = [std::ptr::null(); 3];
            let mut cur = [0usize; 3];
            for t in 0..3 {
                if nb[t] != u32::MAX { let l = &bk.lists[nb[t] as usize]; lp[t] = l.as_ptr(); cur[t] = l.len(); }
            }
            let tail = |t: usize, c: usize| -> i64 { if c > 0 { *lp[t].add(c - 1) as i64 } else { -1 } };
            let (mut v0, mut v1, mut v2) = (tail(0, cur[0]), tail(1, cur[1]), tail(2, cur[2]));
            let (di, dg) = (bk.diag[i], bk.diag.as_ptr());
            let ws = window_start as i64;
            let mut skip_fired = false;
            let mut exhausted = false;
            loop {
                bk.g_idx.clear();
                while bk.g_idx.len() < CHUNK_SIZE {
                    // Largest remaining index across the three buckets.
                    let t = if v0 >= v1 && v0 >= v2 { 0 } else if v1 >= v2 { 1 } else { 2 };
                    let j = [v0, v1, v2][t];
                    if j < ws { exhausted = true; break; }
                    cur[t] -= 1;
                    match t { 0 => v0 = tail(0, cur[0]), 1 => v1 = tail(1, cur[1]), _ => v2 = tail(2, cur[2]) }
                    // Out of band (|diagonal difference| > bw) scores i32::MIN: skip it here.
                    if (*dg.add(j as usize) - di).abs() > bw as i64 { continue; }
                    bk.g_idx.push(j as u32);
                }
                let m = bk.g_idx.len();
                if m == 0 { break; }
                bk.g_ref.clear(); bk.g_qry.clear(); bk.g_span.clear(); bk.g_rid.clear();
                for &j in &bk.g_idx {
                    let j = j as usize;
                    bk.g_ref.push(soa_ref_pos[j]); bk.g_qry.push(soa_query_pos[j]);
                    bk.g_span.push(soa_query_span[j]); bk.g_rid.push(soa_rid_strand[j]);
                }
                bk.g_sc.resize(m, 0);
                compute_chain_scores_batch_avx2(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &bk.g_ref, &bk.g_qry, &bk.g_span, &bk.g_rid,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut bk.g_sc, 0, m,
                );
                for t in 0..m {
                    let sc = bk.g_sc[t];
                    if sc == i32::MIN { continue; }
                    let j = bk.g_idx[t] as usize;
                    let total_sc = sc + scores[j];
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if visited[j] == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                    }
                    let p = predecessors[j];
                    if p >= 0 { visited[p as usize] = i as i32; }
                }
                if skip_fired || exhausted { break; }
            }

            if best_anchor_idx < 0 || a[i].x - a[best_anchor_idx as usize].x > max_dist_x as u64 {
                let mut best = i32::MIN; best_anchor_idx = -1;
                for j in (window_start..i).rev() {
                    let sj = *scores.get_unchecked(j);
                    let better = best < sj;
                    best_anchor_idx = if better { j as i64 } else { best_anchor_idx };
                    best = if better { sj } else { best };
                }
            }

            if best_anchor_idx >= 0 && best_anchor_idx < end_j {
                let sc = compute_chain_score(&a[i], &a[best_anchor_idx as usize],
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip);
                if sc != i32::MIN && best_score < sc + scores[best_anchor_idx as usize] {
                    best_score = sc + scores[best_anchor_idx as usize];
                    best_predecessor = best_anchor_idx;
                }
            }
        } else {
            // Chunked reverse scan with early termination. The skip test usually stops the scan
            // within a few dozen candidates, so the first chunk is smaller (32: -1% CPU vs 64).
            const CHUNK_SIZE: usize = 64;
            const FIRST_CHUNK: usize = 32;
            let mut chunk_end = i;
            let mut skip_fired = false;
            let mut csize = FIRST_CHUNK;
            while chunk_end > window_start {
                let chunk_start = if chunk_end >= window_start + csize { chunk_end - csize } else { window_start };
                csize = CHUNK_SIZE;

                // Bit b set: candidate chunk_start + b scored (invalid ones, mostly off the band,
                // are no-ops below), so only those are visited, in the same descending order.
                let mut valid = chain_scores_masked_avx2(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i], soa_diag[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand, &soa_diag,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, chunk_start, chunk_end,
                );
                // SAFETY: j < i < n and predecessors[j] < j, all within the n-sized buffers.
                let (scb, scr, vis, pre) = (soa_sc_buf.as_ptr(), scores.as_ptr(), visited.as_mut_ptr(), predecessors.as_ptr());
                while valid != 0 {
                    let b = 63 - valid.leading_zeros() as usize;
                    valid &= !(1u64 << b);
                    let j = chunk_start + b;
                    let sc = *scb.add(j);
                    let total_sc = sc + *scr.add(j);
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if *vis.add(j) == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                    }
                    let p = *pre.add(j);
                    if p >= 0 { *vis.add(p as usize) = i as i32; }
                }

                if skip_fired { break; }
                chunk_end = chunk_start;
            }

            if best_anchor_idx < 0 || a[i].x - a[best_anchor_idx as usize].x > max_dist_x as u64 {
                let mut best = i32::MIN; best_anchor_idx = -1;
                for j in (window_start..i).rev() {
                    let sj = *scores.get_unchecked(j);
                    let better = best < sj;
                    best_anchor_idx = if better { j as i64 } else { best_anchor_idx };
                    best = if better { sj } else { best };
                }
            }

            if best_anchor_idx >= 0 && best_anchor_idx < end_j {
                let sc = compute_chain_score(&a[i], &a[best_anchor_idx as usize],
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip);
                if sc != i32::MIN && best_score < sc + scores[best_anchor_idx as usize] {
                    best_score = sc + scores[best_anchor_idx as usize];
                    best_predecessor = best_anchor_idx;
                }
            }
        }

        scores[i] = best_score;
        predecessors[i] = best_predecessor as i32;
        peak_scores[i] = if best_predecessor >= 0
            && peak_scores[best_predecessor as usize] > best_score
        {
            peak_scores[best_predecessor as usize]
        } else {
            best_score
        };

        if best_anchor_idx < 0
            || (a[i].x - a[best_anchor_idx as usize].x <= max_dist_x as u64
                && scores[best_anchor_idx as usize] < scores[i])
        {
            best_anchor_idx = i as i64;
        }

        if global_max_score < best_score {
            global_max_score = best_score;
        }
        if use_buckets { bk.lists[bk.id[i] as usize].push(i as u32); }
    }
    ctx.diag_buckets = bk;


    // --- Backtrack and compact (identical to scalar) ---
    let (u, n_u, n_v) = chain_backtrack(
        n,
        &scores,
        &predecessors,
        &mut peak_scores,
        &mut visited,
        &mut bt_candidates,
        opt.min_cnt,
        opt.min_chain_score,
        real_max_drop,
    );

    if n_u == 0 {
        ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates);
        ctx.scores = scores;
        ctx.peak_scores = peak_scores;
        ctx.visited = visited;
        ctx.soa_ref_pos = soa_ref_pos;
        ctx.soa_query_pos = soa_query_pos;
        ctx.soa_query_span = soa_query_span;
        ctx.soa_ref_id_strand = soa_rid_strand;
        ctx.soa_scores_buf = soa_sc_buf;
        ctx.soa_diag = soa_diag;
        return (Vec::new(), Vec::new());
    }

    // compact_a logic — Step 1: Write chain anchors to b[] in forward order
    let mut b: Vec<Minimizer> = Vec::with_capacity(n_v);
    let mut k = 0usize;
    for &u_val in &u[..n_u] {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        let k0 = k;
        for j in 0..ni {
            let idx = peak_scores[k0 + (ni - j - 1)] as usize;
            b.push(a[idx]);
            k += 1;
        }
    }

    // Step 2: Sort chains by target position of first anchor
    let mut w: Vec<Minimizer> = Vec::with_capacity(n_u);
    let mut k_pos = 0usize;
    for (i, &u_val) in u[..n_u].iter().enumerate() {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        w.push(Minimizer {
            x: b[k_pos].x,
            y: ((k_pos as u64) << 32) | (i as u64),
        });
        k_pos += ni;
    }
    radix_sort_128x(&mut w);

    // Step 3: Reorder u[] and anchors according to sorted order
    let mut u2: Vec<u64> = Vec::with_capacity(n_u);
    let mut b2: Vec<Minimizer> = Vec::with_capacity(n_v);
    for &w_val in &w[..n_u] {
        let j = (w_val.y & 0xFFFFFFFF) as usize;
        let offset = (w_val.y >> 32) as usize;
        let ni = (u[j] & 0xFFFFFFFF) as usize;
        u2.push(u[j]);
        for idx in 0..ni {
            b2.push(b[offset + idx]);
        }
    }

    // Return buffers to context for reuse
    ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates);
    ctx.scores = scores;
    ctx.peak_scores = peak_scores;
    ctx.visited = visited;
    ctx.soa_ref_pos = soa_ref_pos;
    ctx.soa_query_pos = soa_query_pos;
    ctx.soa_query_span = soa_query_span;
    ctx.soa_ref_id_strand = soa_rid_strand;
    ctx.soa_scores_buf = soa_sc_buf;
    ctx.soa_diag = soa_diag;

    (u2, b2)
}}

// =============================================================================
// SSE (x86_64) implementation — 4-wide, for CPUs without AVX2
// =============================================================================

/// Vectorized fast_log2 for 4 floats (SSE2).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[inline]
unsafe fn simd_fast_log2_sse(x: __m128) -> __m128 {
    let z = _mm_castps_si128(x);
    let exp_raw = _mm_and_si128(_mm_srli_epi32(z, 23), _mm_set1_epi32(255));
    let exp = _mm_sub_epi32(exp_raw, _mm_set1_epi32(128));
    let log_2 = _mm_cvtepi32_ps(exp);
    let mant_bits = _mm_or_si128(
        _mm_and_si128(z, _mm_set1_epi32(!(255i32 << 23))),
        _mm_set1_epi32(127i32 << 23),
    );
    let z_f = _mm_castsi128_ps(mant_bits);
    let c0 = _mm_set1_ps(-0.34484843f32);
    let c1 = _mm_set1_ps(2.024_665_8f32);
    let c2 = _mm_set1_ps(-0.674_877_6f32);
    let t1 = _mm_mul_ps(c0, z_f);
    let t2 = _mm_add_ps(t1, c1);
    let t3 = _mm_mul_ps(t2, z_f);
    let poly = _mm_add_ps(t3, c2);
    _mm_add_ps(log_2, poly)
}

/// SSE2 helper: absolute value of i32x4 (no _mm_abs_epi32 in SSE2).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[inline]
unsafe fn sse2_abs_epi32(v: __m128i) -> __m128i {
    let neg = _mm_sub_epi32(_mm_setzero_si128(), v);
    let mask = _mm_cmpgt_epi32(v, _mm_setzero_si128());
    _mm_or_si128(_mm_and_si128(mask, v), _mm_andnot_si128(mask, neg))
}

/// SSE2 helper: min of i32x4 (no _mm_min_epi32 in SSE2, added in SSE4.1).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[inline]
unsafe fn sse2_min_epi32(a: __m128i, b: __m128i) -> __m128i {
    let mask = _mm_cmpgt_epi32(a, b); // mask = a > b
    _mm_or_si128(_mm_and_si128(mask, b), _mm_andnot_si128(mask, a))
}

/// Batch-compute chain scores for 4 predecessors at a time (SSE2).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[allow(clippy::too_many_arguments)]
unsafe fn compute_chain_scores_batch_sse(
    qi: i32, ri: i32, rid_strand_i: u32,
    soa_ref_pos: &[i32], soa_query_pos: &[i32],
    soa_query_span: &[i32], soa_rid_strand: &[u32],
    max_dist_x: i32, max_dist_y: i32, bw: i32,
    chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], start: usize, end: usize,
) { unsafe {
    if start >= end { return; }

    let qi_v = _mm_set1_epi32(qi);
    let ri_v = _mm_set1_epi32(ri);
    let rid_i_v = _mm_set1_epi32(rid_strand_i as i32);
    let max_dx_v = _mm_set1_epi32(max_dist_x);
    let max_dy_v = _mm_set1_epi32(max_dist_y);
    let bw_v = _mm_set1_epi32(bw);
    let zero_v = _mm_setzero_si128();
    let one_v = _mm_set1_epi32(1);
    let min_val_v = _mm_set1_epi32(i32::MIN);
    let pen_gap_v = _mm_set1_ps(chn_pen_gap);
    let pen_skip_v = _mm_set1_ps(chn_pen_skip);
    let half_v = _mm_set1_ps(0.5f32);

    let aligned_end = start + ((end - start) / 4) * 4;
    let mut j = start;

    while j < aligned_end {
        let qj_v = _mm_loadu_si128(soa_query_pos.as_ptr().add(j) as *const __m128i);
        let rj_v = _mm_loadu_si128(soa_ref_pos.as_ptr().add(j) as *const __m128i);
        let span_v = _mm_loadu_si128(soa_query_span.as_ptr().add(j) as *const __m128i);
        let rid_j_v = _mm_loadu_si128(soa_rid_strand.as_ptr().add(j) as *const __m128i);

        let qdiff = _mm_sub_epi32(qi_v, qj_v);
        let rdiff = _mm_sub_epi32(ri_v, rj_v);

        // Build validity mask (SSE2: vector masks, same as AVX2 style)
        let mut valid = _mm_cmpgt_epi32(qdiff, zero_v);
        valid = _mm_andnot_si128(_mm_cmpgt_epi32(qdiff, max_dx_v), valid);
        valid = _mm_and_si128(valid, _mm_cmpeq_epi32(rid_i_v, rid_j_v));
        valid = _mm_andnot_si128(_mm_cmpeq_epi32(rdiff, zero_v), valid);
        valid = _mm_andnot_si128(_mm_cmpgt_epi32(qdiff, max_dy_v), valid);

        let diff = _mm_sub_epi32(rdiff, qdiff);
        let gap_w = sse2_abs_epi32(diff);
        valid = _mm_andnot_si128(_mm_cmpgt_epi32(gap_w, bw_v), valid);

        let min_diff = sse2_min_epi32(rdiff, qdiff);
        let mut score = sse2_min_epi32(span_v, min_diff);

        // Penalty
        let gw_gt_zero = _mm_cmpgt_epi32(gap_w, zero_v);
        let md_gt_span = _mm_cmpgt_epi32(min_diff, span_v);
        let needs_pen = _mm_or_si128(gw_gt_zero, md_gt_span);

        let gap_w_f = _mm_cvtepi32_ps(gap_w);
        let min_diff_f = _mm_cvtepi32_ps(min_diff);
        let lin_pen = _mm_add_ps(_mm_mul_ps(pen_gap_v, gap_w_f), _mm_mul_ps(pen_skip_v, min_diff_f));

        let gw_plus1_f = _mm_cvtepi32_ps(_mm_add_epi32(gap_w, one_v));
        let log2_val = simd_fast_log2_sse(gw_plus1_f);
        let log_pen = _mm_and_ps(_mm_castsi128_ps(gw_gt_zero), log2_val);

        let total_pen = _mm_add_ps(lin_pen, _mm_mul_ps(half_v, log_pen));
        let pen_i32 = _mm_cvttps_epi32(total_pen);

        let pen_masked = _mm_and_si128(needs_pen, pen_i32);
        score = _mm_sub_epi32(score, pen_masked);

        // Blend: valid ? score : MIN — SSE2 doesn't have _mm_blendv_epi8,
        // so use and/andnot/or
        score = _mm_or_si128(_mm_and_si128(valid, score), _mm_andnot_si128(valid, min_val_v));

        _mm_storeu_si128(sc_buf.as_mut_ptr().add(j) as *mut __m128i, score);
        j += 4;
    }

    // Scalar remainder
    while j < end {
        let qj = soa_query_pos[j];
        let rj = soa_ref_pos[j];
        let span_j = soa_query_span[j];
        let rid_j = soa_rid_strand[j];
        let query_diff = qi.wrapping_sub(qj);
        if query_diff <= 0 || query_diff > max_dist_x {
            sc_buf[j] = i32::MIN; j += 1; continue;
        }
        let ref_diff = ri.wrapping_sub(rj);
        if rid_strand_i == rid_j && (ref_diff == 0 || query_diff > max_dist_y) {
            sc_buf[j] = i32::MIN; j += 1; continue;
        }
        let gap_width = (ref_diff - query_diff).abs();
        if rid_strand_i == rid_j && gap_width > bw {
            sc_buf[j] = i32::MIN; j += 1; continue;
        }
        let min_diff = ref_diff.min(query_diff);
        let mut sc = span_j.min(min_diff);
        if gap_width > 0 || min_diff > span_j {
            let lin_pen = chn_pen_gap * (gap_width as f32) + chn_pen_skip * (min_diff as f32);
            let log_pen = if gap_width >= 1 { fast_log2((gap_width + 1) as f32) } else { 0.0f32 };
            sc -= (lin_pen + 0.5f32 * log_pen) as i32;
        }
        sc_buf[j] = sc;
        j += 1;
    }
}}

/// SSE2 SIMD-optimized chaining — 4-wide, for x86_64 CPUs without AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
pub(crate) unsafe fn chain_anchors_sse(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) { unsafe {
    let n = a.len();
    if n == 0 { return (Vec::new(), Vec::new()); }

    let bw = opt.bandwidth;
    let mut max_dist_x = max_dist_x;
    let mut max_dist_y = max_dist_y;
    if max_dist_x < bw { max_dist_x = bw; }
    if max_dist_y < bw { max_dist_y = bw; }

    let real_max_drop = bw;
    let no_chain_skip = opt.max_chain_skip <= 0;

    let mut soa_ref_pos = std::mem::take(&mut ctx.soa_ref_pos);
    let mut soa_query_pos = std::mem::take(&mut ctx.soa_query_pos);
    let mut soa_query_span = std::mem::take(&mut ctx.soa_query_span);
    let mut soa_rid_strand = std::mem::take(&mut ctx.soa_ref_id_strand);
    let mut soa_sc_buf = std::mem::take(&mut ctx.soa_scores_buf);

    soa_ref_pos.resize(n, 0); soa_query_pos.resize(n, 0);
    soa_query_span.resize(n, 0); soa_rid_strand.resize(n, 0);
    soa_sc_buf.resize(n, 0);

    for j in 0..n {
        soa_ref_pos[j] = a[j].ref_pos();
        soa_query_pos[j] = a[j].query_pos();
        soa_query_span[j] = a[j].query_span();
        soa_rid_strand[j] = (a[j].x >> 32) as u32;
    }

    let mut predecessors = std::mem::take(&mut ctx.predecessors); let mut bt_candidates = std::mem::take(&mut ctx.bt_candidates);
    let mut scores = std::mem::take(&mut ctx.scores);
    let mut peak_scores = std::mem::take(&mut ctx.peak_scores);
    let mut visited = std::mem::take(&mut ctx.visited);
    predecessors.resize(n, 0i32); scores.resize(n, 0i32);
    peak_scores.resize(n, 0i32); visited.clear(); visited.resize(n, 0i32);

    let mut global_max_score = 0;
    let mut window_start: usize = 0;
    let mut best_anchor_idx: i64 = -1;

    for i in 0..n {
        let mut best_predecessor: i64 = -1;
        let mut best_score = a[i].query_span();
        let mut skip_count: i32 = 0;

        while window_start < i
            && (a[i].ref_id_strand() != a[window_start].ref_id_strand()
                || (a[i].x > a[window_start].x + max_dist_x as u64))
        { window_start += 1; }

        if !no_chain_skip && (i - window_start) > opt.max_chain_iter as usize {
            window_start = i - opt.max_chain_iter as usize;
        }

        let mut end_j: i64 = if window_start > 0 { window_start as i64 - 1 } else { -1 };

        if no_chain_skip {
            // No early termination possible — compute all scores, then scan
            if i > window_start {
                compute_chain_scores_batch_sse(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, window_start, i,
                );
            }
            for j in window_start..i {
                let sc = soa_sc_buf[j];
                if sc == i32::MIN { continue; }
                let total_sc = sc + scores[j];
                if total_sc > best_score { best_score = total_sc; best_predecessor = j as i64; }
            }
        } else {
            // Chunked reverse scan: compute batch scores for each chunk, then
            // scan with max_chain_skip early termination. Only compute the next
            // chunk if skip hasn't fired — avoids wasting SIMD on predecessors
            // that will never be examined.
            const CHUNK_SIZE: usize = 64;
            let mut chunk_end = i;
            let mut skip_fired = false;
            while chunk_end > window_start {
                let chunk_start = if chunk_end >= window_start + CHUNK_SIZE {
                    chunk_end - CHUNK_SIZE
                } else {
                    window_start
                };

                compute_chain_scores_batch_sse(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, chunk_start, chunk_end,
                );

                for j in (chunk_start..chunk_end).rev() {
                    let sc = soa_sc_buf[j];
                    if sc == i32::MIN { continue; }
                    let total_sc = sc + scores[j];
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if visited[j] == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip {
                            end_j = j as i64;
                            skip_fired = true;
                            break;
                        }
                    }
                    if predecessors[j] >= 0 { visited[predecessors[j] as usize] = i as i32; }
                }

                if skip_fired { break; }
                chunk_end = chunk_start;
            }

            if best_anchor_idx < 0 || a[i].x - a[best_anchor_idx as usize].x > max_dist_x as u64 {
                let mut best = i32::MIN; best_anchor_idx = -1;
                for j in (window_start..i).rev() {
                    if best < scores[j] { best = scores[j]; best_anchor_idx = j as i64; }
                }
            }

            if best_anchor_idx >= 0 && best_anchor_idx < end_j {
                let sc = compute_chain_score(&a[i], &a[best_anchor_idx as usize],
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip);
                if sc != i32::MIN && best_score < sc + scores[best_anchor_idx as usize] {
                    best_score = sc + scores[best_anchor_idx as usize];
                    best_predecessor = best_anchor_idx;
                }
            }
        }

        scores[i] = best_score;
        predecessors[i] = best_predecessor as i32;
        peak_scores[i] = if best_predecessor >= 0
            && peak_scores[best_predecessor as usize] > best_score
        { peak_scores[best_predecessor as usize] } else { best_score };

        if best_anchor_idx < 0
            || (a[i].x - a[best_anchor_idx as usize].x <= max_dist_x as u64
                && scores[best_anchor_idx as usize] < scores[i])
        { best_anchor_idx = i as i64; }

        if global_max_score < best_score { global_max_score = best_score; }
    }

    let (u, n_u, n_v) = chain_backtrack(n, &scores, &predecessors, &mut peak_scores,
        &mut visited, &mut bt_candidates, opt.min_cnt, opt.min_chain_score, real_max_drop);

    if n_u == 0 {
        ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates); ctx.scores = scores;
        ctx.peak_scores = peak_scores; ctx.visited = visited;
        ctx.soa_ref_pos = soa_ref_pos; ctx.soa_query_pos = soa_query_pos;
        ctx.soa_query_span = soa_query_span; ctx.soa_ref_id_strand = soa_rid_strand;
        ctx.soa_scores_buf = soa_sc_buf;
        return (Vec::new(), Vec::new());
    }

    let mut b: Vec<Minimizer> = Vec::with_capacity(n_v);
    let mut k = 0usize;
    for &u_val in &u[..n_u] {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        let k0 = k;
        for j in 0..ni { b.push(a[peak_scores[k0 + (ni - j - 1)] as usize]); k += 1; }
    }
    let mut w: Vec<Minimizer> = Vec::with_capacity(n_u);
    let mut k_pos = 0usize;
    for (i, &u_val) in u[..n_u].iter().enumerate() {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        w.push(Minimizer { x: b[k_pos].x, y: ((k_pos as u64) << 32) | (i as u64) });
        k_pos += ni;
    }
    radix_sort_128x(&mut w);
    let mut u2: Vec<u64> = Vec::with_capacity(n_u);
    let mut b2: Vec<Minimizer> = Vec::with_capacity(n_v);
    for &w_val in &w[..n_u] {
        let j = (w_val.y & 0xFFFFFFFF) as usize;
        let offset = (w_val.y >> 32) as usize;
        let ni = (u[j] & 0xFFFFFFFF) as usize;
        u2.push(u[j]);
        for idx in 0..ni { b2.push(b[offset + idx]); }
    }

    ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates); ctx.scores = scores;
    ctx.peak_scores = peak_scores; ctx.visited = visited;
    ctx.soa_ref_pos = soa_ref_pos; ctx.soa_query_pos = soa_query_pos;
    ctx.soa_query_span = soa_query_span; ctx.soa_ref_id_strand = soa_rid_strand;
    ctx.soa_scores_buf = soa_sc_buf;
    (u2, b2)
}}

// =============================================================================
// AVX-512 (x86_64) implementation — 16-wide, mirrors AVX2 logic
// =============================================================================

/// Vectorized fast_log2 for 16 floats (AVX-512).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
#[inline]
unsafe fn simd_fast_log2_avx512(x: __m512) -> __m512 {
    let z = _mm512_castps_si512(x);
    let exp_raw = _mm512_and_si512(_mm512_srli_epi32(z, 23), _mm512_set1_epi32(255));
    let exp = _mm512_sub_epi32(exp_raw, _mm512_set1_epi32(128));
    let log_2 = _mm512_cvtepi32_ps(exp);
    let mant_bits = _mm512_or_si512(
        _mm512_and_si512(z, _mm512_set1_epi32(!(255i32 << 23))),
        _mm512_set1_epi32(127i32 << 23),
    );
    let z_f = _mm512_castsi512_ps(mant_bits);
    let c0 = _mm512_set1_ps(-0.34484843f32);
    let c1 = _mm512_set1_ps(2.024_665_8f32);
    let c2 = _mm512_set1_ps(-0.674_877_6f32);
    let t1 = _mm512_mul_ps(c0, z_f);
    let t2 = _mm512_add_ps(t1, c1);
    let t3 = _mm512_mul_ps(t2, z_f);
    let poly = _mm512_add_ps(t3, c2);
    _mm512_add_ps(log_2, poly)
}

/// Batch-compute chain scores for 16 predecessors at a time (AVX-512).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f", enable = "avx512bw")]
#[allow(clippy::too_many_arguments)]
unsafe fn compute_chain_scores_batch_avx512(
    qi: i32, ri: i32, rid_strand_i: u32,
    soa_ref_pos: &[i32], soa_query_pos: &[i32],
    soa_query_span: &[i32], soa_rid_strand: &[u32],
    max_dist_x: i32, max_dist_y: i32, bw: i32,
    chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], start: usize, end: usize,
) { unsafe {
    if start >= end { return; }

    let qi_v = _mm512_set1_epi32(qi);
    let ri_v = _mm512_set1_epi32(ri);
    let rid_i_v = _mm512_set1_epi32(rid_strand_i as i32);
    let max_dx_v = _mm512_set1_epi32(max_dist_x);
    let max_dy_v = _mm512_set1_epi32(max_dist_y);
    let bw_v = _mm512_set1_epi32(bw);
    let zero_v = _mm512_setzero_si512();
    let one_v = _mm512_set1_epi32(1);
    let min_val_v = _mm512_set1_epi32(i32::MIN);
    let pen_gap_v = _mm512_set1_ps(chn_pen_gap);
    let pen_skip_v = _mm512_set1_ps(chn_pen_skip);
    let half_v = _mm512_set1_ps(0.5f32);

    // Process 16 predecessors at a time
    let aligned_end = start + ((end - start) / 16) * 16;
    let mut j = start;

    while j < aligned_end {
        let qj_v = _mm512_loadu_si512(soa_query_pos.as_ptr().add(j) as *const __m512i);
        let rj_v = _mm512_loadu_si512(soa_ref_pos.as_ptr().add(j) as *const __m512i);
        let span_v = _mm512_loadu_si512(soa_query_span.as_ptr().add(j) as *const __m512i);
        let rid_j_v = _mm512_loadu_si512(soa_rid_strand.as_ptr().add(j) as *const __m512i);

        let qdiff = _mm512_sub_epi32(qi_v, qj_v);
        let rdiff = _mm512_sub_epi32(ri_v, rj_v);

        // AVX-512 comparisons return __mmask16 instead of vector masks
        let mut valid: __mmask16 = _mm512_cmpgt_epi32_mask(qdiff, zero_v);
        valid &= !_mm512_cmpgt_epi32_mask(qdiff, max_dx_v);
        valid &= _mm512_cmpeq_epi32_mask(rid_i_v, rid_j_v);
        valid &= !_mm512_cmpeq_epi32_mask(rdiff, zero_v);
        valid &= !_mm512_cmpgt_epi32_mask(qdiff, max_dy_v);

        let diff = _mm512_sub_epi32(rdiff, qdiff);
        let gap_w = _mm512_abs_epi32(diff);
        valid &= !_mm512_cmpgt_epi32_mask(gap_w, bw_v);

        if valid == 0 {
            _mm512_storeu_si512(sc_buf.as_mut_ptr().add(j) as *mut __m512i, min_val_v);
            j += 16;
            continue;
        }

        let min_diff = _mm512_min_epi32(rdiff, qdiff);
        let mut score = _mm512_min_epi32(span_v, min_diff);

        // Penalty computation
        let needs_pen: __mmask16 = _mm512_cmpgt_epi32_mask(gap_w, zero_v)
            | _mm512_cmpgt_epi32_mask(min_diff, span_v);

        let gap_w_f = _mm512_cvtepi32_ps(gap_w);
        let min_diff_f = _mm512_cvtepi32_ps(min_diff);
        let lin_pen = _mm512_add_ps(
            _mm512_mul_ps(pen_gap_v, gap_w_f),
            _mm512_mul_ps(pen_skip_v, min_diff_f),
        );

        let gw_plus1_f = _mm512_cvtepi32_ps(_mm512_add_epi32(gap_w, one_v));
        let log2_val = simd_fast_log2_avx512(gw_plus1_f);
        // Zero out log_pen where gap_width == 0
        let gw_gt_zero = _mm512_cmpgt_epi32_mask(gap_w, zero_v);
        let log_pen = _mm512_maskz_mov_ps(gw_gt_zero, log2_val);

        let total_pen = _mm512_add_ps(lin_pen, _mm512_mul_ps(half_v, log_pen));
        let pen_i32 = _mm512_cvttps_epi32(total_pen);

        // Apply penalty only where needs_pen
        let pen_masked = _mm512_maskz_mov_epi32(needs_pen, pen_i32);
        score = _mm512_sub_epi32(score, pen_masked);

        // Mask invalid scores to i32::MIN
        score = _mm512_mask_blend_epi32(valid, min_val_v, score);

        _mm512_storeu_si512(sc_buf.as_mut_ptr().add(j) as *mut __m512i, score);
        j += 16;
    }

    // Tails go through the AVX2 kernel, so its eight-lane/scalar boundaries are reproduced.
    if j < end {
        compute_chain_scores_batch_avx2(qi, ri, rid_strand_i,
            soa_ref_pos, soa_query_pos, soa_query_span, soa_rid_strand,
            max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, sc_buf, j, end);
    }
}}

/// `compute_chain_scores_batch_avx512` for the bucket candidates `idx`, gathered straight from
/// the source arrays instead of copied into per-chunk field arrays. All candidates share anchor
/// i's reference and strand (the anchor window guarantees it), and `diag[k]` is
/// ref[k] - query[k] modulo 2^32, so qdiff + (diag_i - diag_k) is exactly the wrapping reference
/// difference. With `UNIFORM`, every span equals `uniform_span` and its gather is skipped.
/// Tails go through the AVX2 kernel, as in the contiguous kernels.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn chain_scores_bucket_direct_avx512<const UNIFORM: bool>(
    qi: i32, ri: i32, rid: u32, qry: &[i32], span: &[i32], diag: &[i32], idx: &[u32],
    max_dist_x: i32, max_dist_y: i32, bw: i32, chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], uniform_span: i32,
) { unsafe {
    debug_assert_eq!(qry.len(), span.len());
    debug_assert_eq!(qry.len(), diag.len());
    debug_assert!(sc_buf.len() >= idx.len());
    debug_assert!(idx.iter().all(|&k| (k as usize) < qry.len() && k <= i32::MAX as u32));
    let diag_i = ri.wrapping_sub(qi);
    let qi_v = _mm512_set1_epi32(qi);
    let diag_i_v = _mm512_set1_epi32(diag_i);
    let max_dx_v = _mm512_set1_epi32(max_dist_x);
    let max_dy_v = _mm512_set1_epi32(max_dist_y);
    let bw_v = _mm512_set1_epi32(bw);
    let zero_v = _mm512_setzero_si512();
    let one_v = _mm512_set1_epi32(1);
    let min_val_v = _mm512_set1_epi32(i32::MIN);
    let pen_gap_v = _mm512_set1_ps(chn_pen_gap);
    let pen_skip_v = _mm512_set1_ps(chn_pen_skip);
    let half_v = _mm512_set1_ps(0.5f32);
    let aligned_end = idx.len() / 16 * 16;
    let mut j = 0;
    while j < aligned_end {
        let ix = _mm512_loadu_si512(idx.as_ptr().add(j).cast());
        let qj = _mm512_i32gather_epi32::<4>(ix, qry.as_ptr());
        let dj = _mm512_i32gather_epi32::<4>(ix, diag.as_ptr());
        let span_v = if UNIFORM {
            _mm512_set1_epi32(uniform_span)
        } else {
            _mm512_i32gather_epi32::<4>(ix, span.as_ptr())
        };
        let qdiff = _mm512_sub_epi32(qi_v, qj);
        let delta = _mm512_sub_epi32(diag_i_v, dj);
        let rdiff = _mm512_add_epi32(qdiff, delta);
        let gap_w = _mm512_abs_epi32(delta);
        let mut valid = _mm512_cmpgt_epi32_mask(qdiff, zero_v);
        valid &= !_mm512_cmpgt_epi32_mask(qdiff, max_dx_v);
        valid &= !_mm512_cmpgt_epi32_mask(qdiff, max_dy_v);
        valid &= !_mm512_cmpeq_epi32_mask(rdiff, zero_v);
        valid &= !_mm512_cmpgt_epi32_mask(gap_w, bw_v);
        if valid == 0 {
            _mm512_storeu_si512(sc_buf.as_mut_ptr().add(j).cast(), min_val_v);
            j += 16;
            continue;
        }
        let min_diff = _mm512_min_epi32(rdiff, qdiff);
        let mut score = _mm512_min_epi32(span_v, min_diff);
        let gw_gt_zero = _mm512_cmpgt_epi32_mask(gap_w, zero_v);
        let needs_pen = gw_gt_zero | _mm512_cmpgt_epi32_mask(min_diff, span_v);
        let lin_pen = _mm512_add_ps(_mm512_mul_ps(pen_gap_v, _mm512_cvtepi32_ps(gap_w)),
                                    _mm512_mul_ps(pen_skip_v, _mm512_cvtepi32_ps(min_diff)));
        let log2_val = simd_fast_log2_avx512(_mm512_cvtepi32_ps(_mm512_add_epi32(gap_w, one_v)));
        let log_pen = _mm512_maskz_mov_ps(gw_gt_zero, log2_val);
        let total_pen = _mm512_add_ps(lin_pen, _mm512_mul_ps(half_v, log_pen));
        let pen_i32 = _mm512_cvttps_epi32(total_pen);
        score = _mm512_sub_epi32(score, _mm512_maskz_mov_epi32(needs_pen, pen_i32));
        score = _mm512_mask_blend_epi32(valid, min_val_v, score);
        _mm512_storeu_si512(sc_buf.as_mut_ptr().add(j).cast(), score);
        j += 16;
    }
    if j < idx.len() {
        // Rebuild the tail's fields (reference positions from the diagonals) for the AVX2 kernel.
        let m = idx.len() - j;
        let mut tail_ref = [0i32; 16];
        let mut tail_qry = [0i32; 16];
        let mut tail_span = [0i32; 16];
        let tail_rid = [rid; 16];
        for t in 0..m {
            let k = *idx.get_unchecked(j + t) as usize;
            let qk = *qry.get_unchecked(k);
            let delta = diag_i.wrapping_sub(*diag.get_unchecked(k));
            tail_qry[t] = qk;
            tail_ref[t] = ri.wrapping_sub(qi.wrapping_sub(qk).wrapping_add(delta));
            tail_span[t] = if UNIFORM { uniform_span } else { *span.get_unchecked(k) };
        }
        compute_chain_scores_batch_avx2(qi, ri, rid,
            &tail_ref[..m], &tail_qry[..m], &tail_span[..m], &tail_rid[..m],
            max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, &mut sc_buf[j..], 0, m);
    }
}}

/// The span shared by all anchors, or -1 if they differ (anchor spans are 8-bit, never negative).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn bucket_uniform_span(span: &[i32]) -> i32 { unsafe {
    let Some(&first) = span.first() else { return -1 };
    let first_v = _mm512_set1_epi32(first);
    let mut at = 0;
    while at + 16 <= span.len() {
        let lanes = _mm512_loadu_si512(span.as_ptr().add(at).cast());
        if _mm512_cmpeq_epi32_mask(first_v, lanes) != u16::MAX {
            return -1;
        }
        at += 16;
    }
    if span[at..].iter().all(|&s| s == first) { first } else { -1 }
}}

/// Bulk step of the three-way bucket merge. If the last sixteen entries before `cur` in `list`
/// all exceed `floor` (the larger of the other lists' heads and window_start - 1), the merge
/// would take them consecutively, so they are band-filtered in one step and the cursor moves past
/// them. The compress-store keeps lane order, so both the output and the cursor match the scalar
/// merge. Returns false, doing nothing, if the step does not apply.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,popcnt")]
#[inline]
unsafe fn bucket_filter_run16(
    list: *const u32, cur: &mut usize, floor: i64, diags: &[i64], di: i64, bw: i32, out: &mut Vec<u32>,
) -> bool { unsafe {
    if *cur < 16 || out.len() > 48 {
        return false;
    }
    if *list.add(*cur - 16) as i64 <= floor {
        return false;
    }
    debug_assert!(out.capacity() >= 64);
    let reverse = _mm512_setr_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let ix = _mm512_permutexvar_epi32(reverse, _mm512_loadu_si512(list.add(*cur - 16).cast()));
    let low = _mm512_i32gather_epi64::<8>(_mm512_castsi512_si256(ix), diags.as_ptr());
    let high = _mm512_i32gather_epi64::<8>(_mm512_extracti64x4_epi64::<1>(ix), diags.as_ptr());
    let di_v = _mm512_set1_epi64(di);
    let bw_v = _mm512_set1_epi64(bw as i64);
    let m0 = !_mm512_cmpgt_epi64_mask(_mm512_abs_epi64(_mm512_sub_epi64(low, di_v)), bw_v);
    let m1 = !_mm512_cmpgt_epi64_mask(_mm512_abs_epi64(_mm512_sub_epi64(high, di_v)), bw_v);
    let mask = m0 as u16 | ((m1 as u16) << 8);
    let n = out.len();
    _mm512_mask_compressstoreu_epi32(out.as_mut_ptr().add(n).cast(), mask, ix);
    out.set_len(n + mask.count_ones() as usize);
    *cur -= 16;
    true
}}

/// `bucket_filter_run16` on the i32 (wrapping) diagonals: one 16-lane gather instead of two
/// 8-lane i64 gathers. Exact for 0 <= bw <= 2^30: list entries lie in the buckets next to
/// anchor i's, each at most 2^30 wide, so their true diagonal difference from anchor i fits in
/// i32 (even when the diagonals themselves do not) and the wrapping difference equals it.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,popcnt")]
#[inline]
unsafe fn bucket_filter_run16_diag32(
    list: *const u32, cur: &mut usize, floor: i64, diags: &[i32], di: i32, bw: i32, out: &mut Vec<u32>,
) -> bool { unsafe {
    if *cur < 16 || out.len() > 48 {
        return false;
    }
    if *list.add(*cur - 16) as i64 <= floor {
        return false;
    }
    debug_assert!(out.capacity() >= 64);
    let reverse = _mm512_setr_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let ix = _mm512_permutexvar_epi32(reverse, _mm512_loadu_si512(list.add(*cur - 16).cast()));
    let dj = _mm512_i32gather_epi32::<4>(ix, diags.as_ptr());
    let distance = _mm512_abs_epi32(_mm512_sub_epi32(dj, _mm512_set1_epi32(di)));
    let mask = !_mm512_cmpgt_epi32_mask(distance, _mm512_set1_epi32(bw));
    let n = out.len();
    _mm512_mask_compressstoreu_epi32(out.as_mut_ptr().add(n).cast(), mask, ix);
    out.set_len(n + mask.count_ones() as usize);
    *cur -= 16;
    true
}}

/// The (up to) sixteen largest of the first `n` entries of an ascending list, in descending
/// order; missing lanes are -1, below every index.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
#[inline]
unsafe fn bucket_head16(list: *const u32, n: usize) -> __m512i { unsafe {
    if n == 0 {
        return _mm512_set1_epi32(-1);
    }
    let count = n.min(16);
    let mask = if count == 16 { u16::MAX } else { (1u16 << count) - 1 };
    let raw = _mm512_mask_loadu_epi32(_mm512_set1_epi32(-1), mask, list.add(n - count).cast());
    let lanes = _mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    let order = _mm512_sub_epi32(_mm512_set1_epi32(count as i32 - 1), lanes);
    _mm512_mask_permutexvar_epi32(_mm512_set1_epi32(-1), mask, order, raw)
}}

/// The sixteen largest lanes of two descending vectors, descending (bitonic merge).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
#[inline]
unsafe fn merge_desc16(a: __m512i, b: __m512i) -> __m512i {
    let reverse = _mm512_setr_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    let mut v = _mm512_max_epi32(a, _mm512_permutexvar_epi32(reverse, b));
    // Compare-exchange each lane with its partner `ix`: lanes in `mask` keep the maximum.
    macro_rules! step {
        ($ix:expr, $mask:expr) => {{
            let t = _mm512_permutexvar_epi32($ix, v);
            v = _mm512_mask_max_epi32(_mm512_min_epi32(v, t), $mask, v, t);
        }};
    }
    step!(_mm512_setr_epi32(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7), 0x00ff);
    step!(_mm512_setr_epi32(4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11), 0x0f0f);
    step!(_mm512_setr_epi32(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13), 0x3333);
    step!(_mm512_setr_epi32(1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14), 0x5555);
    v
}

/// Sixteen steps of the three-way bucket merge at once: takes the sixteen largest remaining
/// entries across the (disjoint, ascending) lists, moves the cursors past them, and appends
/// those at or after `window_start` and in band to `out`, in descending order. Requires at least
/// sixteen remaining entries in total, entries addressable as i32, and 0 <= bw <= 2^30 (see
/// `bucket_filter_run16_diag32`).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,popcnt")]
#[inline]
unsafe fn bucket_merge_filter16_diag32(
    lists: [*const u32; 3], cur: &mut [usize; 3], window_start: usize,
    diags: &[i32], di: i32, bw: i32, out: &mut Vec<u32>,
) { unsafe {
    debug_assert!(cur.iter().sum::<usize>() >= 16 && out.capacity() >= out.len() + 16);
    let a = bucket_head16(lists[0], cur[0]);
    let b = bucket_head16(lists[1], cur[1]);
    let c = bucket_head16(lists[2], cur[2]);
    let ix = merge_desc16(merge_desc16(a, b), c);
    // Lane 15 is the smallest entry taken.
    let cutoff = _mm512_set1_epi32(_mm_extract_epi32::<3>(_mm512_extracti32x4_epi32::<3>(ix)));
    cur[0] -= _mm512_cmpge_epi32_mask(a, cutoff).count_ones() as usize;
    cur[1] -= _mm512_cmpge_epi32_mask(b, cutoff).count_ones() as usize;
    cur[2] -= _mm512_cmpge_epi32_mask(c, cutoff).count_ones() as usize;
    let dj = _mm512_i32gather_epi32::<4>(ix, diags.as_ptr());
    let delta = _mm512_abs_epi32(_mm512_sub_epi32(dj, _mm512_set1_epi32(di)));
    let mask = _mm512_cmple_epi32_mask(delta, _mm512_set1_epi32(bw))
        & _mm512_cmpge_epi32_mask(ix, _mm512_set1_epi32(window_start as i32));
    let n = out.len();
    _mm512_mask_compressstoreu_epi32(out.as_mut_ptr().add(n).cast(), mask, ix);
    out.set_len(n + mask.count_ones() as usize);
}}

/// Masked reverse-scan chunk, with the same lane validity and score arithmetic as AVX2.
/// Whole off-band groups need only the diagonal load; a zero valid mask skips penalty
/// arithmetic and score stores. Only set bits in the returned mask may be read.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f", enable = "avx512bw", enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn chain_scores_masked_avx512(
    qi: i32, ri: i32, rid_strand_i: u32, diag_i: i32,
    soa_ref_pos: &[i32], soa_query_pos: &[i32], soa_query_span: &[i32], soa_rid_strand: &[u32], soa_diag: &[i32],
    max_dist_x: i32, max_dist_y: i32, bw: i32, chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], start: usize, end: usize,
) -> u64 { unsafe {
    debug_assert!(end - start <= 64);
    let diag_i_v = _mm512_set1_epi32(diag_i);
    let qi_v = _mm512_set1_epi32(qi);
    let ri_v = _mm512_set1_epi32(ri);
    let rid_i_v = _mm512_set1_epi32(rid_strand_i as i32);
    let max_dx_v = _mm512_set1_epi32(max_dist_x);
    let max_dy_v = _mm512_set1_epi32(max_dist_y);
    let bw_v = _mm512_set1_epi32(bw);
    let zero_v = _mm512_setzero_si512();
    let one_v = _mm512_set1_epi32(1);
    let min_val_v = _mm512_set1_epi32(i32::MIN);
    let pen_gap_v = _mm512_set1_ps(chn_pen_gap);
    let pen_skip_v = _mm512_set1_ps(chn_pen_skip);
    let half_v = _mm512_set1_ps(0.5f32);
    let mut mask = 0u64;
    let aligned_end = start + ((end - start) / 16) * 16;
    let mut j = start;
    while j < aligned_end {
        let dj = _mm512_loadu_si512(soa_diag.as_ptr().add(j) as *const __m512i);
        // Keep the batch kernel's signed abs semantics, including i32::MIN.
        let gap_w = _mm512_abs_epi32(_mm512_sub_epi32(diag_i_v, dj));
        let gw_gt_bw = _mm512_cmpgt_epi32_mask(gap_w, bw_v);
        if gw_gt_bw == u16::MAX {
            j += 16;
            continue;
        }
        let qj_v = _mm512_loadu_si512(soa_query_pos.as_ptr().add(j) as *const __m512i);
        let rj_v = _mm512_loadu_si512(soa_ref_pos.as_ptr().add(j) as *const __m512i);
        let span_v = _mm512_loadu_si512(soa_query_span.as_ptr().add(j) as *const __m512i);
        let rid_j_v = _mm512_loadu_si512(soa_rid_strand.as_ptr().add(j) as *const __m512i);
        let qdiff = _mm512_sub_epi32(qi_v, qj_v);
        let rdiff = _mm512_sub_epi32(ri_v, rj_v);
        let mut valid = _mm512_cmpgt_epi32_mask(qdiff, zero_v);
        valid &= !_mm512_cmpgt_epi32_mask(qdiff, max_dx_v);
        valid &= _mm512_cmpeq_epi32_mask(rid_i_v, rid_j_v);
        valid &= !_mm512_cmpeq_epi32_mask(rdiff, zero_v);
        valid &= !_mm512_cmpgt_epi32_mask(qdiff, max_dy_v);
        valid &= !gw_gt_bw;
        if valid == 0 {
            j += 16;
            continue;
        }
        let min_diff = _mm512_min_epi32(rdiff, qdiff);
        let mut score = _mm512_min_epi32(span_v, min_diff);
        let gw_gt_zero = _mm512_cmpgt_epi32_mask(gap_w, zero_v);
        let needs_pen = gw_gt_zero | _mm512_cmpgt_epi32_mask(min_diff, span_v);
        let lin_pen = _mm512_add_ps(_mm512_mul_ps(pen_gap_v, _mm512_cvtepi32_ps(gap_w)),
                                    _mm512_mul_ps(pen_skip_v, _mm512_cvtepi32_ps(min_diff)));
        let log2_val = simd_fast_log2_avx512(_mm512_cvtepi32_ps(_mm512_add_epi32(gap_w, one_v)));
        let log_pen = _mm512_maskz_mov_ps(gw_gt_zero, log2_val);
        let total_pen = _mm512_add_ps(lin_pen, _mm512_mul_ps(half_v, log_pen));
        let pen_i32 = _mm512_cvttps_epi32(total_pen);
        score = _mm512_sub_epi32(score, _mm512_maskz_mov_epi32(needs_pen, pen_i32));
        score = _mm512_mask_blend_epi32(valid, min_val_v, score);
        _mm512_storeu_si512(sc_buf.as_mut_ptr().add(j) as *mut __m512i, score);
        mask |= (valid as u64) << (j - start);
        j += 16;
    }
    // Tails go through the AVX2 kernel, so its eight-lane/scalar boundaries are reproduced.
    if j < end {
        mask |= chain_scores_masked_avx2(qi, ri, rid_strand_i, diag_i,
            soa_ref_pos, soa_query_pos, soa_query_span, soa_rid_strand, soa_diag,
            max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, sc_buf, j, end)
            << (j - start);
    }
    mask
}}

/// Outcome of `advance_nonimproving_block16`.
#[cfg(target_arch = "x86_64")]
#[derive(Debug, PartialEq, Eq)]
enum BlockWalk {
    /// Some candidate improves the best score: walk the block one candidate at a time.
    Improving,
    /// The whole block was walked.
    Advanced,
    /// The skip limit fired at this candidate.
    Stopped(usize),
}

/// Walks a block of sixteen candidates `start..start + 16` (those set in `mask`) of the reverse
/// scan at once when none of them improves on `best`; then only the skip count and the
/// `visited` stamps change. A predecessor always precedes its anchor, so a stamp set within the
/// block can only affect candidates visited later in the reverse scan: computing the block's own
/// stamps up front is exact, even if the skip limit stops the walk partway. Only the visited
/// candidates' stamps are committed (excluding the one that fires the stop).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,avx512f,avx512bw,bmi2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
unsafe fn advance_nonimproving_block16(
    scb: *const i32, scores: *const i32, pre: *const i32, visited: *mut i32,
    start: usize, mask: u16, best: i32, skip: &mut i32, limit: i32, stamp: i32,
) -> BlockWalk { unsafe {
    // The caller guarantees sixteen readable lanes, every nonnegative predecessor below its own
    // anchor index, 0 <= skip <= limit and limit <= i32::MAX - 16.
    let sc = _mm512_loadu_si512(scb.add(start).cast());
    let good = mask & _mm512_cmpneq_epi32_mask(sc, _mm512_set1_epi32(i32::MIN));
    if good == 0 {
        return BlockWalk::Advanced;
    }
    let total = _mm512_add_epi32(sc, _mm512_loadu_si512(scores.add(start).cast()));
    if _mm512_mask_reduce_max_epi32(good, total) > best {
        return BlockWalk::Improving;
    }
    let preds = _mm512_loadu_si512(pre.add(start).cast());
    let zero = _mm512_setzero_si512();
    let nonnegative = _mm512_cmpge_epi32_mask(preds, zero);
    let local = _mm512_sub_epi32(preds, _mm512_set1_epi32(start as i32));
    let in_block = good & nonnegative & _mm512_cmpge_epi32_mask(local, zero)
        & _mm512_cmplt_epi32_mask(local, _mm512_set1_epi32(16));
    let local_marks = _mm512_reduce_or_epi32(_mm512_maskz_sllv_epi32(
        in_block, _mm512_set1_epi32(1), local)) as u16;
    let current_marks = _mm512_cmpeq_epi32_mask(
        _mm512_loadu_si512(visited.add(start).cast()), _mm512_set1_epi32(stamp));
    let marked = good & (current_marks | local_marks);
    let count = (marked as u32).count_ones();
    let budget = (limit - *skip) as u32;
    let (processed, outcome) = if count > budget {
        // The (budget + 1)th marked candidate from the high end fires the stop.
        let rank = count - budget - 1;
        let bit = _pdep_u32(1u32 << rank, marked as u32).trailing_zeros() as usize;
        let processed = good & !(((1u32 << (bit + 1)) - 1) as u16);
        *skip = limit + 1;
        (processed, BlockWalk::Stopped(start + bit))
    } else {
        *skip += count as i32;
        (good, BlockWalk::Advanced)
    };
    // Duplicate scatter destinations all receive the same stamp, so their order is immaterial.
    _mm512_mask_i32scatter_epi32::<4>(visited.cast(), processed & nonnegative,
        preds, _mm512_set1_epi32(stamp));
    outcome
}}

/// AVX-512 chaining. Its kernels also use BMI2 and POPCNT, which every AVX-512BW CPU has; a CPU
/// without them, or an input too large for i32 gather indices, uses the AVX2 chainer instead.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx2")]
pub(crate) unsafe fn chain_anchors_avx512(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) {
    unsafe {
        if is_x86_feature_detected!("bmi2")
            && is_x86_feature_detected!("popcnt")
            && a.len() <= i32::MAX as usize
        {
            chain_anchors_avx512_impl(opt, max_dist_x, max_dist_y, a, ctx)
        } else {
            chain_anchors_avx2(opt, max_dist_x, max_dist_y, a, ctx)
        }
    }
}

/// AVX-512 chaining with AVX2's scan structure (diagonal buckets for dense reads, masked
/// contiguous scan otherwise), 16 lanes wide, plus whole-block walks of non-improving
/// candidates. Requires `a.len() <= i32::MAX`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx2,bmi2,popcnt")]
unsafe fn chain_anchors_avx512_impl(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) { unsafe {
    let n = a.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    debug_assert!(n <= i32::MAX as usize);

    let bw = opt.bandwidth;
    let mut max_dist_x = max_dist_x;
    let mut max_dist_y = max_dist_y;
    if max_dist_x < bw {
        max_dist_x = bw;
    }
    if max_dist_y < bw {
        max_dist_y = bw;
    }

    let real_max_drop = bw;
    let no_chain_skip = opt.max_chain_skip <= 0;

    // --- Extract SoA fields from AoS Minimizer array ---
    let mut soa_ref_pos = std::mem::take(&mut ctx.soa_ref_pos);
    let mut soa_query_pos = std::mem::take(&mut ctx.soa_query_pos);
    let mut soa_query_span = std::mem::take(&mut ctx.soa_query_span);
    let mut soa_rid_strand = std::mem::take(&mut ctx.soa_ref_id_strand);
    let mut soa_sc_buf = std::mem::take(&mut ctx.soa_scores_buf);
    let mut soa_diag = std::mem::take(&mut ctx.soa_diag);

    soa_ref_pos.resize(n, 0);
    soa_query_pos.resize(n, 0);
    soa_query_span.resize(n, 0);
    soa_rid_strand.resize(n, 0);
    soa_sc_buf.resize(n, 0);
    if !no_chain_skip {
        soa_diag.resize(n, 0);
    }

    for j in 0..n {
        soa_ref_pos[j] = a[j].ref_pos();
        soa_query_pos[j] = a[j].query_pos();
        soa_query_span[j] = a[j].query_span();
        soa_rid_strand[j] = (a[j].x >> 32) as u32;
        if !no_chain_skip {
            soa_diag[j] = soa_ref_pos[j].wrapping_sub(soa_query_pos[j]);
        }
    }

    // --- DP buffers ---
    let mut predecessors = std::mem::take(&mut ctx.predecessors);
    let mut bt_candidates = std::mem::take(&mut ctx.bt_candidates);
    let mut scores = std::mem::take(&mut ctx.scores);
    let mut peak_scores = std::mem::take(&mut ctx.peak_scores);
    let mut visited = std::mem::take(&mut ctx.visited);
    predecessors.resize(n, 0i32);
    scores.resize(n, 0i32);
    peak_scores.resize(n, 0i32);
    visited.clear();
    visited.resize(n, 0i32);

    let mut global_max_score = 0;
    let mut window_start: usize = 0;
    let mut best_anchor_idx: i64 = -1;

    let mut bk = std::mem::take(&mut ctx.diag_buckets);
    let use_buckets = !no_chain_skip && n >= DIAG_BUCKET_MIN_N && {
        let (lo, hi) = soa_query_pos.iter().fold((i32::MAX, i32::MIN), |(l, h), &q| (l.min(q), h.max(q)));
        // The score kernels use wrapping i32 position differences, which the i64 bucket
        // diagonals match only without a signed-coordinate wrap; otherwise use the contiguous
        // scan. The density arithmetic is widened for the same reason.
        lo >= 0
            && n >= DIAG_BUCKET_MIN_DENSITY * (hi as i64 - lo as i64 + 1) as usize
            && soa_ref_pos.iter().all(|&r| r >= 0)
    };
    // The i32-diagonal bucket filters need 0 <= bw <= 2^30 (see bucket_filter_run16_diag32).
    let narrow_buckets = (0..=(1 << 30)).contains(&bw);
    let mut uniform_span = -1;
    if use_buckets {
        bk.build(&soa_ref_pos, &soa_query_pos, bw);
        // The bulk merge steps append up to 16 entries without growing the Vec.
        bk.g_idx.reserve(64usize.saturating_sub(bk.g_idx.len()));
        uniform_span = bucket_uniform_span(&soa_query_span);
    }
    // advance_nonimproving_block16 needs the skip count to stay 16 below i32::MAX.
    let block_walk = opt.max_chain_skip <= i32::MAX - 16;

    for i in 0..n {
        let mut best_predecessor: i64 = -1;
        let mut best_score = a[i].query_span();
        let mut skip_count: i32 = 0;

        while window_start < i
            && (a[i].ref_id_strand() != a[window_start].ref_id_strand()
                || (a[i].x > a[window_start].x + max_dist_x as u64))
        {
            window_start += 1;
        }

        if !no_chain_skip && (i - window_start) > opt.max_chain_iter as usize {
            window_start = i - opt.max_chain_iter as usize;
        }

        let mut end_j: i64 = if window_start > 0 {
            window_start as i64 - 1
        } else {
            -1
        };

        if no_chain_skip {
            if i > window_start {
                compute_chain_scores_batch_avx512(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, window_start, i,
                );
            }
            for j in window_start..i {
                let sc = soa_sc_buf[j];
                if sc == i32::MIN { continue; }
                let total_sc = sc + scores[j];
                if total_sc > best_score { best_score = total_sc; best_predecessor = j as i64; }
            }
        } else if use_buckets && i - window_start > DIAG_BUCKET_MIN_WINDOW {
            // Merge only the in-band predecessors, in the original descending visit order, and
            // score them in chunks (a skip stop discards the rest of a chunk).
            const CHUNK_SIZE: usize = 32;
            let nb = bk.nbr[bk.id[i] as usize];
            let mut lp: [*const u32; 3] = [std::ptr::null(); 3];
            let mut cur = [0usize; 3];
            for t in 0..3 {
                if nb[t] != u32::MAX {
                    let list = &bk.lists[nb[t] as usize];
                    lp[t] = list.as_ptr();
                    cur[t] = list.len();
                }
            }
            // SAFETY: a nonzero cursor belongs to a present bucket, and its list is not
            // modified until this anchor's scan ends. All stored indices precede i.
            let tail = |t: usize, c: usize| -> i64 {
                if c > 0 { *lp[t].add(c - 1) as i64 } else { -1 }
            };
            let (mut v0, mut v1, mut v2) = (tail(0, cur[0]), tail(1, cur[1]), tail(2, cur[2]));
            let di = bk.diag[i];
            let ws = window_start as i64;
            let mut skip_fired = false;
            let mut exhausted = false;
            loop {
                bk.g_idx.clear();
                while bk.g_idx.len() < CHUNK_SIZE {
                    let t = if v0 >= v1 && v0 >= v2 { 0 } else if v1 >= v2 { 1 } else { 2 };
                    let j = [v0, v1, v2][t];
                    if j < ws { exhausted = true; break; }
                    if bk.g_idx.len() + 16 <= CHUNK_SIZE && cur[t] >= 16 {
                        let other = match t { 0 => v1.max(v2), 1 => v0.max(v2), _ => v0.max(v1) };
                        let floor = other.max(ws - 1);
                        // A run of sixteen needs sixteen distinct indices above floor.
                        let took_run = j - floor >= 16 && if narrow_buckets {
                            bucket_filter_run16_diag32(
                                lp[t], &mut cur[t], floor, &soa_diag, di as i32, bw, &mut bk.g_idx)
                        } else {
                            bucket_filter_run16(lp[t], &mut cur[t], floor, &bk.diag, di, bw, &mut bk.g_idx)
                        };
                        if took_run {
                            match t {
                                0 => v0 = tail(0, cur[0]),
                                1 => v1 = tail(1, cur[1]),
                                _ => v2 = tail(2, cur[2]),
                            }
                            continue;
                        }
                    }
                    if narrow_buckets && bk.g_idx.len() + 16 <= CHUNK_SIZE && cur.iter().sum::<usize>() >= 16 {
                        bucket_merge_filter16_diag32(
                            lp, &mut cur, window_start, &soa_diag, di as i32, bw, &mut bk.g_idx);
                        v0 = tail(0, cur[0]);
                        v1 = tail(1, cur[1]);
                        v2 = tail(2, cur[2]);
                        continue;
                    }
                    cur[t] -= 1;
                    match t {
                        0 => v0 = tail(0, cur[0]),
                        1 => v1 = tail(1, cur[1]),
                        _ => v2 = tail(2, cur[2]),
                    }
                    // Out-of-band pairs score i32::MIN and have no skip/visited effects.
                    if (bk.diag[j as usize] - di).abs() > bw as i64 { continue; }
                    bk.g_idx.push(j as u32);
                }
                let m = bk.g_idx.len();
                if m == 0 { break; }
                bk.g_sc.resize(m, 0);
                if uniform_span >= 0 {
                    chain_scores_bucket_direct_avx512::<true>(
                        soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                        &soa_query_pos, &soa_query_span, &soa_diag, &bk.g_idx,
                        max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip,
                        &mut bk.g_sc, uniform_span,
                    );
                } else {
                    chain_scores_bucket_direct_avx512::<false>(
                        soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                        &soa_query_pos, &soa_query_span, &soa_diag, &bk.g_idx,
                        max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip,
                        &mut bk.g_sc, uniform_span,
                    );
                }
                for t in 0..m {
                    let sc = bk.g_sc[t];
                    if sc == i32::MIN { continue; }
                    let j = bk.g_idx[t] as usize;
                    let total_sc = sc + scores[j];
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if visited[j] == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                    }
                    if predecessors[j] >= 0 { visited[predecessors[j] as usize] = i as i32; }
                }
                if skip_fired || exhausted { break; }
            }
        } else {
            const CHUNK_SIZE: usize = 64;
            const FIRST_CHUNK: usize = 32;
            let mut chunk_end = i;
            let mut skip_fired = false;
            let mut csize = FIRST_CHUNK;
            while chunk_end > window_start {
                let chunk_start = if chunk_end >= window_start + csize { chunk_end - csize } else { window_start };
                csize = CHUNK_SIZE;

                let mut valid = chain_scores_masked_avx512(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i], soa_diag[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand, &soa_diag,
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, chunk_start, chunk_end,
                );

                // SAFETY: every mask bit indexes chunk_start..chunk_end <= i < n;
                // predecessors were assigned only to earlier anchors. All four buffers have
                // length n and do not alias. Match AVX2's candidate walk without bounds checks.
                let (scb, scr, vis, pre) =
                    (soa_sc_buf.as_ptr(), scores.as_ptr(), visited.as_mut_ptr(), predecessors.as_ptr());
                // Walk the candidates in descending order, sixteen-lane group by group.
                while valid != 0 {
                    let group = (63 - valid.leading_zeros() as usize) & !15;
                    let mut bits = (valid >> group) as u16;
                    valid &= !(0xffffu64 << group);
                    let start = chunk_start + group;
                    if block_walk && best_predecessor >= 0 && start + 16 <= chunk_end && bits.count_ones() >= 4 {
                        match advance_nonimproving_block16(scb, scr, pre, vis, start,
                            bits, best_score, &mut skip_count, opt.max_chain_skip, i as i32) {
                            BlockWalk::Advanced => continue,
                            BlockWalk::Stopped(j) => { end_j = j as i64; skip_fired = true; break; }
                            BlockWalk::Improving => {}
                        }
                    }
                    while bits != 0 {
                        let bit = 15 - bits.leading_zeros() as usize;
                        bits &= !(1u16 << bit);
                        let j = start + bit;
                        let sc = *scb.add(j);
                        if sc == i32::MIN { continue; }
                        let total_sc = sc + *scr.add(j);
                        if total_sc > best_score {
                            best_score = total_sc; best_predecessor = j as i64;
                            if skip_count > 0 { skip_count -= 1; }
                        } else if *vis.add(j) == i as i32 {
                            skip_count += 1;
                            if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                        }
                        let p = *pre.add(j);
                        if p >= 0 { *vis.add(p as usize) = i as i32; }
                    }
                    if skip_fired { break; }
                }

                if skip_fired { break; }
                chunk_end = chunk_start;
            }
        }

        if !no_chain_skip {
            if best_anchor_idx < 0 || a[i].x - a[best_anchor_idx as usize].x > max_dist_x as u64 {
                let mut best = i32::MIN; best_anchor_idx = -1;
                for j in (window_start..i).rev() {
                    let sj = *scores.get_unchecked(j);
                    let better = best < sj;
                    best_anchor_idx = if better { j as i64 } else { best_anchor_idx };
                    best = if better { sj } else { best };
                }
            }

            if best_anchor_idx >= 0 && best_anchor_idx < end_j {
                let sc = compute_chain_score(&a[i], &a[best_anchor_idx as usize],
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip);
                if sc != i32::MIN && best_score < sc + scores[best_anchor_idx as usize] {
                    best_score = sc + scores[best_anchor_idx as usize];
                    best_predecessor = best_anchor_idx;
                }
            }
        }

        scores[i] = best_score;
        predecessors[i] = best_predecessor as i32;
        peak_scores[i] = if best_predecessor >= 0
            && peak_scores[best_predecessor as usize] > best_score
        {
            peak_scores[best_predecessor as usize]
        } else {
            best_score
        };

        if best_anchor_idx < 0
            || (a[i].x - a[best_anchor_idx as usize].x <= max_dist_x as u64
                && scores[best_anchor_idx as usize] < scores[i])
        {
            best_anchor_idx = i as i64;
        }

        if global_max_score < best_score {
            global_max_score = best_score;
        }
        if use_buckets { bk.lists[bk.id[i] as usize].push(i as u32); }
    }
    ctx.diag_buckets = bk;

    // --- Backtrack and compact (identical to AVX2/scalar) ---
    let (u, n_u, n_v) = chain_backtrack(
        n,
        &scores,
        &predecessors,
        &mut peak_scores,
        &mut visited,
        &mut bt_candidates,
        opt.min_cnt,
        opt.min_chain_score,
        real_max_drop,
    );

    if n_u == 0 {
        ctx.predecessors = predecessors;
        ctx.bt_candidates = std::mem::take(&mut bt_candidates);
        ctx.scores = scores;
        ctx.peak_scores = peak_scores;
        ctx.visited = visited;
        ctx.soa_ref_pos = soa_ref_pos;
        ctx.soa_query_pos = soa_query_pos;
        ctx.soa_query_span = soa_query_span;
        ctx.soa_ref_id_strand = soa_rid_strand;
        ctx.soa_scores_buf = soa_sc_buf;
        ctx.soa_diag = soa_diag;
        return (Vec::new(), Vec::new());
    }

    let mut b: Vec<Minimizer> = Vec::with_capacity(n_v);
    let mut k = 0usize;
    for &u_val in &u[..n_u] {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        let k0 = k;
        for j in 0..ni {
            let idx = peak_scores[k0 + (ni - j - 1)] as usize;
            b.push(a[idx]);
            k += 1;
        }
    }

    let mut w: Vec<Minimizer> = Vec::with_capacity(n_u);
    let mut k_pos = 0usize;
    for (i, &u_val) in u[..n_u].iter().enumerate() {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        w.push(Minimizer {
            x: b[k_pos].x,
            y: ((k_pos as u64) << 32) | (i as u64),
        });
        k_pos += ni;
    }
    radix_sort_128x(&mut w);

    let mut u2: Vec<u64> = Vec::with_capacity(n_u);
    let mut b2: Vec<Minimizer> = Vec::with_capacity(n_v);
    for &w_val in &w[..n_u] {
        let j = (w_val.y & 0xFFFFFFFF) as usize;
        let offset = (w_val.y >> 32) as usize;
        let ni = (u[j] & 0xFFFFFFFF) as usize;
        u2.push(u[j]);
        for idx in 0..ni {
            b2.push(b[offset + idx]);
        }
    }

    ctx.predecessors = predecessors;
    ctx.bt_candidates = std::mem::take(&mut bt_candidates);
    ctx.scores = scores;
    ctx.peak_scores = peak_scores;
    ctx.visited = visited;
    ctx.soa_ref_pos = soa_ref_pos;
    ctx.soa_query_pos = soa_query_pos;
    ctx.soa_query_span = soa_query_span;
    ctx.soa_ref_id_strand = soa_rid_strand;
    ctx.soa_scores_buf = soa_sc_buf;
    ctx.soa_diag = soa_diag;

    (u2, b2)
}}

// =============================================================================
// NEON (aarch64) implementation — 4-wide, mirrors AVX2 logic exactly
// =============================================================================

/// Vectorized fast_log2 for 4 floats (NEON).
/// Matches the scalar `fast_log2` exactly (same polynomial, same rounding).
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[inline]
unsafe fn simd_fast_log2_neon(x: float32x4_t) -> float32x4_t {
    let z = vreinterpretq_s32_f32(x);

    // exponent = ((z >> 23) & 255) - 128
    let exp_raw = vandq_s32(vshrq_n_s32(z, 23), vdupq_n_s32(255));
    let exp = vsubq_s32(exp_raw, vdupq_n_s32(128));
    let log_2 = vcvtq_f32_s32(exp);

    // mantissa: clear exponent bits, set exponent to 127 (1.0 bias)
    let mant_bits = vorrq_s32(
        vandq_s32(z, vdupq_n_s32(!(255i32 << 23))),
        vdupq_n_s32(127i32 << 23),
    );
    let z_f = vreinterpretq_f32_s32(mant_bits);

    // Horner's: ((-0.34484843 * z_f) + 2.0246658) * z_f + (-0.6748776)
    // MUST NOT use fused multiply-add — use explicit mul+add to match scalar rounding.
    let c0 = vdupq_n_f32(-0.34484843f32);
    let c1 = vdupq_n_f32(2.024_665_8f32);
    let c2 = vdupq_n_f32(-0.674_877_6f32);

    let t1 = vmulq_f32(c0, z_f);       // c0 * z_f
    let t2 = vaddq_f32(t1, c1);         // c0 * z_f + c1
    let t3 = vmulq_f32(t2, z_f);        // (c0 * z_f + c1) * z_f
    let poly = vaddq_f32(t3, c2);       // ... + c2

    vaddq_f32(log_2, poly)
}

/// Batch-compute chain scores for predecessors j in [start, end) against anchor i.
/// Writes scores to sc_buf[start..end]. Invalid pairs get i32::MIN.
/// 4-wide NEON version — mirrors `compute_chain_scores_batch_avx2` exactly.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[inline]
unsafe fn compute_chain_scores_batch_neon(
    qi: i32,
    ri: i32,
    rid_strand_i: u32,
    soa_ref_pos: &[i32],
    soa_query_pos: &[i32],
    soa_query_span: &[i32],
    soa_rid_strand: &[u32],
    max_dist_x: i32,
    max_dist_y: i32,
    bw: i32,
    chn_pen_gap: f32,
    chn_pen_skip: f32,
    sc_buf: &mut [i32],
    start: usize,
    end: usize,
) {
    unsafe {
        compute_chain_scores_batch_neon_impl::<false>(
            qi, ri, rid_strand_i, soa_ref_pos, soa_query_pos, soa_query_span,
            soa_rid_strand, max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip,
            sc_buf, start, end, 0, &[],
        );
    }
}

/// NEON counterpart of the AVX2 masked scan: skip off-band groups before loading
/// the remaining fields, and return candidates to visit in the original order.
/// Only scores whose bits are set may be read; skipped groups remain untouched.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[inline]
unsafe fn chain_scores_masked_neon(
    qi: i32, ri: i32, rid_strand_i: u32, diag_i: i32,
    soa_ref_pos: &[i32], soa_query_pos: &[i32], soa_query_span: &[i32],
    soa_rid_strand: &[u32], soa_diag: &[i32],
    max_dist_x: i32, max_dist_y: i32, bw: i32, chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], start: usize, end: usize,
) -> u64 {
    unsafe {
        compute_chain_scores_batch_neon_impl::<true>(
            qi, ri, rid_strand_i, soa_ref_pos, soa_query_pos, soa_query_span,
            soa_rid_strand, max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip,
            sc_buf, start, end, diag_i, soa_diag,
        )
    }
}

// Share the score arithmetic so masked and contiguous scans cannot drift.
// The const flag removes mask/diagonal work from the existing contiguous kernel.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
#[inline]
unsafe fn compute_chain_scores_batch_neon_impl<const MASKED: bool>(
    qi: i32, ri: i32, rid_strand_i: u32,
    soa_ref_pos: &[i32], soa_query_pos: &[i32], soa_query_span: &[i32],
    soa_rid_strand: &[u32],
    max_dist_x: i32, max_dist_y: i32, bw: i32, chn_pen_gap: f32, chn_pen_skip: f32,
    sc_buf: &mut [i32], start: usize, end: usize, diag_i: i32, soa_diag: &[i32],
) -> u64 { unsafe {
    let mut mask = 0u64;
    if start >= end {
        return mask;
    }
    debug_assert!(!MASKED || end - start <= 64);

    let qi_v = vdupq_n_s32(qi);
    let ri_v = vdupq_n_s32(ri);
    let rid_i_v = vdupq_n_s32(rid_strand_i as i32);
    let max_dx_v = vdupq_n_s32(max_dist_x);
    let max_dy_v = vdupq_n_s32(max_dist_y);
    let bw_v = vdupq_n_s32(bw);
    let zero_v = vdupq_n_s32(0);
    let one_v = vdupq_n_s32(1);
    let min_val_v = vdupq_n_s32(i32::MIN);
    let pen_gap_v = vdupq_n_f32(chn_pen_gap);
    let pen_skip_v = vdupq_n_f32(chn_pen_skip);
    let half_v = vdupq_n_f32(0.5f32);
    let diag_i_v = vdupq_n_s32(diag_i);
    let lane_bits = vld1q_u32([1u32, 2, 4, 8].as_ptr());

    // Process 4 predecessors at a time
    let aligned_end = start + ((end - start) / 4) * 4;
    let mut j = start;

    while j < aligned_end {
        let band_gap = if MASKED {
            // Wrapping i32 diagonals match the score kernel even across the signed
            // coordinate boundary. abs(i32::MIN) stays negative, as in the old kernel.
            let dj = vld1q_s32(soa_diag.as_ptr().add(j));
            let gap = vabsq_s32(vsubq_s32(diag_i_v, dj));
            if vminvq_u32(vcgtq_s32(gap, bw_v)) == u32::MAX {
                j += 4;
                continue;
            }
            gap
        } else {
            zero_v
        };
        // Load 4 predecessors' SoA data
        let qj_v = vld1q_s32(soa_query_pos.as_ptr().add(j));
        let rj_v = vld1q_s32(soa_ref_pos.as_ptr().add(j));
        let span_v = vld1q_s32(soa_query_span.as_ptr().add(j));
        let rid_j_v = vld1q_s32(soa_rid_strand.as_ptr().add(j) as *const i32);

        // query_diff = qi - qj
        let qdiff = vsubq_s32(qi_v, qj_v);
        // ref_diff = ri - rj
        let rdiff = vsubq_s32(ri_v, rj_v);

        // --- Build validity mask ---
        // 1. query_diff > 0
        let mut valid = vcgtq_s32(qdiff, zero_v);

        // 2. query_diff <= max_dist_x  (NOT (query_diff > max_dist_x))
        let qdiff_gt_max_dx = vcgtq_s32(qdiff, max_dx_v);
        valid = vbicq_u32(valid, qdiff_gt_max_dx);

        // 3. Same ref_id_strand
        let same_rid = vceqq_s32(rid_i_v, rid_j_v);
        valid = vandq_u32(valid, same_rid);

        // 4. ref_diff != 0
        let rdiff_zero = vceqq_s32(rdiff, zero_v);
        valid = vbicq_u32(valid, rdiff_zero);

        // 5. query_diff <= max_dist_y
        let qdiff_gt_max_dy = vcgtq_s32(qdiff, max_dy_v);
        valid = vbicq_u32(valid, qdiff_gt_max_dy);

        // 6. gap_width = abs(ref_diff - query_diff); gap_width <= bw
        let gap_w = if MASKED { band_gap } else { vabsq_s32(vsubq_s32(rdiff, qdiff)) };
        let gw_gt_bw = vcgtq_s32(gap_w, bw_v);
        valid = vbicq_u32(valid, gw_gt_bw);
        let bits = if MASKED {
            let bits = vaddvq_u32(vandq_u32(valid, lane_bits)) as u64;
            if bits == 0 {
                j += 4;
                continue;
            }
            bits
        } else { 0 };

        // --- Score computation ---
        let min_diff = vminq_s32(rdiff, qdiff);
        let mut score = vminq_s32(span_v, min_diff);

        // --- Penalty computation ---
        let gw_gt_zero = vcgtq_s32(gap_w, zero_v);
        let md_gt_span = vcgtq_s32(min_diff, span_v);
        let needs_pen = vorrq_u32(gw_gt_zero, md_gt_span);

        // lin_pen = chn_pen_gap * gap_width + chn_pen_skip * min_diff
        let gap_w_f = vcvtq_f32_s32(gap_w);
        let min_diff_f = vcvtq_f32_s32(min_diff);
        let mul1 = vmulq_f32(pen_gap_v, gap_w_f);
        let mul2 = vmulq_f32(pen_skip_v, min_diff_f);
        let lin_pen = vaddq_f32(mul1, mul2);

        // log_pen = fast_log2(gap_width + 1) if gap_width >= 1, else 0
        let gw_plus1 = vaddq_s32(gap_w, one_v);
        let gw_plus1_f = vcvtq_f32_s32(gw_plus1);
        let log2_val = simd_fast_log2_neon(gw_plus1_f);
        // Mask: zero out log_pen where gap_width == 0
        let log_pen = vreinterpretq_f32_u32(vandq_u32(gw_gt_zero, vreinterpretq_u32_f32(log2_val)));

        // total_pen = lin_pen + 0.5 * log_pen
        let half_log = vmulq_f32(half_v, log_pen);
        let total_pen = vaddq_f32(lin_pen, half_log);

        // Convert to i32 (truncation toward zero, matches Rust `as i32`)
        let pen_i32 = vcvtq_s32_f32(total_pen);

        // Apply penalty only where needs_pen is set
        let pen_masked = vandq_s32(vreinterpretq_s32_u32(needs_pen), pen_i32);
        score = vsubq_s32(score, pen_masked);

        // Mask invalid scores to i32::MIN
        score = vbslq_s32(valid, score, min_val_v);

        // Store 4 scores
        vst1q_s32(sc_buf.as_mut_ptr().add(j), score);
        if MASKED { mask |= bits << (j - start); }

        j += 4;
    }

    // Scalar remainder
    while j < end {
        let qj = soa_query_pos[j];
        let rj = soa_ref_pos[j];
        let span_j = soa_query_span[j];
        let rid_j = soa_rid_strand[j];

        let query_diff = qi.wrapping_sub(qj);
        if query_diff <= 0 || query_diff > max_dist_x {
            sc_buf[j] = i32::MIN;
            j += 1;
            continue;
        }

        let ref_diff = ri.wrapping_sub(rj);
        if rid_strand_i == rid_j && (ref_diff == 0 || query_diff > max_dist_y) {
            sc_buf[j] = i32::MIN;
            j += 1;
            continue;
        }

        let gap_width = (ref_diff - query_diff).abs();
        if rid_strand_i == rid_j && gap_width > bw {
            sc_buf[j] = i32::MIN;
            j += 1;
            continue;
        }

        let min_diff = ref_diff.min(query_diff);
        let mut sc = span_j.min(min_diff);

        if gap_width > 0 || min_diff > span_j {
            let lin_pen = chn_pen_gap * (gap_width as f32) + chn_pen_skip * (min_diff as f32);
            let log_pen = if gap_width >= 1 {
                fast_log2((gap_width + 1) as f32)
            } else {
                0.0f32
            };
            sc -= (lin_pen + 0.5f32 * log_pen) as i32;
        }

        sc_buf[j] = sc;
        if MASKED && sc != i32::MIN { mask |= 1u64 << (j - start); }
        j += 1;
    }
    mask
}}

/// NEON SIMD-optimized chaining for aarch64.
/// Phase 1: Batch-compute scores with NEON (4-wide). Phase 2: Scalar skip/visited logic.
/// Produces bit-exact output matching `chain_anchors_scalar`.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub(crate) unsafe fn chain_anchors_neon(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) {
    unsafe { chain_anchors_neon_impl::<true, true>(opt, max_dist_x, max_dist_y, a, ctx) }
}

// Const switches let differential tests disable either optimization while sharing the
// score arithmetic and backtracking. Production enables both at compile time.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn chain_anchors_neon_impl<const USE_DIAG_BUCKETS: bool, const USE_MASKED_SCAN: bool>(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) { unsafe {
    let n = a.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    let bw = opt.bandwidth;
    let mut max_dist_x = max_dist_x;
    let mut max_dist_y = max_dist_y;
    if max_dist_x < bw {
        max_dist_x = bw;
    }
    if max_dist_y < bw {
        max_dist_y = bw;
    }

    let real_max_drop = bw;
    let no_chain_skip = opt.max_chain_skip <= 0;

    // --- Extract SoA fields from AoS Minimizer array ---
    let mut soa_ref_pos = std::mem::take(&mut ctx.soa_ref_pos);
    let mut soa_query_pos = std::mem::take(&mut ctx.soa_query_pos);
    let mut soa_query_span = std::mem::take(&mut ctx.soa_query_span);
    let mut soa_rid_strand = std::mem::take(&mut ctx.soa_ref_id_strand);
    let mut soa_sc_buf = std::mem::take(&mut ctx.soa_scores_buf);
    let mut soa_diag = std::mem::take(&mut ctx.soa_diag);

    soa_ref_pos.resize(n, 0);
    soa_query_pos.resize(n, 0);
    soa_query_span.resize(n, 0);
    soa_rid_strand.resize(n, 0);
    soa_sc_buf.resize(n, 0);
    if USE_MASKED_SCAN && !no_chain_skip { soa_diag.resize(n, 0); }

    for j in 0..n {
        soa_ref_pos[j] = a[j].ref_pos();
        soa_query_pos[j] = a[j].query_pos();
        soa_query_span[j] = a[j].query_span();
        soa_rid_strand[j] = (a[j].x >> 32) as u32;
        if USE_MASKED_SCAN && !no_chain_skip {
            soa_diag[j] = soa_ref_pos[j].wrapping_sub(soa_query_pos[j]);
        }
    }

    // --- DP buffers ---
    let mut predecessors = std::mem::take(&mut ctx.predecessors); let mut bt_candidates = std::mem::take(&mut ctx.bt_candidates);
    let mut scores = std::mem::take(&mut ctx.scores);
    let mut peak_scores = std::mem::take(&mut ctx.peak_scores);
    let mut visited = std::mem::take(&mut ctx.visited);
    predecessors.resize(n, 0i32);
    scores.resize(n, 0i32);
    peak_scores.resize(n, 0i32);
    visited.clear();
    visited.resize(n, 0i32);

    let mut global_max_score = 0;
    let mut window_start: usize = 0;
    let mut best_anchor_idx: i64 = -1;

    let mut bk = std::mem::take(&mut ctx.diag_buckets);
    let use_buckets = USE_DIAG_BUCKETS && !no_chain_skip && n >= DIAG_BUCKET_MIN_N && {
        let (lo, hi) = soa_query_pos.iter().fold((i32::MAX, i32::MIN), |(l, h), &q| (l.min(q), h.max(q)));
        // The score kernel uses wrapping i32 position differences. The i64 bucket
        // diagonals are equivalent only without a signed-coordinate wrap; retain the
        // contiguous scan for larger positions and widen the density arithmetic too.
        lo >= 0
            && n >= DIAG_BUCKET_MIN_DENSITY * (hi as i64 - lo as i64 + 1) as usize
            && soa_ref_pos.iter().all(|&r| r >= 0)
    };
    if use_buckets { bk.build(&soa_ref_pos, &soa_query_pos, bw); }

    for i in 0..n {
        let mut best_predecessor: i64 = -1;
        let mut best_score = a[i].query_span();
        let mut skip_count: i32 = 0;

        while window_start < i
            && (a[i].ref_id_strand() != a[window_start].ref_id_strand()
                || (a[i].x > a[window_start].x + max_dist_x as u64))
        {
            window_start += 1;
        }

        if !no_chain_skip && (i - window_start) > opt.max_chain_iter as usize {
            window_start = i - opt.max_chain_iter as usize;
        }

        let mut end_j: i64 = if window_start > 0 {
            window_start as i64 - 1
        } else {
            -1
        };

        if no_chain_skip {
            if i > window_start {
                compute_chain_scores_batch_neon(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, window_start, i,
                );
            }
            for j in window_start..i {
                let sc = soa_sc_buf[j];
                if sc == i32::MIN { continue; }
                let total_sc = sc + scores[j];
                if total_sc > best_score { best_score = total_sc; best_predecessor = j as i64; }
            }
        } else if use_buckets && i - window_start > DIAG_BUCKET_MIN_WINDOW {
            // Merge only the in-band predecessors, in the original descending visit order.
            const CHUNK_SIZE: usize = 64;
            let nb = bk.nbr[bk.id[i] as usize];
            let mut lp: [*const u32; 3] = [std::ptr::null(); 3];
            let mut cur = [0usize; 3];
            for t in 0..3 {
                if nb[t] != u32::MAX {
                    let list = &bk.lists[nb[t] as usize];
                    lp[t] = list.as_ptr();
                    cur[t] = list.len();
                }
            }
            // SAFETY: a nonzero cursor belongs to a present bucket, and its list is not
            // modified until this anchor's scan ends. All stored indices precede i.
            let tail = |t: usize, c: usize| -> i64 {
                if c > 0 { *lp[t].add(c - 1) as i64 } else { -1 }
            };
            let (mut v0, mut v1, mut v2) = (tail(0, cur[0]), tail(1, cur[1]), tail(2, cur[2]));
            let di = bk.diag[i];
            let ws = window_start as i64;
            let mut skip_fired = false;
            let mut exhausted = false;
            loop {
                bk.g_idx.clear();
                while bk.g_idx.len() < CHUNK_SIZE {
                    let t = if v0 >= v1 && v0 >= v2 { 0 } else if v1 >= v2 { 1 } else { 2 };
                    let j = [v0, v1, v2][t];
                    if j < ws { exhausted = true; break; }
                    cur[t] -= 1;
                    match t {
                        0 => v0 = tail(0, cur[0]),
                        1 => v1 = tail(1, cur[1]),
                        _ => v2 = tail(2, cur[2]),
                    }
                    // Out-of-band pairs score i32::MIN and have no skip/visited effects.
                    if (bk.diag[j as usize] - di).abs() > bw as i64 { continue; }
                    bk.g_idx.push(j as u32);
                }
                let m = bk.g_idx.len();
                if m == 0 { break; }
                bk.g_ref.clear(); bk.g_qry.clear(); bk.g_span.clear(); bk.g_rid.clear();
                for &j in &bk.g_idx {
                    let j = j as usize;
                    bk.g_ref.push(soa_ref_pos[j]); bk.g_qry.push(soa_query_pos[j]);
                    bk.g_span.push(soa_query_span[j]); bk.g_rid.push(soa_rid_strand[j]);
                }
                bk.g_sc.resize(m, 0);
                compute_chain_scores_batch_neon(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &bk.g_ref, &bk.g_qry, &bk.g_span, &bk.g_rid,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut bk.g_sc, 0, m,
                );
                for t in 0..m {
                    let sc = bk.g_sc[t];
                    if sc == i32::MIN { continue; }
                    let j = bk.g_idx[t] as usize;
                    let total_sc = sc + scores[j];
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if visited[j] == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                    }
                    if predecessors[j] >= 0 { visited[predecessors[j] as usize] = i as i32; }
                }
                if skip_fired || exhausted { break; }
            }
        } else {
            const CHUNK_SIZE: usize = 64;
            const FIRST_CHUNK: usize = 32;
            let mut chunk_end = i;
            let mut skip_fired = false;
            let mut csize = if USE_MASKED_SCAN { FIRST_CHUNK } else { CHUNK_SIZE };
            while chunk_end > window_start {
                let chunk_start = if chunk_end >= window_start + csize { chunk_end - csize } else { window_start };
                csize = CHUNK_SIZE;

                let mut valid = if USE_MASKED_SCAN {
                    chain_scores_masked_neon(
                        soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i], soa_diag[i],
                        &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand, &soa_diag,
                        max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip,
                        &mut soa_sc_buf, chunk_start, chunk_end,
                    )
                } else {
                    compute_chain_scores_batch_neon(
                        soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                        &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                        max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip,
                        &mut soa_sc_buf, chunk_start, chunk_end,
                    );
                    (chunk_start..chunk_end).fold(0u64, |m, j| {
                        m | (((soa_sc_buf[j] != i32::MIN) as u64) << (j - chunk_start))
                    })
                };

                while valid != 0 {
                    let bit = 63 - valid.leading_zeros() as usize;
                    valid &= !(1u64 << bit);
                    let j = chunk_start + bit;
                    let sc = soa_sc_buf[j];
                    if sc == i32::MIN { continue; }
                    let total_sc = sc + scores[j];
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if visited[j] == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                    }
                    if predecessors[j] >= 0 { visited[predecessors[j] as usize] = i as i32; }
                }

                if skip_fired { break; }
                chunk_end = chunk_start;
            }
        }

        if !no_chain_skip {
            if best_anchor_idx < 0 || a[i].x - a[best_anchor_idx as usize].x > max_dist_x as u64 {
                let mut best = i32::MIN; best_anchor_idx = -1;
                for j in (window_start..i).rev() {
                    if best < scores[j] { best = scores[j]; best_anchor_idx = j as i64; }
                }
            }

            if best_anchor_idx >= 0 && best_anchor_idx < end_j {
                let sc = compute_chain_score(&a[i], &a[best_anchor_idx as usize],
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip);
                if sc != i32::MIN && best_score < sc + scores[best_anchor_idx as usize] {
                    best_score = sc + scores[best_anchor_idx as usize];
                    best_predecessor = best_anchor_idx;
                }
            }
        }

        scores[i] = best_score;
        predecessors[i] = best_predecessor as i32;
        peak_scores[i] = if best_predecessor >= 0
            && peak_scores[best_predecessor as usize] > best_score
        {
            peak_scores[best_predecessor as usize]
        } else {
            best_score
        };

        if best_anchor_idx < 0
            || (a[i].x - a[best_anchor_idx as usize].x <= max_dist_x as u64
                && scores[best_anchor_idx as usize] < scores[i])
        {
            best_anchor_idx = i as i64;
        }

        if global_max_score < best_score {
            global_max_score = best_score;
        }
        if use_buckets { bk.lists[bk.id[i] as usize].push(i as u32); }
    }
    ctx.diag_buckets = bk;

    // --- Backtrack and compact (identical to AVX2/scalar) ---
    let (u, n_u, n_v) = chain_backtrack(
        n,
        &scores,
        &predecessors,
        &mut peak_scores,
        &mut visited,
        &mut bt_candidates,
        opt.min_cnt,
        opt.min_chain_score,
        real_max_drop,
    );

    if n_u == 0 {
        ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates);
        ctx.scores = scores;
        ctx.peak_scores = peak_scores;
        ctx.visited = visited;
        ctx.soa_ref_pos = soa_ref_pos;
        ctx.soa_query_pos = soa_query_pos;
        ctx.soa_query_span = soa_query_span;
        ctx.soa_ref_id_strand = soa_rid_strand;
        ctx.soa_scores_buf = soa_sc_buf;
        ctx.soa_diag = soa_diag;
        return (Vec::new(), Vec::new());
    }

    let mut b: Vec<Minimizer> = Vec::with_capacity(n_v);
    let mut k = 0usize;
    for &u_val in &u[..n_u] {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        let k0 = k;
        for j in 0..ni {
            let idx = peak_scores[k0 + (ni - j - 1)] as usize;
            b.push(a[idx]);
            k += 1;
        }
    }

    let mut w: Vec<Minimizer> = Vec::with_capacity(n_u);
    let mut k_pos = 0usize;
    for (i, &u_val) in u[..n_u].iter().enumerate() {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        w.push(Minimizer {
            x: b[k_pos].x,
            y: ((k_pos as u64) << 32) | (i as u64),
        });
        k_pos += ni;
    }
    radix_sort_128x(&mut w);

    let mut u2: Vec<u64> = Vec::with_capacity(n_u);
    let mut b2: Vec<Minimizer> = Vec::with_capacity(n_v);
    for &w_val in &w[..n_u] {
        let j = (w_val.y & 0xFFFFFFFF) as usize;
        let offset = (w_val.y >> 32) as usize;
        let ni = (u[j] & 0xFFFFFFFF) as usize;
        u2.push(u[j]);
        for idx in 0..ni {
            b2.push(b[offset + idx]);
        }
    }

    ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates);
    ctx.scores = scores;
    ctx.peak_scores = peak_scores;
    ctx.visited = visited;
    ctx.soa_ref_pos = soa_ref_pos;
    ctx.soa_query_pos = soa_query_pos;
    ctx.soa_query_span = soa_query_span;
    ctx.soa_ref_id_strand = soa_rid_strand;
    ctx.soa_scores_buf = soa_sc_buf;
    ctx.soa_diag = soa_diag;

    (u2, b2)
}}

// =============================================================================
// WASM SIMD128 implementation — 4-wide, mirrors NEON logic exactly
// =============================================================================

/// Vectorized fast_log2 for 4 floats (WASM SIMD128).
#[cfg(target_arch = "wasm32")]
#[target_feature(enable = "simd128")]
#[inline]
unsafe fn simd_fast_log2_wasm(x: v128) -> v128 {
    let z = x; // reinterpret as i32x4 via bit ops

    // exponent = ((z >> 23) & 255) - 128
    let exp_raw = v128_and(i32x4_shr(z, 23), i32x4_splat(255));
    let exp = i32x4_sub(exp_raw, i32x4_splat(128));
    let log_2 = f32x4_convert_i32x4(exp);

    // mantissa: clear exponent bits, set exponent to 127
    let mant_bits = v128_or(
        v128_and(z, i32x4_splat(!(255i32 << 23))),
        i32x4_splat(127i32 << 23),
    );
    let z_f = mant_bits; // reinterpret as f32x4

    // Horner's: ((-0.34484843 * z_f) + 2.0246658) * z_f + (-0.6748776)
    let c0 = f32x4_splat(-0.34484843f32);
    let c1 = f32x4_splat(2.024_665_8f32);
    let c2 = f32x4_splat(-0.674_877_6f32);

    let t1 = f32x4_mul(c0, z_f);
    let t2 = f32x4_add(t1, c1);
    let t3 = f32x4_mul(t2, z_f);
    let poly = f32x4_add(t3, c2);

    f32x4_add(log_2, poly)
}

/// Batch-compute chain scores for 4 predecessors at a time (WASM SIMD128).
#[cfg(target_arch = "wasm32")]
#[target_feature(enable = "simd128")]
#[inline]
unsafe fn compute_chain_scores_batch_wasm(
    qi: i32,
    ri: i32,
    rid_strand_i: u32,
    soa_ref_pos: &[i32],
    soa_query_pos: &[i32],
    soa_query_span: &[i32],
    soa_rid_strand: &[u32],
    max_dist_x: i32,
    max_dist_y: i32,
    bw: i32,
    chn_pen_gap: f32,
    chn_pen_skip: f32,
    sc_buf: &mut [i32],
    start: usize,
    end: usize,
) { unsafe {
    if start >= end { return; }

    let qi_v = i32x4_splat(qi);
    let ri_v = i32x4_splat(ri);
    let rid_i_v = i32x4_splat(rid_strand_i as i32);
    let max_dx_v = i32x4_splat(max_dist_x);
    let max_dy_v = i32x4_splat(max_dist_y);
    let bw_v = i32x4_splat(bw);
    let zero_v = i32x4_splat(0);
    let one_v = i32x4_splat(1);
    let min_val_v = i32x4_splat(i32::MIN);
    let pen_gap_v = f32x4_splat(chn_pen_gap);
    let pen_skip_v = f32x4_splat(chn_pen_skip);
    let half_v = f32x4_splat(0.5f32);

    let aligned_end = start + ((end - start) / 4) * 4;
    let mut j = start;

    while j < aligned_end {
        let qj_v = v128_load(soa_query_pos.as_ptr().add(j) as *const v128);
        let rj_v = v128_load(soa_ref_pos.as_ptr().add(j) as *const v128);
        let span_v = v128_load(soa_query_span.as_ptr().add(j) as *const v128);
        let rid_j_v = v128_load(soa_rid_strand.as_ptr().add(j) as *const v128);

        let qdiff = i32x4_sub(qi_v, qj_v);
        let rdiff = i32x4_sub(ri_v, rj_v);

        // Build validity mask
        let mut valid = i32x4_gt(qdiff, zero_v);
        let qdiff_gt_max_dx = i32x4_gt(qdiff, max_dx_v);
        valid = v128_andnot(valid, qdiff_gt_max_dx);  // valid & ~gt_max_dx
        let same_rid = i32x4_eq(rid_i_v, rid_j_v);
        valid = v128_and(valid, same_rid);
        let rdiff_zero = i32x4_eq(rdiff, zero_v);
        valid = v128_andnot(valid, rdiff_zero);
        let qdiff_gt_max_dy = i32x4_gt(qdiff, max_dy_v);
        valid = v128_andnot(valid, qdiff_gt_max_dy);

        let diff = i32x4_sub(rdiff, qdiff);
        let gap_w = i32x4_abs(diff);
        let gw_gt_bw = i32x4_gt(gap_w, bw_v);
        valid = v128_andnot(valid, gw_gt_bw);

        // Score computation
        let min_diff = i32x4_min(rdiff, qdiff);
        let mut score = i32x4_min(span_v, min_diff);

        // Penalty
        let gw_gt_zero = i32x4_gt(gap_w, zero_v);
        let md_gt_span = i32x4_gt(min_diff, span_v);
        let needs_pen = v128_or(gw_gt_zero, md_gt_span);

        let gap_w_f = f32x4_convert_i32x4(gap_w);
        let min_diff_f = f32x4_convert_i32x4(min_diff);
        let mul1 = f32x4_mul(pen_gap_v, gap_w_f);
        let mul2 = f32x4_mul(pen_skip_v, min_diff_f);
        let lin_pen = f32x4_add(mul1, mul2);

        let gw_plus1 = i32x4_add(gap_w, one_v);
        let gw_plus1_f = f32x4_convert_i32x4(gw_plus1);
        let log2_val = simd_fast_log2_wasm(gw_plus1_f);
        let log_pen = v128_and(gw_gt_zero, log2_val);

        let half_log = f32x4_mul(half_v, log_pen);
        let total_pen = f32x4_add(lin_pen, half_log);

        let pen_i32 = i32x4_trunc_sat_f32x4(total_pen);

        let pen_masked = v128_and(needs_pen, pen_i32);
        score = i32x4_sub(score, pen_masked);

        // Mask invalid to i32::MIN
        score = v128_bitselect(score, min_val_v, valid);

        v128_store(sc_buf.as_mut_ptr().add(j) as *mut v128, score);
        j += 4;
    }

    // Scalar remainder
    while j < end {
        let qj = soa_query_pos[j];
        let rj = soa_ref_pos[j];
        let span_j = soa_query_span[j];
        let rid_j = soa_rid_strand[j];

        let query_diff = qi.wrapping_sub(qj);
        if query_diff <= 0 || query_diff > max_dist_x {
            sc_buf[j] = i32::MIN; j += 1; continue;
        }
        let ref_diff = ri.wrapping_sub(rj);
        if rid_strand_i == rid_j && (ref_diff == 0 || query_diff > max_dist_y) {
            sc_buf[j] = i32::MIN; j += 1; continue;
        }
        let gap_width = (ref_diff - query_diff).abs();
        if rid_strand_i == rid_j && gap_width > bw {
            sc_buf[j] = i32::MIN; j += 1; continue;
        }
        let min_diff = ref_diff.min(query_diff);
        let mut sc = span_j.min(min_diff);
        if gap_width > 0 || min_diff > span_j {
            let lin_pen = chn_pen_gap * (gap_width as f32) + chn_pen_skip * (min_diff as f32);
            let log_pen = if gap_width >= 1 { fast_log2((gap_width + 1) as f32) } else { 0.0f32 };
            sc -= (lin_pen + 0.5f32 * log_pen) as i32;
        }
        sc_buf[j] = sc;
        j += 1;
    }
}}

/// WASM SIMD128 chaining — 4-wide, same structure as NEON/AVX2.
#[cfg(target_arch = "wasm32")]
#[target_feature(enable = "simd128")]
pub(crate) unsafe fn chain_anchors_wasm(
    opt: &ChainingParams,
    max_dist_x: i32,
    max_dist_y: i32,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) { unsafe {
    let n = a.len();
    if n == 0 { return (Vec::new(), Vec::new()); }

    let bw = opt.bandwidth;
    let mut max_dist_x = max_dist_x;
    let mut max_dist_y = max_dist_y;
    if max_dist_x < bw { max_dist_x = bw; }
    if max_dist_y < bw { max_dist_y = bw; }

    let real_max_drop = bw;
    let no_chain_skip = opt.max_chain_skip <= 0;

    let mut soa_ref_pos = std::mem::take(&mut ctx.soa_ref_pos);
    let mut soa_query_pos = std::mem::take(&mut ctx.soa_query_pos);
    let mut soa_query_span = std::mem::take(&mut ctx.soa_query_span);
    let mut soa_rid_strand = std::mem::take(&mut ctx.soa_ref_id_strand);
    let mut soa_sc_buf = std::mem::take(&mut ctx.soa_scores_buf);

    soa_ref_pos.resize(n, 0);
    soa_query_pos.resize(n, 0);
    soa_query_span.resize(n, 0);
    soa_rid_strand.resize(n, 0);
    soa_sc_buf.resize(n, 0);

    for j in 0..n {
        soa_ref_pos[j] = a[j].ref_pos();
        soa_query_pos[j] = a[j].query_pos();
        soa_query_span[j] = a[j].query_span();
        soa_rid_strand[j] = (a[j].x >> 32) as u32;
    }

    let mut predecessors = std::mem::take(&mut ctx.predecessors); let mut bt_candidates = std::mem::take(&mut ctx.bt_candidates);
    let mut scores = std::mem::take(&mut ctx.scores);
    let mut peak_scores = std::mem::take(&mut ctx.peak_scores);
    let mut visited = std::mem::take(&mut ctx.visited);
    predecessors.resize(n, 0i32);
    scores.resize(n, 0i32);
    peak_scores.resize(n, 0i32);
    visited.clear();
    visited.resize(n, 0i32);

    let mut global_max_score = 0;
    let mut window_start: usize = 0;
    let mut best_anchor_idx: i64 = -1;

    for i in 0..n {
        let mut best_predecessor: i64 = -1;
        let mut best_score = a[i].query_span();
        let mut skip_count: i32 = 0;

        while window_start < i
            && (a[i].ref_id_strand() != a[window_start].ref_id_strand()
                || (a[i].x > a[window_start].x + max_dist_x as u64))
        { window_start += 1; }

        if !no_chain_skip && (i - window_start) > opt.max_chain_iter as usize {
            window_start = i - opt.max_chain_iter as usize;
        }

        let mut end_j: i64 = if window_start > 0 { window_start as i64 - 1 } else { -1 };

        if no_chain_skip {
            if i > window_start {
                compute_chain_scores_batch_wasm(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, window_start, i,
                );
            }
            for j in window_start..i {
                let sc = soa_sc_buf[j];
                if sc == i32::MIN { continue; }
                let total_sc = sc + scores[j];
                if total_sc > best_score { best_score = total_sc; best_predecessor = j as i64; }
            }
        } else {
            const CHUNK_SIZE: usize = 64;
            let mut chunk_end = i;
            let mut skip_fired = false;
            while chunk_end > window_start {
                let chunk_start = if chunk_end >= window_start + CHUNK_SIZE { chunk_end - CHUNK_SIZE } else { window_start };

                compute_chain_scores_batch_wasm(
                    soa_query_pos[i], soa_ref_pos[i], soa_rid_strand[i],
                    &soa_ref_pos, &soa_query_pos, &soa_query_span, &soa_rid_strand,
                    max_dist_x, max_dist_y, bw,
                    opt.chn_pen_gap, opt.chn_pen_skip,
                    &mut soa_sc_buf, chunk_start, chunk_end,
                );

                for j in (chunk_start..chunk_end).rev() {
                    let sc = soa_sc_buf[j];
                    if sc == i32::MIN { continue; }
                    let total_sc = sc + scores[j];
                    if total_sc > best_score {
                        best_score = total_sc; best_predecessor = j as i64;
                        if skip_count > 0 { skip_count -= 1; }
                    } else if visited[j] == i as i32 {
                        skip_count += 1;
                        if skip_count > opt.max_chain_skip { end_j = j as i64; skip_fired = true; break; }
                    }
                    if predecessors[j] >= 0 { visited[predecessors[j] as usize] = i as i32; }
                }

                if skip_fired { break; }
                chunk_end = chunk_start;
            }

            if best_anchor_idx < 0 || a[i].x - a[best_anchor_idx as usize].x > max_dist_x as u64 {
                let mut best = i32::MIN; best_anchor_idx = -1;
                for j in (window_start..i).rev() {
                    if best < scores[j] { best = scores[j]; best_anchor_idx = j as i64; }
                }
            }

            if best_anchor_idx >= 0 && best_anchor_idx < end_j {
                let sc = compute_chain_score(&a[i], &a[best_anchor_idx as usize],
                    max_dist_x, max_dist_y, bw, opt.chn_pen_gap, opt.chn_pen_skip);
                if sc != i32::MIN && best_score < sc + scores[best_anchor_idx as usize] {
                    best_score = sc + scores[best_anchor_idx as usize];
                    best_predecessor = best_anchor_idx;
                }
            }
        }

        scores[i] = best_score;
        predecessors[i] = best_predecessor as i32;
        peak_scores[i] = if best_predecessor >= 0
            && peak_scores[best_predecessor as usize] > best_score
        { peak_scores[best_predecessor as usize] } else { best_score };

        if best_anchor_idx < 0
            || (a[i].x - a[best_anchor_idx as usize].x <= max_dist_x as u64
                && scores[best_anchor_idx as usize] < scores[i])
        { best_anchor_idx = i as i64; }

        if global_max_score < best_score { global_max_score = best_score; }
    }

    let (u, n_u, n_v) = chain_backtrack(
        n, &scores, &predecessors, &mut peak_scores, &mut visited,
        &mut bt_candidates,
        opt.min_cnt, opt.min_chain_score, real_max_drop,
    );

    if n_u == 0 {
        ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates); ctx.scores = scores;
        ctx.peak_scores = peak_scores; ctx.visited = visited;
        ctx.soa_ref_pos = soa_ref_pos; ctx.soa_query_pos = soa_query_pos;
        ctx.soa_query_span = soa_query_span; ctx.soa_ref_id_strand = soa_rid_strand;
        ctx.soa_scores_buf = soa_sc_buf;
        return (Vec::new(), Vec::new());
    }

    let mut b: Vec<Minimizer> = Vec::with_capacity(n_v);
    let mut k = 0usize;
    for &u_val in &u[..n_u] {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        let k0 = k;
        for j in 0..ni {
            let idx = peak_scores[k0 + (ni - j - 1)] as usize;
            b.push(a[idx]);
            k += 1;
        }
    }

    let mut w: Vec<Minimizer> = Vec::with_capacity(n_u);
    let mut k_pos = 0usize;
    for (i, &u_val) in u[..n_u].iter().enumerate() {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        w.push(Minimizer { x: b[k_pos].x, y: ((k_pos as u64) << 32) | (i as u64) });
        k_pos += ni;
    }
    radix_sort_128x(&mut w);

    let mut u2: Vec<u64> = Vec::with_capacity(n_u);
    let mut b2: Vec<Minimizer> = Vec::with_capacity(n_v);
    for &w_val in &w[..n_u] {
        let j = (w_val.y & 0xFFFFFFFF) as usize;
        let offset = (w_val.y >> 32) as usize;
        let ni = (u[j] & 0xFFFFFFFF) as usize;
        u2.push(u[j]);
        for idx in 0..ni { b2.push(b[offset + idx]); }
    }

    ctx.predecessors = predecessors; ctx.bt_candidates = std::mem::take(&mut bt_candidates); ctx.scores = scores;
    ctx.peak_scores = peak_scores; ctx.visited = visited;
    ctx.soa_ref_pos = soa_ref_pos; ctx.soa_query_pos = soa_query_pos;
    ctx.soa_query_span = soa_query_span; ctx.soa_ref_id_strand = soa_rid_strand;
    ctx.soa_scores_buf = soa_sc_buf;

    (u2, b2)
}}

#[cfg(all(test, any(target_arch = "aarch64", target_arch = "x86_64")))]
#[path = "chain_simd_masked_tests.rs"]
mod masked_tests;

#[cfg(all(test, target_arch = "x86_64"))]
#[path = "chain_simd_avx512_tests.rs"]
mod avx512_tests;

#[cfg(all(test, target_arch = "aarch64"))]
#[path = "chain_simd_neon_tests.rs"]
mod neon_tests;

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;
    use crate::align::chain::{chain_anchors_scalar, fast_log2};
    use crate::align::map::ChainingParams;

    #[test]
    fn test_simd_log2_vs_scalar() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let test_values: Vec<f32> = vec![
            1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 100.0, 1000.0, 0.5, 0.001, 1.5, 7.7, 50.3, 999.9,
            2.5, 123.456,
        ];

        unsafe {
            // Process 8 at a time
            for chunk in test_values.chunks(8) {
                let mut input = [1.0f32; 8];
                for (i, &v) in chunk.iter().enumerate() {
                    input[i] = v;
                }
                let v = _mm256_loadu_ps(input.as_ptr());
                let result = simd_fast_log2_avx2(v);
                let mut output = [0.0f32; 8];
                _mm256_storeu_ps(output.as_mut_ptr(), result);

                for (i, &v) in chunk.iter().enumerate() {
                    let scalar = fast_log2(v);
                    assert_eq!(
                        output[i].to_bits(),
                        scalar.to_bits(),
                        "log2 mismatch for {}: simd={} scalar={}",
                        v,
                        output[i],
                        scalar
                    );
                }
            }
        }
    }

    #[test]
    fn test_chain_simd_vs_scalar_synthetic() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        let k = 15;
        let span_mask = (k as u64) << 32;
        let mut anchors_simd: Vec<Minimizer> = Vec::new();

        // Chain 1: Linear perfect match
        let chain1 = vec![
            Minimizer { x: 100, y: span_mask | 100 },
            Minimizer { x: 120, y: span_mask | 120 },
            Minimizer { x: 150, y: span_mask | 150 },
        ];
        // Chain 2: Isolated
        let chain2 = vec![Minimizer { x: 500, y: span_mask | 500 }];
        // Chain 3: Indel
        let chain3 = vec![
            Minimizer { x: 1000, y: span_mask | 1000 },
            Minimizer { x: 1030, y: span_mask | 1020 },
        ];

        // Add enough anchors to exceed the 32-anchor threshold
        for i in 0..30 {
            anchors_simd.push(Minimizer {
                x: 2000 + i * 20,
                y: span_mask | (2000 + i * 20),
            });
        }

        for a in &chain1 {
            anchors_simd.push(*a);
        }
        for a in &chain2 {
            anchors_simd.push(*a);
        }
        for a in &chain3 {
            anchors_simd.push(*a);
        }
        let mut anchors_scalar: Vec<Minimizer> = anchors_simd.clone();

        // Sort by x (ref position) as required
        anchors_simd.sort_by_key(|a| a.x);
        anchors_scalar.sort_by_key(|a| a.x);

        let opt = ChainingParams {
            min_cnt: 1,
            min_chain_score: 10,
            max_gap: 500,
            max_gap_ref: -1,
            max_dist_x: 500,
            max_dist_y: 500,
            bandwidth: 500,
            bandwidth_long: 500,
            max_chain_skip: 25,
            max_chain_iter: 50,
            chn_pen_gap: 0.5,
            chn_pen_skip: 0.5,
            chain_gap_scale: 0.8,
            rmq_rescue_size: 1000,
            rmq_rescue_ratio: 0.1,
            rmq_inner_dist: 1000,
            rmq_size_cap: 100000,
            seg_rmq: false, rescue_join: false, rescue_max_density: 0,
        };

        let mut ctx_simd = ChainingBuffers::new();
        let mut ctx_scalar = ChainingBuffers::new();

        let (u_simd, chains_simd) = unsafe {
            chain_anchors_avx2(&opt, 500, 500, &mut anchors_simd, &mut ctx_simd)
        };
        let (u_scalar, chains_scalar) =
            chain_anchors_scalar(&opt, 500, 500, &mut anchors_scalar, &mut ctx_scalar);

        assert_eq!(u_simd.len(), u_scalar.len(), "chain count mismatch");
        for (i, (us, uc)) in u_simd.iter().zip(u_scalar.iter()).enumerate() {
            assert_eq!(us, uc, "chain {} u mismatch: simd={} scalar={}", i, us, uc);
        }
        assert_eq!(
            chains_simd.len(),
            chains_scalar.len(),
            "anchor count mismatch"
        );
        for (i, (as_, ac)) in chains_simd.iter().zip(chains_scalar.iter()).enumerate() {
            assert_eq!(
                as_.x, ac.x,
                "anchor {} x mismatch: simd={} scalar={}",
                i, as_.x, ac.x
            );
            assert_eq!(
                as_.y, ac.y,
                "anchor {} y mismatch: simd={} scalar={}",
                i, as_.y, ac.y
            );
        }
    }

    #[test]
    fn test_chain_simd_larger_window() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        // Create a long collinear chain with 200 anchors to stress SIMD
        let k = 15;
        let span_mask = (k as u64) << 32;
        let mut anchors_simd: Vec<Minimizer> = Vec::new();

        for i in 0..200 {
            anchors_simd.push(Minimizer {
                x: (i * 10) as u64,
                y: span_mask | ((i * 10) as u64),
            });
        }
        let mut anchors_scalar: Vec<Minimizer> = anchors_simd.clone();

        let opt = ChainingParams {
            min_cnt: 1,
            min_chain_score: 10,
            max_gap: 5000,
            max_gap_ref: -1,
            max_dist_x: 5000,
            max_dist_y: 5000,
            bandwidth: 5000,
            bandwidth_long: 5000,
            max_chain_skip: 25,
            max_chain_iter: 5000,
            chn_pen_gap: 0.5,
            chn_pen_skip: 0.5,
            chain_gap_scale: 0.8,
            rmq_rescue_size: 1000,
            rmq_rescue_ratio: 0.1,
            rmq_inner_dist: 1000,
            rmq_size_cap: 100000,
            seg_rmq: false, rescue_join: false, rescue_max_density: 0,
        };

        let mut ctx_simd = ChainingBuffers::new();
        let mut ctx_scalar = ChainingBuffers::new();

        let (u_simd, chains_simd) = unsafe {
            chain_anchors_avx2(&opt, 5000, 5000, &mut anchors_simd, &mut ctx_simd)
        };
        let (u_scalar, chains_scalar) =
            chain_anchors_scalar(&opt, 5000, 5000, &mut anchors_scalar, &mut ctx_scalar);

        assert_eq!(u_simd.len(), u_scalar.len(), "chain count mismatch");
        for (i, (us, uc)) in u_simd.iter().zip(u_scalar.iter()).enumerate() {
            assert_eq!(us, uc, "chain {} u mismatch", i);
        }
        assert_eq!(chains_simd.len(), chains_scalar.len(), "anchor count mismatch");
        for (i, (as_, ac)) in chains_simd.iter().zip(chains_scalar.iter()).enumerate() {
            assert_eq!(as_.x, ac.x, "anchor {} x mismatch", i);
            assert_eq!(as_.y, ac.y, "anchor {} y mismatch", i);
        }
    }

    /// Dense tandem-repeat-like anchors (> DIAG_BUCKET_MIN_N, windows > DIAG_BUCKET_MIN_WINDOW):
    /// the diagonal-bucket scan must give exactly the scalar chainer's output.
    #[test]
    fn test_chain_simd_diag_buckets_dense() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }
        let k = 19u64;
        let span_mask = k << 32;
        let mut seed: u64 = 0x2545_f491_4f6c_dd1d;
        let mut rnd = |m: u64| { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; seed % m };
        let mut anchors: Vec<Minimizer> = Vec::new();
        // Read of 20 kb over a tandem array (unit 37 bp, 400 copies) plus a 900 bp-shifted copy
        // and random hits; each query minimizer hits many shifted copies.
        let base = 1_000_000u64;
        let mut q = 0u64;
        while q < 20_000 {
            for _ in 0..30u64 {
                let shift = rnd(400) * 37;
                let x = base + q + shift;
                anchors.push(Minimizer { x, y: span_mask | q });
            }
            if rnd(4) == 0 { anchors.push(Minimizer { x: base + q + 900, y: span_mask | q }); }
            if rnd(3) == 0 { anchors.push(Minimizer { x: base + rnd(40_000), y: span_mask | q }); }
            q += 7 + rnd(6);
        }
        anchors.sort_by_key(|a| (a.x, a.y));
        anchors.dedup();
        assert!(anchors.len() > 4 * DIAG_BUCKET_MIN_N);
        for (bw, skip) in [(500, 25), (100, 10), (5000, 25)] {
            let opt = ChainingParams {
                min_cnt: 3, min_chain_score: 40, max_gap: 5000, max_gap_ref: -1,
                max_dist_x: 5000, max_dist_y: 5000, bandwidth: bw, bandwidth_long: 20000,
                max_chain_skip: skip, max_chain_iter: 5000, chn_pen_gap: 0.152, chn_pen_skip: 0.0,
                chain_gap_scale: 0.8, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1, rmq_inner_dist: 1000,
                rmq_size_cap: 100000, seg_rmq: false, rescue_join: false, rescue_max_density: 0,
            };
            let (mut a1, mut a2) = (anchors.clone(), anchors.clone());
            let (mut c1, mut c2) = (ChainingBuffers::new(), ChainingBuffers::new());
            let (u1, ch1) = unsafe { chain_anchors_avx2(&opt, 5000, 5000, &mut a1, &mut c1) };
            let (u2, ch2) = chain_anchors_scalar(&opt, 5000, 5000, &mut a2, &mut c2);
            assert_eq!(u1, u2, "bw {bw}: chains differ");
            assert!(ch1.iter().zip(&ch2).all(|(x, y)| x.x == y.x && x.y == y.y) && ch1.len() == ch2.len(), "bw {bw}: anchors differ");
        }
    }

    #[test]
    fn test_chain_no_skip_fast_path() {
        // Tests the --max-chain-skip=0 fast path (no skip heuristic).
        // This produces DIFFERENT output from the skip path, so we just verify
        // it runs without panics and produces valid chains.
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        let k = 15;
        let span_mask = (k as u64) << 32;

        // Create anchors with gaps that would trigger skip heuristic
        let mut anchors: Vec<Minimizer> = Vec::new();
        // Main chain: evenly spaced
        for i in 0..50 {
            anchors.push(Minimizer {
                x: (i * 20) as u64,
                y: span_mask | ((i * 20) as u64),
            });
        }
        // Noise anchors interleaved (different ref_id to be filtered)
        for i in 0..20 {
            anchors.push(Minimizer {
                x: (1u64 << 33) | ((i * 30) as u64), // different ref_id
                y: span_mask | ((i * 30 + 5) as u64),
            });
        }
        anchors.sort_by_key(|a| a.x);

        let opt_skip = ChainingParams {
            min_cnt: 1,
            min_chain_score: 10,
            max_gap: 5000,
            max_gap_ref: -1,
            max_dist_x: 5000,
            max_dist_y: 5000,
            bandwidth: 5000,
            bandwidth_long: 5000,
            max_chain_skip: 25,
            max_chain_iter: 5000,
            chn_pen_gap: 0.5,
            chn_pen_skip: 0.5,
            chain_gap_scale: 0.8,
            rmq_rescue_size: 1000,
            rmq_rescue_ratio: 0.1,
            rmq_inner_dist: 1000,
            rmq_size_cap: 100000,
            seg_rmq: false, rescue_join: false, rescue_max_density: 0,
        };

        let mut opt_noskip = opt_skip.clone();
        opt_noskip.max_chain_skip = 0;

        let mut anchors_skip = anchors.clone();
        let mut anchors_noskip = anchors;

        let mut ctx_skip = ChainingBuffers::new();
        let mut ctx_noskip = ChainingBuffers::new();

        let (u_skip, chains_skip) = unsafe {
            chain_anchors_avx2(&opt_skip, 5000, 5000, &mut anchors_skip, &mut ctx_skip)
        };
        let (u_noskip, chains_noskip) = unsafe {
            chain_anchors_avx2(&opt_noskip, 5000, 5000, &mut anchors_noskip, &mut ctx_noskip)
        };

        // Both should produce valid output (non-empty for this input)
        assert!(!u_skip.is_empty(), "skip path should produce chains");
        assert!(!u_noskip.is_empty(), "no-skip path should produce chains");

        // Verify chain scores are positive
        for (i, u) in u_skip.iter().enumerate() {
            let score = (*u >> 32) as i32;
            assert!(score > 0, "skip chain {} has non-positive score {}", i, score);
        }
        for (i, u) in u_noskip.iter().enumerate() {
            let score = (*u >> 32) as i32;
            assert!(score > 0, "no-skip chain {} has non-positive score {}", i, score);
        }

        // Verify anchor counts are consistent
        let skip_anchor_count: usize = u_skip.iter().map(|u| (u & 0xFFFFFFFF) as usize).sum();
        assert_eq!(skip_anchor_count, chains_skip.len(),
            "skip: u anchor count doesn't match chains length");
        let noskip_anchor_count: usize = u_noskip.iter().map(|u| (u & 0xFFFFFFFF) as usize).sum();
        assert_eq!(noskip_anchor_count, chains_noskip.len(),
            "no-skip: u anchor count doesn't match chains length");
    }
}

#[cfg(test)]
mod diag_bucket_tests {
    use super::*;

    #[test]
    fn integer_hasher_buckets_match_std_hashmap() {
        let mut seed = 0x3891_abd7_44ae_0f11u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        // One DiagBuckets for all cases: rebuilding must also reset reused buffers.
        let mut bk = DiagBuckets::default();
        for case in 0..700 {
            let n = [0, 1, 2, 31, 1023, 1024, 2048, 6000][case % 8];
            let q: Vec<_> = (0..n).map(|_| rnd() as i32).collect();
            let r: Vec<_> = (0..n).map(|_| if case % 3 == 0 { 17 } else { rnd() as i32 }).collect();
            let bw = [-1, 0, 1, 2, 15, 512, 20000, i32::MAX][case % 8];
            bk.build(&r, &q, bw);

            // Reference: bucket ids in first-seen order, from a std HashMap.
            let shift = if bw > 1 { 32 - ((bw - 1) as u32).leading_zeros() } else { 0 };
            let diag: Vec<i64> = r.iter().zip(&q).map(|(&r, &q)| r as i64 - q as i64).collect();
            let mut map = std::collections::HashMap::new();
            let mut keys = Vec::new();
            let id: Vec<u32> = diag.iter().map(|&d| {
                *map.entry(d >> shift).or_insert_with(|| {
                    keys.push(d >> shift);
                    keys.len() as u32 - 1
                })
            }).collect();
            let get = |k: i64| map.get(&k).copied().unwrap_or(u32::MAX);
            let nbr: Vec<_> = keys.iter().map(|&k| [get(k - 1), get(k), get(k + 1)]).collect();

            assert_eq!(bk.diag, diag, "diag case={case}");
            assert_eq!(bk.id, id, "id case={case}");
            assert_eq!(bk.keys, keys, "keys case={case}");
            assert_eq!(bk.nbr, nbr, "nbr case={case}");
            assert!(bk.lists[..keys.len()].iter().all(Vec::is_empty), "lists case={case}");
            for l in &mut bk.lists {
                l.push(17);
            }
        }
    }
}
