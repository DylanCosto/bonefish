//! AVX2 dual-affine banded DP with the semantics of `extend_dual_affine_scalar` (bit-identical
//! results), for the presets that do not need minimap2-identical output.
//!
//! The default AVX2 kernel reproduces minimap2's SSE ksw2 behaviour exactly: bands rounded out
//! to 16-cell blocks, cells outside the band computed from stale scores, a separate score pass
//! per anti-diagonal. That costs work on every anti-diagonal, and with narrow bands (a few
//! vectors per anti-diagonal) the fixed work dominates. This kernel follows the scalar
//! reference instead: each anti-diagonal covers exactly [st0, en0], cells right of the band keep
//! their initial (gap-open) values, and the scores are computed inside the main loop.
//!
//! Selected with `set_clean_dual(true)`; `extend_dual_affine` then routes AVX2 calls here.

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

use std::sync::atomic::{AtomicBool, Ordering};

use super::common::*;

static CLEAN_DUAL: AtomicBool = AtomicBool::new(false);

/// Route dual-affine DP to the scalar-semantics AVX2 kernel (process-wide).
pub fn set_clean_dual(on: bool) {
    CLEAN_DUAL.store(on, Ordering::Relaxed);
}

#[inline]
pub fn clean_dual() -> bool {
    CLEAN_DUAL.load(Ordering::Relaxed) || *CLEAN_DUAL_ENV
}

/// BONEFISH_EXP_CLEAN_DP=1 selects the kernel too (experiments).
static CLEAN_DUAL_ENV: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("BONEFISH_EXP_CLEAN_DP").is_ok_and(|v| v != "0"));

/// AVX2 kernel with `extend_dual_affine_scalar` semantics. Same arguments as the other dual
/// kernels; generic (non-uniform mismatch) scoring is not handled here (callers route it to
/// the scalar kernel).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn extend_dual_affine_avx2_clean(
    qseq: &[u8], tseq: &[u8], alphabet_size: i8, score_matrix: &[i8],
    gap_open: i8, gap_extend: i8, gap_open2: i8, gap_extend2: i8,
    bandwidth: i32, z_drop: i32, end_bonus: i32, flags: i32, adaptive: i32, band_shift: i32,
    result: &mut DpResult,
) { unsafe {
    // One specialised copy per flag combination: with the flags as constants the inner loop has
    // no per-vector branches and its masks stay in vector registers.
    let m1 = alphabet_size.wrapping_sub(1) as u8;
    // Branch-free max reductions (vectorised), not an early-exit scan.
    let has_n = qseq.iter().fold(0u8, |m, &c| m.max(c)).max(tseq.iter().fold(0u8, |m, &c| m.max(c))) >= m1;
    let approx_max = (flags & APPROX_MAX) != 0;
    let with_cigar = (flags & SCORE_ONLY) == 0;
    let right_align = (flags & RIGHT_ALIGN) != 0;
    if approx_max && adaptive <= 0 && (0..=31).contains(&bandwidth) && !*REG1_OFF {
        macro_rules! go1 {
            ($c:literal, $r:literal, $n:literal) => {
                reg1_impl::<$c, $r, $n>(qseq, tseq, alphabet_size, score_matrix, gap_open, gap_extend, gap_open2, gap_extend2,
                                        bandwidth, z_drop, end_bonus, flags, band_shift, result)
            };
        }
        match (with_cigar, right_align, has_n) {
            (false, false, false) => go1!(false, false, false),
            (false, false, true) => go1!(false, false, true),
            (false, true, false) => go1!(false, true, false),
            (false, true, true) => go1!(false, true, true),
            (true, false, false) => go1!(true, false, false),
            (true, false, true) => go1!(true, false, true),
            (true, true, false) => go1!(true, true, false),
            (true, true, true) => go1!(true, true, true),
        }
        return;
    }
    macro_rules! go {
        ($a:literal, $c:literal, $r:literal, $n:literal) => {
            clean_impl::<$a, $c, $r, $n>(qseq, tseq, alphabet_size, score_matrix, gap_open, gap_extend, gap_open2, gap_extend2,
                                         bandwidth, z_drop, end_bonus, flags, adaptive, band_shift, result)
        };
    }
    match (approx_max, with_cigar, right_align, has_n) {
        (false, false, false, false) => go!(false, false, false, false),
        (false, false, false, true) => go!(false, false, false, true),
        (false, false, true, false) => go!(false, false, true, false),
        (false, false, true, true) => go!(false, false, true, true),
        (false, true, false, false) => go!(false, true, false, false),
        (false, true, false, true) => go!(false, true, false, true),
        (false, true, true, false) => go!(false, true, true, false),
        (false, true, true, true) => go!(false, true, true, true),
        (true, false, false, false) => go!(true, false, false, false),
        (true, false, false, true) => go!(true, false, false, true),
        (true, false, true, false) => go!(true, false, true, false),
        (true, false, true, true) => go!(true, false, true, true),
        (true, true, false, false) => go!(true, true, false, false),
        (true, true, false, true) => go!(true, true, false, true),
        (true, true, true, false) => go!(true, true, true, false),
        (true, true, true, true) => go!(true, true, true, true),
    }
}}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn clean_impl<const APPROX: bool, const CIGAR: bool, const RIGHT: bool, const HAS_N: bool>(
    qseq: &[u8],
    tseq: &[u8],
    alphabet_size: i8,
    score_matrix: &[i8],
    gap_open: i8,
    gap_extend: i8,
    gap_open2: i8,
    gap_extend2: i8,
    bandwidth: i32,
    z_drop: i32,
    end_bonus: i32,
    flags: i32,
    adaptive: i32,
    band_shift: i32,
    result: &mut DpResult,
) { unsafe {
    let query_len = qseq.len();
    let target_len = tseq.len();
    let approx_max = APPROX;
    let with_cigar = CIGAR;
    let right_align = RIGHT;

    if alphabet_size <= 1 || query_len == 0 || target_len == 0 {
        return;
    }

    let (gap_open, gap_extend, gap_open2, gap_extend2) =
        if (gap_open2 as i32 + gap_extend2 as i32) < (gap_open as i32 + gap_extend as i32) {
            (gap_open2, gap_extend2, gap_open, gap_extend)
        } else {
            (gap_open, gap_extend, gap_open2, gap_extend2)
        };
    let (qi, ei, q2i, e2i) = (gap_open as i32, gap_extend as i32, gap_open2 as i32, gap_extend2 as i32);

    let mut long_thres: i32 = if ei != e2i { (q2i - qi) / (ei - e2i) - 1 } else { 0 };
    if (q2i + e2i + long_thres * e2i) > (qi + ei + long_thres * ei) {
        long_thres += 1;
    }
    let long_diff: i8 = (long_thres * (ei - e2i) - (q2i - qi) - e2i) as i8;

    let neg_qe = (-(qi + ei)) as i8;
    let neg_q2e2 = (-(q2i + e2i)) as i8;

    let zero_ = _mm256_setzero_si256();
    let q_ = _mm256_set1_epi8(gap_open);
    let q2_ = _mm256_set1_epi8(gap_open2);
    let qe_ = _mm256_set1_epi8((qi + ei) as i8);
    let qe2_ = _mm256_set1_epi8((q2i + e2i) as i8);
    let neg_qe_ = _mm256_set1_epi8(neg_qe);
    let neg_q2e2_ = _mm256_set1_epi8(neg_q2e2);
    let sc_mch_ = _mm256_set1_epi8(score_matrix[0]);
    let sc_mis_ = _mm256_set1_epi8(score_matrix[1]);
    let n_sc = score_matrix[(alphabet_size as usize * alphabet_size as usize) - 1];
    let sc_n_ = _mm256_set1_epi8(if n_sc == 0 { -gap_extend2 } else { n_sc });
    let m1_ = _mm256_set1_epi8(alphabet_size - 1);
    let sc_mch_mis_ = _mm256_set1_epi8(score_matrix[0].wrapping_sub(score_matrix[1]));
    let has_n = HAS_N;
    let flag1_ = _mm256_set1_epi8(1);
    let flag2_ = _mm256_set1_epi8(2);
    let flag3_ = _mm256_set1_epi8(3);
    let flag4_ = _mm256_set1_epi8(4);
    let flag8_ = _mm256_set1_epi8(0x08);
    let flag16_ = _mm256_set1_epi8(0x10);
    let flag32_ = _mm256_set1_epi8(0x20);
    let flag64_ = _mm256_set1_epi8(0x40);

    let wl = if bandwidth < 0 { target_len.max(query_len) as i32 } else { bandwidth };
    // Adaptive band: a window of `adaptive` cells per anti-diagonal that moves one step towards
    // the anti-diagonal's best cell (exact max tracking) or the greedy path cell (approximate
    // tracking) each time, instead of the static band |t - q| <= wl. Its start and end never
    // decrease, so the scalar boundary semantics carry over unchanged; the matrix bound
    // st0 >= r - qlen + 1 pulls it onto the end corner for global alignment.
    let adaptive = if adaptive > 0 { adaptive.min(wl.max(1)) } else { 0 };
    let band_w = if adaptive > 0 { adaptive } else { wl };
    let mut a_st: i32 = 0;
    let mut a_step: i32 = 0;

    // Layout: six state arrays and the target, each tlen + 64 bytes (vector loads and the
    // post-loop reset run up to 32 bytes past en0); the reversed query, qlen + 64 bytes; then
    // traceback rows and band offsets.
    let tpad = (target_len + 64).next_multiple_of(32);
    let qpad = (query_len + 64).next_multiple_of(32);
    let n_col = query_len.min(target_len).min((band_w + 1) as usize).div_ceil(32) + 1;
    let stride = n_col * 32;
    let n_diag = query_len + target_len - 1;
    let state_bytes = 7 * tpad + qpad;
    let (p_off, off_off, off_end_off, total) = if with_cigar {
        let p_off = state_bytes;
        let off_off = (p_off + n_diag * stride + 31) & !31;
        let off_end_off = (off_off + n_diag * 4 + 31) & !31;
        (p_off, off_off, off_end_off, off_end_off + n_diag * 4)
    } else {
        (0, 0, 0, state_bytes)
    };
    let mem = AlignedMemory::new(total, 32);
    let base = mem.as_ptr();
    let u = base;
    let v = u.add(tpad);
    let x = v.add(tpad);
    let y = x.add(tpad);
    let x2 = y.add(tpad);
    let y2 = x2.add(tpad);
    let sf = y2.add(tpad);
    let qr = sf.add(tpad);
    for a in [u, v, x, y] {
        std::ptr::write_bytes(a, neg_qe as u8, tpad);
    }
    std::ptr::write_bytes(x2, neg_q2e2 as u8, tpad);
    std::ptr::write_bytes(y2, neg_q2e2 as u8, tpad);
    std::ptr::copy_nonoverlapping(tseq.as_ptr(), sf, target_len);
    std::ptr::write_bytes(sf.add(target_len), 0, tpad - target_len);
    for t in 0..query_len {
        *qr.add(t) = qseq[query_len - 1 - t];
    }
    std::ptr::write_bytes(qr.add(query_len), 0, qpad - query_len);
    let (p, band_off, band_off_end) = if with_cigar {
        (base.add(p_off), base.add(off_off) as *mut i32, base.add(off_end_off) as *mut i32)
    } else {
        (std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut())
    };

    let (h_vec, h_ptr) = alloc_h_array(approx_max, tpad / 32, 32);
    let _ = &h_vec;
    // See the default kernel's `packed_h` for the range argument.
    let packed_h = query_len + target_len <= 20000 && wl >= 32;

    init_dp_result(result);

    let mut last_st: i32 = -1;
    let mut last_en: i32 = -1;
    let mut h0: i32 = 0;
    let mut last_h0_t: i32 = 0;
    let qe_scalar = qi + ei;

    for r in 0..n_diag as i32 {
        let mut st0 = 0i32;
        let mut en0 = target_len as i32 - 1;
        if st0 < r - query_len as i32 + 1 { st0 = r - query_len as i32 + 1; }
        if en0 > r { en0 = r; }
        if adaptive > 0 {
            a_st = (a_st + a_step).max(st0);
            st0 = a_st;
            en0 = en0.min(st0 + adaptive - 1);
        } else {
            // Static band: diagonals t - q in [band_shift - wl, band_shift + wl].
            if st0 < (r - wl + band_shift + 1) >> 1 { st0 = (r - wl + band_shift + 1) >> 1; }
            if en0 > (r + wl + band_shift) >> 1 { en0 = (r + wl + band_shift) >> 1; }
        }
        if st0 > en0 {
            result.zdropped = 1;
            break;
        }

        let (x1, x21, v1): (i8, i8, i8) = if st0 > 0 {
            if st0 > last_st && st0 - 1 <= last_en {
                let i = (st0 - 1) as usize;
                (*x.add(i) as i8, *x2.add(i) as i8, *v.add(i) as i8)
            } else {
                (neg_qe, neg_q2e2, neg_qe)
            }
        } else {
            let v1 = if r == 0 { neg_qe } else if r < long_thres { -gap_extend } else if r == long_thres { long_diff } else { -gap_extend2 };
            (neg_qe, neg_q2e2, v1)
        };
        if en0 >= r {
            let i = r as usize;
            *y.add(i) = neg_qe as u8;
            *y2.add(i) = neg_q2e2 as u8;
            *u.add(i) = (if r == 0 { neg_qe } else if r < long_thres { -gap_extend } else if r == long_thres { long_diff } else { -gap_extend2 }) as u8;
        }

        // Previous raw vectors (byte 31 carries into the next shift); the first carries the
        // left boundary.
        // The shift only reads byte 31 of the previous vector: a broadcast is enough.
        let mut px = _mm256_set1_epi8(x1);
        let mut px2 = _mm256_set1_epi8(x21);
        let mut pv = _mm256_set1_epi8(v1);
        let qrr = qr.offset(query_len as isize - 1 - r as isize);
        let st = st0 as usize;
        let en = en0 as usize;
        let prow = if with_cigar {
            *band_off.add(r as usize) = st0;
            *band_off_end.add(r as usize) = en0;
            p.add(r as usize * stride).wrapping_sub(st)
        } else {
            std::ptr::null_mut()
        };

        let mut t = st;
        while t <= en {
            let sq = _mm256_loadu_si256(sf.add(t) as *const __m256i);
            let qv = _mm256_loadu_si256(qrr.add(t) as *const __m256i);
            let mut z = _mm256_add_epi8(sc_mis_, _mm256_and_si256(_mm256_cmpeq_epi8(sq, qv), sc_mch_mis_));
            if has_n {
                // nt4 codes are 0..=m1 (m1 = N): either base is N iff their max is m1.
                z = _mm256_blendv_epi8(z, sc_n_, _mm256_cmpeq_epi8(_mm256_max_epu8(sq, qv), m1_));
            }

            let xt = _mm256_loadu_si256(x.add(t) as *const __m256i);
            let xt1 = super::common::avx2_shl1_prev(xt, px);
            px = xt;
            let vt = _mm256_loadu_si256(v.add(t) as *const __m256i);
            let vt1 = super::common::avx2_shl1_prev(vt, pv);
            pv = vt;
            let x2t = _mm256_loadu_si256(x2.add(t) as *const __m256i);
            let x2t1 = super::common::avx2_shl1_prev(x2t, px2);
            px2 = x2t;
            let ut = _mm256_loadu_si256(u.add(t) as *const __m256i);
            let a = _mm256_add_epi8(xt1, vt1);
            let b = _mm256_add_epi8(_mm256_loadu_si256(y.add(t) as *const __m256i), ut);
            let a2 = _mm256_add_epi8(x2t1, vt1);
            let b2 = _mm256_add_epi8(_mm256_loadu_si256(y2.add(t) as *const __m256i), ut);

            if !with_cigar {
                z = _mm256_max_epi8(z, a);
                z = _mm256_max_epi8(z, b);
                z = _mm256_max_epi8(z, a2);
                z = _mm256_max_epi8(z, b2);
                z = _mm256_min_epi8(z, sc_mch_);
                _mm256_storeu_si256(u.add(t) as *mut __m256i, _mm256_sub_epi8(z, vt1));
                _mm256_storeu_si256(v.add(t) as *mut __m256i, _mm256_sub_epi8(z, ut));
                let tmp1 = _mm256_sub_epi8(z, q_);
                let tmp2 = _mm256_sub_epi8(z, q2_);
                _mm256_storeu_si256(x.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(a, tmp1), zero_), qe_));
                _mm256_storeu_si256(y.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(b, tmp1), zero_), qe_));
                _mm256_storeu_si256(x2.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(a2, tmp2), zero_), qe2_));
                _mm256_storeu_si256(y2.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(b2, tmp2), zero_), qe2_));
            } else {
                let mut d;
                if !right_align {
                    let tmp = _mm256_cmpgt_epi8(a, z);
                    d = _mm256_and_si256(tmp, flag1_);
                    z = _mm256_max_epi8(z, a);
                    d = _mm256_blendv_epi8(d, flag2_, _mm256_cmpgt_epi8(b, z));
                    z = _mm256_max_epi8(z, b);
                    d = _mm256_blendv_epi8(d, flag3_, _mm256_cmpgt_epi8(a2, z));
                    z = _mm256_max_epi8(z, a2);
                    d = _mm256_blendv_epi8(d, flag4_, _mm256_cmpgt_epi8(b2, z));
                    z = _mm256_max_epi8(z, b2);
                } else {
                    d = _mm256_andnot_si256(_mm256_cmpgt_epi8(z, a), flag1_);
                    z = _mm256_max_epi8(z, a);
                    d = _mm256_blendv_epi8(flag2_, d, _mm256_cmpgt_epi8(z, b));
                    z = _mm256_max_epi8(z, b);
                    d = _mm256_blendv_epi8(flag3_, d, _mm256_cmpgt_epi8(z, a2));
                    z = _mm256_max_epi8(z, a2);
                    d = _mm256_blendv_epi8(flag4_, d, _mm256_cmpgt_epi8(z, b2));
                    z = _mm256_max_epi8(z, b2);
                }
                z = _mm256_min_epi8(z, sc_mch_);
                _mm256_storeu_si256(u.add(t) as *mut __m256i, _mm256_sub_epi8(z, vt1));
                _mm256_storeu_si256(v.add(t) as *mut __m256i, _mm256_sub_epi8(z, ut));
                let tmp1 = _mm256_sub_epi8(z, q_);
                let tmp2 = _mm256_sub_epi8(z, q2_);
                let an = _mm256_sub_epi8(a, tmp1);
                let bn = _mm256_sub_epi8(b, tmp1);
                let a2n = _mm256_sub_epi8(a2, tmp2);
                let b2n = _mm256_sub_epi8(b2, tmp2);
                // Left alignment extends a gap when the new value is > 0, right alignment when >= 0.
                let (ma, mb, ma2, mb2) = if !right_align {
                    (_mm256_cmpgt_epi8(an, zero_), _mm256_cmpgt_epi8(bn, zero_),
                     _mm256_cmpgt_epi8(a2n, zero_), _mm256_cmpgt_epi8(b2n, zero_))
                } else {
                    (_mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, an), _mm256_set1_epi8(-1)),
                     _mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, bn), _mm256_set1_epi8(-1)),
                     _mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, a2n), _mm256_set1_epi8(-1)),
                     _mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, b2n), _mm256_set1_epi8(-1)))
                };
                _mm256_storeu_si256(x.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_and_si256(ma, an), qe_));
                _mm256_storeu_si256(y.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_and_si256(mb, bn), qe_));
                _mm256_storeu_si256(x2.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_and_si256(ma2, a2n), qe2_));
                _mm256_storeu_si256(y2.add(t) as *mut __m256i, _mm256_sub_epi8(_mm256_and_si256(mb2, b2n), qe2_));
                d = _mm256_or_si256(d, _mm256_and_si256(ma, flag8_));
                d = _mm256_or_si256(d, _mm256_and_si256(mb, flag16_));
                d = _mm256_or_si256(d, _mm256_and_si256(ma2, flag32_));
                d = _mm256_or_si256(d, _mm256_and_si256(mb2, flag64_));
                _mm256_storeu_si256(prow.add(t) as *mut __m256i, d);
            }
            t += 32;
        }
        // The last vector may have written up to 31 cells right of en0; cells right of the band
        // must keep their initial values (the scalar kernel never touches them).
        let es = en + 1;
        // Nothing was written past en0 when the band is a whole number of vectors.
        if (en - st + 1) % 32 != 0 {
            _mm256_storeu_si256(u.add(es) as *mut __m256i, neg_qe_);
            _mm256_storeu_si256(v.add(es) as *mut __m256i, neg_qe_);
            _mm256_storeu_si256(x.add(es) as *mut __m256i, neg_qe_);
            _mm256_storeu_si256(y.add(es) as *mut __m256i, neg_qe_);
            _mm256_storeu_si256(x2.add(es) as *mut __m256i, neg_q2e2_);
            _mm256_storeu_si256(y2.add(es) as *mut __m256i, neg_q2e2_);
        }

        if !approx_max {
            let mut max_h: i32;
            let mut max_t: i32;
            if r > 0 {
                let h_en0 = if en0 > 0 {
                    *h_ptr.add(en - 1) + *u.add(en) as i8 as i32
                } else {
                    *h_ptr.add(en) + *v.add(en) as i8 as i32
                };
                *h_ptr.add(en) = h_en0;
                max_h = h_en0;
                max_t = en0;
                let en1 = st0 + (en0 - st0) / 4 * 4;
                let en8 = st0 + (en1 - st0) / 8 * 8;
                let mut tt = st0;
                if packed_h {
                    // Packed (H, lane, position) keys: see the default kernel.
                    let mut low = _mm256_setr_epi32(32767, 32767 - (1 << 13), 32767 - (2 << 13), 32767 - (3 << 13),
                                                    32766, 32766 - (1 << 13), 32766 - (2 << 13), 32766 - (3 << 13));
                    let two = _mm256_set1_epi32(2);
                    let mut mk = _mm256_set1_epi32(i32::MIN);
                    while tt < en8 {
                        let hp = h_ptr.add(tt as usize);
                        let vv = _mm256_cvtepi8_epi32(_mm_loadl_epi64(v.add(tt as usize) as *const __m128i));
                        let h1 = _mm256_add_epi32(_mm256_loadu_si256(hp as *const __m256i), vv);
                        _mm256_storeu_si256(hp as *mut __m256i, h1);
                        mk = _mm256_max_epi32(mk, _mm256_or_si256(_mm256_slli_epi32(h1, 15), low));
                        low = _mm256_sub_epi32(low, two);
                        tt += 8;
                    }
                    let mut m4 = _mm_max_epi32(_mm256_castsi256_si128(mk), _mm256_extracti128_si256(mk, 1));
                    if tt < en1 {
                        let hp = h_ptr.add(tt as usize);
                        let vv = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(std::ptr::read_unaligned(v.add(tt as usize) as *const i32)));
                        let h1 = _mm_add_epi32(_mm_loadu_si128(hp as *const __m128i), vv);
                        _mm_storeu_si128(hp as *mut __m128i, h1);
                        m4 = _mm_max_epi32(m4, _mm_or_si128(_mm_slli_epi32(h1, 15), _mm256_castsi256_si128(low)));
                        tt += 4;
                    }
                    m4 = _mm_max_epi32(m4, _mm_shuffle_epi32(m4, 0b01_00_11_10));
                    m4 = _mm_max_epi32(m4, _mm_shuffle_epi32(m4, 0b10_11_00_01));
                    let key = _mm_cvtsi128_si32(m4);
                    if key != i32::MIN && (key >> 15) > max_h {
                        let rank = 32767 - (key & 0x7fff);
                        max_h = key >> 15;
                        max_t = st0 + (rank & 0x1fff) * 4 + (rank >> 13);
                    }
                } else {
                    let mut lane_h = [max_h; 4];
                    let mut lane_t = [max_t; 4];
                    while tt < en1 {
                        for i in 0..4i32 {
                            let pos = (tt + i) as usize;
                            *h_ptr.add(pos) += *v.add(pos) as i8 as i32;
                            if *h_ptr.add(pos) > lane_h[i as usize] {
                                lane_h[i as usize] = *h_ptr.add(pos);
                                lane_t[i as usize] = tt + i;
                            }
                        }
                        tt += 4;
                    }
                    for i in 0..4 {
                        if max_h < lane_h[i] { max_h = lane_h[i]; max_t = lane_t[i]; }
                    }
                }
                while tt < en0 {
                    *h_ptr.add(tt as usize) += *v.add(tt as usize) as i8 as i32;
                    if *h_ptr.add(tt as usize) > max_h {
                        max_h = *h_ptr.add(tt as usize);
                        max_t = tt;
                    }
                    tt += 1;
                }
            } else {
                *h_ptr = *v as i8 as i32 - qe_scalar;
                max_h = *h_ptr;
                max_t = 0;
            }
            if en0 == target_len as i32 - 1 && *h_ptr.add(en) > result.max_target_end_score {
                result.max_target_end_score = *h_ptr.add(en);
                result.max_target_end_query_pos = r - en0;
            }
            if r - st0 == query_len as i32 - 1 && *h_ptr.add(st) > result.max_query_end_score {
                result.max_query_end_score = *h_ptr.add(st);
                result.max_query_end_target_pos = st0;
            }
            if max_h > result.max {
                result.max = max_h;
                result.max_score_target_pos = max_t;
                result.max_score_query_pos = r - max_t;
            } else if max_t >= result.max_score_target_pos && (r - max_t) >= result.max_score_query_pos {
                let tl = max_t - result.max_score_target_pos;
                let ql = (r - max_t) - result.max_score_query_pos;
                let l = if tl > ql { tl - ql } else { ql - tl };
                if z_drop >= 0 && (result.max - max_h) > z_drop + l * e2i {
                    result.zdropped = 1;
                    break;
                }
            }
            if adaptive > 0 {
                // Move towards this anti-diagonal's best cell once the window is full: the next
                // start advances iff that cell lies right of the window's centre.
                a_step = (en0 - st0 + 1 >= adaptive && 2 * max_t > st0 + en0) as i32;
            }
            if r == query_len as i32 + target_len as i32 - 2 && en0 == target_len as i32 - 1 {
                result.score = *h_ptr.add(target_len - 1);
            }
        } else {
            if r > 0 {
                if last_h0_t >= st0 && last_h0_t <= en0 && last_h0_t + 1 >= st0 && last_h0_t < en0 {
                    // which neighbour the greedy path takes is data-dependent: select, don't branch
                    let d0 = *v.add(last_h0_t as usize) as i8 as i32;
                    let d1 = *u.add((last_h0_t + 1) as usize) as i8 as i32;
                    let right = d0 <= d1;
                    h0 += if right { d1 } else { d0 };
                    last_h0_t += right as i32;
                } else if last_h0_t >= st0 && last_h0_t <= en0 {
                    h0 += *v.add(last_h0_t as usize) as i8 as i32;
                } else {
                    last_h0_t += 1;
                    h0 += *u.add(last_h0_t as usize) as i8 as i32;
                }
            } else {
                h0 = *v as i8 as i32 - qe_scalar;
                last_h0_t = 0;
            }
            if adaptive > 0 {
                a_step = (en0 - st0 + 1 >= adaptive && 2 * last_h0_t > st0 + en0) as i32;
            }
            if (flags & APPROX_DROP) != 0 {
                if h0 > result.max {
                    result.max = h0;
                    result.max_score_target_pos = last_h0_t;
                    result.max_score_query_pos = r - last_h0_t;
                } else if last_h0_t >= result.max_score_target_pos && (r - last_h0_t) >= result.max_score_query_pos {
                    let tl = last_h0_t - result.max_score_target_pos;
                    let ql = (r - last_h0_t) - result.max_score_query_pos;
                    let l = if tl > ql { tl - ql } else { ql - tl };
                    if z_drop >= 0 && (result.max - h0) > z_drop + l * e2i {
                        result.zdropped = 1;
                        break;
                    }
                }
            }
            if r == query_len as i32 + target_len as i32 - 2 && en0 == target_len as i32 - 1 {
                result.score = h0;
            }
        }

        last_st = st0;
        last_en = en0;
    }

    if approx_max && result.score == NEG_INF {
        result.score = result.max;
    }
    if with_cigar {
        traceback_dual_affine(result, query_len, target_len, end_bonus, flags, n_col, 32, p, band_off, band_off_end);
    }
}}

/// BONEFISH_REG1_OFF=1: use the memory-array kernel for single-vector bands too (comparisons).
static REG1_OFF: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("BONEFISH_REG1_OFF").is_ok_and(|v| v != "0"));

/// Bytes 1..=31 of `v` moved down one lane, `fill` byte 0 into lane 31.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn shr1_next(v: __m256i, fill: __m256i) -> __m256i {
    _mm256_alignr_epi8(_mm256_permute2x128_si256(v, fill, 0x21), v, 1)
}

/// `clean_impl` for approximate-max tracking with a static band of at most 32 cells per
/// anti-diagonal (bandwidth <= 31), the common short gap fill. The six difference vectors live in
/// registers indexed from the anti-diagonal's band start instead of in target-length arrays:
/// each anti-diagonal reads the previous one's values at the same target positions, so when the
/// band start advances by one the registers move down one lane, the incoming top lane taking the
/// value the arrays hold beyond any band (their initial value). Same recurrences, boundary values
/// and post-band resets as `clean_impl`, so the results are identical; it avoids the array
/// round trip (stores re-read one byte off, which defeats store forwarding).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
unsafe fn reg1_impl<const CIGAR: bool, const RIGHT: bool, const HAS_N: bool>(
    qseq: &[u8], tseq: &[u8], alphabet_size: i8, score_matrix: &[i8],
    gap_open: i8, gap_extend: i8, gap_open2: i8, gap_extend2: i8,
    bandwidth: i32, z_drop: i32, end_bonus: i32, flags: i32, band_shift: i32,
    result: &mut DpResult,
) { unsafe {
    let query_len = qseq.len();
    let target_len = tseq.len();
    let with_cigar = CIGAR;
    let right_align = RIGHT;
    if alphabet_size <= 1 || query_len == 0 || target_len == 0 {
        return;
    }
    let (gap_open, gap_extend, gap_open2, gap_extend2) =
        if (gap_open2 as i32 + gap_extend2 as i32) < (gap_open as i32 + gap_extend as i32) {
            (gap_open2, gap_extend2, gap_open, gap_extend)
        } else {
            (gap_open, gap_extend, gap_open2, gap_extend2)
        };
    let (qi, ei, q2i, e2i) = (gap_open as i32, gap_extend as i32, gap_open2 as i32, gap_extend2 as i32);
    let mut long_thres: i32 = if ei != e2i { (q2i - qi) / (ei - e2i) - 1 } else { 0 };
    if (q2i + e2i + long_thres * e2i) > (qi + ei + long_thres * ei) {
        long_thres += 1;
    }
    let long_diff: i8 = (long_thres * (ei - e2i) - (q2i - qi) - e2i) as i8;
    let neg_qe = (-(qi + ei)) as i8;
    let neg_q2e2 = (-(q2i + e2i)) as i8;

    let zero_ = _mm256_setzero_si256();
    let q_ = _mm256_set1_epi8(gap_open);
    let q2_ = _mm256_set1_epi8(gap_open2);
    let qe_ = _mm256_set1_epi8((qi + ei) as i8);
    let qe2_ = _mm256_set1_epi8((q2i + e2i) as i8);
    let neg_qe_ = _mm256_set1_epi8(neg_qe);
    let neg_q2e2_ = _mm256_set1_epi8(neg_q2e2);
    let sc_mch_ = _mm256_set1_epi8(score_matrix[0]);
    let sc_mis_ = _mm256_set1_epi8(score_matrix[1]);
    let n_sc = score_matrix[(alphabet_size as usize * alphabet_size as usize) - 1];
    let sc_n_ = _mm256_set1_epi8(if n_sc == 0 { -gap_extend2 } else { n_sc });
    let m1_ = _mm256_set1_epi8(alphabet_size - 1);
    let sc_mch_mis_ = _mm256_set1_epi8(score_matrix[0].wrapping_sub(score_matrix[1]));
    let flag1_ = _mm256_set1_epi8(1);
    let flag2_ = _mm256_set1_epi8(2);
    let flag3_ = _mm256_set1_epi8(3);
    let flag4_ = _mm256_set1_epi8(4);
    let flag8_ = _mm256_set1_epi8(0x08);
    let flag16_ = _mm256_set1_epi8(0x10);
    let flag32_ = _mm256_set1_epi8(0x20);
    let flag64_ = _mm256_set1_epi8(0x40);
    let lane_ = _mm256_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
    let wl = bandwidth;

    // Target and reversed query, padded for 32-byte loads; traceback rows and band offsets.
    let tpad = (target_len + 64).next_multiple_of(32);
    let qpad = (query_len + 64).next_multiple_of(32);
    let n_col = query_len.min(target_len).min((wl + 1) as usize).div_ceil(32) + 1;
    let stride = n_col * 32;
    let n_diag = query_len + target_len - 1;
    let state_bytes = tpad + qpad;
    let (p_off, off_off, off_end_off, total) = if with_cigar {
        let p_off = state_bytes;
        let off_off = (p_off + n_diag * stride + 31) & !31;
        let off_end_off = (off_off + n_diag * 4 + 31) & !31;
        (p_off, off_off, off_end_off, off_end_off + n_diag * 4)
    } else {
        (0, 0, 0, state_bytes)
    };
    let mem = AlignedMemory::new(total, 32);
    let base = mem.as_ptr();
    let sf = base;
    let qr = sf.add(tpad);
    std::ptr::copy_nonoverlapping(tseq.as_ptr(), sf, target_len);
    std::ptr::write_bytes(sf.add(target_len), 0, tpad - target_len);
    for t in 0..query_len {
        *qr.add(t) = qseq[query_len - 1 - t];
    }
    std::ptr::write_bytes(qr.add(query_len), 0, qpad - query_len);
    let (p, band_off, band_off_end) = if with_cigar {
        (base.add(p_off), base.add(off_off) as *mut i32, base.add(off_end_off) as *mut i32)
    } else {
        (std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut())
    };

    init_dp_result(result);

    // Registers: lane j holds the array value at target position `lo + j`.
    let (mut ur, mut vr, mut xr, mut yr) = (neg_qe_, neg_qe_, neg_qe_, neg_qe_);
    let (mut x2r, mut y2r) = (neg_q2e2_, neg_q2e2_);
    let mut lo: i32 = 0;
    let mut last_st: i32 = -1;
    let mut last_en: i32 = -1;
    let mut h0: i32 = 0;
    let mut last_h0_t: i32 = 0;
    let qe_scalar = qi + ei;
    #[repr(align(32))]
    struct Lanes([i8; 32]);
    let mut ub = Lanes([0; 32]);
    let mut vb = Lanes([0; 32]);
    let top_u = |r: i32| -> i8 { if r == 0 { neg_qe } else if r < long_thres { -gap_extend } else if r == long_thres { long_diff } else { -gap_extend2 } };

    for r in 0..n_diag as i32 {
        let mut st0 = 0i32;
        let mut en0 = target_len as i32 - 1;
        if st0 < r - query_len as i32 + 1 { st0 = r - query_len as i32 + 1; }
        if en0 > r { en0 = r; }
        if st0 < (r - wl + band_shift + 1) >> 1 { st0 = (r - wl + band_shift + 1) >> 1; }
        if en0 > (r + wl + band_shift) >> 1 { en0 = (r + wl + band_shift) >> 1; }
        if st0 > en0 {
            result.zdropped = 1;
            break;
        }
        debug_assert!(en0 - st0 < 32);

        // Left boundary: the previous anti-diagonal's value at st0 - 1 (lane st0 - 1 - lo).
        let (x1, x21, v1): (i8, i8, i8) = if st0 > 0 {
            if st0 > last_st && st0 - 1 <= last_en {
                let l = (st0 - 1 - lo) as usize;
                _mm256_store_si256(ub.0.as_mut_ptr() as *mut __m256i, xr);
                let xv = ub.0[l];
                _mm256_store_si256(ub.0.as_mut_ptr() as *mut __m256i, x2r);
                let x2v = ub.0[l];
                _mm256_store_si256(ub.0.as_mut_ptr() as *mut __m256i, vr);
                (xv, x2v, ub.0[l])
            } else {
                (neg_qe, neg_q2e2, neg_qe)
            }
        } else {
            (neg_qe, neg_q2e2, top_u(r))
        };
        // Re-index the registers from st0 (the band start advances by at most one).
        while lo < st0 {
            ur = shr1_next(ur, neg_qe_);
            vr = shr1_next(vr, neg_qe_);
            xr = shr1_next(xr, neg_qe_);
            yr = shr1_next(yr, neg_qe_);
            x2r = shr1_next(x2r, neg_q2e2_);
            y2r = shr1_next(y2r, neg_q2e2_);
            lo += 1;
        }
        if en0 >= r {
            // New top cell at target position r.
            let m = _mm256_cmpeq_epi8(lane_, _mm256_set1_epi8((r - st0) as i8));
            yr = _mm256_blendv_epi8(yr, neg_qe_, m);
            y2r = _mm256_blendv_epi8(y2r, neg_q2e2_, m);
            ur = _mm256_blendv_epi8(ur, _mm256_set1_epi8(top_u(r)), m);
        }

        let qrr = qr.offset(query_len as isize - 1 - r as isize);
        let st = st0 as usize;
        let en = en0 as usize;
        let prow = if with_cigar {
            *band_off.add(r as usize) = st0;
            *band_off_end.add(r as usize) = en0;
            p.add(r as usize * stride).wrapping_sub(st)
        } else {
            std::ptr::null_mut()
        };

        let sq = _mm256_loadu_si256(sf.add(st) as *const __m256i);
        let qv = _mm256_loadu_si256(qrr.add(st) as *const __m256i);
        let mut z = _mm256_add_epi8(sc_mis_, _mm256_and_si256(_mm256_cmpeq_epi8(sq, qv), sc_mch_mis_));
        if HAS_N {
            z = _mm256_blendv_epi8(z, sc_n_, _mm256_cmpeq_epi8(_mm256_max_epu8(sq, qv), m1_));
        }
        let xt1 = super::common::avx2_shl1_prev(xr, _mm256_set1_epi8(x1));
        let vt1 = super::common::avx2_shl1_prev(vr, _mm256_set1_epi8(v1));
        let x2t1 = super::common::avx2_shl1_prev(x2r, _mm256_set1_epi8(x21));
        let ut = ur;
        let a = _mm256_add_epi8(xt1, vt1);
        let b = _mm256_add_epi8(yr, ut);
        let a2 = _mm256_add_epi8(x2t1, vt1);
        let b2 = _mm256_add_epi8(y2r, ut);
        if !with_cigar {
            z = _mm256_max_epi8(z, a);
            z = _mm256_max_epi8(z, b);
            z = _mm256_max_epi8(z, a2);
            z = _mm256_max_epi8(z, b2);
            z = _mm256_min_epi8(z, sc_mch_);
            ur = _mm256_sub_epi8(z, vt1);
            vr = _mm256_sub_epi8(z, ut);
            let tmp1 = _mm256_sub_epi8(z, q_);
            let tmp2 = _mm256_sub_epi8(z, q2_);
            xr = _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(a, tmp1), zero_), qe_);
            yr = _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(b, tmp1), zero_), qe_);
            x2r = _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(a2, tmp2), zero_), qe2_);
            y2r = _mm256_sub_epi8(_mm256_max_epi8(_mm256_sub_epi8(b2, tmp2), zero_), qe2_);
        } else {
            let mut d;
            if !right_align {
                let tmp = _mm256_cmpgt_epi8(a, z);
                d = _mm256_and_si256(tmp, flag1_);
                z = _mm256_max_epi8(z, a);
                d = _mm256_blendv_epi8(d, flag2_, _mm256_cmpgt_epi8(b, z));
                z = _mm256_max_epi8(z, b);
                d = _mm256_blendv_epi8(d, flag3_, _mm256_cmpgt_epi8(a2, z));
                z = _mm256_max_epi8(z, a2);
                d = _mm256_blendv_epi8(d, flag4_, _mm256_cmpgt_epi8(b2, z));
                z = _mm256_max_epi8(z, b2);
            } else {
                d = _mm256_andnot_si256(_mm256_cmpgt_epi8(z, a), flag1_);
                z = _mm256_max_epi8(z, a);
                d = _mm256_blendv_epi8(flag2_, d, _mm256_cmpgt_epi8(z, b));
                z = _mm256_max_epi8(z, b);
                d = _mm256_blendv_epi8(flag3_, d, _mm256_cmpgt_epi8(z, a2));
                z = _mm256_max_epi8(z, a2);
                d = _mm256_blendv_epi8(flag4_, d, _mm256_cmpgt_epi8(z, b2));
                z = _mm256_max_epi8(z, b2);
            }
            z = _mm256_min_epi8(z, sc_mch_);
            ur = _mm256_sub_epi8(z, vt1);
            vr = _mm256_sub_epi8(z, ut);
            let tmp1 = _mm256_sub_epi8(z, q_);
            let tmp2 = _mm256_sub_epi8(z, q2_);
            let an = _mm256_sub_epi8(a, tmp1);
            let bn = _mm256_sub_epi8(b, tmp1);
            let a2n = _mm256_sub_epi8(a2, tmp2);
            let b2n = _mm256_sub_epi8(b2, tmp2);
            let (ma, mb, ma2, mb2) = if !right_align {
                (_mm256_cmpgt_epi8(an, zero_), _mm256_cmpgt_epi8(bn, zero_),
                 _mm256_cmpgt_epi8(a2n, zero_), _mm256_cmpgt_epi8(b2n, zero_))
            } else {
                (_mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, an), _mm256_set1_epi8(-1)),
                 _mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, bn), _mm256_set1_epi8(-1)),
                 _mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, a2n), _mm256_set1_epi8(-1)),
                 _mm256_andnot_si256(_mm256_cmpgt_epi8(zero_, b2n), _mm256_set1_epi8(-1)))
            };
            xr = _mm256_sub_epi8(_mm256_and_si256(ma, an), qe_);
            yr = _mm256_sub_epi8(_mm256_and_si256(mb, bn), qe_);
            x2r = _mm256_sub_epi8(_mm256_and_si256(ma2, a2n), qe2_);
            y2r = _mm256_sub_epi8(_mm256_and_si256(mb2, b2n), qe2_);
            d = _mm256_or_si256(d, _mm256_and_si256(ma, flag8_));
            d = _mm256_or_si256(d, _mm256_and_si256(mb, flag16_));
            d = _mm256_or_si256(d, _mm256_and_si256(ma2, flag32_));
            d = _mm256_or_si256(d, _mm256_and_si256(mb2, flag64_));
            _mm256_storeu_si256(prow.add(st) as *mut __m256i, d);
        }
        // Lanes right of the band keep the arrays' initial values (the post-band reset).
        if (en - st + 1) % 32 != 0 {
            let out = _mm256_cmpgt_epi8(lane_, _mm256_set1_epi8((en0 - st0) as i8));
            ur = _mm256_blendv_epi8(ur, neg_qe_, out);
            vr = _mm256_blendv_epi8(vr, neg_qe_, out);
            xr = _mm256_blendv_epi8(xr, neg_qe_, out);
            yr = _mm256_blendv_epi8(yr, neg_qe_, out);
            x2r = _mm256_blendv_epi8(x2r, neg_q2e2_, out);
            y2r = _mm256_blendv_epi8(y2r, neg_q2e2_, out);
        }

        // Greedy-path score (approximate max), as clean_impl.
        _mm256_store_si256(ub.0.as_mut_ptr() as *mut __m256i, ur);
        _mm256_store_si256(vb.0.as_mut_ptr() as *mut __m256i, vr);
        let uat = |t: i32| ub.0[(t - st0) as usize] as i32;
        let vat = |t: i32| vb.0[(t - st0) as usize] as i32;
        if r > 0 {
            if last_h0_t >= st0 && last_h0_t <= en0 && last_h0_t + 1 >= st0 && last_h0_t < en0 {
                let d0 = vat(last_h0_t);
                let d1 = uat(last_h0_t + 1);
                let right = d0 <= d1;
                h0 += if right { d1 } else { d0 };
                last_h0_t += right as i32;
            } else if last_h0_t >= st0 && last_h0_t <= en0 {
                h0 += vat(last_h0_t);
            } else {
                last_h0_t += 1;
                h0 += uat(last_h0_t);
            }
        } else {
            h0 = vat(0) - qe_scalar;
            last_h0_t = 0;
        }
        if (flags & APPROX_DROP) != 0 {
            if h0 > result.max {
                result.max = h0;
                result.max_score_target_pos = last_h0_t;
                result.max_score_query_pos = r - last_h0_t;
            } else if last_h0_t >= result.max_score_target_pos && (r - last_h0_t) >= result.max_score_query_pos {
                let tl = last_h0_t - result.max_score_target_pos;
                let ql = (r - last_h0_t) - result.max_score_query_pos;
                let l = if tl > ql { tl - ql } else { ql - tl };
                if z_drop >= 0 && (result.max - h0) > z_drop + l * e2i {
                    result.zdropped = 1;
                    break;
                }
            }
        }
        if r == query_len as i32 + target_len as i32 - 2 && en0 == target_len as i32 - 1 {
            result.score = h0;
        }
        last_st = st0;
        last_en = en0;
    }

    if result.score == NEG_INF {
        result.score = result.max;
    }
    if with_cigar {
        traceback_dual_affine(result, query_len, target_len, end_bonus, flags, n_col, 32, p, band_off, band_off_end);
    }
}}

#[cfg(all(test, target_arch = "x86_64"))]
mod reg1_tests {
    use super::*;

    /// The register kernel must reproduce clean_impl exactly on every eligible call.
    #[test]
    fn reg1_matches_clean_impl() {
        if !is_x86_feature_detected!("avx2") { return; }
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = move || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let mat = { let mut m = [0i8; 25]; for i in 0..5 { for j in 0..5 { m[i * 5 + j] = if i == 4 || j == 4 { -1 } else if i == j { 2 } else { -4 }; } } m };
        let mut n_cases = 0;
        for case in 0..20000 {
            let tl = 1 + (rnd() % 400) as usize;
            let mut t: Vec<u8> = (0..tl).map(|_| (rnd() & 3) as u8).collect();
            let div = rnd() % 30; // edits per 100 bases
            let mut q = Vec::new();
            for &c in &t {
                let r = rnd() % 100;
                if r < div / 3 { q.push((rnd() & 3) as u8); }          // substitution
                else if r < div / 2 { }                                 // deletion
                else if r < div { q.push(c); q.push((rnd() & 3) as u8); } // insertion
                else { q.push(c); }
            }
            if q.is_empty() { q.push(0); }
            if case % 7 == 0 { let i = (rnd() as usize) % q.len(); q[i] = 4; }
            if case % 11 == 0 { let i = (rnd() as usize) % t.len(); t[i] = 4; }
            let bw = (rnd() % 32) as i32;
            let shift = if case % 5 == 0 { (rnd() % 9) as i32 - 4 } else { 0 };
            let zdrop = [400, 100, 30, -1][(rnd() % 4) as usize];
            let end_bonus = if case % 3 == 0 { 5 } else { -1 };
            let mut flags = APPROX_MAX;
            if rnd() & 1 == 0 { flags |= RIGHT_ALIGN; }
            if rnd() % 4 == 0 { flags |= SCORE_ONLY; }
            if rnd() & 1 == 0 { flags |= APPROX_DROP; }
            let has_n = q.iter().chain(t.iter()).any(|&c| c >= 4);
            let (cig, right) = (flags & SCORE_ONLY == 0, flags & RIGHT_ALIGN != 0);
            let mut a = DpResult::default();
            let mut b = DpResult::default();
            unsafe {
                macro_rules! both {
                    ($c:literal, $r:literal, $n:literal) => {{
                        clean_impl::<true, $c, $r, $n>(&q, &t, 5, &mat, 4, 2, 24, 1, bw, zdrop, end_bonus, flags, 0, shift, &mut a);
                        reg1_impl::<$c, $r, $n>(&q, &t, 5, &mat, 4, 2, 24, 1, bw, zdrop, end_bonus, flags, shift, &mut b);
                    }};
                }
                match (cig, right, has_n) {
                    (false, false, false) => both!(false, false, false),
                    (false, false, true) => both!(false, false, true),
                    (false, true, false) => both!(false, true, false),
                    (false, true, true) => both!(false, true, true),
                    (true, false, false) => both!(true, false, false),
                    (true, false, true) => both!(true, false, true),
                    (true, true, false) => both!(true, true, false),
                    (true, true, true) => both!(true, true, true),
                }
            }
            let key = |r: &DpResult| (r.max, r.max_score_query_pos, r.max_score_target_pos, r.max_query_end_score,
                r.max_query_end_target_pos, r.max_target_end_score, r.max_target_end_query_pos, r.score, r.reach_end,
                r.zdropped, r.cigar.clone());
            assert_eq!(key(&a), key(&b), "case {case}: qlen {} tlen {} bw {bw} shift {shift} zdrop {zdrop} flags {flags:#x}", q.len(), t.len());
            n_cases += 1;
        }
        assert_eq!(n_cases, 20000);
    }
}
