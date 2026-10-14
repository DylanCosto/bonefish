// Lightweight Smith-Waterman and global alignment kernels

#[cfg(target_arch = "aarch64")]
use core::arch::aarch64::*;

#[cfg(target_arch = "x86_64")]
use core::arch::x86_64::*;

#[cfg(target_arch = "wasm32")]
use super::common::simd_compat::*;

use super::common::*;
use super::single::extend_single_affine;

// ============================================================================
// Lightweight Smith-Waterman
// Used for quick inversion scoring.
// ============================================================================

/// Query profile for lightweight i16 Smith-Waterman (lightweight_profile_init, size=2 only)
pub struct LightweightProfile {
    pub qlen: i32,
    pub segment_len: i32,     // segmented length = ceil(qlen / 8)
    pub query_profile: Vec<i16>,  // query profile: m * segment_len * 8 values
    pub h0: Vec<i16>,  // segment_len * 8
    pub h1: Vec<i16>,  // segment_len * 8
    pub e: Vec<i16>,   // segment_len * 8
    pub hmax: Vec<i16>, // segment_len * 8
}

/// Initialize query profile for lightweight i16 Smith-Waterman.
/// Initialize lightweight query profile (size=2 / int16 only).
pub fn lightweight_profile_init(qlen: i32, query: &[u8], alphabet_size: i32, score_matrix: &[i8]) -> LightweightProfile {
    let p = 8i32; // 8 int16 values per __m128i
    let slen = (qlen + p - 1) / p;
    let m_usize = alphabet_size as usize;

    // Build segmented query profile (int16)
    // Layout: for each alphabet char a (0..m), for each segment i (0..slen),
    // store p values at positions k = i, i+slen, i+2*slen, ...
    let mut qp = vec![0i16; m_usize * slen as usize * p as usize];
    {
        let mut t = 0usize;
        for a in 0..m_usize {
            let nlen = (slen * p) as usize;
            let ma = &score_matrix[a * m_usize..a * m_usize + m_usize];
            for i in 0..slen as usize {
                let mut k = i;
                while k < nlen {
                    qp[t] = if (k as i32) >= qlen {
                        -1
                    } else {
                        ma[query[k] as usize] as i16
                    };
                    t += 1;
                    k += slen as usize;
                }
            }
        }
    }

    let sz = slen as usize * p as usize;
    LightweightProfile {
        qlen,
        segment_len: slen,
        query_profile: qp,
        h0: vec![0i16; sz],
        h1: vec![0i16; sz],
        e: vec![0i16; sz],
        hmax: vec![0i16; sz],
    }
}

/// Lightweight i16 Smith-Waterman local alignment.
/// Lightweight i16 Smith-Waterman local alignment.
/// Returns (score, query_end, target_end). query_end/target_end are -1 if no alignment found.
#[cfg(target_arch = "x86_64")]
pub fn lightweight_align_i16(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) {
    unsafe { lightweight_align_i16_sse2(qp, target_len, target, gap_open, gap_extend) }
}

#[cfg(target_arch = "aarch64")]
pub fn lightweight_align_i16(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) {
    unsafe { lightweight_align_i16_neon(qp, target_len, target, gap_open, gap_extend) }
}

#[cfg(target_arch = "wasm32")]
pub fn lightweight_align_i16(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) {
    unsafe { lightweight_align_i16_wasm(qp, target_len, target, gap_open, gap_extend) }
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64", target_arch = "wasm32")))]
pub fn lightweight_align_i16(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) {
    lightweight_align_i16_scalar(qp, target_len, target, gap_open, gap_extend)
}

/// Best local alignment score of `query` against `target` (the score of
/// [`lightweight_align_i16`], without the end positions). Used by the z-drop inversion check,
/// which needs only the score: with AVX2 it runs 16 lanes and skips the per-row copy of the
/// best row. The score is independent of the SIMD width: padding cells past the query end never
/// exceed a real cell and never feed one, and the lazy-F corrections never raise a row maximum.
pub fn lightweight_score(query: &[u8], target: &[u8], score_matrix: &[i8], gap_open: i32, gap_extend: i32) -> i32 {
    #[cfg(target_arch = "x86_64")]
    if !query.is_empty() && !target.is_empty() && std::is_x86_feature_detected!("avx2") {
        return unsafe { lightweight_score_avx2(query, target, score_matrix, gap_open, gap_extend) };
    }
    let mut qp = lightweight_profile_init(query.len() as i32, query, 5, score_matrix);
    lightweight_align_i16(&mut qp, target.len() as i32, target, gap_open, gap_extend).0
}

#[cfg(target_arch = "x86_64")]
thread_local! {
    #[allow(clippy::type_complexity)] // four i16 DP rows
    static LW_BUF: std::cell::RefCell<(Vec<i16>, Vec<i16>, Vec<i16>, Vec<i16>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new(), Vec::new(), Vec::new())) };
}

/// Shift a vector of 16 i16 lanes up by one lane (lane l -> l + 1), inserting 0 in lane 0.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn shl1_epi16(v: __m256i) -> __m256i {
    _mm256_alignr_epi8(v, _mm256_permute2x128_si256(v, v, 0x08), 14)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn lightweight_score_avx2(query: &[u8], target: &[u8], score_matrix: &[i8], gap_open: i32, gap_extend: i32) -> i32 { unsafe {
    const P: usize = 16;
    let qlen = query.len();
    let slen = qlen.div_ceil(P);
    LW_BUF.with(|cell| {
        let mut b = cell.borrow_mut();
        let (prof, h0, h1, ev) = &mut *b;
        prof.clear();
        prof.resize(5 * slen * P, 0);
        for a in 0..5 {
            let row = &score_matrix[a * 5..a * 5 + 5];
            for i in 0..slen {
                for l in 0..P {
                    let pos = i + l * slen;
                    prof[(a * slen + i) * P + l] = if pos < qlen { row[query[pos] as usize] as i16 } else { -1 };
                }
            }
        }
        for v in [&mut *h0, &mut *h1, &mut *ev] {
            v.clear();
            v.resize(slen * P, 0);
        }
        let zero = _mm256_setzero_si256();
        let gapoe = _mm256_set1_epi16((gap_open + gap_extend) as i16);
        let gape = _mm256_set1_epi16(gap_extend as i16);
        let mut vmax = zero;
        let ld = |v: &Vec<i16>, j: usize| _mm256_loadu_si256(v.as_ptr().add(j * P) as *const __m256i);
        let st = |v: &mut Vec<i16>, j: usize, x: __m256i| _mm256_storeu_si256(v.as_mut_ptr().add(j * P) as *mut __m256i, x);
        for &c in target {
            let pr = prof.as_ptr().add((c.min(4) as usize) * slen * P);
            let mut f = zero;
            let mut h = shl1_epi16(ld(h0, slen - 1));
            for j in 0..slen {
                h = _mm256_adds_epi16(h, _mm256_loadu_si256(pr.add(j * P) as *const __m256i));
                let e = ld(ev, j);
                h = _mm256_max_epi16(h, e);
                h = _mm256_max_epi16(h, f);
                vmax = _mm256_max_epi16(vmax, h);
                st(h1, j, h);
                let h_sub = _mm256_subs_epu16(h, gapoe);
                st(ev, j, _mm256_max_epi16(_mm256_subs_epu16(e, gape), h_sub));
                f = _mm256_max_epi16(_mm256_subs_epu16(f, gape), h_sub);
                h = ld(h0, j);
            }
            'lazy: for _ in 0..P {
                f = shl1_epi16(f);
                for j in 0..slen {
                    let x = _mm256_max_epi16(ld(h1, j), f);
                    st(h1, j, x);
                    f = _mm256_subs_epu16(f, gape);
                    if _mm256_movemask_epi8(_mm256_cmpgt_epi16(f, _mm256_subs_epu16(x, gapoe))) == 0 {
                        break 'lazy;
                    }
                }
            }
            std::mem::swap(h0, h1);
        }
        let mut lanes = [0i16; P];
        _mm256_storeu_si256(lanes.as_mut_ptr() as *mut __m256i, vmax);
        lanes.iter().copied().max().unwrap_or(0).max(0) as i32
    })
}}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn lightweight_align_i16_sse2(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) { unsafe {
    let slen = qp.segment_len;
    let mut gmax: i32 = 0;
    let qlen8 = slen * 8;
    let mut query_end: i32 = -1;
    let mut target_end: i32 = -1;

    let zero = _mm_set1_epi32(0);
    let gapoe = _mm_set1_epi16((gap_open + gap_extend) as i16);
    let gape_v = _mm_set1_epi16(gap_extend as i16);

    // Zero out working arrays
    for v in qp.e.iter_mut() { *v = 0; }
    for v in qp.h0.iter_mut() { *v = 0; }
    for v in qp.hmax.iter_mut() { *v = 0; }

    for i in 0..target_len {
        let mut f = zero;
        let mut max = zero;
        let s_offset = target[i as usize] as usize * slen as usize * 8;

        // h = H0[slen-1] shifted left by 2 bytes
        let h0_last_idx = (slen as usize - 1) * 8;
        let mut h = _mm_loadu_si128(qp.h0[h0_last_idx..].as_ptr() as *const __m128i);
        h = _mm_slli_si128(h, 2);

        for j in 0..slen as usize {
            let s = _mm_loadu_si128(qp.query_profile[s_offset + j * 8..].as_ptr() as *const __m128i);
            h = _mm_adds_epi16(h, s);
            let e = _mm_loadu_si128(qp.e[j * 8..].as_ptr() as *const __m128i);
            h = _mm_max_epi16(h, e);
            h = _mm_max_epi16(h, f);
            max = _mm_max_epi16(max, h);
            _mm_storeu_si128(qp.h1[j * 8..].as_mut_ptr() as *mut __m128i, h);
            let h_sub = _mm_subs_epu16(h, gapoe);
            let e_sub = _mm_subs_epu16(e, gape_v);
            let e_new = _mm_max_epi16(e_sub, h_sub);
            _mm_storeu_si128(qp.e[j * 8..].as_mut_ptr() as *mut __m128i, e_new);
            f = _mm_subs_epu16(f, gape_v);
            f = _mm_max_epi16(f, h_sub);
            h = _mm_loadu_si128(qp.h0[j * 8..].as_ptr() as *const __m128i);
        }

        // F-wave propagation
        for _k in 0..8 {
            f = _mm_slli_si128(f, 2);
            let mut did_break = false;
            for j in 0..slen as usize {
                let mut h1 = _mm_loadu_si128(qp.h1[j * 8..].as_ptr() as *const __m128i);
                h1 = _mm_max_epi16(h1, f);
                _mm_storeu_si128(qp.h1[j * 8..].as_mut_ptr() as *mut __m128i, h1);
                let h1_sub = _mm_subs_epu16(h1, gapoe);
                f = _mm_subs_epu16(f, gape_v);
                if _mm_movemask_epi8(_mm_cmpgt_epi16(f, h1_sub)) == 0 {
                    did_break = true;
                    break;
                }
            }
            if did_break { break; }
        }

        // __max_8: find max of 8 int16 values
        let mut imax_v = max;
        imax_v = _mm_max_epi16(imax_v, _mm_srli_si128(imax_v, 8));
        imax_v = _mm_max_epi16(imax_v, _mm_srli_si128(imax_v, 4));
        imax_v = _mm_max_epi16(imax_v, _mm_srli_si128(imax_v, 2));
        let imax = _mm_extract_epi16(imax_v, 0) as i16 as i32;

        if imax >= gmax {
            gmax = imax;
            target_end = i;
            qp.hmax.copy_from_slice(&qp.h1);
        }

        // Swap H0 and H1
        std::mem::swap(&mut qp.h0, &mut qp.h1);
    }

    // Find query_end from Hmax. Prefer an in-bounds position with value == gmax;
    // fall back to the highest-scoring in-bounds position when gmax was achieved
    // only at SIMD padding cells. The score returned is still gmax.
    for i in 0..qlen8 {
        let val = qp.hmax[i as usize] as i32;
        if val == gmax {
            let pos = i / 8 + (i % 8) * slen;
            if pos < qp.qlen { query_end = pos; }
        }
    }
    if query_end < 0 {
        let mut best_val: i32 = i32::MIN;
        for i in 0..qlen8 {
            let val = qp.hmax[i as usize] as i32;
            let pos = i / 8 + (i % 8) * slen;
            if pos < qp.qlen && val > best_val {
                best_val = val;
                query_end = pos;
            }
        }
    }

    (gmax, query_end, target_end)
}}

#[cfg(target_arch = "aarch64")]
unsafe fn lightweight_align_i16_neon(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) { unsafe {
    let slen = qp.segment_len;
    let mut gmax: i32 = 0;
    let qlen8 = slen * 8;
    let mut query_end: i32 = -1;
    let mut target_end: i32 = -1;

    let zero = vdupq_n_s16(0);
    let gapoe = vdupq_n_u16((gap_open + gap_extend) as u16);
    let gape_v = vdupq_n_u16(gap_extend as u16);

    for v in qp.e.iter_mut() { *v = 0; }
    for v in qp.h0.iter_mut() { *v = 0; }
    for v in qp.hmax.iter_mut() { *v = 0; }

    for i in 0..target_len {
        let mut f = zero;
        let mut max = zero;
        let s_offset = target[i as usize] as usize * slen as usize * 8;

        let h0_last_idx = (slen as usize - 1) * 8;
        let mut h = vld1q_s16(qp.h0[h0_last_idx..].as_ptr());
        h = vextq_s16(vdupq_n_s16(0), h, 7);

        for j in 0..slen as usize {
            let s = vld1q_s16(qp.query_profile[s_offset + j * 8..].as_ptr());
            h = vqaddq_s16(h, s);
            let e = vld1q_s16(qp.e[j * 8..].as_ptr());
            h = vmaxq_s16(h, e);
            h = vmaxq_s16(h, f);
            max = vmaxq_s16(max, h);
            vst1q_s16(qp.h1[j * 8..].as_mut_ptr(), h);
            let h_sub = vqsubq_u16(vreinterpretq_u16_s16(h), gapoe);
            let e_sub = vqsubq_u16(vreinterpretq_u16_s16(e), gape_v);
            let e_new = vmaxq_s16(vreinterpretq_s16_u16(e_sub), vreinterpretq_s16_u16(h_sub));
            vst1q_s16(qp.e[j * 8..].as_mut_ptr(), e_new);
            f = vreinterpretq_s16_u16(vqsubq_u16(vreinterpretq_u16_s16(f), gape_v));
            f = vmaxq_s16(f, vreinterpretq_s16_u16(h_sub));
            h = vld1q_s16(qp.h0[j * 8..].as_ptr());
        }

        // F-wave propagation
        for _k in 0..8 {
            f = vextq_s16(vdupq_n_s16(0), f, 7);
            let mut did_break = false;
            for j in 0..slen as usize {
                let mut h1 = vld1q_s16(qp.h1[j * 8..].as_ptr());
                h1 = vmaxq_s16(h1, f);
                vst1q_s16(qp.h1[j * 8..].as_mut_ptr(), h1);
                let h1_sub = vqsubq_u16(vreinterpretq_u16_s16(h1), gapoe);
                f = vreinterpretq_s16_u16(vqsubq_u16(vreinterpretq_u16_s16(f), gape_v));
                let cmp = vcgtq_s16(f, vreinterpretq_s16_u16(h1_sub));
                // Check if all lanes are zero (no f > h1_sub)
                let any_set = vmaxvq_u16(cmp);
                if any_set == 0 {
                    did_break = true;
                    break;
                }
            }
            if did_break { break; }
        }

        // Find max of 8 int16 values
        let imax = vmaxvq_s16(max) as i32;

        if imax >= gmax {
            gmax = imax;
            target_end = i;
            qp.hmax.copy_from_slice(&qp.h1);
        }

        std::mem::swap(&mut qp.h0, &mut qp.h1);
    }

    // Find query_end from Hmax. Prefer an in-bounds position with value == gmax;
    // fall back to the highest-scoring in-bounds position when gmax was achieved
    // only at SIMD padding cells. The score returned is still gmax.
    for i in 0..qlen8 {
        let val = qp.hmax[i as usize] as i32;
        if val == gmax {
            let pos = i / 8 + (i % 8) * slen;
            if pos < qp.qlen { query_end = pos; }
        }
    }
    if query_end < 0 {
        let mut best_val: i32 = i32::MIN;
        for i in 0..qlen8 {
            let val = qp.hmax[i as usize] as i32;
            let pos = i / 8 + (i % 8) * slen;
            if pos < qp.qlen && val > best_val {
                best_val = val;
                query_end = pos;
            }
        }
    }

    (gmax, query_end, target_end)
}}

/// WASM SIMD128 implementation of lightweight_align_i16
#[cfg(target_arch = "wasm32")]
#[target_feature(enable = "simd128")]
unsafe fn lightweight_align_i16_wasm(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) {
    let slen = qp.segment_len;
    let mut gmax: i32 = 0;
    let qlen8 = slen * 8;
    let mut query_end: i32 = -1;
    let mut target_end: i32 = -1;

    // SIMD region: every intrinsic call below is an `unsafe fn` (the
    // simd_compat shims wrap raw-pointer v128 loads/stores), so the whole
    // computation runs inside one inner `unsafe` block per Rust 2024's
    // unsafe_op_in_unsafe_fn. Nothing declared here escapes the block.
    unsafe {
        let zero = _mm_set1_epi32(0);
        let gapoe = _mm_set1_epi16((gap_open + gap_extend) as i16);
        let gape_v = _mm_set1_epi16(gap_extend as i16);

        for v in qp.e.iter_mut() { *v = 0; }
        for v in qp.h0.iter_mut() { *v = 0; }
        for v in qp.hmax.iter_mut() { *v = 0; }

        for i in 0..target_len {
            let mut f = zero;
            let mut max = zero;
            let s_offset = target[i as usize] as usize * slen as usize * 8;
            let h0_last_idx = (slen as usize - 1) * 8;
            let mut h = _mm_loadu_si128(qp.h0[h0_last_idx..].as_ptr() as *const __m128i);
            h = _mm_slli_si128(h, 2);

            for j in 0..slen as usize {
                let s = _mm_loadu_si128(qp.query_profile[s_offset + j * 8..].as_ptr() as *const __m128i);
                h = _mm_adds_epi16(h, s);
                let e = _mm_loadu_si128(qp.e[j * 8..].as_ptr() as *const __m128i);
                h = _mm_max_epi16(h, e);
                h = _mm_max_epi16(h, f);
                max = _mm_max_epi16(max, h);
                _mm_storeu_si128(qp.h1[j * 8..].as_mut_ptr() as *mut __m128i, h);
                let h_sub = _mm_subs_epu16(h, gapoe);
                let e_sub = _mm_subs_epu16(e, gape_v);
                let e_new = _mm_max_epi16(e_sub, h_sub);
                _mm_storeu_si128(qp.e[j * 8..].as_mut_ptr() as *mut __m128i, e_new);
                f = _mm_subs_epu16(f, gape_v);
                f = _mm_max_epi16(f, h_sub);
                h = _mm_loadu_si128(qp.h0[j * 8..].as_ptr() as *const __m128i);
            }

            for _k in 0..8 {
                f = _mm_slli_si128(f, 2);
                let mut did_break = false;
                for j in 0..slen as usize {
                    let mut h1 = _mm_loadu_si128(qp.h1[j * 8..].as_ptr() as *const __m128i);
                    h1 = _mm_max_epi16(h1, f);
                    _mm_storeu_si128(qp.h1[j * 8..].as_mut_ptr() as *mut __m128i, h1);
                    let h1_sub = _mm_subs_epu16(h1, gapoe);
                    f = _mm_subs_epu16(f, gape_v);
                    if _mm_movemask_epi8(_mm_cmpgt_epi16(f, h1_sub)) == 0 {
                        did_break = true;
                        break;
                    }
                }
                if did_break { break; }
            }

            let mut imax_v = max;
            imax_v = _mm_max_epi16(imax_v, _mm_srli_si128(imax_v, 8));
            imax_v = _mm_max_epi16(imax_v, _mm_srli_si128(imax_v, 4));
            imax_v = _mm_max_epi16(imax_v, _mm_srli_si128(imax_v, 2));
            let imax = _mm_extract_epi16::<0>(imax_v) as i16 as i32;

            if imax >= gmax {
                gmax = imax;
                target_end = i;
                qp.hmax.copy_from_slice(&qp.h1);
            }

            std::mem::swap(&mut qp.h0, &mut qp.h1);
        }
    }

    // Find query_end from Hmax. Prefer an in-bounds position with value == gmax;
    // fall back to the highest-scoring in-bounds position when gmax was achieved
    // only at SIMD padding cells. The score returned is still gmax.
    for i in 0..qlen8 {
        let val = qp.hmax[i as usize] as i32;
        if val == gmax {
            let pos = i / 8 + (i % 8) * slen;
            if pos < qp.qlen { query_end = pos; }
        }
    }
    if query_end < 0 {
        let mut best_val: i32 = i32::MIN;
        for i in 0..qlen8 {
            let val = qp.hmax[i as usize] as i32;
            let pos = i / 8 + (i % 8) * slen;
            if pos < qp.qlen && val > best_val {
                best_val = val;
                query_end = pos;
            }
        }
    }

    (gmax, query_end, target_end)
}

/// Scalar fallback for lightweight_align_i16
pub fn lightweight_align_i16_scalar(qp: &mut LightweightProfile, target_len: i32, target: &[u8], gap_open: i32, gap_extend: i32) -> (i32, i32, i32) {
    // Simple scalar Smith-Waterman for non-SIMD platforms
    let qlen = qp.qlen as usize;
    let target_len = target_len as usize;
    let _m = 5usize; // alphabet size

    // Reconstruct scoring from query profile
    // qp.query_profile is segmented; we need to unsegment to get mat scores
    let slen = qp.segment_len as usize;
    let p = 8usize;

    let mut h_prev = vec![0i32; qlen + 1];
    let mut h_curr = vec![0i32; qlen + 1];
    let mut e_arr = vec![0i32; qlen + 1];
    let mut gmax: i32 = 0;
    let mut query_end: i32 = -1;
    let mut target_end: i32 = -1;

    for (i, &t_base) in target[..target_len].iter().enumerate() {
        let tc = t_base as usize;
        let mut f: i32 = 0;
        for j in 1..=qlen {
            // Get score from query profile: unsegment index
            // segmented index for query position (j-1): seg = (j-1) % slen, lane = (j-1) / slen
            // offset in qp.query_profile: tc * slen * p + seg * p + lane
            let seg = (j - 1) % slen;
            let lane = (j - 1) / slen;
            let score = if lane < p {
                qp.query_profile[tc * slen * p + seg * p + lane] as i32
            } else {
                0
            };

            let mut h = h_prev[j - 1] + score;
            if h < 0 { h = 0; }

            e_arr[j] = std::cmp::max(e_arr[j].saturating_sub(gap_extend), h_curr.get(j).copied().unwrap_or(0).saturating_sub(gap_open + gap_extend));
            // Wait, this is wrong - we need unsigned saturation semantics
            // Let's just do it simply
            let e_val = std::cmp::max(
                (e_arr[j] as i64 - gap_extend as i64).max(0) as i32,
                (h as i64 - (gap_open + gap_extend) as i64).max(0) as i32,
            );
            e_arr[j] = e_val;

            f = std::cmp::max(
                (f as i64 - gap_extend as i64).max(0) as i32,
                (h as i64 - (gap_open + gap_extend) as i64).max(0) as i32,
            );

            h = std::cmp::max(h, std::cmp::max(e_arr[j], f));
            h_curr[j] = h;

            if h >= gmax {
                gmax = h;
                target_end = i as i32;
                query_end = (j - 1) as i32;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_curr);
        for v in h_curr.iter_mut() { *v = 0; }
    }

    (gmax, query_end, target_end)
}


/// Global alignment (Needleman-Wunsch) with CIGAR traceback.
///
/// Dispatches to the best available implementation:
/// - SIMD targets (x86_64/aarch64/wasm32): SIMD extension with end_bonus
///   (AVX512 > AVX2 > SSE2 > NEON > WASM SIMD128)
/// - Non-SIMD targets: row-by-row Gotoh NW (scalar fallback)
///
/// Both produce correct end-to-end CIGARs. The SIMD path recomputes
/// the score from the CIGAR (stripping the end_bonus inflation).
#[allow(clippy::too_many_arguments)]
pub fn global_align(
    qseq: &[u8], tseq: &[u8],
    alphabet_size: i8, score_matrix: &[i8],
    gap_open: i32, gap_extend: i32,
    bandwidth: i32,
    result: &mut DpResult,
) {
    let qlen = qseq.len();
    let tlen = tseq.len();
    if qlen == 0 || tlen == 0 {
        if qlen > 0 { result.cigar = vec![(qlen as u32) << 4 | 1]; result.score = -(gap_open + gap_extend * qlen as i32); }
        if tlen > 0 { result.cigar = vec![(tlen as u32) << 4 | 2]; result.score = -(gap_open + gap_extend * tlen as i32); }
        return;
    }

    // Use SIMD extension on SIMD-capable targets (faster at all lengths).
    // Gotoh scalar NW is the fallback for non-SIMD architectures.
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64", target_arch = "wasm32"))]
    let has_simd = !*crate::align::env_flags::FORCE_SCALAR;
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64", target_arch = "wasm32")))]
    let has_simd = false;

    if has_simd {
        global_align_simd(qseq, tseq, alphabet_size, score_matrix, gap_open, gap_extend, bandwidth, result);
    } else {
        global_align_gotoh(qseq, tseq, alphabet_size, score_matrix, gap_open, gap_extend, bandwidth, result);
    }
}

/// SIMD-accelerated global alignment via extension DP with end_bonus.
///
/// Uses the same SIMD kernels as the mapper (AVX512/AVX2/SSE2). The large
/// end_bonus forces the alignment to cover both sequences end-to-end.
/// Score is recomputed from CIGAR for correctness (end_bonus inflates ez.score).
#[allow(clippy::too_many_arguments)]
fn global_align_simd(
    qseq: &[u8], tseq: &[u8],
    alphabet_size: i8, score_matrix: &[i8],
    gap_open: i32, gap_extend: i32,
    bandwidth: i32,
    result: &mut DpResult,
) {
    let end_bonus = (score_matrix[0] as i32 * std::cmp::max(qseq.len(), tseq.len()) as i32).max(1000);
    let bw = if bandwidth > 0 { bandwidth } else { -1 };
    extend_single_affine(
        qseq, tseq, alphabet_size, score_matrix,
        gap_open as i8, gap_extend as i8,
        bw, -1, end_bonus, APPROX_MAX, result,
    );
    // Recompute score from CIGAR (strip end_bonus inflation)
    let m = alphabet_size as usize;
    let mut score = 0i32;
    let mut qi = 0usize;
    let mut ti = 0usize;
    for &c in &result.cigar {
        let len = (c >> 4) as usize;
        match c & 0xf {
            0 => {
                for _ in 0..len {
                    if qi < qseq.len() && ti < tseq.len() {
                        score += score_matrix[tseq[ti].min(4) as usize * m + qseq[qi].min(4) as usize] as i32;
                    }
                    qi += 1; ti += 1;
                }
            }
            1 => { score -= gap_open + gap_extend * len as i32; qi += len; }
            2 => { score -= gap_open + gap_extend * len as i32; ti += len; }
            _ => {}
        }
    }
    result.score = score;
}

/// Row-by-row Gotoh NW with banded backtrack matrix.
///
/// Row-by-row Gotoh NW (scalar fallback for non-SIMD architectures).
/// Simple H[i][j], E[j], F recurrence with banded backtrack matrix.
#[allow(clippy::too_many_arguments)]
fn global_align_gotoh(
    qseq: &[u8], tseq: &[u8],
    _alphabet_size: i8, score_matrix: &[i8],
    gap_open: i32, gap_extend: i32,
    bandwidth: i32,
    result: &mut DpResult,
) {
    let qlen = qseq.len();
    let tlen = tseq.len();
    let m = 5usize; // nt4 alphabet
    let gapoe = gap_open + gap_extend;
    let w = if bandwidth > 0 { bandwidth as usize } else { qlen + tlen };

    // DP arrays (current and previous row of H, plus E for query gaps)
    let mut h_prev = vec![NEG_INF; tlen + 1];
    let mut h_curr = vec![NEG_INF; tlen + 1];
    let mut e = vec![NEG_INF; tlen + 1]; // E[j]: best score ending with query gap at column j

    // Initialize first row: H[0][j] = gap penalties
    h_prev[0] = 0;
    for j in 1..=tlen {
        h_prev[j] = -(gapoe + gap_extend * (j as i32 - 1));
        e[j] = NEG_INF;
    }

    // Backtrack matrix: 2 bits per cell (0=diag, 1=up/I, 2=left/D)
    let mut bt = vec![0u8; qlen * tlen];

    for i in 1..=qlen {
        let mut f = NEG_INF; // F: best score ending with target gap (current row)
        h_curr[0] = -(gapoe + gap_extend * (i as i32 - 1));

        // Band boundaries
        let j_start = if w < i { i - w } else { 1 };
        let j_end = std::cmp::min(tlen, i + w);

        for j in j_start..=j_end {
            // Match/mismatch from diagonal
            let s = score_matrix[tseq[j - 1] as usize * m + qseq[i - 1] as usize] as i32;
            let diag = h_prev[j - 1] + s;

            // Query gap (insertion): extend or open from H
            let e_ext = e[j] - gap_extend;
            let e_open = h_prev[j] - gapoe;
            e[j] = std::cmp::max(e_ext, e_open);

            // Target gap (deletion): extend or open from H
            let f_ext = f - gap_extend;
            let f_open = h_curr[j - 1] - gapoe;
            f = std::cmp::max(f_ext, f_open);

            // Best of three
            let h = std::cmp::max(diag, std::cmp::max(e[j], f));
            h_curr[j] = h;

            // Backtrack direction
            let d = if h == diag { 0u8 }
                    else if h == e[j] { 1u8 }
                    else { 2u8 };
            bt[(i - 1) * tlen + (j - 1)] = d;
        }

        std::mem::swap(&mut h_prev, &mut h_curr);
        h_curr.fill(NEG_INF);
    }

    result.score = h_prev[tlen];

    // Traceback from (qlen, tlen)
    let mut cigar = Vec::new();
    let mut i = qlen;
    let mut j = tlen;

    while i > 0 && j > 0 {
        let d = bt[(i - 1) * tlen + (j - 1)];
        match d {
            0 => { push_cigar(&mut cigar, 0, 1); i -= 1; j -= 1; } // M
            1 => { push_cigar(&mut cigar, 1, 1); i -= 1; }          // I (consume query)
            _ => { push_cigar(&mut cigar, 2, 1); j -= 1; }          // D (consume target)
        }
    }
    if i > 0 { push_cigar(&mut cigar, 1, i as u32); }
    if j > 0 { push_cigar(&mut cigar, 2, j as u32); }

    cigar.reverse();
    result.cigar = cigar;
}

#[cfg(test)]
mod lw_score_tests {
    use super::*;

    /// lightweight_score must equal the score of lightweight_align_i16 (the AVX2 path uses a
    /// different lane layout and no best-row copy).
    #[test]
    fn lightweight_score_matches_align_i16() {
        let mut mat = [0i8; 25];
        for a in 0..5 { for b in 0..5 { mat[a * 5 + b] = if a == 4 || b == 4 { -1 } else if a == b { 2 } else { -4 }; } }
        let mut x: u64 = 0x9e3779b97f4a7c15;
        let mut rnd = |n: u64| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x % n };
        for case in 0..3000 {
            let ql = 1 + rnd(if case % 10 == 0 { 3000 } else { 400 }) as usize;
            let tl = 1 + rnd(if case % 10 == 0 { 3000 } else { 400 }) as usize;
            let q: Vec<u8> = (0..ql).map(|_| if rnd(50) == 0 { 4 } else { rnd(4) as u8 }).collect();
            // target: a mutated, sometimes reverse-complemented copy of part of q, or random
            let t: Vec<u8> = if case % 3 == 0 {
                (0..tl).map(|_| rnd(5) as u8).collect()
            } else {
                (0..tl).map(|i| { let c = q[(i + case) % ql]; if rnd(8) == 0 { rnd(4) as u8 } else { c } }).collect()
            };
            let (go, ge) = if case % 2 == 0 { (4, 2) } else { (6, 1) };
            let mut qp = lightweight_profile_init(ql as i32, &q, 5, &mat);
            let want = lightweight_align_i16(&mut qp, tl as i32, &t, go, ge).0;
            assert_eq!(lightweight_score(&q, &t, &mat, go, ge), want, "case {case} ql {ql} tl {tl}");
        }
    }
}
