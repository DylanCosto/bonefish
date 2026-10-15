//! 16-lane port of the clean dual-affine kernel (`dual_clean.rs`) for targets without AVX2:
//! NEON on aarch64 and SSE4.1 on x86_64. Same recurrences, band (static, shifted or adaptive),
//! boundary values and post-band resets as `extend_dual_affine_avx2_clean`, with the same
//! traceback layout, so its results are identical; only the vector width differs. Without it,
//! those targets fell back to the ksw2-style kernels (rounded bands, stale scores) and ignored
//! the adaptive and shifted bands, so their output differed from x86 AVX2.
//!
//! The kernel is written once over [`V16`]; the SSE4.1 instance lets x86 tests check the shared
//! code against the AVX2 kernel, which leaves only the one-to-one NEON intrinsic mapping
//! platform specific.

use super::common::*;

/// The byte-vector operations the kernel needs (16 lanes of i8, masks all-ones or zero).
pub(super) trait V16: Copy {
    unsafe fn load(p: *const u8) -> Self;
    unsafe fn store(self, p: *mut u8);
    unsafe fn splat(x: i8) -> Self;
    /// Wrapping lane-wise addition.
    unsafe fn add(self, o: Self) -> Self;
    /// Wrapping lane-wise subtraction.
    unsafe fn sub(self, o: Self) -> Self;
    unsafe fn max(self, o: Self) -> Self;
    unsafe fn min(self, o: Self) -> Self;
    /// Unsigned maximum.
    unsafe fn max_u(self, o: Self) -> Self;
    unsafe fn eq(self, o: Self) -> Self;
    /// Signed `self > o`.
    unsafe fn gt(self, o: Self) -> Self;
    unsafe fn and(self, o: Self) -> Self;
    unsafe fn or(self, o: Self) -> Self;
    /// `!self & o`.
    unsafe fn andnot(self, o: Self) -> Self;
    /// Lanes of `o` where `m` is set, else `self`.
    unsafe fn blend(self, o: Self, m: Self) -> Self;
    /// Lanes moved up one: `[prev[15], self[0], ..., self[14]]`.
    unsafe fn shl1_prev(self, prev: Self) -> Self;
}

#[cfg(target_arch = "x86_64")]
mod sse {
    use core::arch::x86_64::*;

    #[derive(Clone, Copy)]
    pub(in super::super) struct Sse(__m128i);

    impl super::V16 for Sse {
        #[inline(always)] unsafe fn load(p: *const u8) -> Self { unsafe { Sse(_mm_loadu_si128(p as *const __m128i)) } }
        #[inline(always)] unsafe fn store(self, p: *mut u8) { unsafe { _mm_storeu_si128(p as *mut __m128i, self.0) } }
        #[inline(always)] unsafe fn splat(x: i8) -> Self { unsafe { Sse(_mm_set1_epi8(x)) } }
        #[inline(always)] unsafe fn add(self, o: Self) -> Self { unsafe { Sse(_mm_add_epi8(self.0, o.0)) } }
        #[inline(always)] unsafe fn sub(self, o: Self) -> Self { unsafe { Sse(_mm_sub_epi8(self.0, o.0)) } }
        #[inline(always)] unsafe fn max(self, o: Self) -> Self { unsafe { Sse(_mm_max_epi8(self.0, o.0)) } }
        #[inline(always)] unsafe fn min(self, o: Self) -> Self { unsafe { Sse(_mm_min_epi8(self.0, o.0)) } }
        #[inline(always)] unsafe fn max_u(self, o: Self) -> Self { unsafe { Sse(_mm_max_epu8(self.0, o.0)) } }
        #[inline(always)] unsafe fn eq(self, o: Self) -> Self { unsafe { Sse(_mm_cmpeq_epi8(self.0, o.0)) } }
        #[inline(always)] unsafe fn gt(self, o: Self) -> Self { unsafe { Sse(_mm_cmpgt_epi8(self.0, o.0)) } }
        #[inline(always)] unsafe fn and(self, o: Self) -> Self { unsafe { Sse(_mm_and_si128(self.0, o.0)) } }
        #[inline(always)] unsafe fn or(self, o: Self) -> Self { unsafe { Sse(_mm_or_si128(self.0, o.0)) } }
        #[inline(always)] unsafe fn andnot(self, o: Self) -> Self { unsafe { Sse(_mm_andnot_si128(self.0, o.0)) } }
        #[inline(always)] unsafe fn blend(self, o: Self, m: Self) -> Self { unsafe { Sse(_mm_blendv_epi8(self.0, o.0, m.0)) } }
        #[inline(always)] unsafe fn shl1_prev(self, prev: Self) -> Self { unsafe { Sse(_mm_alignr_epi8(self.0, prev.0, 15)) } }
    }
}

#[cfg(target_arch = "aarch64")]
mod neon {
    use core::arch::aarch64::*;

    #[derive(Clone, Copy)]
    pub(in super::super) struct Neon(int8x16_t);

    impl super::V16 for Neon {
        #[inline(always)] unsafe fn load(p: *const u8) -> Self { unsafe { Neon(vld1q_s8(p as *const i8)) } }
        #[inline(always)] unsafe fn store(self, p: *mut u8) { unsafe { vst1q_s8(p as *mut i8, self.0) } }
        #[inline(always)] unsafe fn splat(x: i8) -> Self { unsafe { Neon(vdupq_n_s8(x)) } }
        #[inline(always)] unsafe fn add(self, o: Self) -> Self { unsafe { Neon(vaddq_s8(self.0, o.0)) } }
        #[inline(always)] unsafe fn sub(self, o: Self) -> Self { unsafe { Neon(vsubq_s8(self.0, o.0)) } }
        #[inline(always)] unsafe fn max(self, o: Self) -> Self { unsafe { Neon(vmaxq_s8(self.0, o.0)) } }
        #[inline(always)] unsafe fn min(self, o: Self) -> Self { unsafe { Neon(vminq_s8(self.0, o.0)) } }
        #[inline(always)] unsafe fn max_u(self, o: Self) -> Self {
            unsafe { Neon(vreinterpretq_s8_u8(vmaxq_u8(vreinterpretq_u8_s8(self.0), vreinterpretq_u8_s8(o.0)))) }
        }
        #[inline(always)] unsafe fn eq(self, o: Self) -> Self { unsafe { Neon(vreinterpretq_s8_u8(vceqq_s8(self.0, o.0))) } }
        #[inline(always)] unsafe fn gt(self, o: Self) -> Self { unsafe { Neon(vreinterpretq_s8_u8(vcgtq_s8(self.0, o.0))) } }
        #[inline(always)] unsafe fn and(self, o: Self) -> Self { unsafe { Neon(vandq_s8(self.0, o.0)) } }
        #[inline(always)] unsafe fn or(self, o: Self) -> Self { unsafe { Neon(vorrq_s8(self.0, o.0)) } }
        // vbicq(a, b) = a & !b
        #[inline(always)] unsafe fn andnot(self, o: Self) -> Self { unsafe { Neon(vbicq_s8(o.0, self.0)) } }
        // vbslq(m, a, b) = m ? a : b (bitwise; masks are whole lanes)
        #[inline(always)] unsafe fn blend(self, o: Self, m: Self) -> Self {
            unsafe { Neon(vbslq_s8(vreinterpretq_u8_s8(m.0), o.0, self.0)) }
        }
        // vextq(a, b, 15) = [a[15], b[0], ..., b[14]]
        #[inline(always)] unsafe fn shl1_prev(self, prev: Self) -> Self { unsafe { Neon(vextq_s8::<15>(prev.0, self.0)) } }
    }
}

/// Dispatch over the flag combinations, as the AVX2 entry point does.
macro_rules! clean16_entry {
    ($v:ty, $qseq:ident, $tseq:ident, $alphabet_size:ident, $score_matrix:ident, $gap_open:ident, $gap_extend:ident,
     $gap_open2:ident, $gap_extend2:ident, $bandwidth:ident, $z_drop:ident, $end_bonus:ident, $flags:ident,
     $adaptive:ident, $band_shift:ident, $result:ident) => {{
        let m1 = $alphabet_size.wrapping_sub(1) as u8;
        let has_n = $qseq.iter().fold(0u8, |m, &c| m.max(c)).max($tseq.iter().fold(0u8, |m, &c| m.max(c))) >= m1;
        let approx_max = ($flags & APPROX_MAX) != 0;
        let with_cigar = ($flags & SCORE_ONLY) == 0;
        let right_align = ($flags & RIGHT_ALIGN) != 0;
        macro_rules! go {
            ($a:literal, $c:literal, $r:literal, $n:literal) => {
                clean16_impl::<$v, $a, $c, $r, $n>($qseq, $tseq, $alphabet_size, $score_matrix, $gap_open, $gap_extend,
                    $gap_open2, $gap_extend2, $bandwidth, $z_drop, $end_bonus, $flags, $adaptive, $band_shift, $result)
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
    }};
}

/// SSE4.1 instance of the 16-lane clean kernel (x86_64 without AVX2; also the x86 test vehicle
/// for the shared code). Same arguments as `extend_dual_affine_avx2_clean`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.1")]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn extend_dual_affine_sse41_clean(
    qseq: &[u8], tseq: &[u8], alphabet_size: i8, score_matrix: &[i8],
    gap_open: i8, gap_extend: i8, gap_open2: i8, gap_extend2: i8,
    bandwidth: i32, z_drop: i32, end_bonus: i32, flags: i32, adaptive: i32, band_shift: i32,
    result: &mut DpResult,
) { unsafe {
    clean16_entry!(sse::Sse, qseq, tseq, alphabet_size, score_matrix, gap_open, gap_extend, gap_open2, gap_extend2,
                   bandwidth, z_drop, end_bonus, flags, adaptive, band_shift, result)
}}

/// NEON instance of the 16-lane clean kernel. Same arguments as `extend_dual_affine_avx2_clean`.
#[cfg(target_arch = "aarch64")]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn extend_dual_affine_neon_clean(
    qseq: &[u8], tseq: &[u8], alphabet_size: i8, score_matrix: &[i8],
    gap_open: i8, gap_extend: i8, gap_open2: i8, gap_extend2: i8,
    bandwidth: i32, z_drop: i32, end_bonus: i32, flags: i32, adaptive: i32, band_shift: i32,
    result: &mut DpResult,
) { unsafe {
    clean16_entry!(neon::Neon, qseq, tseq, alphabet_size, score_matrix, gap_open, gap_extend, gap_open2, gap_extend2,
                   bandwidth, z_drop, end_bonus, flags, adaptive, band_shift, result)
}}

/// `clean_impl` of `dual_clean.rs` with 16-lane vectors. Memory layout, padding and traceback
/// stride are the AVX2 kernel's (multiples of 32), so the traceback is shared unchanged.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
#[allow(clippy::manual_is_multiple_of)] // kernel left as written; `% W` matches the vector width
unsafe fn clean16_impl<V: V16, const APPROX: bool, const CIGAR: bool, const RIGHT: bool, const HAS_N: bool>(
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
    const W: usize = 16;
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

    let zero_ = V::splat(0);
    let ones_ = V::splat(-1);
    let q_ = V::splat(gap_open);
    let q2_ = V::splat(gap_open2);
    let qe_ = V::splat((qi + ei) as i8);
    let qe2_ = V::splat((q2i + e2i) as i8);
    let neg_qe_ = V::splat(neg_qe);
    let neg_q2e2_ = V::splat(neg_q2e2);
    let sc_mch_ = V::splat(score_matrix[0]);
    let sc_mis_ = V::splat(score_matrix[1]);
    let n_sc = score_matrix[(alphabet_size as usize * alphabet_size as usize) - 1];
    let sc_n_ = V::splat(if n_sc == 0 { -gap_extend2 } else { n_sc });
    let m1_ = V::splat(alphabet_size - 1);
    let sc_mch_mis_ = V::splat(score_matrix[0].wrapping_sub(score_matrix[1]));
    let has_n = HAS_N;
    let flag1_ = V::splat(1);
    let flag2_ = V::splat(2);
    let flag3_ = V::splat(3);
    let flag4_ = V::splat(4);
    let flag8_ = V::splat(0x08);
    let flag16_ = V::splat(0x10);
    let flag32_ = V::splat(0x20);
    let flag64_ = V::splat(0x40);

    let wl = if bandwidth < 0 { target_len.max(query_len) as i32 } else { bandwidth };
    let adaptive = if adaptive > 0 { adaptive.min(wl.max(1)) } else { 0 };
    let band_w = if adaptive > 0 { adaptive } else { wl };
    let mut a_st: i32 = 0;
    let mut a_step: i32 = 0;

    // Same layout as the AVX2 kernel (the 16-lane loop needs less padding than it provides).
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

        // The shift only reads lane 15 of the previous vector: a broadcast carries the left
        // boundary into the first one.
        let mut px = V::splat(x1);
        let mut px2 = V::splat(x21);
        let mut pv = V::splat(v1);
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
            let sq = V::load(sf.add(t));
            let qv = V::load(qrr.add(t));
            let mut z = sc_mis_.add(sq.eq(qv).and(sc_mch_mis_));
            if has_n {
                z = z.blend(sc_n_, sq.max_u(qv).eq(m1_));
            }

            let xt = V::load(x.add(t));
            let xt1 = xt.shl1_prev(px);
            px = xt;
            let vt = V::load(v.add(t));
            let vt1 = vt.shl1_prev(pv);
            pv = vt;
            let x2t = V::load(x2.add(t));
            let x2t1 = x2t.shl1_prev(px2);
            px2 = x2t;
            let ut = V::load(u.add(t));
            let a = xt1.add(vt1);
            let b = V::load(y.add(t)).add(ut);
            let a2 = x2t1.add(vt1);
            let b2 = V::load(y2.add(t)).add(ut);

            if !with_cigar {
                z = z.max(a).max(b).max(a2).max(b2).min(sc_mch_);
                z.sub(vt1).store(u.add(t));
                z.sub(ut).store(v.add(t));
                let tmp1 = z.sub(q_);
                let tmp2 = z.sub(q2_);
                a.sub(tmp1).max(zero_).sub(qe_).store(x.add(t));
                b.sub(tmp1).max(zero_).sub(qe_).store(y.add(t));
                a2.sub(tmp2).max(zero_).sub(qe2_).store(x2.add(t));
                b2.sub(tmp2).max(zero_).sub(qe2_).store(y2.add(t));
            } else {
                let mut d;
                if !right_align {
                    d = a.gt(z).and(flag1_);
                    z = z.max(a);
                    d = d.blend(flag2_, b.gt(z));
                    z = z.max(b);
                    d = d.blend(flag3_, a2.gt(z));
                    z = z.max(a2);
                    d = d.blend(flag4_, b2.gt(z));
                    z = z.max(b2);
                } else {
                    d = z.gt(a).andnot(flag1_);
                    z = z.max(a);
                    d = flag2_.blend(d, z.gt(b));
                    z = z.max(b);
                    d = flag3_.blend(d, z.gt(a2));
                    z = z.max(a2);
                    d = flag4_.blend(d, z.gt(b2));
                    z = z.max(b2);
                }
                z = z.min(sc_mch_);
                z.sub(vt1).store(u.add(t));
                z.sub(ut).store(v.add(t));
                let tmp1 = z.sub(q_);
                let tmp2 = z.sub(q2_);
                let an = a.sub(tmp1);
                let bn = b.sub(tmp1);
                let a2n = a2.sub(tmp2);
                let b2n = b2.sub(tmp2);
                // Left alignment extends a gap when the new value is > 0, right alignment when >= 0.
                let (ma, mb, ma2, mb2) = if !right_align {
                    (an.gt(zero_), bn.gt(zero_), a2n.gt(zero_), b2n.gt(zero_))
                } else {
                    (zero_.gt(an).andnot(ones_), zero_.gt(bn).andnot(ones_),
                     zero_.gt(a2n).andnot(ones_), zero_.gt(b2n).andnot(ones_))
                };
                ma.and(an).sub(qe_).store(x.add(t));
                mb.and(bn).sub(qe_).store(y.add(t));
                ma2.and(a2n).sub(qe2_).store(x2.add(t));
                mb2.and(b2n).sub(qe2_).store(y2.add(t));
                d = d.or(ma.and(flag8_));
                d = d.or(mb.and(flag16_));
                d = d.or(ma2.and(flag32_));
                d = d.or(mb2.and(flag64_));
                d.store(prow.add(t));
            }
            t += W;
        }
        // The last vector may have written up to 15 cells right of en0; cells right of the band
        // must keep their initial values (the scalar kernel never touches them).
        let es = en + 1;
        if (en - st + 1) % W != 0 {
            neg_qe_.store(u.add(es));
            neg_qe_.store(v.add(es));
            neg_qe_.store(x.add(es));
            neg_qe_.store(y.add(es));
            neg_q2e2_.store(x2.add(es));
            neg_q2e2_.store(y2.add(es));
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
                let mut tt = st0;
                // Four interleaved lanes, as the AVX2 kernel's packed and unpacked reductions.
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
                a_step = (en0 - st0 + 1 >= adaptive && 2 * max_t > st0 + en0) as i32;
            }
            if r == query_len as i32 + target_len as i32 - 2 && en0 == target_len as i32 - 1 {
                result.score = *h_ptr.add(target_len - 1);
            }
        } else {
            if r > 0 {
                if last_h0_t >= st0 && last_h0_t <= en0 && last_h0_t + 1 >= st0 && last_h0_t < en0 {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// FNV-1a digest of every case's results from the platform's clean kernel (AVX2 on x86_64,
    /// NEON on aarch64): equal on every platform when the kernels agree. Recompute with
    /// `cargo test --release clean16 -- --nocapture` after an intended kernel change.
    const GOLDEN: u64 = 0x4fcd_48c4_d927_0cf3;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
            self.0 = x;
            x.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn below(&mut self, n: u64) -> u64 { self.next() % n }
    }

    struct Case { q: Vec<u8>, t: Vec<u8>, mat: [i8; 25], bw: i32, zdrop: i32, end_bonus: i32, flags: i32, adaptive: i32, shift: i32 }

    /// Read-like pairs (substitutions, short and long indels, some Ns) under every flag
    /// combination, static, shifted and adaptive bands, band widths around the vector sizes.
    fn cases() -> Vec<Case> {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut out = Vec::new();
        for i in 0..3000u64 {
            let len = match i % 10 { 0 => 1 + rng.below(8) as usize, 9 => 300 + rng.below(1200) as usize, _ => 10 + rng.below(300) as usize };
            let n_rate = if i % 7 == 0 { 20 } else { 0 };
            let t: Vec<u8> = (0..len).map(|_| if rng.below(1000) < n_rate { 4 } else { rng.below(4) as u8 }).collect();
            let long = if i % 5 < 2 { 60 } else { 4 };
            let mut q = Vec::with_capacity(len + 64);
            let mut k = 0;
            while k < t.len() {
                let r = rng.below(1000);
                if r < 40 { q.push(rng.below(4) as u8); k += 1; }
                else if r < 55 { k += 1 + rng.below(long) as usize; }
                else if r < 70 { for _ in 0..1 + rng.below(long) { q.push(rng.below(4) as u8); } }
                else { q.push(t[k]); k += 1; }
            }
            if q.is_empty() { q.push(1); }
            let n_sc = if i % 11 == 0 { 0 } else { -1 };
            let mut mat = [0i8; 25];
            for a in 0..5 { for b in 0..5 { mat[a * 5 + b] = if a == 4 || b == 4 { n_sc } else if a == b { 2 } else { -4 }; } }
            let bw = [-1, 0, 3, 15, 16, 17, 31, 32, 33, 47, 64, 100, 200, 751][rng.below(14) as usize];
            let adaptive = if i % 3 == 0 { [8, 16, 33, 64][rng.below(4) as usize] } else { 0 };
            let shift = if adaptive == 0 && i % 4 == 1 { rng.below(41) as i32 - 20 } else { 0 };
            let mut flags = 0;
            if rng.below(2) == 0 { flags |= APPROX_MAX; }
            if rng.below(3) == 0 { flags |= SCORE_ONLY; }
            if rng.below(3) == 0 { flags |= RIGHT_ALIGN; }
            if rng.below(4) == 0 { flags |= EXTENSION_ONLY; }
            if flags & APPROX_MAX != 0 && rng.below(2) == 0 { flags |= APPROX_DROP; }
            let zdrop = [-1, 50, 100, 400][rng.below(4) as usize];
            let end_bonus = if rng.below(4) == 0 { 10 } else { -1 };
            out.push(Case { q, t, mat, bw, zdrop, end_bonus, flags, adaptive, shift });
        }
        out
    }

    type Kernel = unsafe fn(&[u8], &[u8], i8, &[i8], i8, i8, i8, i8, i32, i32, i32, i32, i32, i32, &mut DpResult);

    fn run(k: Kernel, c: &Case) -> DpResult {
        let mut r = DpResult::default();
        unsafe { k(&c.q, &c.t, 5, &c.mat, 4, 2, 24, 1, c.bw, c.zdrop, c.end_bonus, c.flags, c.adaptive, c.shift, &mut r); }
        r
    }

    fn fields(r: &DpResult) -> ([i32; 10], &[u32]) {
        ([r.max, r.max_score_query_pos, r.max_score_target_pos, r.max_query_end_score, r.max_query_end_target_pos,
          r.max_target_end_score, r.max_target_end_query_pos, r.score, r.reach_end, r.zdropped], &r.cigar)
    }

    fn digest(rs: &[DpResult]) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut put = |v: u32| for b in v.to_le_bytes() { h ^= b as u64; h = h.wrapping_mul(0x0100_0000_01b3); };
        for r in rs {
            let (f, cig) = fields(r);
            for v in f { put(v as u32); }
            put(cig.len() as u32);
            for &c in cig { put(c); }
        }
        h
    }

    #[allow(unreachable_code)]
    fn platform_kernel() -> Option<Kernel> {
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("avx2") { return Some(super::super::dual_clean::extend_dual_affine_avx2_clean); }
        #[cfg(target_arch = "aarch64")]
        return Some(extend_dual_affine_neon_clean);
        None
    }

    /// Static bands: the clean kernels equal the scalar reference.
    #[test]
    fn clean16_static_band_equals_scalar() {
        let mut kernels: Vec<(&str, Kernel)> = Vec::new();
        #[cfg(target_arch = "x86_64")]
        if is_x86_feature_detected!("sse4.1") { kernels.push(("sse4.1", extend_dual_affine_sse41_clean)); }
        if let Some(k) = platform_kernel() { kernels.push(("platform", k)); }
        let mut n = 0;
        for (i, c) in cases().iter().enumerate() {
            if c.adaptive != 0 || c.shift != 0 { continue; }
            let mut want = DpResult::default();
            super::super::extend_dual_affine_scalar(&c.q, &c.t, 5, &c.mat, 4, 2, 24, 1, c.bw, c.zdrop, c.end_bonus, c.flags, &mut want);
            for (name, k) in &kernels {
                let got = run(*k, c);
                assert_eq!(fields(&got), fields(&want), "{name} vs scalar, case {i} (qlen {} tlen {} bw {} flags {:#x})", c.q.len(), c.t.len(), c.bw, c.flags);
            }
            n += 1;
        }
        assert!(n > 1000);
    }

    /// The SSE4.1 instance of the 16-lane kernel equals the AVX2 kernel on every case,
    /// including adaptive and shifted bands (which the scalar reference does not implement).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn clean16_sse41_equals_avx2() {
        if !is_x86_feature_detected!("avx2") || !is_x86_feature_detected!("sse4.1") { return; }
        for (i, c) in cases().iter().enumerate() {
            let want = run(super::super::dual_clean::extend_dual_affine_avx2_clean, c);
            let got = run(extend_dual_affine_sse41_clean, c);
            assert_eq!(fields(&got), fields(&want), "case {i} (qlen {} tlen {} bw {} adaptive {} shift {} flags {:#x})",
                c.q.len(), c.t.len(), c.bw, c.adaptive, c.shift, c.flags);
        }
    }

    /// Cross-platform check: the platform's clean kernel reproduces the recorded digest.
    #[test]
    fn clean16_golden_digest() {
        let Some(k) = platform_kernel() else { return };
        let rs: Vec<DpResult> = cases().iter().map(|c| run(k, c)).collect();
        let d = digest(&rs);
        println!("clean kernel digest: {d:#018x}");
        assert_eq!(d, GOLDEN, "clean kernel results differ from the recorded x86 AVX2 digest");
    }
}
