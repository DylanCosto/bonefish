//! Register-resident clean DP for static or adaptive bands of 33..=64 active cells: the
//! 64-lane counterpart of `dual_clean::reg1_impl`, with `clean_impl`'s results. The six
//! difference vectors stay in registers indexed from the band start. AVX-512VBMI makes each
//! full-register byte shift a single permute, and a one-cell band advance cancels the left
//! shift of x/v/x2, so an anti-diagonal needs three shifts. `dual_clean512` dispatches here;
//! CPUs without VBMI use its array kernel.
use super::common::*;
use core::arch::x86_64::*;

#[target_feature(enable = "avx2,avx512f,avx512bw,avx512vbmi")]
#[inline]
unsafe fn shl1_prev_vbmi(v: __m512i, prev: __m512i) -> __m512i {
    unsafe {
        const INDEX: [u8; 64] = [
            63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
            85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
            105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
            122, 123, 124, 125, 126,
        ];
        _mm512_permutex2var_epi8(prev, _mm512_loadu_si512(INDEX.as_ptr().cast()), v)
    }
}
#[target_feature(enable = "avx2,avx512f,avx512bw,avx512vbmi")]
#[inline]
unsafe fn shr1_next_vbmi(v: __m512i, fill: __m512i) -> __m512i {
    unsafe {
        const INDEX: [u8; 64] = [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
            25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46,
            47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
        ];
        _mm512_permutex2var_epi8(v, _mm512_loadu_si512(INDEX.as_ptr().cast()), fill)
    }
}

#[target_feature(enable = "avx2,avx512f,avx512bw,avx512vbmi")]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn reg64_impl<
    const APPROX: bool,
    const CIGAR: bool,
    const RIGHT: bool,
    const HAS_N: bool,
>(
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
) {
    unsafe {
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
        let (qi, ei, q2i, e2i) = (
            gap_open as i32,
            gap_extend as i32,
            gap_open2 as i32,
            gap_extend2 as i32,
        );
        let mut long_thres: i32 = if ei != e2i {
            (q2i - qi) / (ei - e2i) - 1
        } else {
            0
        };
        if (q2i + e2i + long_thres * e2i) > (qi + ei + long_thres * ei) {
            long_thres += 1;
        }
        let long_diff: i8 = (long_thres * (ei - e2i) - (q2i - qi) - e2i) as i8;
        let neg_qe = (-(qi + ei)) as i8;
        let neg_q2e2 = (-(q2i + e2i)) as i8;

        let zero_ = _mm512_setzero_si512();
        let q_ = _mm512_set1_epi8(gap_open);
        let q2_ = _mm512_set1_epi8(gap_open2);
        let qe_ = _mm512_set1_epi8((qi + ei) as i8);
        let qe2_ = _mm512_set1_epi8((q2i + e2i) as i8);
        let neg_qe_ = _mm512_set1_epi8(neg_qe);
        let neg_q2e2_ = _mm512_set1_epi8(neg_q2e2);
        let sc_mch_ = _mm512_set1_epi8(score_matrix[0]);
        let sc_mis_ = _mm512_set1_epi8(score_matrix[1]);
        let n_sc = score_matrix[(alphabet_size as usize * alphabet_size as usize) - 1];
        let sc_n_ = _mm512_set1_epi8(if n_sc == 0 { -gap_extend2 } else { n_sc });
        let m1_ = _mm512_set1_epi8(alphabet_size - 1);
        let sc_mch_mis_ = _mm512_set1_epi8(score_matrix[0].wrapping_sub(score_matrix[1]));
        let flag1_ = _mm512_set1_epi8(1);
        let flag2_ = _mm512_set1_epi8(2);
        let flag3_ = _mm512_set1_epi8(3);
        let flag4_ = _mm512_set1_epi8(4);
        let flag8_ = _mm512_set1_epi8(0x08);
        let flag16_ = _mm512_set1_epi8(0x10);
        let flag32_ = _mm512_set1_epi8(0x20);
        let flag64_ = _mm512_set1_epi8(0x40);
        let wl = if bandwidth < 0 {
            query_len.max(target_len) as i32
        } else {
            bandwidth
        };
        let adaptive = if adaptive > 0 {
            adaptive.min(wl.max(1))
        } else {
            0
        };
        let mut a_st: i32 = 0;
        let mut a_step: i32 = 0;

        // Target and reversed query, padded for 64-byte loads; traceback rows and band offsets.
        let tpad = (target_len + 128).next_multiple_of(64);
        let qpad = (query_len + 128).next_multiple_of(64);
        let n_col = super::dual_clean::clean_band_cells(query_len, target_len, bandwidth, adaptive)
            .div_ceil(64);
        let stride = n_col * 64;
        let n_diag = query_len + target_len - 1;
        let state_bytes = tpad + qpad;
        let (p_off, off_off, off_end_off, total) = if with_cigar {
            let p_off = state_bytes;
            let off_off = (p_off + n_diag * stride + 63) & !63;
            let off_end_off = (off_off + n_diag * 4 + 63) & !63;
            (p_off, off_off, off_end_off, off_end_off + n_diag * 4)
        } else {
            (0, 0, 0, state_bytes)
        };
        let mem = AlignedMemory::new(total, 64);
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
            (
                base.add(p_off),
                base.add(off_off) as *mut i32,
                base.add(off_end_off) as *mut i32,
            )
        } else {
            (
                std::ptr::null_mut(),
                std::ptr::null_mut(),
                std::ptr::null_mut(),
            )
        };

        let (h_vec, h_ptr) = alloc_h_array(approx_max, tpad / 64, 64);
        let _ = &h_vec;
        let packed_h = query_len + target_len <= 20000 && wl >= 32;
        init_dp_result(result);

        // Registers: lane j holds the array value at target position `lo + j`.
        let (mut ur, mut vr, mut xr, mut yr) = (neg_qe_, neg_qe_, neg_qe_, neg_qe_);
        let (mut x2r, mut y2r) = (neg_q2e2_, neg_q2e2_);
        let mut lo: i32 = 0;
        let mut h0: i32 = 0;
        let mut last_h0_t: i32 = 0;
        let qe_scalar = qi + ei;
        #[repr(align(64))]
        struct Lanes([i8; 64]);
        let mut ub = Lanes([0; 64]);
        let mut vb = Lanes([0; 64]);
        let top_u = |r: i32| -> i8 {
            if r == 0 {
                neg_qe
            } else if r < long_thres {
                -gap_extend
            } else if r == long_thres {
                long_diff
            } else {
                -gap_extend2
            }
        };

        for r in 0..n_diag as i32 {
            let mut st0 = 0i32;
            let mut en0 = target_len as i32 - 1;
            if st0 < r - query_len as i32 + 1 {
                st0 = r - query_len as i32 + 1;
            }
            if en0 > r {
                en0 = r;
            }
            if adaptive > 0 {
                a_st = (a_st + a_step).max(st0);
                st0 = a_st;
                en0 = en0.min(st0 + adaptive - 1);
            } else {
                if st0 < (r - wl + band_shift + 1) >> 1 {
                    st0 = (r - wl + band_shift + 1) >> 1;
                }
                if en0 > (r + wl + band_shift) >> 1 {
                    en0 = (r + wl + band_shift) >> 1;
                }
            }
            if st0 > en0 {
                result.zdropped = 1;
                break;
            }
            debug_assert!(en0 - st0 < 64);

            // A one-cell band advance cancels the recurrence's left shift for x/v/x2.
            // Keep those old vectors directly; only the same-position u/y/y2 inputs move right.
            // Static and adaptive band starts advance by at most one cell per anti-diagonal.
            let (xt1, vt1, x2t1) = if lo < st0 {
                debug_assert_eq!(st0, lo + 1);
                ur = shr1_next_vbmi(ur, neg_qe_);
                yr = shr1_next_vbmi(yr, neg_qe_);
                y2r = shr1_next_vbmi(y2r, neg_q2e2_);
                lo = st0;
                (xr, vr, x2r)
            } else {
                let v1 = if st0 == 0 { top_u(r) } else { neg_qe };
                (
                    shl1_prev_vbmi(xr, neg_qe_),
                    shl1_prev_vbmi(vr, _mm512_set1_epi8(v1)),
                    shl1_prev_vbmi(x2r, neg_q2e2_),
                )
            };
            if en0 >= r {
                // New top cell at target position r.
                let m = 1u64 << (r - st0);
                yr = _mm512_mask_mov_epi8(yr, m, neg_qe_);
                y2r = _mm512_mask_mov_epi8(y2r, m, neg_q2e2_);
                ur = _mm512_mask_mov_epi8(ur, m, _mm512_set1_epi8(top_u(r)));
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

            let sq = _mm512_loadu_si512(sf.add(st) as *const __m512i);
            let qv = _mm512_loadu_si512(qrr.add(st) as *const __m512i);
            let mut z = _mm512_add_epi8(
                sc_mis_,
                _mm512_maskz_mov_epi8(_mm512_cmpeq_epi8_mask(sq, qv), sc_mch_mis_),
            );
            if HAS_N {
                z = _mm512_mask_blend_epi8(
                    _mm512_cmpeq_epi8_mask(_mm512_max_epu8(sq, qv), m1_),
                    z,
                    sc_n_,
                );
            }
            let ut = ur;
            let a = _mm512_add_epi8(xt1, vt1);
            let b = _mm512_add_epi8(yr, ut);
            let a2 = _mm512_add_epi8(x2t1, vt1);
            let b2 = _mm512_add_epi8(y2r, ut);
            if !with_cigar {
                z = _mm512_max_epi8(z, a);
                z = _mm512_max_epi8(z, b);
                z = _mm512_max_epi8(z, a2);
                z = _mm512_max_epi8(z, b2);
                z = _mm512_min_epi8(z, sc_mch_);
                ur = _mm512_sub_epi8(z, vt1);
                vr = _mm512_sub_epi8(z, ut);
                let tmp1 = _mm512_sub_epi8(z, q_);
                let tmp2 = _mm512_sub_epi8(z, q2_);
                xr = _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(a, tmp1), zero_), qe_);
                yr = _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(b, tmp1), zero_), qe_);
                x2r = _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(a2, tmp2), zero_), qe2_);
                y2r = _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(b2, tmp2), zero_), qe2_);
            } else {
                let mut d;
                if !right_align {
                    let tmp = _mm512_cmpgt_epi8_mask(a, z);
                    d = _mm512_maskz_mov_epi8(tmp, flag1_);
                    z = _mm512_max_epi8(z, a);
                    d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b, z), d, flag2_);
                    z = _mm512_max_epi8(z, b);
                    d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(a2, z), d, flag3_);
                    z = _mm512_max_epi8(z, a2);
                    d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(b2, z), d, flag4_);
                    z = _mm512_max_epi8(z, b2);
                } else {
                    d = _mm512_maskz_mov_epi8(!_mm512_cmpgt_epi8_mask(z, a), flag1_);
                    z = _mm512_max_epi8(z, a);
                    d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b), flag2_, d);
                    z = _mm512_max_epi8(z, b);
                    d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, a2), flag3_, d);
                    z = _mm512_max_epi8(z, a2);
                    d = _mm512_mask_blend_epi8(_mm512_cmpgt_epi8_mask(z, b2), flag4_, d);
                    z = _mm512_max_epi8(z, b2);
                }
                z = _mm512_min_epi8(z, sc_mch_);
                ur = _mm512_sub_epi8(z, vt1);
                vr = _mm512_sub_epi8(z, ut);
                let tmp1 = _mm512_sub_epi8(z, q_);
                let tmp2 = _mm512_sub_epi8(z, q2_);
                let an = _mm512_sub_epi8(a, tmp1);
                let bn = _mm512_sub_epi8(b, tmp1);
                let a2n = _mm512_sub_epi8(a2, tmp2);
                let b2n = _mm512_sub_epi8(b2, tmp2);
                let (ma, mb, ma2, mb2) = if !right_align {
                    (
                        _mm512_cmpgt_epi8_mask(an, zero_),
                        _mm512_cmpgt_epi8_mask(bn, zero_),
                        _mm512_cmpgt_epi8_mask(a2n, zero_),
                        _mm512_cmpgt_epi8_mask(b2n, zero_),
                    )
                } else {
                    (
                        !_mm512_cmpgt_epi8_mask(zero_, an),
                        !_mm512_cmpgt_epi8_mask(zero_, bn),
                        !_mm512_cmpgt_epi8_mask(zero_, a2n),
                        !_mm512_cmpgt_epi8_mask(zero_, b2n),
                    )
                };
                xr = _mm512_sub_epi8(_mm512_maskz_mov_epi8(ma, an), qe_);
                yr = _mm512_sub_epi8(_mm512_maskz_mov_epi8(mb, bn), qe_);
                x2r = _mm512_sub_epi8(_mm512_maskz_mov_epi8(ma2, a2n), qe2_);
                y2r = _mm512_sub_epi8(_mm512_maskz_mov_epi8(mb2, b2n), qe2_);
                // Direction uses bits 0..2 and continuation bits are disjoint:
                // masked addition is exactly OR, with no possible byte carry.
                d = _mm512_mask_add_epi8(d, ma, d, flag8_);
                d = _mm512_mask_add_epi8(d, mb, d, flag16_);
                d = _mm512_mask_add_epi8(d, ma2, d, flag32_);
                d = _mm512_mask_add_epi8(d, mb2, d, flag64_);
                _mm512_storeu_si512(prow.add(st) as *mut __m512i, d);
            }
            // Inactive lanes preserve the scalar kernel's initial boundary values.
            let active = u64::MAX >> (64 - (en - st + 1));
            ur = _mm512_mask_mov_epi8(neg_qe_, active, ur);
            vr = _mm512_mask_mov_epi8(neg_qe_, active, vr);
            xr = _mm512_mask_mov_epi8(neg_qe_, active, xr);
            yr = _mm512_mask_mov_epi8(neg_qe_, active, yr);
            x2r = _mm512_mask_mov_epi8(neg_q2e2_, active, x2r);
            y2r = _mm512_mask_mov_epi8(neg_q2e2_, active, y2r);

            // Exact-max reduction keeps the array kernel's four-lane tie priority. The
            // difference states stay in registers; only u/v need spill slots for score tracking.
            _mm512_store_si512(ub.0.as_mut_ptr().cast(), ur);
            _mm512_store_si512(vb.0.as_mut_ptr().cast(), vr);
            let u = ub.0.as_ptr() as *const u8;
            let v = vb.0.as_ptr() as *const u8;
            if !approx_max {
                let mut max_h: i32;
                let mut max_t: i32;
                if r > 0 {
                    let h_en0 = if en0 > 0 {
                        *h_ptr.add(en - 1) + *u.add((en) - st) as i8 as i32
                    } else {
                        *h_ptr.add(en) + *v.add((en) - st) as i8 as i32
                    };
                    *h_ptr.add(en) = h_en0;
                    max_h = h_en0;
                    max_t = en0;
                    let en1 = st0 + (en0 - st0) / 4 * 4;
                    let en8 = st0 + (en1 - st0) / 8 * 8;
                    let mut tt = st0;
                    if packed_h {
                        // Packed (H, lane, position) keys, as in `clean_impl`.
                        let mut low = _mm256_setr_epi32(
                            32767,
                            32767 - (1 << 13),
                            32767 - (2 << 13),
                            32767 - (3 << 13),
                            32766,
                            32766 - (1 << 13),
                            32766 - (2 << 13),
                            32766 - (3 << 13),
                        );
                        let two = _mm256_set1_epi32(2);
                        let mut mk = _mm256_set1_epi32(i32::MIN);
                        while tt < en8 {
                            let hp = h_ptr.add(tt as usize);
                            let vv = _mm256_cvtepi8_epi32(_mm_loadl_epi64(
                                v.add((tt as usize) - st) as *const __m128i,
                            ));
                            let h1 = _mm256_add_epi32(_mm256_loadu_si256(hp as *const __m256i), vv);
                            _mm256_storeu_si256(hp as *mut __m256i, h1);
                            mk = _mm256_max_epi32(
                                mk,
                                _mm256_or_si256(_mm256_slli_epi32(h1, 15), low),
                            );
                            low = _mm256_sub_epi32(low, two);
                            tt += 8;
                        }
                        let mut m4 = _mm_max_epi32(
                            _mm256_castsi256_si128(mk),
                            _mm256_extracti128_si256(mk, 1),
                        );
                        if tt < en1 {
                            let hp = h_ptr.add(tt as usize);
                            let vv = _mm_cvtepi8_epi32(_mm_cvtsi32_si128(
                                std::ptr::read_unaligned(v.add((tt as usize) - st) as *const i32),
                            ));
                            let h1 = _mm_add_epi32(_mm_loadu_si128(hp as *const __m128i), vv);
                            _mm_storeu_si128(hp as *mut __m128i, h1);
                            m4 = _mm_max_epi32(
                                m4,
                                _mm_or_si128(_mm_slli_epi32(h1, 15), _mm256_castsi256_si128(low)),
                            );
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
                                *h_ptr.add(pos) += *v.add((pos) - st) as i8 as i32;
                                if *h_ptr.add(pos) > lane_h[i as usize] {
                                    lane_h[i as usize] = *h_ptr.add(pos);
                                    lane_t[i as usize] = tt + i;
                                }
                            }
                            tt += 4;
                        }
                        for i in 0..4 {
                            if max_h < lane_h[i] {
                                max_h = lane_h[i];
                                max_t = lane_t[i];
                            }
                        }
                    }
                    while tt < en0 {
                        *h_ptr.add(tt as usize) += *v.add((tt as usize) - st) as i8 as i32;
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
                } else if max_t >= result.max_score_target_pos
                    && (r - max_t) >= result.max_score_query_pos
                {
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
                    if last_h0_t >= st0
                        && last_h0_t <= en0
                        && last_h0_t + 1 >= st0
                        && last_h0_t < en0
                    {
                        // which neighbour the greedy path takes is data-dependent: select, don't branch
                        let d0 = *v.add((last_h0_t as usize) - st) as i8 as i32;
                        let d1 = *u.add(((last_h0_t + 1) as usize) - st) as i8 as i32;
                        let right = d0 <= d1;
                        h0 += if right { d1 } else { d0 };
                        last_h0_t += right as i32;
                    } else if last_h0_t >= st0 && last_h0_t <= en0 {
                        h0 += *v.add((last_h0_t as usize) - st) as i8 as i32;
                    } else {
                        last_h0_t += 1;
                        h0 += *u.add((last_h0_t as usize) - st) as i8 as i32;
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
                    } else if last_h0_t >= result.max_score_target_pos
                        && (r - last_h0_t) >= result.max_score_query_pos
                    {
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
        }

        if approx_max && result.score == NEG_INF {
            result.score = result.max;
        }
        if with_cigar {
            traceback_dual_affine(
                result,
                query_len,
                target_len,
                end_bonus,
                flags,
                n_col,
                64,
                p,
                band_off,
                band_off_end,
            );
        }
    }
}
