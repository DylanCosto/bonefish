//! 64-lane AVX-512 clean dual-affine DP: the semantics, band handling and tie-breaking of
//! `dual_clean::clean_impl` (bit-identical results) for bands wider than one 32-byte vector.
//! `dual_clean::extend_dual_affine_clean` routes calls here when `use_avx512()` allows; bands
//! of 33..=64 cells use the register-resident kernel (`dual_clean_reg64`) on CPUs with
//! AVX-512VBMI.
//!
//! Beyond the vector width, the array kernel differs from `clean_impl` only in ways that leave
//! results unchanged:
//! - x, v and x2 are double-buffered so their one-cell-shifted inputs are plain loads;
//! - bands of `TRACE_PLANES_MIN_CELLS` or more store traceback as bit planes (see `common.rs`),
//!   and traceback matrices of `TRACE_STREAM_MIN_BYTES` or more use non-temporal stores;
//! - with the default scoring (`shifted_clip_eligible`) the gap updates use a shorter form;
//! - the exact-max score pass keeps `clean_impl`'s packed (score, lane, position) keys, so ties
//!   resolve the same way regardless of vector width.

use super::common::*;
use core::arch::x86_64::*;

// The macros below are used only inside `avx2,avx512f,avx512bw` target-feature functions.

/// Five-way max of the score and the four gap candidates, returning the new score and the four
/// bit-plane masks (`common.rs` layout). Ties resolve as in the byte-direction code: towards the
/// earlier candidate for left alignment, the later one for right alignment.
macro_rules! trace_plane_argmax {
    ($score:expr,$a:expr,$b:expr,$a2:expr,$b2:expr,$right:expr) => {{
        let (s, a, b, a2, b2) = ($score, $a, $b, $a2, $b2);
        let short = _mm512_max_epi8(a, b);
        let long = _mm512_max_epi8(a2, b2);
        let gaps = _mm512_max_epi8(short, long);
        let z = _mm512_max_epi8(s, gaps);
        let short_b = if $right {
            !_mm512_cmpgt_epi8_mask(a, b)
        } else {
            _mm512_cmpgt_epi8_mask(b, a)
        };
        let long_b = if $right {
            !_mm512_cmpgt_epi8_mask(a2, b2)
        } else {
            _mm512_cmpgt_epi8_mask(b2, a2)
        };
        let long_wins = if $right {
            !_mm512_cmpgt_epi8_mask(short, long)
        } else {
            _mm512_cmpgt_epi8_mask(long, short)
        };
        let gap_wins = if $right {
            !_mm512_cmpgt_epi8_mask(s, gaps)
        } else {
            _mm512_cmpgt_epi8_mask(gaps, s)
        };
        (z, gap_wins, long_wins, short_b, long_b)
    }};
}

/// Stores four mask registers as adjacent u64 planes at `$ptr`, with non-temporal stores when
/// `$nt`. Not `pure`/`nomem`/`readonly`: these asm blocks write memory. The scratch register is
/// an output, so it cannot alias `ptr`.
macro_rules! store_four_trace_planes {
    ($ptr:expr, $nt:expr, $m0:expr, $m1:expr, $m2:expr, $m3:expr) => {{
        let address = $ptr;
        if $nt {
            let _scratch: u64;
            core::arch::asm!(
                "kmovq {tmp}, {m0}", "movnti qword ptr [{ptr}], {tmp}",
                "kmovq {tmp}, {m1}", "movnti qword ptr [{ptr}+8], {tmp}",
                "kmovq {tmp}, {m2}", "movnti qword ptr [{ptr}+16], {tmp}",
                "kmovq {tmp}, {m3}", "movnti qword ptr [{ptr}+24], {tmp}",
                ptr = in(reg) address, tmp = out(reg) _scratch,
                m0 = in(kreg) $m0, m1 = in(kreg) $m1, m2 = in(kreg) $m2, m3 = in(kreg) $m3,
                options(nostack, preserves_flags),
            );
        } else {
            core::arch::asm!(
                "kmovq qword ptr [{ptr}], {m0}",
                "kmovq qword ptr [{ptr}+8], {m1}",
                "kmovq qword ptr [{ptr}+16], {m2}",
                "kmovq qword ptr [{ptr}+24], {m3}",
                ptr = in(reg) address,
                m0 = in(kreg) $m0, m1 = in(kreg) $m1, m2 = in(kreg) $m2, m3 = in(kreg) $m3,
                options(nostack, preserves_flags),
            );
        }
    }};
}

/// `$d` with `$f8`/`$f16`/`$f32`/`$f64` added in the lanes of masks `$ma`/`$mb`/`$ma2`/`$mb2`.
/// Written as asm because LLVM turns some of these disjoint-bit additions into OR plus a masked
/// move; this keeps exactly four masked additions, all register operands.
macro_rules! append_continuation_bits_native {
    ($d:expr, $ma:expr, $mb:expr, $ma2:expr, $mb2:expr, $f8:expr, $f16:expr, $f32:expr, $f64:expr) => {{
        let mut value = $d;
        core::arch::asm!(
            "vpaddb {value}{{{ma}}}, {value}, {f8}",
            "vpaddb {value}{{{mb}}}, {value}, {f16}",
            "vpaddb {value}{{{ma2}}}, {value}, {f32}",
            "vpaddb {value}{{{mb2}}}, {value}, {f64}",
            value = inout(zmm_reg) value,
            ma = in(kreg) $ma, mb = in(kreg) $mb, ma2 = in(kreg) $ma2, mb2 = in(kreg) $mb2,
            f8 = in(zmm_reg) $f8, f16 = in(zmm_reg) $f16, f32 = in(zmm_reg) $f32, f64 = in(zmm_reg) $f64,
            options(pure, nomem, nostack, preserves_flags),
        );
        value
    }};
}

/// One 16-cell block of the packed exact-max pass: adds the 16 score deltas at `v` to the H
/// scores at `h` and folds the (H << 15 | `low`) keys into `mk`. `low` encodes the lane and
/// position so the maximum key resolves ties in the 4-lane order of `clean_impl`.
#[target_feature(enable = "avx2,avx512f,avx512bw")]
#[inline]
unsafe fn update_h16_packed(h: *mut i32, v: *const u8, low: __m512i, mk: __m512i) -> __m512i {
    unsafe {
        let vv = _mm512_cvtepi8_epi32(_mm_loadu_si128(v.cast()));
        let next = _mm512_add_epi32(_mm512_loadu_si512(h.cast()), vv);
        _mm512_storeu_si512(h.cast(), next);
        _mm512_max_epi32(mk, _mm512_or_si512(_mm512_slli_epi32(next, 15), low))
    }
}

/// Whether the shortened gap arithmetic (`FAST`) is valid: the default scoring (match 2,
/// mismatch -4, N score in -4..=2, gaps 4,2,24,1 in either order), a static band and lengths
/// far from i32 overflow. Other scoring keeps the general wrapping-byte expressions.
#[allow(clippy::too_many_arguments)]
#[inline]
fn shifted_clip_eligible(
    ql: usize,
    tl: usize,
    alphabet: i8,
    mat: &[i8],
    q: i8,
    e: i8,
    q2: i8,
    e2: i8,
    bw: i32,
    ad: i32,
    shift: i32,
) -> bool {
    if alphabet != 5
        || mat.len() < 25
        || mat[0] != 2
        || mat[1] != -4
        || !(-4..=2).contains(&mat[24])
        || ad > 0
    {
        return false;
    }
    if !((q, e, q2, e2) == (4, 2, 24, 1) || (q, e, q2, e2) == (24, 1, 4, 2)) {
        return false;
    }
    if ql > i32::MAX as usize || tl > i32::MAX as usize {
        return false;
    }
    let wl = if bw < 0 { ql.max(tl) as i64 } else { bw as i64 };
    ql as i64 + tl as i64 + wl + (shift as i64).abs() + 2 <= i32::MAX as i64
}

#[cfg(test)]
mod tests;

/// Bands of at least this many cells store traceback as bit planes; narrower bands (one or two
/// vectors per anti-diagonal) keep direction bytes.
const TRACE_PLANES_MIN_CELLS: usize = 128;

/// Traceback matrices of at least this many bytes are written with non-temporal stores, which
/// pay off only for large matrices. Rows and vector offsets are multiples of 64 bytes, so each
/// store writes a whole cache line.
const TRACE_STREAM_MIN_BYTES: usize = 1 << 20;

/// Bands of `REG64_MIN_CELLS..=64` cells use `dual_clean_reg64` when the CPU has AVX-512VBMI
/// (narrower bands are routed to the 256-bit kernels before reaching this module).
const REG64_MIN_CELLS: usize = 33;

#[target_feature(enable = "avx2,avx512f,avx512bw")]
#[allow(clippy::too_many_arguments)]
pub(super) unsafe fn extend_dual_affine_avx512_clean(
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
        // One specialised copy per flag combination: with the flags as constants the inner loop has
        // no per-vector branches and its masks stay in vector registers.
        let m1 = alphabet_size.wrapping_sub(1) as u8;
        // Branch-free max reductions (vectorised), not an early-exit scan.
        let has_n = qseq
            .iter()
            .fold(0u8, |m, &c| m.max(c))
            .max(tseq.iter().fold(0u8, |m, &c| m.max(c)))
            >= m1;
        let approx_max = (flags & APPROX_MAX) != 0;
        let with_cigar = (flags & SCORE_ONLY) == 0;
        let right_align = (flags & RIGHT_ALIGN) != 0;
        let cells =
            super::dual_clean::clean_band_cells(qseq.len(), tseq.len(), bandwidth, adaptive);
        let use_registers =
            (REG64_MIN_CELLS..=64).contains(&cells) && is_x86_feature_detected!("avx512vbmi");
        macro_rules! go {
            ($a:literal, $c:literal, $r:literal, $n:literal) => {
                if use_registers {
                    super::dual_clean_reg64::reg64_impl::<$a, $c, $r, $n>(
                        qseq,
                        tseq,
                        alphabet_size,
                        score_matrix,
                        gap_open,
                        gap_extend,
                        gap_open2,
                        gap_extend2,
                        bandwidth,
                        z_drop,
                        end_bonus,
                        flags,
                        adaptive,
                        band_shift,
                        result,
                    )
                } else {
                    clean_impl512::<$a, $c, $r, $n>(
                        qseq,
                        tseq,
                        alphabet_size,
                        score_matrix,
                        gap_open,
                        gap_extend,
                        gap_open2,
                        gap_extend2,
                        bandwidth,
                        z_drop,
                        end_bonus,
                        flags,
                        adaptive,
                        band_shift,
                        result,
                    )
                }
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
    }
}

/// The 64-lane array kernel: picks the traceback layout, store mode and gap arithmetic.
#[target_feature(enable = "avx2,avx512f,avx512bw")]
#[allow(clippy::too_many_arguments)]
unsafe fn clean_impl512<
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
        if alphabet_size <= 1 || qseq.is_empty() || tseq.is_empty() {
            return;
        }
        let cells =
            super::dual_clean::clean_band_cells(qseq.len(), tseq.len(), bandwidth, adaptive);
        let planes = CIGAR && cells >= TRACE_PLANES_MIN_CELLS;
        let stride = cells.div_ceil(64) * 64;
        let n_diag = qseq.len() + tseq.len() - 1;
        let stream = CIGAR && n_diag.saturating_mul(stride) >= TRACE_STREAM_MIN_BYTES;
        let fast = CIGAR
            && shifted_clip_eligible(
                qseq.len(),
                tseq.len(),
                alphabet_size,
                score_matrix,
                gap_open,
                gap_extend,
                gap_open2,
                gap_extend2,
                bandwidth,
                adaptive,
                band_shift,
            );
        macro_rules! go {
            ($p:literal, $s:literal, $f:literal) => {
                clean_impl512_kernel::<APPROX, CIGAR, RIGHT, HAS_N, $p, $s, $f, false>(
                    qseq,
                    tseq,
                    alphabet_size,
                    score_matrix,
                    gap_open,
                    gap_extend,
                    gap_open2,
                    gap_extend2,
                    bandwidth,
                    z_drop,
                    end_bonus,
                    flags,
                    adaptive,
                    band_shift,
                    result,
                )
            };
        }
        match (planes, stream, fast) {
            (false, false, false) => go!(false, false, false),
            (false, false, true) => go!(false, false, true),
            (false, true, false) => go!(false, true, false),
            (false, true, true) => go!(false, true, true),
            (true, false, false) => go!(true, false, false),
            (true, false, true) => go!(true, false, true),
            (true, true, false) => go!(true, true, false),
            (true, true, true) => go!(true, true, true),
        }
    }
}

/// Array kernel body. `PLANES`: bit-plane rather than byte traceback; `STREAM`: non-temporal
/// traceback stores; `FAST`: shortened gap arithmetic (requires `shifted_clip_eligible`);
/// `CHECK`: test-only assertions that the `FAST` form agrees with the general one on every lane.
#[target_feature(enable = "avx2,avx512f,avx512bw")]
#[allow(clippy::too_many_arguments)]
unsafe fn clean_impl512_kernel<
    const APPROX: bool,
    const CIGAR: bool,
    const RIGHT: bool,
    const HAS_N: bool,
    const PLANES: bool,
    const STREAM: bool,
    const FAST: bool,
    const CHECK: bool,
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
        let e_ = _mm512_set1_epi8(gap_extend);
        let e2_ = _mm512_set1_epi8(gap_extend2);
        let neg_qe_ = _mm512_set1_epi8(neg_qe);
        let neg_q2e2_ = _mm512_set1_epi8(neg_q2e2);
        let sc_mch_ = _mm512_set1_epi8(score_matrix[0]);
        let sc_mis_ = _mm512_set1_epi8(score_matrix[1]);
        let n_sc = score_matrix[(alphabet_size as usize * alphabet_size as usize) - 1];
        let sc_n_ = _mm512_set1_epi8(if n_sc == 0 { -gap_extend2 } else { n_sc });
        let m1_ = _mm512_set1_epi8(alphabet_size - 1);
        let has_n = HAS_N;
        let flag1_ = _mm512_set1_epi8(1);
        let flag2_ = _mm512_set1_epi8(2);
        let flag3_ = _mm512_set1_epi8(3);
        let flag4_ = _mm512_set1_epi8(4);
        let flag8_ = _mm512_set1_epi8(0x08);
        let flag16_ = _mm512_set1_epi8(0x10);
        let flag32_ = _mm512_set1_epi8(0x20);
        let flag64_ = _mm512_set1_epi8(0x40);

        let wl = if bandwidth < 0 {
            target_len.max(query_len) as i32
        } else {
            bandwidth
        };
        // Adaptive band: a window of `adaptive` cells per anti-diagonal that moves one step towards
        // the anti-diagonal's best cell (exact max tracking) or the greedy path cell (approximate
        // tracking) each time, instead of the static band |t - q| <= wl. Its start and end never
        // decrease, so the scalar boundary semantics carry over unchanged; the matrix bound
        // st0 >= r - qlen + 1 pulls it onto the end corner for global alignment.
        let adaptive = if adaptive > 0 {
            adaptive.min(wl.max(1))
        } else {
            0
        };
        let mut a_st: i32 = 0;
        let mut a_step: i32 = 0;

        // Layout: nine state arrays (u, v, x, y, x2, y2, and the next rows of v, x, x2) and the
        // target, each tlen + 128 bytes (loads and the post-loop reset run up to 64 bytes past en0);
        // the reversed query, qlen + 128 bytes; then traceback rows and band offsets.
        let tpad = (target_len + 128).next_multiple_of(64);
        let qpad = (query_len + 128).next_multiple_of(64);
        let n_col = super::dual_clean::clean_band_cells(query_len, target_len, bandwidth, adaptive)
            .div_ceil(64);
        let stride = n_col * 64;
        let n_diag = query_len + target_len - 1;
        let stream_trace = STREAM;
        let state_bytes = 10 * tpad + qpad;
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
        let u = base;
        let mut v = u.add(tpad);
        let mut x = v.add(tpad);
        let y = x.add(tpad);
        let mut x2 = y.add(tpad);
        let y2 = x2.add(tpad);
        // x, v and x2 are read one cell to the left of where they are written: keep the previous
        // anti-diagonal's rows in separate buffers (swapped after each row) so the shifted inputs
        // are plain unaligned loads instead of cross-lane byte shifts.
        let mut v_next = y2.add(tpad);
        let mut x_next = v_next.add(tpad);
        let mut x2_next = x_next.add(tpad);
        let sf = x2_next.add(tpad);
        let qr = sf.add(tpad);
        for a in [u, v, x, y, v_next, x_next] {
            std::ptr::write_bytes(a, neg_qe as u8, tpad);
        }
        for a in [x2, y2, x2_next] {
            std::ptr::write_bytes(a, neg_q2e2 as u8, tpad);
        }
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
        // See `dual_clean::clean_impl`'s `packed_h` for the range argument.
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
                // Static band: diagonals t - q in [band_shift - wl, band_shift + wl].
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

            let (x1, x21, v1): (i8, i8, i8) = if st0 > 0 {
                if st0 > last_st && st0 - 1 <= last_en {
                    let i = (st0 - 1) as usize;
                    (*x.add(i) as i8, *x2.add(i) as i8, *v.add(i) as i8)
                } else {
                    (neg_qe, neg_q2e2, neg_qe)
                }
            } else {
                let v1 = if r == 0 {
                    neg_qe
                } else if r < long_thres {
                    -gap_extend
                } else if r == long_thres {
                    long_diff
                } else {
                    -gap_extend2
                };
                (neg_qe, neg_q2e2, v1)
            };
            if en0 >= r {
                let i = r as usize;
                *y.add(i) = neg_qe as u8;
                *y2.add(i) = neg_q2e2 as u8;
                *u.add(i) = (if r == 0 {
                    neg_qe
                } else if r < long_thres {
                    -gap_extend
                } else if r == long_thres {
                    long_diff
                } else {
                    -gap_extend2
                }) as u8;
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

            let mut t = st;
            let mut first = true;
            while t <= en {
                let sq = _mm512_loadu_si512(sf.add(t) as *const __m512i);
                let qv = _mm512_loadu_si512(qrr.add(t) as *const __m512i);
                let mut z =
                    _mm512_mask_blend_epi8(_mm512_cmpeq_epi8_mask(sq, qv), sc_mis_, sc_mch_);
                if has_n {
                    // nt4 codes are 0..=m1 (m1 = N): either base is N iff their max is m1.
                    z = _mm512_mask_blend_epi8(
                        _mm512_cmpeq_epi8_mask(_mm512_max_epu8(sq, qv), m1_),
                        z,
                        sc_n_,
                    );
                }

                // Cell t - 1 of the previous anti-diagonal (x, v, x2 are not written during the row).
                let (xt1, vt1, x2t1) = if first {
                    // Lane 0 takes the left boundary value; the byte before st0 is not part of
                    // the band (and may lie outside the array when st0 == 0).
                    (
                        _mm512_mask_loadu_epi8(
                            _mm512_set1_epi8(x1),
                            !1u64,
                            x.add(t).wrapping_sub(1).cast(),
                        ),
                        _mm512_mask_loadu_epi8(
                            _mm512_set1_epi8(v1),
                            !1u64,
                            v.add(t).wrapping_sub(1).cast(),
                        ),
                        _mm512_mask_loadu_epi8(
                            _mm512_set1_epi8(x21),
                            !1u64,
                            x2.add(t).wrapping_sub(1).cast(),
                        ),
                    )
                } else {
                    (
                        _mm512_loadu_si512(x.add(t - 1).cast()),
                        _mm512_loadu_si512(v.add(t - 1).cast()),
                        _mm512_loadu_si512(x2.add(t - 1).cast()),
                    )
                };
                let ut = _mm512_loadu_si512(u.add(t) as *const __m512i);
                let a = _mm512_add_epi8(xt1, vt1);
                let b = _mm512_add_epi8(_mm512_loadu_si512(y.add(t) as *const __m512i), ut);
                let a2 = _mm512_add_epi8(x2t1, vt1);
                let b2 = _mm512_add_epi8(_mm512_loadu_si512(y2.add(t) as *const __m512i), ut);

                if !with_cigar {
                    z = _mm512_max_epi8(z, a);
                    z = _mm512_max_epi8(z, b);
                    z = _mm512_max_epi8(z, a2);
                    z = _mm512_max_epi8(z, b2);
                    z = _mm512_min_epi8(z, sc_mch_);
                    _mm512_storeu_si512(u.add(t) as *mut __m512i, _mm512_sub_epi8(z, vt1));
                    _mm512_storeu_si512(v_next.add(t) as *mut __m512i, _mm512_sub_epi8(z, ut));
                    let tmp1 = _mm512_sub_epi8(z, q_);
                    let tmp2 = _mm512_sub_epi8(z, q2_);
                    _mm512_storeu_si512(
                        x_next.add(t) as *mut __m512i,
                        _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(a, tmp1), zero_), qe_),
                    );
                    _mm512_storeu_si512(
                        y.add(t) as *mut __m512i,
                        _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(b, tmp1), zero_), qe_),
                    );
                    _mm512_storeu_si512(
                        x2_next.add(t) as *mut __m512i,
                        _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(a2, tmp2), zero_), qe2_),
                    );
                    _mm512_storeu_si512(
                        y2.add(t) as *mut __m512i,
                        _mm512_sub_epi8(_mm512_max_epi8(_mm512_sub_epi8(b2, tmp2), zero_), qe2_),
                    );
                } else {
                    // Direction bytes (byte layout); the bit-plane layout stores masks instead.
                    let mut d = zero_;
                    if PLANES {
                        let (next, m0, m1, m2, m3) = trace_plane_argmax!(z, a, b, a2, b2, RIGHT);
                        z = next;
                        store_four_trace_planes!(prow.add(t), stream_trace, m0, m1, m2, m3);
                    } else if !right_align {
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
                    _mm512_storeu_si512(u.add(t) as *mut __m512i, _mm512_sub_epi8(z, vt1));
                    _mm512_storeu_si512(v_next.add(t) as *mut __m512i, _mm512_sub_epi8(z, ut));
                    // FAST: x = max(a - z - e, -(q + e)), the same value as max(a - z + q, 0) - (q + e)
                    // without the final subtraction; `shifted_clip_eligible` bounds the operands so
                    // neither form wraps. The continuation thresholds shift by -(q + e) to match.
                    let tmp1 = if FAST {
                        _mm512_add_epi8(z, e_)
                    } else {
                        _mm512_sub_epi8(z, q_)
                    };
                    let tmp2 = if FAST {
                        _mm512_add_epi8(z, e2_)
                    } else {
                        _mm512_sub_epi8(z, q2_)
                    };
                    let an = _mm512_sub_epi8(a, tmp1);
                    let bn = _mm512_sub_epi8(b, tmp1);
                    let a2n = _mm512_sub_epi8(a2, tmp2);
                    let b2n = _mm512_sub_epi8(b2, tmp2);
                    let threshold = if FAST { neg_qe_ } else { zero_ };
                    let threshold2 = if FAST { neg_q2e2_ } else { zero_ };
                    // Left alignment extends a gap when the new value is > 0, right alignment when >= 0.
                    let (ma, mb, ma2, mb2) = if !right_align {
                        (
                            _mm512_cmpgt_epi8_mask(an, threshold),
                            _mm512_cmpgt_epi8_mask(bn, threshold),
                            _mm512_cmpgt_epi8_mask(a2n, threshold2),
                            _mm512_cmpgt_epi8_mask(b2n, threshold2),
                        )
                    } else {
                        (
                            !_mm512_cmpgt_epi8_mask(threshold, an),
                            !_mm512_cmpgt_epi8_mask(threshold, bn),
                            !_mm512_cmpgt_epi8_mask(threshold2, a2n),
                            !_mm512_cmpgt_epi8_mask(threshold2, b2n),
                        )
                    };
                    _mm512_storeu_si512(
                        x_next.add(t) as *mut __m512i,
                        if FAST {
                            _mm512_max_epi8(an, neg_qe_)
                        } else {
                            _mm512_sub_epi8(_mm512_max_epi8(an, zero_), qe_)
                        },
                    );
                    _mm512_storeu_si512(
                        y.add(t) as *mut __m512i,
                        if FAST {
                            _mm512_max_epi8(bn, neg_qe_)
                        } else {
                            _mm512_sub_epi8(_mm512_max_epi8(bn, zero_), qe_)
                        },
                    );
                    _mm512_storeu_si512(
                        x2_next.add(t) as *mut __m512i,
                        if FAST {
                            _mm512_max_epi8(a2n, neg_q2e2_)
                        } else {
                            _mm512_sub_epi8(_mm512_max_epi8(a2n, zero_), qe2_)
                        },
                    );
                    _mm512_storeu_si512(
                        y2.add(t) as *mut __m512i,
                        if FAST {
                            _mm512_max_epi8(b2n, neg_q2e2_)
                        } else {
                            _mm512_sub_epi8(_mm512_max_epi8(b2n, zero_), qe2_)
                        },
                    );
                    // Test-only (CHECK): on every active lane, the operands stay in the bounds the
                    // FAST form relies on, and its values and masks equal the general form's.
                    if CHECK && FAST {
                        let active = u64::MAX >> 63usize.saturating_sub(en - t);
                        let old_tmp = _mm512_sub_epi8(z, q_);
                        let old_tmp2 = _mm512_sub_epi8(z, q2_);
                        let old_values = [
                            _mm512_sub_epi8(a, old_tmp),
                            _mm512_sub_epi8(b, old_tmp),
                            _mm512_sub_epi8(a2, old_tmp2),
                            _mm512_sub_epi8(b2, old_tmp2),
                        ];
                        let shifted = [an, bn, a2n, b2n];
                        let masks = [ma, mb, ma2, mb2];
                        for k in 0..4 {
                            let low = _mm512_set1_epi8(if k < 2 { -10 } else { -9 });
                            let high = _mm512_set1_epi8(if k < 2 { 4 } else { 24 });
                            let bad = _mm512_cmpgt_epi8_mask(low, old_values[k])
                                | _mm512_cmpgt_epi8_mask(old_values[k], high);
                            assert_eq!(bad & active, 0, "active transformed-gap bound k{k}");
                            let penalty = if k < 2 { qe_ } else { qe2_ };
                            let floor = if k < 2 { neg_qe_ } else { neg_q2e2_ };
                            let original =
                                _mm512_sub_epi8(_mm512_max_epi8(old_values[k], zero_), penalty);
                            let shortened = _mm512_max_epi8(shifted[k], floor);
                            assert_eq!(
                                (!_mm512_cmpeq_epi8_mask(original, shortened)) & active,
                                0,
                                "active clipped gap k{k}"
                            );
                            let mask = if RIGHT {
                                !_mm512_cmpgt_epi8_mask(zero_, old_values[k])
                            } else {
                                _mm512_cmpgt_epi8_mask(old_values[k], zero_)
                            };
                            assert_eq!(
                                (mask ^ masks[k]) & active,
                                0,
                                "active continuation tie k{k}"
                            );
                        }
                    }
                    if PLANES {
                        store_four_trace_planes!(prow.add(t + 32), stream_trace, ma, mb, ma2, mb2);
                    } else {
                        // Direction occupies bits 0..=2 and each continuation bit is added once, so
                        // the masked additions cannot carry: they equal OR, one instruction each.
                        d = append_continuation_bits_native!(
                            d, ma, mb, ma2, mb2, flag8_, flag16_, flag32_, flag64_
                        );
                        if stream_trace {
                            // p and stride are 64-byte aligned and t - st advances by 64.
                            _mm512_stream_si512(prow.add(t).cast(), d);
                        } else {
                            _mm512_storeu_si512(prow.add(t).cast(), d);
                        }
                    }
                }
                t += 64;
                first = false;
            }
            std::mem::swap(&mut v, &mut v_next);
            std::mem::swap(&mut x, &mut x_next);
            std::mem::swap(&mut x2, &mut x2_next);
            // The last vector may have written up to 63 cells right of en0; cells right of the band
            // must keep their initial values (the scalar kernel never touches them).
            let es = en + 1;
            // Nothing was written past en0 when the band is a whole number of vectors.
            if !(en - st + 1).is_multiple_of(64) {
                _mm512_storeu_si512(u.add(es) as *mut __m512i, neg_qe_);
                _mm512_storeu_si512(v.add(es) as *mut __m512i, neg_qe_);
                _mm512_storeu_si512(x.add(es) as *mut __m512i, neg_qe_);
                _mm512_storeu_si512(y.add(es) as *mut __m512i, neg_qe_);
                _mm512_storeu_si512(x2.add(es) as *mut __m512i, neg_q2e2_);
                _mm512_storeu_si512(y2.add(es) as *mut __m512i, neg_q2e2_);
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
                        // Packed (H, lane, position) keys as in `clean_impl`, 16 at a time over whole
                        // 16-cell blocks (en0 itself was handled above), then the same 8/4/scalar tails.
                        let en16 = st0 + (en1 - st0) / 16 * 16;
                        let mut low16 = _mm512_setr_epi32(
                            32767, 24575, 16383, 8191, 32766, 24574, 16382, 8190, 32765, 24573,
                            16381, 8189, 32764, 24572, 16380, 8188,
                        );
                        let four = _mm512_set1_epi32(4);
                        let mut mk16 = _mm512_set1_epi32(i32::MIN);
                        while tt < en16 {
                            mk16 = update_h16_packed(
                                h_ptr.add(tt as usize),
                                v.add(tt as usize),
                                low16,
                                mk16,
                            );
                            low16 = _mm512_sub_epi32(low16, four);
                            tt += 16;
                        }
                        let mut low = _mm512_castsi512_si256(low16);
                        let two = _mm256_set1_epi32(2);
                        let mut mk = _mm256_max_epi32(
                            _mm512_castsi512_si256(mk16),
                            _mm512_extracti64x4_epi64::<1>(mk16),
                        );
                        while tt < en8 {
                            let hp = h_ptr.add(tt as usize);
                            let vv = _mm256_cvtepi8_epi32(_mm_loadl_epi64(
                                v.add(tt as usize) as *const __m128i
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
                                std::ptr::read_unaligned(v.add(tt as usize) as *const i32),
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
                                *h_ptr.add(pos) += *v.add(pos) as i8 as i32;
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
                    // Move towards this anti-diagonal's best cell once the window is full: the next
                    // start advances iff that cell lies right of the window's centre.
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

            last_st = st0;
            last_en = en0;
        }

        // Non-temporal stores must be fenced before the traceback reads them or AlignedMemory
        // returns the buffer to its cache, including after early z-drop/band-exhaustion exits.
        if stream_trace {
            _mm_sfence();
        }
        if approx_max && result.score == NEG_INF {
            result.score = result.max;
        }
        if with_cigar {
            if PLANES {
                traceback_dual_affine_bitplanes64(
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
            } else {
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
}
