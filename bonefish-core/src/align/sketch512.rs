//! AVX-512 minimizer sketching: the sequence is split into eight intervals, sketched in parallel
//! as eight 64-bit lanes that each run the scalar sketcher's emission state machine.
//!
//! Each lane first warms up on the k+w-1 bases before its interval. A minimizer belongs to the
//! lane in whose interval it is emitted (not the one holding its seed), so concatenating the
//! lanes' output gives exactly the scalar output. Only odd k is supported: odd k-mers cannot be
//! their own reverse complement, so the symmetric k-mers the scalar sketcher skips never occur.
use super::Minimizer;
use core::arch::x86_64::*;
use core::mem::transmute;
use std::cell::RefCell;

/// Per-thread buffers reused across calls; all but `nt` and `output` hold O(w) rows.
#[derive(Default)]
struct Scratch {
    nt: Vec<u8>,
    xs: Vec<[u64; 8]>,
    ys: Vec<[u64; 8]>,
    output: [Vec<Minimizer>; 8],
    block_x: Vec<[u64; 8]>,
    block_y: Vec<[u64; 8]>,
    block_count: Vec<[u64; 8]>,
    block_active: Vec<u8>,
    block_own: Vec<u8>,
    suffix_x: Vec<[u64; 8]>,
    suffix_y: Vec<[u64; 8]>,
    suffix_pos: Vec<[u64; 8]>,
    suffix_dup: Vec<u8>,
}
thread_local! { static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default()); }

#[inline(always)]
fn emit(mut lanes: u8, x: [u64; 8], y: [u64; 8], output: &mut [Vec<Minimizer>; 8]) {
    while lanes != 0 {
        let lane = lanes.trailing_zeros() as usize;
        lanes &= lanes - 1;
        output[lane].push(Minimizer {
            x: x[lane],
            y: y[lane],
        });
    }
}

/// nt4-encodes `seq` into `nt` (U as T), followed by 64 padding bytes of 4 (as for an ambiguous base).
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn encode(seq: &[u8], nt: &mut Vec<u8>) {
    unsafe {
        nt.resize(seq.len() + 64, 4);
        let lower = _mm512_set1_epi8(0x20);
        let four = _mm512_set1_epi8(4);
        let mut i = 0;
        while i + 64 <= seq.len() {
            let v = _mm512_or_si512(_mm512_loadu_si512(seq.as_ptr().add(i).cast()), lower);
            let mut c = four;
            c = _mm512_mask_mov_epi8(
                c,
                _mm512_cmpeq_epi8_mask(v, _mm512_set1_epi8(b'a' as i8)),
                _mm512_setzero_si512(),
            );
            c = _mm512_mask_mov_epi8(
                c,
                _mm512_cmpeq_epi8_mask(v, _mm512_set1_epi8(b'c' as i8)),
                _mm512_set1_epi8(1),
            );
            c = _mm512_mask_mov_epi8(
                c,
                _mm512_cmpeq_epi8_mask(v, _mm512_set1_epi8(b'g' as i8)),
                _mm512_set1_epi8(2),
            );
            c = _mm512_mask_mov_epi8(
                c,
                _mm512_cmpeq_epi8_mask(v, _mm512_set1_epi8(b't' as i8))
                    | _mm512_cmpeq_epi8_mask(v, _mm512_set1_epi8(b'u' as i8)),
                _mm512_set1_epi8(3),
            );
            _mm512_storeu_si512(nt.as_mut_ptr().add(i).cast(), c);
            i += 64;
        }
        for j in i..seq.len() {
            nt[j] = super::encode_base(seq[j]);
        }
        nt[seq.len()..].fill(4);
    }
}

/// Appends the minimizers of `seq[..len]` to `p`, identical to `sketch_sequence_scalar`.
/// The caller ensures AVX-512BW is available; requires odd `k <= 28` and `0 < w < 256`.
#[target_feature(enable = "avx512f,avx512bw")]
pub(super) unsafe fn sketch(
    seq: &[u8],
    len: usize,
    w: usize,
    k: usize,
    rid: usize,
    p: &mut Vec<Minimizer>,
) {
    unsafe {
        if len == 0 {
            return;
        }
        assert!(len <= seq.len() && w > 0 && w < 256 && k > 0 && k <= 28 && k % 2 == 1);
        SCRATCH.with(|cell| {
            let mut scratch = cell.borrow_mut();
            encode(&seq[..len], &mut scratch.nt);
            let Scratch {
                nt,
                xs,
                ys,
                output,
                block_x,
                block_y,
                block_count,
                block_active,
                block_own,
                suffix_x,
                suffix_y,
                suffix_pos,
                suffix_dup,
            } = &mut *scratch;
            block_x.resize(w, [0; 8]);
            block_y.resize(w, [0; 8]);
            block_count.resize(w, [0; 8]);
            block_active.resize(w, 0);
            block_own.resize(w, 0);
            suffix_x.resize(w, [u64::MAX; 8]);
            suffix_x.fill([u64::MAX; 8]);
            suffix_y.resize(w, [u64::MAX; 8]);
            suffix_y.fill([u64::MAX; 8]);
            suffix_pos.resize(w, [0; 8]);
            suffix_dup.resize(w, 0);
            suffix_dup.fill(0);
            xs.resize(w, [0; 8]);
            xs.fill([0; 8]);
            ys.resize(w, [0; 8]);
            ys.fill([0; 8]);
            for out in output.iter_mut() {
                out.clear();
            }
            let size = len.div_ceil(8);
            let overlap = k + w - 1;
            let owned: [u64; 8] = std::array::from_fn(|lane| (lane * size).min(len) as u64);
            let starts: [u64; 8] =
                std::array::from_fn(|lane| owned[lane].saturating_sub(overlap as u64));
            let ends: [u64; 8] = std::array::from_fn(|lane| ((lane + 1) * size).min(len) as u64);
            let owner = _mm512_loadu_si512(owned.as_ptr().cast());
            let end = _mm512_loadu_si512(ends.as_ptr().cast());
            let mut pos = _mm512_loadu_si512(starts.as_ptr().cast());
            let zero = _mm512_setzero_si512();
            let one = _mm512_set1_epi64(1);
            let max = _mm512_set1_epi64(-1);
            let kvec = _mm512_set1_epi64(k as i64);
            let keymask = _mm512_set1_epi64(((1u64 << (2 * k)) - 1) as i64);
            let first_window = _mm512_set1_epi64((w + k - 1) as i64);
            let later_window = _mm512_set1_epi64((w + k) as i64);
            let shift = _mm512_set1_epi64((2 * (k - 1)) as i64);
            let ridvec = _mm512_set1_epi64(((rid as u64) << 32) as i64);
            let mut fw = zero;
            let mut rc = zero;
            let mut count = zero;
            let mut minx = max;
            let mut miny = max;
            let mut minpos = zero;
            let mut step = 0usize;
            let mut input = zero;
            loop {
                // Hash one window-sized block before evaluating its sliding minima.
                // This keeps scratch O(w), independent of read length.
                let mut block_len = 0;
                for row in 0..w {
                    let active = _mm512_cmplt_epu64_mask(pos, end);
                    if active == 0 {
                        break;
                    }
                    // Each gathered word holds eight successive nt4 bases. A lane
                    // stops at its interval end; the nt buffer has 64 padding bytes.
                    if step & 7 == 0 {
                        input =
                            _mm512_mask_i64gather_epi64::<1>(zero, active, pos, nt.as_ptr().cast());
                    }
                    let c = _mm512_and_si512(input, _mm512_set1_epi64(255));
                    let valid = active & _mm512_cmplt_epu64_mask(c, _mm512_set1_epi64(4));
                    let next_fw =
                        _mm512_and_si512(_mm512_or_si512(_mm512_slli_epi64::<2>(fw), c), keymask);
                    let next_rc = _mm512_or_si512(
                        _mm512_srli_epi64::<2>(rc),
                        _mm512_sllv_epi64(_mm512_xor_si512(c, _mm512_set1_epi64(3)), shift),
                    );
                    fw = _mm512_mask_mov_epi64(fw, valid, next_fw);
                    rc = _mm512_mask_mov_epi64(rc, valid, next_rc);
                    count = _mm512_maskz_mov_epi64(valid, _mm512_add_epi64(count, one));
                    let strand = _mm512_cmpge_epu64_mask(fw, rc);
                    let mut bits = _mm512_mask_blend_epi64(strand, fw, rc);
                    bits = _mm512_and_si512(
                        _mm512_add_epi64(
                            _mm512_xor_si512(bits, max),
                            _mm512_slli_epi64::<21>(bits),
                        ),
                        keymask,
                    );
                    bits = _mm512_xor_si512(bits, _mm512_srli_epi64::<24>(bits));
                    bits = _mm512_and_si512(
                        _mm512_add_epi64(
                            _mm512_add_epi64(bits, _mm512_slli_epi64::<3>(bits)),
                            _mm512_slli_epi64::<8>(bits),
                        ),
                        keymask,
                    );
                    bits = _mm512_xor_si512(bits, _mm512_srli_epi64::<14>(bits));
                    bits = _mm512_and_si512(
                        _mm512_add_epi64(
                            _mm512_add_epi64(bits, _mm512_slli_epi64::<2>(bits)),
                            _mm512_slli_epi64::<4>(bits),
                        ),
                        keymask,
                    );
                    bits = _mm512_xor_si512(bits, _mm512_srli_epi64::<28>(bits));
                    bits = _mm512_and_si512(
                        _mm512_add_epi64(bits, _mm512_slli_epi64::<31>(bits)),
                        keymask,
                    );
                    let good = valid & _mm512_cmpge_epu64_mask(count, kvec);
                    let infox = _mm512_mask_mov_epi64(
                        max,
                        good,
                        _mm512_or_si512(
                            _mm512_slli_epi64::<8>(bits),
                            _mm512_min_epu64(count, kvec),
                        ),
                    );
                    let infoy = _mm512_mask_mov_epi64(
                        max,
                        good,
                        _mm512_or_si512(
                            ridvec,
                            _mm512_or_si512(
                                _mm512_slli_epi64::<1>(pos),
                                _mm512_maskz_mov_epi64(strand, one),
                            ),
                        ),
                    );
                    block_x[row] = transmute::<__m512i, [u64; 8]>(infox);
                    block_y[row] = transmute::<__m512i, [u64; 8]>(infoy);
                    block_count[row] = transmute::<__m512i, [u64; 8]>(count);
                    block_active[row] = active;
                    block_own[row] = active & _mm512_cmpge_epu64_mask(pos, owner);
                    input = _mm512_srli_epi64::<8>(input);
                    step += 1;
                    pos = _mm512_add_epi64(pos, one);
                    block_len += 1;
                }
                if block_len == 0 {
                    break;
                }
                let mut prefix_x = max;
                let mut prefix_y = max;
                let mut prefix_pos = zero;
                let mut prefix_dup = 0u8;
                for bp in 0..block_len {
                    let active = block_active[bp];
                    let infox = _mm512_loadu_si512(block_x[bp].as_ptr().cast());
                    let infoy = _mm512_loadu_si512(block_y[bp].as_ptr().cast());
                    let count = _mm512_loadu_si512(block_count[bp].as_ptr().cast());
                    let less = _mm512_cmplt_epu64_mask(infox, prefix_x);
                    let equal = _mm512_cmpeq_epi64_mask(infox, prefix_x);
                    prefix_dup = (prefix_dup & !less) | equal;
                    let take = less | equal;
                    prefix_x = _mm512_mask_mov_epi64(prefix_x, take, infox);
                    prefix_y = _mm512_mask_mov_epi64(prefix_y, take, infoy);
                    prefix_pos =
                        _mm512_mask_mov_epi64(prefix_pos, take, _mm512_set1_epi64(bp as i64));
                    let (sx, sy, sp, sd) = if bp + 1 < w {
                        (
                            _mm512_loadu_si512(suffix_x[bp + 1].as_ptr().cast()),
                            _mm512_loadu_si512(suffix_y[bp + 1].as_ptr().cast()),
                            _mm512_loadu_si512(suffix_pos[bp + 1].as_ptr().cast()),
                            suffix_dup[bp + 1],
                        )
                    } else {
                        (max, max, zero, 0)
                    };
                    let pless = _mm512_cmplt_epu64_mask(prefix_x, sx);
                    let equal = _mm512_cmpeq_epi64_mask(prefix_x, sx);
                    let take = pless | equal;
                    let window_x = _mm512_mask_blend_epi64(take, sx, prefix_x);
                    let window_y = _mm512_mask_blend_epi64(take, sy, prefix_y);
                    let window_pos = _mm512_mask_blend_epi64(take, sp, prefix_pos);
                    let window_dup = (prefix_dup & pless) | (sd & !take) | equal;
                    xs[bp] = transmute::<__m512i, [u64; 8]>(infox);
                    ys[bp] = transmute::<__m512i, [u64; 8]>(infoy);
                    let own = block_own[bp];
                    let old_good = _mm512_cmpneq_epi64_mask(minx, max);
                    let first = own & old_good & _mm512_cmpeq_epi64_mask(count, first_window);
                    if first != 0 {
                        for row in (bp + 1..w).chain(0..bp) {
                            let x = _mm512_loadu_si512(xs[row].as_ptr().cast());
                            let y = _mm512_loadu_si512(ys[row].as_ptr().cast());
                            let lanes = first
                                & _mm512_cmpeq_epi64_mask(x, minx)
                                & _mm512_cmpneq_epi64_mask(y, miny);
                            emit(lanes, xs[row], ys[row], output);
                        }
                    }
                    let new_min = active & _mm512_cmple_epu64_mask(infox, minx);
                    let expired = active
                        & !new_min
                        & _mm512_cmpeq_epi64_mask(minpos, _mm512_set1_epi64(bp as i64));
                    let emit_new =
                        own & new_min & old_good & _mm512_cmpge_epu64_mask(count, later_window);
                    let emit_exp =
                        own & expired & old_good & _mm512_cmpge_epu64_mask(count, first_window);
                    emit(
                        emit_new | emit_exp,
                        transmute::<__m512i, [u64; 8]>(minx),
                        transmute::<__m512i, [u64; 8]>(miny),
                        output,
                    );
                    minx = _mm512_mask_mov_epi64(minx, new_min, infox);
                    miny = _mm512_mask_mov_epi64(miny, new_min, infoy);
                    minpos = _mm512_mask_mov_epi64(minpos, new_min, _mm512_set1_epi64(bp as i64));
                    if expired != 0 {
                        minx = _mm512_mask_mov_epi64(minx, expired, window_x);
                        miny = _mm512_mask_mov_epi64(miny, expired, window_y);
                        minpos = _mm512_mask_mov_epi64(minpos, expired, window_pos);
                        let ready = own
                            & expired
                            & window_dup
                            & _mm512_cmpneq_epi64_mask(minx, max)
                            & _mm512_cmpge_epu64_mask(count, first_window);
                        if ready & window_dup != 0 {
                            for row in (bp + 1..w).chain(0..=bp) {
                                let rx = _mm512_loadu_si512(xs[row].as_ptr().cast());
                                let ry = _mm512_loadu_si512(ys[row].as_ptr().cast());
                                emit(
                                    ready
                                        & _mm512_cmpeq_epi64_mask(rx, minx)
                                        & _mm512_cmpneq_epi64_mask(ry, miny),
                                    xs[row],
                                    ys[row],
                                    output,
                                );
                            }
                        }
                    }
                }
                if block_len < w {
                    break;
                }
                // Reverse suffix minima of this complete block. Strict comparison
                // retains the rightmost occurrence; equal values set a duplicate bit.
                let mut sx = max;
                let mut sy = max;
                let mut sp = zero;
                let mut dup = 0u8;
                for row in (0..w).rev() {
                    let x = _mm512_loadu_si512(block_x[row].as_ptr().cast());
                    let y = _mm512_loadu_si512(block_y[row].as_ptr().cast());
                    let less = _mm512_cmplt_epu64_mask(x, sx);
                    let equal = _mm512_cmpeq_epi64_mask(x, sx);
                    dup = (dup & !less) | equal;
                    sx = _mm512_mask_mov_epi64(sx, less, x);
                    sy = _mm512_mask_mov_epi64(sy, less, y);
                    sp = _mm512_mask_mov_epi64(sp, less, _mm512_set1_epi64(row as i64));
                    suffix_x[row] = transmute::<__m512i, [u64; 8]>(sx);
                    suffix_y[row] = transmute::<__m512i, [u64; 8]>(sy);
                    suffix_pos[row] = transmute::<__m512i, [u64; 8]>(sp);
                    suffix_dup[row] = dup;
                }
            }
            // Only the final stream owns the end-of-sequence flush.
            let x = transmute::<__m512i, [u64; 8]>(minx);
            let y = transmute::<__m512i, [u64; 8]>(miny);
            // Empty trailing ownership intervals must not add duplicate warm-up results.
            let last = (0..8).rev().find(|&lane| owned[lane] < ends[lane]).unwrap();
            if x[last] != u64::MAX {
                output[last].push(Minimizer {
                    x: x[last],
                    y: y[last],
                });
            }
            for out in output.iter() {
                p.extend_from_slice(out);
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rng(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }
    #[test]
    fn avx512_minimizers_match_scalar_and_avx2_all_fields() {
        if !crate::align::dp::simd_test_available("avx512bw", is_x86_feature_detected!("avx512bw"))
        {
            return;
        }
        let mut seed = 0xa3c5_99e1_2d37_b041;
        for case in 0..2200 {
            let k = [1, 3, 5, 9, 15, 19, 23, 27][case % 8];
            let w = [1, 2, 5, 10, 19, 32, 63, 127, 255][case % 9];
            let len = if case < 80 {
                case
            } else if case % 7 == 0 {
                8 * (k + w) + case % 17
            } else {
                rng(&mut seed) as usize % 6000
            };
            let alphabet = if case % 3 == 0 {
                b"ACGTacgtNnx!?".as_slice()
            } else {
                b"ACGT".as_slice()
            };
            let seq: Vec<_> = (0..len)
                .map(|i| {
                    if case % 11 == 0 {
                        b"ATATATAC"[i % 8]
                    } else if case % 13 == 0 {
                        b'A'
                    } else {
                        alphabet[rng(&mut seed) as usize % alphabet.len()]
                    }
                })
                .collect();
            let rid = [0, 1, 12345, u32::MAX as usize][case % 4];
            let mut expected = Vec::new();
            super::super::sketch_sequence_scalar(&seq, len, w, k, rid, &mut expected);
            let sentinel = Minimizer { x: 17, y: 23 };
            let mut actual = vec![sentinel];
            unsafe {
                sketch(&seq, len, w, k, rid, &mut actual);
            }
            assert_eq!(&actual[1..], expected, "case={case} len={len} k={k} w={w}");
            assert_eq!(actual[0], sentinel);
            if is_x86_feature_detected!("avx2") {
                let mut avx2 = Vec::new();
                unsafe {
                    super::super::sketch_chunked_avx2(&seq, len, w, k, rid, &mut avx2);
                }
                assert_eq!(expected, avx2, "AVX2 reference case={case}");
            }
        }
    }
}
