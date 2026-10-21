//! Differential tests: the SIMD chainers' spliced (CDNA) instantiations against the scalar
//! spliced chainer (`chain_anchors_cdna`'s reference), on randomized transcript-like anchors.
use super::*;
use crate::align::chain::{chain_anchors_cdna, chain_anchors_scalar_impl, chain_score};

fn random(seed: &mut u64, modulus: u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed % modulus
}

/// The rna/cdna presets' chaining parameters: reference gap and bandwidth 200 kb, query gap 2 kb.
fn spliced_params() -> ChainingParams {
    ChainingParams {
        min_cnt: 3, min_chain_score: 40, max_gap: 2000, max_gap_ref: 200_000,
        max_dist_x: 200_000, max_dist_y: 2000, bandwidth: 200_000, bandwidth_long: 200_000,
        max_chain_skip: 25, max_chain_iter: 5000, chn_pen_gap: 0.112, chn_pen_skip: 0.0,
        chain_gap_scale: 0.8, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1,
        rmq_inner_dist: 1000, rmq_size_cap: 100000, seg_rmq: false,
        rescue_join: false, rescue_max_density: 0,
    }
}

fn anchor(rid_strand: u64, r: u64, span: u64, q: u64) -> Minimizer {
    Minimizer { x: (rid_strand << 32) | r, y: (span << 32) | q }
}

/// Reads of `n_reads` transcripts with exons separated by introns of up to ~250 kb (some beyond
/// the 200 kb reference gap limit), each hitting a few paralogous copies (on several targets and
/// both strands) with indels in the query, plus tandem-repeat hits and random hits; sorted by
/// reference position as the mapper sorts them.
fn transcript_anchors(seed: &mut u64, n_reads: usize) -> Vec<Minimizer> {
    let mut anchors = Vec::new();
    for _ in 0..n_reads {
        let qlen = 200 + random(seed, 4000);
        let copies = 1 + random(seed, 4);
        for _ in 0..copies {
            let rid_strand = random(seed, 6);
            let mut r = 1_000 + random(seed, 50_000_000);
            let mut q = random(seed, 40);
            while q < qlen {
                let exon = 30 + random(seed, 400);
                let end = (q + exon).min(qlen);
                while q < end {
                    if random(seed, 3) > 0 {
                        anchors.push(anchor(rid_strand, r, 10 + random(seed, 10), q));
                    }
                    let step = 1 + random(seed, 12);
                    q += step;
                    r += match random(seed, 40) {
                        0 => step + random(seed, 30),
                        1 => step.saturating_sub(random(seed, 30)),
                        _ => step,
                    };
                }
                r += match random(seed, 8) {
                    0 => 150_000 + random(seed, 100_000),
                    1 | 2 => random(seed, 2_000),
                    _ => 80 + random(seed, 20_000),
                };
                if random(seed, 10) == 0 {
                    q += random(seed, 2_500);
                }
            }
        }
        // Tandem repeat: many shifted copies of each query position.
        if random(seed, 3) == 0 {
            let base = random(seed, 50_000_000);
            let rid_strand = random(seed, 6);
            let mut q = 0;
            while q < qlen.min(1500) {
                for _ in 0..8 {
                    anchors.push(anchor(rid_strand, base + q + random(seed, 60) * 41, 14, q));
                }
                q += 5 + random(seed, 6);
            }
        }
        for _ in 0..random(seed, 200) {
            anchors.push(anchor(random(seed, 6), random(seed, 60_000_000), 10 + random(seed, 10), random(seed, qlen)));
        }
    }
    anchors.sort_by_key(|a| (a.x, a.y));
    anchors
}

fn chain_with(
    chainer: impl Fn(&ChainingParams, i32, i32, &mut [Minimizer], &mut ChainingBuffers) -> (Vec<u64>, Vec<Minimizer>),
    anchors: &[Minimizer], opt: &ChainingParams, ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>, Vec<i32>, Vec<i32>) {
    let mut input = anchors.to_vec();
    let (u, chains) = chainer(opt, opt.max_dist_x, opt.max_dist_y, &mut input, ctx);
    let n = anchors.len();
    (u, chains, ctx.scores[..n].to_vec(), ctx.predecessors[..n].to_vec())
}

/// Chains, DP scores and predecessor ties of every available SIMD chainer (reusing one context
/// across inputs, as the mapper does) against the scalar spliced chainer.
fn assert_simd_matches_scalar(anchors: &[Minimizer], opt: &ChainingParams, ctxs: &mut [ChainingBuffers]) {
    let reference = chain_with(chain_anchors_scalar_impl::<true>, anchors, opt, &mut ChainingBuffers::new());
    let check = |name: &str, got: (Vec<u64>, Vec<Minimizer>, Vec<i32>, Vec<i32>)| {
        let what = format!("{name}, n {}, skip {}, iter {}", anchors.len(), opt.max_chain_skip, opt.max_chain_iter);
        assert_eq!(got.0, reference.0, "chains: {what}");
        assert_eq!(got.1, reference.1, "chain anchors: {what}");
        assert_eq!(got.2, reference.2, "DP scores: {what}");
        assert_eq!(got.3, reference.3, "DP predecessors: {what}");
    };
    let mut checked = 0;
    #[cfg(target_arch = "x86_64")]
    {
        check("SSE2", chain_with(|o, x, y, a, c| unsafe { chain_anchors_sse::<true>(o, x, y, a, c) }, anchors, opt, &mut ctxs[0]));
        checked += 1;
        if crate::align::dp::simd_test_available("avx2", is_x86_feature_detected!("avx2")) {
            check("AVX2", chain_with(|o, x, y, a, c| unsafe { chain_anchors_avx2::<true>(o, x, y, a, c) }, anchors, opt, &mut ctxs[1]));
            checked += 1;
        }
        let avx512 = is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw")
            && is_x86_feature_detected!("avx2") && is_x86_feature_detected!("bmi2") && is_x86_feature_detected!("popcnt");
        if crate::align::dp::simd_test_available("avx512bw", avx512) {
            check("AVX-512", chain_with(|o, x, y, a, c| unsafe { chain_anchors_avx512_impl::<true>(o, x, y, a, c) }, anchors, opt, &mut ctxs[2]));
            checked += 1;
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        check("NEON", chain_with(|o, x, y, a, c| unsafe { chain_anchors_neon_impl::<true, true, true>(o, x, y, a, c) }, anchors, opt, &mut ctxs[0]));
        check("NEON contiguous", chain_with(|o, x, y, a, c| unsafe { chain_anchors_neon_impl::<false, false, true>(o, x, y, a, c) }, anchors, opt, &mut ctxs[1]));
        checked += 2;
    }
    // The public entry point, with the runtime dispatch.
    check("dispatch", chain_with(chain_anchors_cdna, anchors, opt, &mut ctxs[3]));
    assert!(checked > 0);
}

#[test]
fn test_cdna_simd_matches_scalar_randomized() {
    let mut ctxs: Vec<_> = (0..4).map(|_| ChainingBuffers::new()).collect();
    for case in 0..48u64 {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(case + 1);
        let anchors = transcript_anchors(&mut seed, 1 + case as usize % 6);
        let mut opt = spliced_params();
        opt.max_chain_skip = [25, 25, 1, 0, -1, 5][case as usize % 6];
        opt.max_chain_iter = [5000, 5000, 50, 7][case as usize % 4];
        opt.chn_pen_gap = [0.112, 0.12, 0.5, 0.0][case as usize % 4];
        opt.chn_pen_skip = [0.0, 0.0, 0.05][case as usize % 3];
        opt.max_dist_y = [2000, 2000, 300, 200_000][case as usize % 4];
        opt.min_chain_score = [40, 10][case as usize % 2];
        opt.min_cnt = [3, 1][case as usize % 2];
        assert_simd_matches_scalar(&anchors, &opt, &mut ctxs);
    }
}

#[test]
fn test_cdna_simd_small_and_boundary_inputs() {
    let mut ctxs: Vec<_> = (0..4).map(|_| ChainingBuffers::new()).collect();
    let opt = spliced_params();
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let anchors = transcript_anchors(&mut seed, 3);
    // Prefixes: below the dispatch threshold and the vector widths, and the first chunk sizes.
    for n in [0, 1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 200] {
        assert_simd_matches_scalar(&anchors[..n.min(anchors.len())], &opt, &mut ctxs);
    }
    // Gaps at the limits: reference gaps of bandwidth - 1 ..= bandwidth + 1 (with equal and
    // unequal query gaps), query gaps at max_dist_y, and anchors on adjacent targets and strands.
    let mut edge = Vec::new();
    for rid_strand in [0u64, 1, 2, 7] {
        let mut r = 10_000u64;
        let mut q = 0u64;
        for k in 0..60u64 {
            edge.push(anchor(rid_strand, r, 14, q));
            let (dr, dq) = match k % 6 {
                0 => (200_000 + k % 3, 5),
                1 => (199_999, 0),
                2 => (2000 + k % 3, 1999 + k % 3),
                3 => (17, 2000),
                4 => (3, 2001),
                _ => (25, 25),
            };
            r += dr;
            q += dq;
        }
    }
    edge.sort_by_key(|a| (a.x, a.y));
    // Query positions falling along the reference: every later anchor would be a valid
    // predecessor, so a scan reaching past anchor i (a vector lane above it) changes the chains.
    let mut inverted: Vec<_> = (0..40u64).map(|k| anchor(0, 1_000 + 7 * k, 14, 5_000 - 9 * k)).collect();
    inverted.extend((0..40u64).map(|k| anchor(0, 1_003 + 11 * k, 15, 4_000 + 10 * k)));
    inverted.sort_by_key(|a| (a.x, a.y));
    for skip in [25, 0] {
        let mut o = opt.clone();
        o.max_chain_skip = skip;
        o.min_cnt = 1;
        o.min_chain_score = 10;
        assert_simd_matches_scalar(&edge, &o, &mut ctxs);
        for n in [8, 9, 12, 16, 17, 33, 80] {
            assert_simd_matches_scalar(&inverted[..n], &o, &mut ctxs);
        }
    }
}

/// Spliced gap scores of the SIMD kernels, lane by lane, against the scalar `chain_score::<true>`
/// (reference gaps up to 250 kb, query gaps across the 2 kb limit, both gap orders).
#[test]
fn test_cdna_kernel_scores_match_scalar_chain_score() {
    let mut seed = 0x8765_12ab_3456_78efu64;
    for case in 0..400 {
        let n = 80;
        let qi = 300_000 + random(&mut seed, 1000) as i32;
        let ri = 900_000 + random(&mut seed, 1000) as i32;
        let q: Vec<i32> = (0..n).map(|_| qi - random(&mut seed, 2600) as i32 + 10).collect();
        let r: Vec<i32> = (0..n).map(|k| match k % 3 {
            0 => ri - random(&mut seed, 250_000) as i32,
            1 => ri - (qi - q[k]) - random(&mut seed, 40) as i32 + 20,
            _ => ri - random(&mut seed, 3000) as i32 + 10,
        }).collect();
        let span: Vec<i32> = (0..n).map(|_| 10 + random(&mut seed, 246) as i32).collect();
        let rid = vec![3u32; n];
        let (max_dist_x, max_dist_y, bw) = (200_000, [2000, 300][case % 2], [200_000, 5000][case / 2 % 2]);
        let (gap, skip) = ([0.112f32, 0.12, 0.5][case % 3], [0.0f32, 0.05][case % 2]);
        let start = case % 9;
        let end = n - case % 5;
        let ai = Minimizer { x: (3u64 << 32) | ri as u32 as u64, y: (14u64 << 32) | qi as u32 as u64 };
        let expected: Vec<i32> = (0..n).map(|k| {
            let aj = Minimizer { x: (3u64 << 32) | r[k] as u32 as u64, y: ((span[k] as u64) << 32) | q[k] as u32 as u64 };
            chain_score::<true>(&ai, &aj, max_dist_x, max_dist_y, bw, gap, skip)
        }).collect();
        let mut got = vec![1234567; n];
        #[cfg(target_arch = "x86_64")]
        {
            unsafe { compute_chain_scores_batch_sse::<true>(qi, ri, 3, &r, &q, &span, &rid,
                max_dist_x, max_dist_y, bw, gap, skip, &mut got, start, end) };
            assert_eq!(got[start..end], expected[start..end], "SSE2 case {case}");
        }
        #[cfg(target_arch = "x86_64")]
        if crate::align::dp::simd_test_available("avx2", is_x86_feature_detected!("avx2")) {
            unsafe { compute_chain_scores_batch_avx2::<true>(qi, ri, 3, &r, &q, &span, &rid,
                max_dist_x, max_dist_y, bw, gap, skip, &mut got, start, end) };
            assert_eq!(got[start..end], expected[start..end], "AVX2 case {case}");
            for chunk_end in (start + 1..=end).rev().step_by(64) {
                let chunk_start = chunk_end.saturating_sub(64).max(start);
                got.fill(1234567);
                let mask = unsafe { chain_scores_spliced_avx2(qi, ri, &r, &q, &span,
                    max_dist_x, max_dist_y, bw, gap, skip, &mut got, chunk_start, chunk_end) };
                for k in chunk_start..chunk_end {
                    let valid = mask & (1u64 << (k - chunk_start)) != 0;
                    assert_eq!(valid, expected[k] != i32::MIN, "AVX2 spliced mask case {case} k {k}");
                    if valid { assert_eq!(got[k], expected[k], "AVX2 spliced case {case} k {k}"); }
                }
            }
        }
        #[cfg(target_arch = "x86_64")]
        if crate::align::dp::simd_test_available("avx512bw",
            is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw"))
        {
            unsafe { compute_chain_scores_batch_avx512::<true>(qi, ri, 3, &r, &q, &span, &rid,
                max_dist_x, max_dist_y, bw, gap, skip, &mut got, start, end) };
            assert_eq!(got[start..end], expected[start..end], "AVX-512 case {case}");
            for chunk_end in (start + 1..=end).rev().step_by(64) {
                let chunk_start = chunk_end.saturating_sub(64).max(start);
                got.fill(1234567);
                let mask = unsafe { chain_scores_spliced_avx512(qi, ri, &r, &q, &span,
                    max_dist_x, max_dist_y, bw, gap, skip, &mut got, chunk_start, chunk_end) };
                for k in chunk_start..chunk_end {
                    let valid = mask & (1u64 << (k - chunk_start)) != 0;
                    assert_eq!(valid, expected[k] != i32::MIN, "AVX-512 spliced mask case {case} k {k}");
                    if valid { assert_eq!(got[k], expected[k], "AVX-512 spliced case {case} k {k}"); }
                }
            }
        }
        #[cfg(target_arch = "aarch64")]
        {
            unsafe { compute_chain_scores_batch_neon::<true>(qi, ri, 3, &r, &q, &span, &rid,
                max_dist_x, max_dist_y, bw, gap, skip, &mut got, start, end) };
            assert_eq!(got[start..end], expected[start..end], "NEON case {case}");
            for chunk_end in (start + 1..=end).rev().step_by(64) {
                let chunk_start = chunk_end.saturating_sub(64).max(start);
                got.fill(1234567);
                let mask = unsafe { chain_scores_spliced_neon(qi, ri, &r, &q, &span,
                    max_dist_x, max_dist_y, bw, gap, skip, &mut got, chunk_start, chunk_end) };
                for k in chunk_start..chunk_end {
                    let valid = mask & (1u64 << (k - chunk_start)) != 0;
                    assert_eq!(valid, expected[k] != i32::MIN, "NEON spliced mask case {case} k {k}");
                    if valid { assert_eq!(got[k], expected[k], "NEON spliced case {case} k {k}"); }
                }
            }
        }
    }
}
