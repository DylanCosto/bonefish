//! Shared masked-score conformance fixtures for native AVX2, AVX-512 and NEON backends.
use super::*;

fn available() -> bool {
    #[cfg(target_arch = "x86_64")]
    { crate::align::dp::simd_test_available("avx2", is_x86_feature_detected!("avx2")) }
    #[cfg(target_arch = "aarch64")]
    { true }
}

/// Masked kernels against the contiguous kernels, genomic and spliced (CDNA) scoring; the
/// spliced chainer's query-gap limit (`max_dist_y`) is below its bandwidth, so it is checked with
/// a tight one too.
#[allow(clippy::too_many_arguments)]
fn check(qi: i32, ri: i32, r: &[i32], q: &[i32], rid: &[u32], start: usize, end: usize, bw: i32) {
    check_mode::<false>(qi, ri, r, q, rid, start, end, bw, 5000);
    check_mode::<true>(qi, ri, r, q, rid, start, end, bw, 5000);
    check_mode::<true>(qi, ri, r, q, rid, start, end, bw, 150);
}

#[allow(clippy::too_many_arguments)]
fn check_mode<const CDNA: bool>(qi: i32, ri: i32, r: &[i32], q: &[i32], rid: &[u32], start: usize, end: usize,
    bw: i32, max_dist_y: i32) {
    let spans: Vec<_> = (0..r.len()).map(|j| [1, 15, 19, 255][j % 4]).collect();
    let diagonals: Vec<_> = r.iter().zip(q).map(|(&r, &q)| r.wrapping_sub(q)).collect();
    let mut full = vec![1234567; r.len()];
    let mut masked = full.clone();
    let mask;
    unsafe {
        #[cfg(target_arch = "x86_64")]
        {
            compute_chain_scores_batch_avx2::<CDNA>(qi, ri, 0, r, q, &spans, rid,
                5000, max_dist_y, bw, 0.152, 0.05, &mut full, start, end);
            mask = chain_scores_masked_avx2::<CDNA>(qi, ri, 0, ri.wrapping_sub(qi),
                r, q, &spans, rid, &diagonals, 5000, max_dist_y, bw, 0.152, 0.05,
                &mut masked, start, end);
        }
        #[cfg(target_arch = "aarch64")]
        {
            compute_chain_scores_batch_neon::<CDNA>(qi, ri, 0, r, q, &spans, rid,
                5000, max_dist_y, bw, 0.152, 0.05, &mut full, start, end);
            mask = chain_scores_masked_neon::<CDNA>(qi, ri, 0, ri.wrapping_sub(qi),
                r, q, &spans, rid, &diagonals, 5000, max_dist_y, bw, 0.152, 0.05,
                &mut masked, start, end);
        }
    }
    let expected = (start..end).fold(0u64, |m, j| m | (((full[j] != i32::MIN) as u64) << (j - start)));
    let what = format!("cdna={CDNA} start={start} end={end} bandwidth={bw} max_dist_y={max_dist_y}");
    assert_eq!(mask, expected, "mask {what}");
    for j in start..end {
        if mask & (1u64 << (j - start)) != 0 {
            assert_eq!(masked[j], full[j], "score at {j}, {what}");
        }
    }
    for j in (0..start).chain(end..r.len()) {
        assert_eq!(masked[j], 1234567, "wrote outside the score chunk");
    }
    // The spliced kernels may store scores outside the chunk; only the set bits are read. They
    // assume that every candidate is on anchor i's target and strand (a chaining window holds no
    // other), so only those candidates are compared.
    let same_target = (start..end).fold(0u64, |m, j| if rid[j] == 0 { m | (1u64 << (j - start)) } else { m });
    if CDNA {
        let mut spliced = vec![1234567; r.len()];
        let spliced_mask;
        unsafe {
            #[cfg(target_arch = "x86_64")]
            {
                spliced_mask = chain_scores_spliced_avx2(qi, ri, r, q, &spans,
                    5000, max_dist_y, bw, 0.152, 0.05, &mut spliced, start, end);
            }
            #[cfg(target_arch = "aarch64")]
            {
                spliced_mask = chain_scores_spliced_neon(qi, ri, r, q, &spans,
                    5000, max_dist_y, bw, 0.152, 0.05, &mut spliced, start, end);
            }
        }
        assert_eq!(spliced_mask & same_target, expected & same_target, "spliced mask {what}");
        for j in start..end {
            if spliced_mask & same_target & (1u64 << (j - start)) != 0 {
                assert_eq!(spliced[j], full[j], "spliced score at {j}, {what}");
            }
        }
    }
    #[cfg(target_arch = "x86_64")]
    if crate::align::dp::simd_test_available("avx512bw",
        is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw"))
    {
        let mut wide_full = vec![1234567; r.len()];
        let mut wide_masked = wide_full.clone();
        let wide_mask = unsafe {
            compute_chain_scores_batch_avx512::<CDNA>(qi, ri, 0, r, q, &spans, rid,
                5000, max_dist_y, bw, 0.152, 0.05, &mut wide_full, start, end);
            chain_scores_masked_avx512::<CDNA>(qi, ri, 0, ri.wrapping_sub(qi),
                r, q, &spans, rid, &diagonals, 5000, max_dist_y, bw, 0.152, 0.05,
                &mut wide_masked, start, end)
        };
        assert_eq!(wide_full, full, "AVX-512 full scores, {what}");
        assert_eq!(wide_mask, expected, "AVX-512 validity mask, {what}");
        for j in start..end {
            if wide_mask & (1u64 << (j - start)) != 0 {
                assert_eq!(wide_masked[j], full[j], "AVX-512 score at {j}, {what}");
            }
        }
        for j in (0..start).chain(end..r.len()) {
            assert_eq!(wide_masked[j], 1234567, "AVX-512 wrote outside the score chunk");
        }
        if CDNA {
            let mut spliced = vec![1234567; r.len()];
            let spliced_mask = unsafe {
                chain_scores_spliced_avx512(qi, ri, r, q, &spans,
                    5000, max_dist_y, bw, 0.152, 0.05, &mut spliced, start, end)
            };
            assert_eq!(spliced_mask & same_target, expected & same_target, "AVX-512 spliced mask {what}");
            for j in start..end {
                if spliced_mask & same_target & (1u64 << (j - start)) != 0 {
                    assert_eq!(spliced[j], full[j], "AVX-512 spliced score at {j}, {what}");
                }
            }
        }
    }
}

#[test]
fn test_masked_scores_all_chunk_lengths_and_offsets() {
    if !available() { return; }
    let (qi, ri) = (10000, 80000);
    let q: Vec<_> = (0..96i32).map(|j| if j % 11 == 0 { qi + 10 } else { qi - 10 - j * 3 }).collect();
    let r: Vec<_> = q.iter().enumerate().map(|(j, &qj)| {
        let delta = if (j / 8) % 3 == 0 { 10000 } else { (j % 21) as i32 - 10 };
        ri - (qi - qj) + delta
    }).collect();
    let rid = vec![0; q.len()];
    for start in [0, 1, 3, 8] {
        for len in 0..=64 {
            for bw in [0, 1, 19, 100, 500] {
                check(qi, ri, &r, &q, &rid, start, start + len, bw);
            }
        }
    }
}

#[test]
fn test_masked_scores_mixed_targets_in_vector_groups() {
    if !available() { return; }
    let q: Vec<_> = (0..64i32).map(|j| 10000 - j - 1).collect();
    let r: Vec<_> = (0..64i32).map(|j| 80000 - j - 1).collect();
    let rid: Vec<_> = (0..64).map(|j| if j % 3 == 0 { 0x80000001 } else { 0 }).collect();
    for start in [0, 1, 7] {
        for len in 0..=64-start {
            check(10000, 80000, &r, &q, &rid, start, start+len, 100);
        }
    }
}

#[test]
fn test_masked_scores_signed_coordinate_crossings() {
    if !available() { return; }
    for (qi, ri) in [(5000i32, i32::MIN + 32), (i32::MIN + 32, 50000)] {
        let q: Vec<_> = (0..80i32).map(|j| qi.wrapping_sub(j + 1)).collect();
        let r: Vec<_> = (0..80i32).map(|j| ri.wrapping_sub(j + 1)).collect();
        let rid = vec![0; q.len()];
        for len in 0..=64 {
            check(qi, ri, &r, &q, &rid, 1, 1 + len, 100);
        }
    }
}
