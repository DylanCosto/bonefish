//! Shared masked-score conformance fixtures for native AVX2, AVX-512 and NEON backends.
use super::*;

fn available() -> bool {
    #[cfg(target_arch = "x86_64")]
    { crate::align::dp::simd_test_available("avx2", is_x86_feature_detected!("avx2")) }
    #[cfg(target_arch = "aarch64")]
    { true }
}

#[allow(clippy::too_many_arguments)]
fn check(qi: i32, ri: i32, r: &[i32], q: &[i32], rid: &[u32], start: usize, end: usize, bw: i32) {
    let spans: Vec<_> = (0..r.len()).map(|j| [1, 15, 19, 255][j % 4]).collect();
    let diagonals: Vec<_> = r.iter().zip(q).map(|(&r, &q)| r.wrapping_sub(q)).collect();
    let mut full = vec![1234567; r.len()];
    let mut masked = full.clone();
    let mask;
    unsafe {
        #[cfg(target_arch = "x86_64")]
        {
            compute_chain_scores_batch_avx2(qi, ri, 0, r, q, &spans, rid,
                5000, 5000, bw, 0.152, 0.05, &mut full, start, end);
            mask = chain_scores_masked_avx2(qi, ri, 0, ri.wrapping_sub(qi),
                r, q, &spans, rid, &diagonals, 5000, 5000, bw, 0.152, 0.05,
                &mut masked, start, end);
        }
        #[cfg(target_arch = "aarch64")]
        {
            compute_chain_scores_batch_neon(qi, ri, 0, r, q, &spans, rid,
                5000, 5000, bw, 0.152, 0.05, &mut full, start, end);
            mask = chain_scores_masked_neon(qi, ri, 0, ri.wrapping_sub(qi),
                r, q, &spans, rid, &diagonals, 5000, 5000, bw, 0.152, 0.05,
                &mut masked, start, end);
        }
    }
    let expected = (start..end).fold(0u64, |m, j| m | (((full[j] != i32::MIN) as u64) << (j - start)));
    assert_eq!(mask, expected, "mask start={start} end={end} bandwidth={bw}");
    for j in start..end {
        if mask & (1u64 << (j - start)) != 0 {
            assert_eq!(masked[j], full[j], "score at {j}, bandwidth {bw}");
        }
    }
    for j in (0..start).chain(end..r.len()) {
        assert_eq!(masked[j], 1234567, "wrote outside the score chunk");
    }
    #[cfg(target_arch = "x86_64")]
    if crate::align::dp::simd_test_available("avx512bw",
        is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw"))
    {
        let mut wide_full = vec![1234567; r.len()];
        let mut wide_masked = wide_full.clone();
        let wide_mask = unsafe {
            compute_chain_scores_batch_avx512(qi, ri, 0, r, q, &spans, rid,
                5000, 5000, bw, 0.152, 0.05, &mut wide_full, start, end);
            chain_scores_masked_avx512(qi, ri, 0, ri.wrapping_sub(qi),
                r, q, &spans, rid, &diagonals, 5000, 5000, bw, 0.152, 0.05,
                &mut wide_masked, start, end)
        };
        assert_eq!(wide_full, full, "AVX-512 full scores, start={start} end={end}");
        assert_eq!(wide_mask, expected, "AVX-512 validity mask, start={start} end={end}");
        for j in start..end {
            if wide_mask & (1u64 << (j - start)) != 0 {
                assert_eq!(wide_masked[j], full[j], "AVX-512 score at {j}");
            }
        }
        for j in (0..start).chain(end..r.len()) {
            assert_eq!(wide_masked[j], 1234567, "AVX-512 wrote outside the score chunk");
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
