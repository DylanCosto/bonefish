//! Differential tests for the NEON diagonal-bucket path and its contiguous fallback.
use super::*;
use crate::align::chain::chain_anchors_scalar;

fn params(bandwidth: i32, max_chain_skip: i32) -> ChainingParams {
    ChainingParams {
        min_cnt: 3, min_chain_score: 40, max_gap: 5000, max_gap_ref: -1,
        max_dist_x: 5000, max_dist_y: 5000, bandwidth, bandwidth_long: 20000,
        max_chain_skip, max_chain_iter: 5000, chn_pen_gap: 0.152, chn_pen_skip: 0.0,
        chain_gap_scale: 0.8, rmq_rescue_size: 1000, rmq_rescue_ratio: 0.1,
        rmq_inner_dist: 1000, rmq_size_cap: 100000, seg_rmq: false,
        rescue_join: false, rescue_max_density: 0,
    }
}

fn random(seed: &mut u64, modulus: u64) -> u64 {
    *seed ^= *seed << 13;
    *seed ^= *seed >> 7;
    *seed ^= *seed << 17;
    *seed % modulus
}

fn assert_matches_contiguous(
    anchors: &[Minimizer],
    opt: &ChainingParams,
    bucket_ctx: &mut ChainingBuffers,
    contiguous_ctx: &mut ChainingBuffers,
) {
    let (mut a1, mut a2) = (anchors.to_vec(), anchors.to_vec());
    let (u1, ch1) = unsafe {
        chain_anchors_neon_impl::<true, true, false>(opt, opt.max_dist_x, opt.max_dist_y, &mut a1, bucket_ctx)
    };
    let (u2, ch2) = unsafe {
        chain_anchors_neon_impl::<false, false, false>(opt, opt.max_dist_x, opt.max_dist_y, &mut a2, contiguous_ctx)
    };
    assert_eq!(u1, u2, "chain scores/counts, bandwidth {}", opt.bandwidth);
    assert_eq!(ch1, ch2, "chain anchor order, bandwidth {}", opt.bandwidth);
    let n = anchors.len();
    assert_eq!(&bucket_ctx.scores[..n], &contiguous_ctx.scores[..n], "DP scores");
    assert_eq!(&bucket_ctx.predecessors[..n], &contiguous_ctx.predecessors[..n], "DP predecessor ties");
    // The scalar chainer always applies the skip heuristic, including when skip <= 0.
    // For positive skip, it is an independent oracle for both NEON scan strategies.
    if opt.max_chain_skip > 0 {
        let mut scalar_ctx = ChainingBuffers::new();
        let mut a3 = anchors.to_vec();
        let (u3, ch3) = chain_anchors_scalar(
            opt, opt.max_dist_x, opt.max_dist_y, &mut a3, &mut scalar_ctx,
        );
        assert_eq!(u1, u3, "scalar chains");
        assert_eq!(ch1, ch3, "scalar anchor order");
        assert_eq!(&bucket_ctx.scores[..n], &scalar_ctx.scores[..n], "scalar DP scores");
        assert_eq!(&bucket_ctx.predecessors[..n], &scalar_ctx.predecessors[..n], "scalar predecessor ties");
    }
}

fn tandem_anchors() -> Vec<Minimizer> {
    let mut seed = 0x2545_f491_4f6c_dd1d;
    let mut anchors = Vec::new();
    let mut q = 0;
    while q < 20_000 {
        for _ in 0..30 {
            let shift = random(&mut seed, 400) * 37;
            anchors.push(Minimizer { x: 1_000_000 + q + shift, y: (19 << 32) | q });
        }
        if random(&mut seed, 4) == 0 {
            anchors.push(Minimizer { x: 1_000_000 + q + 900, y: (19 << 32) | q });
        }
        if random(&mut seed, 3) == 0 {
            anchors.push(Minimizer { x: 1_000_000 + random(&mut seed, 40_000), y: (19 << 32) | q });
        }
        q += 7 + random(&mut seed, 6);
    }
    anchors.sort_by_key(|a| (a.x, a.y));
    anchors.dedup();
    anchors
}

#[test]
fn test_neon_diag_buckets_dense() {
    let anchors = tandem_anchors();
    assert!(anchors.len() > 4 * DIAG_BUCKET_MIN_N);
    let (mut buckets, mut contiguous) = (ChainingBuffers::new(), ChainingBuffers::new());
    for (bw, skip) in [(500, 25), (100, 10), (5000, 25)] {
        assert_matches_contiguous(&anchors, &params(bw, skip), &mut buckets, &mut contiguous);
        assert_eq!(buckets.diag_buckets.id.len(), anchors.len());
        let nb = buckets.diag_buckets.nbr.len();
        assert_eq!(buckets.diag_buckets.lists[..nb].iter().map(Vec::len).sum::<usize>(), anchors.len());
    }
}

#[test]
fn test_neon_diag_buckets_negative_diagonals_and_reuse() {
    let mut anchors = Vec::new();
    // Neighboring buckets on either side of zero and of power-of-two boundaries;
    // three target/strand groups also share bucket lists across target transitions.
    for rid_strand in [0u64, 1, 0x8000_0000] {
        for q in 8000..8256u64 {
            for diagonal in [-129i64, -128, -127, -1, 0, 1, 127, 128, 129, 4000] {
                let r = (q as i64 + diagonal) as u64;
                anchors.push(Minimizer { x: (rid_strand << 32) | r, y: (19 << 32) | q });
            }
        }
    }
    anchors.sort_by_key(|a| (a.x, a.y));
    let (mut buckets, mut contiguous) = (ChainingBuffers::new(), ChainingBuffers::new());
    for bw in [1, 127, 128, 129, 500] {
        assert_matches_contiguous(&anchors, &params(bw, 1), &mut buckets, &mut contiguous);
    }
    // Check restoration even when backtracking returns no chains, then exercise a
    // small read, an empty read, and the dense input again with the same scratch.
    let mut no_chains = params(128, 25);
    no_chains.min_cnt = anchors.len() as i32 + 1;
    assert_matches_contiguous(&anchors, &no_chains, &mut buckets, &mut contiguous);
    assert_matches_contiguous(&anchors[..100], &params(128, 25), &mut buckets, &mut contiguous);
    assert_matches_contiguous(&[], &params(128, 25), &mut buckets, &mut contiguous);
    assert_matches_contiguous(&anchors, &params(128, 25), &mut buckets, &mut contiguous);
}

#[test]
fn test_neon_diag_buckets_window_and_activation_boundaries() {
    let mut anchors = Vec::new();
    // Exactly 1024 anchors over a 512-base query interval: both activation thresholds.
    for q in 0..512u64 {
        for shift in [0, 1000] {
            anchors.push(Minimizer { x: 100_000 + q + shift, y: (19 << 32) | q });
        }
    }
    anchors.sort_by_key(|a| (a.x, a.y));
    let (mut buckets, mut contiguous) = (ChainingBuffers::new(), ChainingBuffers::new());
    for max_iter in [127, 128, 129, 5000] {
        let mut opt = params(100, 10);
        opt.max_chain_iter = max_iter;
        assert_matches_contiguous(&anchors, &opt, &mut buckets, &mut contiguous);
    }
    for skip in [0, -1] {
        assert_matches_contiguous(&anchors, &params(100, skip), &mut buckets, &mut contiguous);
    }
    assert_matches_contiguous(&anchors[..1023], &params(100, 25), &mut buckets, &mut contiguous);
    let sparse: Vec<_> = (0..1100u64).map(|q| Minimizer {
        x: 100_000 + q * 20, y: (19 << 32) | (q * 20),
    }).collect();
    assert_matches_contiguous(&sparse, &params(500, 25), &mut buckets, &mut contiguous);
}

#[test]
fn test_neon_diag_buckets_randomized() {
    let (mut buckets, mut contiguous) = (ChainingBuffers::new(), ChainingBuffers::new());
    for case in 0..24u64 {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(case + 1);
        let mut anchors = Vec::new();
        for q in (0..1200u64).step_by(6) {
            for _ in 0..20 {
                let rid_strand = random(&mut seed, 3);
                let r = 100_000 + q + random(&mut seed, 12000);
                let span = 9 + random(&mut seed, 22);
                anchors.push(Minimizer { x: (rid_strand << 32) | r, y: (span << 32) | q });
            }
        }
        anchors.sort_by_key(|a| (a.x, a.y));
        let mut opt = params([1, 16, 31, 100, 500, 2000][case as usize % 6], [1, 10, 25][case as usize % 3]);
        opt.chn_pen_gap = [0.0, 0.152, 0.5][case as usize % 3];
        opt.chn_pen_skip = [0.0, 0.05][case as usize % 2];
        opt.max_chain_iter = [129, 500, 5000][case as usize % 3];
        assert_matches_contiguous(&anchors, &opt, &mut buckets, &mut contiguous);
    }
}

#[test]
fn test_neon_diag_buckets_signed_coordinate_boundary() {
    for query_boundary in [false, true] {
        let mut anchors = Vec::new();
        for q in 0..64u64 {
            for shift in 0..32u64 {
                let (r, query) = if query_boundary {
                    (100_000 + q + shift, i32::MAX as u64 - 40 + q)
                } else {
                    (i32::MAX as u64 - 40 + q + shift, q)
                };
                anchors.push(Minimizer { x: r, y: (19 << 32) | query });
            }
        }
        anchors.sort_by_key(|a| (a.x, a.y));
        let (mut buckets, mut contiguous) = (ChainingBuffers::new(), ChainingBuffers::new());
        assert_matches_contiguous(&anchors, &params(100, 25), &mut buckets, &mut contiguous);
        assert!(buckets.diag_buckets.id.is_empty(), "signed-coordinate wrap must use the contiguous fallback");
    }
}

#[test]
fn test_neon_masked_scan_sparse_anchors() {
    let (mut buckets, mut contiguous) = (ChainingBuffers::new(), ChainingBuffers::new());
    for n in [31, 32, 128, 1023, 1024, 4096] {
        let mut seed = 0x2545_f491_4f6c_dd1d;
        let mut anchors: Vec<_> = (0..n as u64).map(|j| {
            let q = j * 13;
            let shift = if j % 5 == 0 { random(&mut seed, 20000) } else { 0 };
            Minimizer { x: 100_000 + q + shift, y: (19 << 32) | q }
        }).collect();
        anchors.sort_by_key(|a| (a.x, a.y));
        for (bw, skip) in [(100, 1), (500, 25), (5000, 25), (500, 0)] {
            assert_matches_contiguous(&anchors, &params(bw, skip), &mut buckets, &mut contiguous);
        }
    }
}
