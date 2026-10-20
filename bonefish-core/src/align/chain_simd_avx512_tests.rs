//! Differential tests for the AVX-512 diagonal-bucket path and its contiguous fallback.
use super::*;
use crate::align::chain::chain_anchors_scalar;

#[test]
fn avx512_bucket_uniform_spans_match_scalar() {
    if !available() {return}
    let mut seed=0x6781_1234_12ab_cdefu64;
    for n in 0..258 {
        for pass in 0..8 {
            let first=random(&mut seed,1u64<<32) as i32;
            let mut v=vec![first;n];
            if n>0 && pass%2==0 {let at=random(&mut seed,n as u64) as usize;v[at]=first.wrapping_add(1);}
            let expected=if n>0 && v.iter().all(|&x|x==v[0]) {v[0]} else {-1};
            assert_eq!(unsafe{bucket_uniform_span(&v)},expected,"n={n} pass={pass}");
        }
    }
}

#[test]
fn avx512_bucket_direct_matches_original_scores_and_tail_conversions() {
    if !available() {return}
    let mut seed=0x8765_12ab_3456_78efu64;
    for case in 0..2600 {
        let n=96usize;let m=case%65;
        let full=case%3==0;
        let qi=if full {random(&mut seed,1u64<<32) as i32} else {5000};
        let ri=if full {random(&mut seed,1u64<<32) as i32} else {10000};
        let rid=random(&mut seed,1u64<<32) as u32;
        let q:Vec<_>=(0..n).map(|_|random(&mut seed,if full {1u64<<32} else {12000}) as i32).collect();
        let r:Vec<_>=(0..n).map(|_|random(&mut seed,if full {1u64<<32} else {16000}) as i32).collect();
        let uniform=case%2==0;
        let span:Vec<_>=(0..n).map(|_|if uniform {19} else {random(&mut seed,256) as i32}).collect();
        let diag:Vec<_>=r.iter().zip(&q).map(|(&r,&q)|r.wrapping_sub(q)).collect();
        let idx:Vec<_>=(0..m).map(|_|random(&mut seed,n as u64) as u32).collect();
        let gq:Vec<_>=idx.iter().map(|&k|q[k as usize]).collect();
        let gr:Vec<_>=idx.iter().map(|&k|r[k as usize]).collect();
        let gs:Vec<_>=idx.iter().map(|&k|span[k as usize]).collect();
        let gids=vec![rid;m];
        let dx=[0,1,100,5000,i32::MAX][case%5];
        let dy=[i32::MAX,5000,100,1,0][case%5];
        let bw=[-1,0,50,100,1000,i32::MAX][case%6];
        let gap=[0.0,-0.0,0.152,-0.152,f32::NAN,f32::INFINITY,f32::NEG_INFINITY,f32::MAX][case%8];
        let skip=[0.0,-0.0,0.05,-0.05,f32::NAN,f32::INFINITY][case%6];
        let mut expected=vec![1234567;m+8];let mut actual=expected.clone();
        unsafe {
            compute_chain_scores_batch_avx512(qi,ri,rid,&gr,&gq,&gs,&gids,dx,dy,bw,gap,skip,&mut expected,0,m);
            chain_scores_bucket_direct_avx512::<false>(qi,ri,rid,&q,&span,&diag,&idx,dx,dy,bw,gap,skip,&mut actual,-1);
        }
        assert_eq!(actual,expected,"direct case={case} m={m} full={full}");
        if uniform {
            actual.fill(1234567);
            unsafe {chain_scores_bucket_direct_avx512::<true>(qi,ri,rid,&q,&span,&diag,&idx,dx,dy,bw,gap,skip,&mut actual,19);}
            assert_eq!(actual,expected,"uniform case={case} m={m} full={full}");
        }
    }
}

fn available() -> bool {
    is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw") && is_x86_feature_detected!("avx2")
}

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

#[test]
fn avx512_bucket_gathers_preserve_fields_and_visit_order() {
    if !available() { return; }
    let mut seed = 0x9e37_1234_5678_abcd;
    let refs: Vec<_> = (0..257).map(|_| random(&mut seed, 1u64 << 32) as i32).collect();
    let qry: Vec<_> = (0..257).map(|_| random(&mut seed, 1u64 << 32) as i32).collect();
    let span: Vec<_> = (0..257).map(|_| random(&mut seed, 1u64 << 32) as i32).collect();
    let rid: Vec<_> = (0..257).map(|_| random(&mut seed, 1u64 << 32) as u32).collect();
    let mut bk = DiagBuckets::default();
    // Vary chunk size up and down to exercise reuse, SIMD tails and repeated indices.
    for pass in 0..256 {
        let n = (pass * 37) % 65;
        bk.g_idx = (0..n).map(|_| random(&mut seed, 257) as u32).collect();
        unsafe { gather_bucket_fields_avx512(&mut bk, &refs, &qry, &span, &rid); }
        assert_eq!(bk.g_ref, bk.g_idx.iter().map(|&i| refs[i as usize]).collect::<Vec<_>>());
        assert_eq!(bk.g_qry, bk.g_idx.iter().map(|&i| qry[i as usize]).collect::<Vec<_>>());
        assert_eq!(bk.g_span, bk.g_idx.iter().map(|&i| span[i as usize]).collect::<Vec<_>>());
        assert_eq!(bk.g_rid, bk.g_idx.iter().map(|&i| rid[i as usize]).collect::<Vec<_>>());
    }
}

#[test]
fn avx512_zero_diagonal_penalties_match_avx2() {
    if !available() { return; }
    let q: Vec<_> = (0..96).map(|i| 10000 - i * 7 - 1).collect();
    for crossing in [false, true] {
    let r: Vec<_> = q.iter().map(|&v| {
        let rd = (10000i32 - v).wrapping_add(if crossing { i32::MIN } else { 0 });
        80000i32.wrapping_sub(rd)
    }).collect();
    let span: Vec<_> = (0..96).map(|i| [1, 19, 255][i % 3]).collect();
    let rid: Vec<_> = (0..96).map(|i| if i % 5 == 0 { 1 } else { 0 }).collect();
    let diag: Vec<_> = r.iter().zip(&q).map(|(&r,&q)| r.wrapping_sub(q)).collect();
    for gap in [0.152, -0.152, f32::MAX, f32::INFINITY, f32::NAN] {
        for skip in [0.0, -0.0, 0.05, -0.05] {
            for start in [0, 1, 7] { for len in 0..=64 {
                let mut expected = vec![1234567; 96];
                let mut full = expected.clone();
                let mut masked = expected.clone();
                let bits = unsafe {
                    compute_chain_scores_batch_avx2(10000,80000,0,&r,&q,&span,&rid,5000,5000,100,gap,skip,&mut expected,start,start+len);
                    compute_chain_scores_batch_avx512(10000,80000,0,&r,&q,&span,&rid,5000,5000,100,gap,skip,&mut full,start,start+len);
                    chain_scores_masked_avx512(10000,80000,0,70000,&r,&q,&span,&rid,&diag,5000,5000,100,gap,skip,&mut masked,start,start+len)
                };
                assert_eq!(full, expected, "gap={gap} skip={skip} start={start} len={len}");
                for i in start..start+len {
                    if bits & (1u64 << (i-start)) != 0 { assert_eq!(masked[i],expected[i]); }
                }
            }}
        }
    }
    }
}


#[test]
fn avx512_nonimproving_blocks_match_scalar_visits_and_stops() {
    if !available() || !is_x86_feature_detected!("bmi2") || !is_x86_feature_detected!("popcnt") { return; }
    let mut seed = 0x7a12_3456_89ab_cdefu64;
    for case in 0..12000 {
        let n = 96usize;
        let start = 16 + random(&mut seed, 64) as usize;
        let mask = if case % 5 == 0 { u16::MAX } else { random(&mut seed, 65536) as u16 };
        let best = [0, 20, 1000, i32::MAX][case % 4];
        let limit = [0, 1, 2, 7, 16, 25, 100, i32::MAX - 32][case % 8];
        let initial = random(&mut seed, limit as u64 + 1) as i32;
        let stamp = 999;
        let sc: Vec<_> = (0..n).map(|j| if j % 11 == 0 { i32::MIN } else if case % 4 == 3 {
            random(&mut seed, 1u64 << 32) as i32
        } else { random(&mut seed, 200) as i32 - 100 }).collect();
        let scores: Vec<_> = (0..n).map(|_| if case % 4 == 3 { random(&mut seed, 1u64 << 32) as i32 }
            else { random(&mut seed, 200) as i32 }).collect();
        let pre: Vec<_> = (0..n).map(|j| if j == 0 || case % 5 == 0 { -1 } else if case % 3 == 0 {
            j as i32 - 1
        } else if case % 7 == 0 { 0 } else { random(&mut seed, j as u64) as i32 }).collect();
        let visited: Vec<_> = (0..n).map(|_| if case % 6 == 0 || random(&mut seed, 3) == 0 { stamp } else { -1 }).collect();
        let improving = (0..16).any(|bit| mask & (1 << bit) != 0 && sc[start+bit] != i32::MIN
            && sc[start+bit].wrapping_add(scores[start+bit]) > best);
        let mut actual = visited.clone();
        let mut skip = initial;
        let got = unsafe { advance_nonimproving_block16(sc.as_ptr(),scores.as_ptr(),pre.as_ptr(),actual.as_mut_ptr(),
            start,mask,best,&mut skip,limit,stamp) };
        if improving {
            assert_eq!(got,BlockWalk::Improving,"case {case}");
            assert_eq!(actual,visited,"improving block must have no side effects");assert_eq!(skip,initial);
            continue;
        }
        let mut expected=visited.clone();let mut es=initial;let mut stop=None;
        for bit in (0..16).rev() {
            if mask & (1 << bit)==0 { continue; }
            let j=start+bit;if sc[j]==i32::MIN { continue; }
            if expected[j]==stamp { es+=1;if es>limit { stop=Some(j);break; } }
            if pre[j]>=0 { expected[pre[j] as usize]=stamp; }
        }
        let want=if let Some(j)=stop {BlockWalk::Stopped(j)} else {BlockWalk::Advanced};
        assert_eq!(got,want,"stop case {case} mask={mask:#x} initial={initial} limit={limit}");
        assert_eq!(skip,es,"skip case {case}");assert_eq!(actual,expected,"visited case {case}");
    }
}

fn assert_matches_contiguous(
    anchors: &[Minimizer],
    opt: &ChainingParams,
    bucket_ctx: &mut ChainingBuffers,
    contiguous_ctx: &mut ChainingBuffers,
) {
    let (mut a1, mut a2) = (anchors.to_vec(), anchors.to_vec());
    let (u1, ch1) = unsafe {
        chain_anchors_avx512_impl::<true, true>(opt, opt.max_dist_x, opt.max_dist_y, &mut a1, bucket_ctx)
    };
    let (u2, ch2) = unsafe {
        chain_anchors_avx512_impl::<false, false>(opt, opt.max_dist_x, opt.max_dist_y, &mut a2, contiguous_ctx)
    };
    assert_eq!(u1, u2, "chain scores/counts, bandwidth {}", opt.bandwidth);
    assert_eq!(ch1, ch2, "chain anchor order, bandwidth {}", opt.bandwidth);
    let n = anchors.len();
    if is_x86_feature_detected!("bmi2") && is_x86_feature_detected!("popcnt") {
        for fast in [0,1,2,3,4,5,6] {
            let mut input=anchors.to_vec();let mut c=ChainingBuffers::new();
            let (u,ch)=unsafe {
                if fast==5 {chain_anchors_avx512_bucket_chunk16::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c)}
                else if fast==6 {chain_anchors_avx512_bucket_chunk32::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c)}
                else if fast==4 {chain_anchors_avx512_bucket_diag32::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c)}
                else if fast==3 {chain_anchors_avx512_bucket_bulk::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c)}
                else if fast==2 {chain_anchors_avx512_bucket_direct::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c)}
                else if fast==1 { chain_anchors_avx512_block_walk::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c) }
                else { chain_anchors_avx512_feature_control::<true,true>(opt,opt.max_dist_x,opt.max_dist_y,&mut input,&mut c) }
            };
            assert_eq!(u1,u,"block-walk chain descriptors fast={fast}");assert_eq!(ch1,ch,"block-walk anchor order fast={fast}");
            assert_eq!(&bucket_ctx.scores[..n],&c.scores[..n],"block-walk scores fast={fast}");
            assert_eq!(&bucket_ctx.predecessors[..n],&c.predecessors[..n],"block-walk ties fast={fast}");
            assert_eq!(&bucket_ctx.visited[..n],&c.visited[..n],"block-walk visited stamps fast={fast}");
        }
    }
    assert_eq!(&bucket_ctx.scores[..n], &contiguous_ctx.scores[..n], "DP scores");
    assert_eq!(&bucket_ctx.predecessors[..n], &contiguous_ctx.predecessors[..n], "DP predecessor ties");
    if anchors.iter().all(|a| a.ref_pos() >= 0 && a.query_pos() >= 0) {
        let mut input = anchors.to_vec();
        let mut avx2_ctx = ChainingBuffers::new();
        let (u, ch) = unsafe { chain_anchors_avx2(opt, opt.max_dist_x, opt.max_dist_y, &mut input, &mut avx2_ctx) };
        assert_eq!(u1, u, "AVX2 chain descriptors");
        assert_eq!(ch1, ch, "AVX2 anchor order");
        assert_eq!(&bucket_ctx.scores[..n], &avx2_ctx.scores[..n], "AVX2 DP scores");
        assert_eq!(&bucket_ctx.predecessors[..n], &avx2_ctx.predecessors[..n], "AVX2 predecessor ties");
    }
    // The scalar chainer always applies the skip heuristic, including when skip <= 0.
    // For positive skip, it is an independent oracle for both AVX-512 scan strategies.
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
fn test_avx512_diag_buckets_dense() {
    if !available() { eprintln!("Skipping AVX512 chaining: CPU support unavailable"); return; }
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
fn test_avx512_diag_buckets_negative_diagonals_and_reuse() {
    if !available() { eprintln!("Skipping AVX512 chaining: CPU support unavailable"); return; }
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
fn test_avx512_diag_buckets_window_and_activation_boundaries() {
    if !available() { eprintln!("Skipping AVX512 chaining: CPU support unavailable"); return; }
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
fn test_avx512_diag_buckets_randomized() {
    if !available() { eprintln!("Skipping AVX512 chaining: CPU support unavailable"); return; }
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
fn test_avx512_diag_buckets_signed_coordinate_boundary() {
    if !available() { eprintln!("Skipping AVX512 chaining: CPU support unavailable"); return; }
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
fn test_avx512_masked_scan_sparse_anchors() {
    if !available() { eprintln!("Skipping AVX512 chaining: CPU support unavailable"); return; }
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

#[test]
fn avx512_bucket_merge_runs_preserve_indices_and_cursors() {
    if !available() || !is_x86_feature_detected!("popcnt") {return}
    let mut seed=0x7123_4567_89ab_cdefu64;let mut bulk_calls=0;
    for case in 0..2400 {
        let n=96+random(&mut seed,416) as usize;
        let mut lists:[Vec<u32>;3]=std::array::from_fn(|_|Vec::new());
        let run=1+case%48;
        for j in 0..n {
            let k=if case%7==0 {0} else if case%5==0 {random(&mut seed,3) as usize} else {(j/run)%3};
            lists[k].push(j as u32);
        }
        let initial:[usize;3]=std::array::from_fn(|k|if case%4==0 {random(&mut seed,lists[k].len() as u64+1) as usize} else {lists[k].len()});
        let ws=if case%11==0 {0} else {random(&mut seed,n as u64) as i64};
        let di=(random(&mut seed,1u64<<32) as i32) as i64;
        let bw=[-1,0,17,100,1000,i32::MAX][case%6];
        let diag:Vec<i64>=(0..n).map(|_|if case%9==0 {(random(&mut seed,1u64<<32) as i32) as i64}
            else {di+random(&mut seed,3001) as i64-1500}).collect();
        let mut expected:Vec<u32>=lists.iter().zip(initial).flat_map(|(l,c)|l[..c].iter().copied())
            .filter(|&j|j as i64>=ws && (diag[j as usize]-di).abs()<=bw as i64).collect();
        expected.sort_unstable_by(|a,b|b.cmp(a));expected.truncate(64);
        let cutoff=if expected.len()==64 {*expected.last().unwrap() as i64} else {ws};
        let want_cur:[usize;3]=std::array::from_fn(|k|lists[k][..initial[k]].partition_point(|&j|(j as i64)<cutoff));
        let mut cur=initial;let mut out=Vec::with_capacity(64);
        loop {
            if out.len()==64 {break}
            let head:[i64;3]=std::array::from_fn(|k|if cur[k]>0 {lists[k][cur[k]-1] as i64} else {-1});
            let t=if head[0]>=head[1] && head[0]>=head[2] {0} else if head[1]>=head[2] {1} else {2};
            let j=head[t];if j<ws {break}
            let other=match t {0=>head[1].max(head[2]),1=>head[0].max(head[2]),_=>head[0].max(head[1])};
            let floor=other.max(ws-1);
            if j-floor>=16 && unsafe{bucket_filter_run16(lists[t].as_ptr(),&mut cur[t],floor,&diag,di,bw,&mut out)} {
                bulk_calls+=1;continue
            }
            cur[t]-=1;
            if (diag[j as usize]-di).abs()<=bw as i64 {out.push(j as u32);}
        }
        assert_eq!(out,expected,"merge indices case={case}");
        assert_eq!(cur,want_cur,"merge cursors case={case}");
    }
    assert!(bulk_calls>100,"bulk path was not exercised: {bulk_calls}");
}

#[test]
fn avx512_neighbor_bucket_diag32_matches_full_i64_at_extreme_boundaries() {
    if !available() || !is_x86_feature_detected!("popcnt") {return}
    let mut seed=0x5327_8ab1_12cd_567fu64;
    let bands=[0,1,2,17,31,100,1000,1<<29,(1<<29)+1,1<<30];
    let mut calls=0;
    for case in 0..12000 {
        let bw=bands[case%bands.len()];
        let shift=if bw>1 {32-((bw-1) as u32).leading_zeros()} else {0};
        let width=1i64<<shift;
        let key=[-4294967295i64,-3,-1,0,1,3,4294967295][case%7];
        let di=key*width+if case%3==0 {0} else if case%3==1 {width-1} else {random(&mut seed,width as u64) as i64};
        let offset=[-1,0,1][case%3];
        let full:Vec<i64>=(0..80).map(|j|(key+offset)*width+if j%8==0 {0} else if j%8==1 {width-1} else {random(&mut seed,width as u64) as i64}).collect();
        assert!(full.iter().all(|&d|(d-di).abs()<=i32::MAX as i64));
        let narrow:Vec<i32>=full.iter().map(|&d|d as i32).collect();
        let list:Vec<u32>=(0..80).map(|i|i as u32).collect();
        let initial=case%81;let floor=[-1,0,16,64,78][case%5];let prefix=case%50;
        let mut a=(0..prefix as u32).collect::<Vec<_>>();a.reserve(64-a.len());
        let mut b=a.clone();b.reserve(64-b.len());let mut ca=initial;let mut cb=initial;
        let ra=unsafe{bucket_filter_run16(list.as_ptr(),&mut ca,floor,&full,di,bw,&mut a)};
        let rb=unsafe{bucket_filter_run16_diag32(list.as_ptr(),&mut cb,floor,&narrow,di as i32,bw,&mut b)};
        assert_eq!(ra,rb,"accepted case={case}");assert_eq!(ca,cb,"cursor case={case}");assert_eq!(a,b,"indices case={case} band={bw}");
        calls+=usize::from(rb);
    }
    assert!(calls>1000,"narrow filter not exercised: {calls}");
}

fn merge_random(s:&mut u64)->u64 {*s^=*s<<13;*s^=*s>>7;*s^=*s<<17;*s}
#[test]
fn vector_merge_matches_reverse_merge_and_filter() {
 if !available(){return}
 let mut seed=0x8ab3_8917_fba2_0031u64;
 for case in 0..20000 {
  let n=16+(merge_random(&mut seed)%500) as usize;
  let mut lists=[vec![],vec![],vec![]];
  for i in 0..n {lists[(merge_random(&mut seed)%3) as usize].push(i as u32);}
  let mut cur=[lists[0].len(),lists[1].len(),lists[2].len()];let before=cur;
  let ws=(merge_random(&mut seed)%(n+1) as u64) as usize;
  let di=if case%3==0 {i32::MAX-1000} else {0};
  let bw=(merge_random(&mut seed)%2000) as i32;
  let diags:Vec<_>=(0..n).map(|_|di.wrapping_add((merge_random(&mut seed)%6000) as i32-3000)).collect();
  let mut expected=Vec::with_capacity(64);let mut ec=cur;
  for _ in 0..16 {
   let t=(0..3).filter(|&t|ec[t]>0).max_by_key(|&t|lists[t][ec[t]-1]).unwrap();
   ec[t]-=1;let j=lists[t][ec[t]] as usize;
   if j>=ws && diags[j].wrapping_sub(di).abs()<=bw {expected.push(j as u32)}
  }
  let mut actual=Vec::with_capacity(64);
  unsafe{bucket_merge_filter16_diag32([lists[0].as_ptr(),lists[1].as_ptr(),lists[2].as_ptr()],&mut cur,ws,&diags,di,bw,&mut actual)};
  assert_eq!(actual,expected,"case={case} before={before:?}");assert_eq!(cur,ec,"cursor case={case}");
 }
}
