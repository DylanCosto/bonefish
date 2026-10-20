use super::*;

fn available() -> bool {
    let yes = is_x86_feature_detected!("avx2") && is_x86_feature_detected!("avx512bw");
    if !yes { eprintln!("Skipping AVX512 clean DP tests: CPU support unavailable"); }
    yes
}

fn assert_result(a: &DpResult, b: &DpResult, context: &str) {
    assert_eq!(
        (a.max, a.max_score_query_pos, a.max_score_target_pos, a.max_query_end_score,
         a.max_query_end_target_pos, a.max_target_end_score, a.max_target_end_query_pos,
         a.score, a.reach_end, a.zdropped, a.cigar_len, &a.cigar),
        (b.max, b.max_score_query_pos, b.max_score_target_pos, b.max_query_end_score,
         b.max_query_end_target_pos, b.max_target_end_score, b.max_target_end_query_pos,
         b.score, b.reach_end, b.zdropped, b.cigar_len, &b.cigar),
        "{context}"
    );
}

#[test]
fn avx512_clean_shift_crosses_every_lane_boundary() {
    if !available() { return; }
    let v: [u8; 64] = std::array::from_fn(|i| i as u8);
    let prev: [u8; 64] = std::array::from_fn(|i| (128 + i) as u8);
    let mut got = [0u8; 64];
    unsafe {
        let x = shl1_prev(_mm512_loadu_si512(v.as_ptr().cast()), _mm512_loadu_si512(prev.as_ptr().cast()));
        _mm512_storeu_si512(got.as_mut_ptr().cast(), x);
    }
    let want: [u8; 64] = std::array::from_fn(|i| if i == 0 { prev[63] } else { v[i - 1] });
    assert_eq!(got, want);
}

#[test]
fn avx512_clean_dispatch_keeps_narrow_bands_on_avx2() {
    use super::super::dual_clean::clean_band_cells;
    assert_eq!(clean_band_cells(1000, 1000, 31, 0), 32);
    assert_eq!(clean_band_cells(1000, 1000, 32, 0), 33);
    assert_eq!(clean_band_cells(32, 1000, 500, 0), 32);
    assert_eq!(clean_band_cells(1000, 31, -1, 0), 31);
    assert_eq!(clean_band_cells(1000, 1000, 500, 32), 32);
    assert_eq!(clean_band_cells(1000, 1000, 500, 33), 33);
    assert_eq!(clean_band_cells(1000, 1000, 16, 64), 16);
}

fn matrix(n_score: i8) -> [i8; 25] {
    std::array::from_fn(|i| {
        let (a, b) = (i / 5, i % 5);
        if a == 4 || b == 4 { n_score } else if a == b { 2 } else { -4 }
    })
}

#[allow(clippy::too_many_arguments)]
fn compare(q: &[u8], t: &[u8], mat: &[i8], gaps: (i8, i8, i8, i8), bw: i32,
           adaptive: i32, shift: i32, flags: i32, zdrop: i32, end_bonus: i32, context: &str) {
    let mut a = DpResult::default();
    let mut b = DpResult::default();
    unsafe {
        super::super::dual_clean::extend_dual_affine_avx2_clean(q, t, 5, mat,
            gaps.0, gaps.1, gaps.2, gaps.3, bw, zdrop, end_bonus, flags, adaptive, shift, &mut a);
        extend_dual_affine_avx512_clean(q, t, 5, mat,
            gaps.0, gaps.1, gaps.2, gaps.3, bw, zdrop, end_bonus, flags, adaptive, shift, &mut b);
    }
    assert_result(&a, &b, context);
    // Also verify runtime dispatch, including calls that must stay on AVX2.
    let mut routed = DpResult::default();
    unsafe {
        super::super::dual_clean::extend_dual_affine_clean(q, t, 5, mat,
            gaps.0, gaps.1, gaps.2, gaps.3, bw, zdrop, end_bonus, flags, adaptive, shift, &mut routed);
    }
    assert_result(&a, &routed, context);
}

#[test]
fn avx512_clean_matches_avx2_boundaries_and_ties() {
    if !available() { return; }
    let lengths = [1, 15, 16, 31, 32, 33, 63, 64, 65, 95, 96, 127, 128, 129, 255, 256, 257];
    for &len in &lengths {
        for &bw in &[0, 15, 31, 32, 33, 63, 64, 65, 127, 128, -1] {
            for f in 0..16 {
                let flags = (if f & 1 != 0 { APPROX_MAX | APPROX_DROP } else { 0 })
                    | (if f & 2 != 0 { RIGHT_ALIGN } else { 0 })
                    | (if f & 4 != 0 { SCORE_ONLY } else { 0 })
                    | (if f & 8 != 0 { EXTENSION_ONLY | REV_CIGAR } else { 0 });
                let t = vec![0; len];
                let mut q = vec![0; len + len / 3];
                if len > 32 { q[len / 2] = 1; }
                compare(&q, &t, &matrix(-1), (4, 2, 24, 1), bw, 0, 0, flags, 40, 5,
                    &format!("boundary len={len} bw={bw} flags={flags:#x}"));
            }
        }
    }
}

#[test]
fn avx512_clean_matches_avx2_random_static_shifted_adaptive() {
    if !available() { return; }
    let mut x = 0xBAD5_EED5_1234_6789u64;
    let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
    let widths = [0, 1, 15, 16, 30, 31, 32, 33, 62, 63, 64, 65, 95, 127, 128, 129, 255, -1];
    for case in 0..6000 {
        let len = 1 + (rnd() % 600) as usize;
        let mut t: Vec<u8> = (0..len).map(|_| (rnd() & 3) as u8).collect();
        let mut q = Vec::new();
        for &c in &t {
            match rnd() % 32 {
                0..=2 => q.push((rnd() & 3) as u8),
                3..=4 => {},
                5..=6 => { q.push(c); q.push((rnd() & 3) as u8); },
                _ => q.push(c),
            }
        }
        if q.is_empty() { q.push(0); }
        if case % 7 == 0 { let pos = rnd() as usize % q.len(); q[pos] = 4; }
        if case % 11 == 0 { let pos = rnd() as usize % t.len(); t[pos] = 4; }
        let bw = widths[case % widths.len()];
        let adaptive = if case % 3 == 0 { [1, 16, 31, 32, 33, 63, 64, 65, 128][rnd() as usize % 9] } else { 0 };
        let shift = if adaptive == 0 && case % 4 == 0 { (rnd() % 31) as i32 - 15 } else { 0 };
        let mut flags = 0;
        for flag in [APPROX_MAX, APPROX_DROP, RIGHT_ALIGN, SCORE_ONLY, EXTENSION_ONLY, REV_CIGAR] {
            if rnd() & 1 != 0 { flags |= flag; }
        }
        let zdrop = [-1, 10, 40, 100, 400][rnd() as usize % 5];
        let end_bonus = [-1, 0, 5, 20][rnd() as usize % 4];
        let gaps = [(4, 2, 24, 1), (24, 1, 4, 2), (6, 2, 12, 1), (4, 1, 8, 1)][case % 4];
        let mat = matrix(if case % 5 == 0 { 0 } else { -1 });
        compare(&q, &t, &mat, gaps, bw, adaptive, shift, flags, zdrop, end_bonus,
            &format!("random case={case} qlen={} tlen={} bw={bw} adaptive={adaptive} shift={shift} flags={flags:#x}", q.len(), t.len()));
    }
}

#[test]
fn avx512_register64_long_exact_score_range() {
    if !available() { return; }
    // The combined sequence length exceeds the packed-key guard. Exercise exact tracking
    // through the uncompressed H-score fallback with a register-resident difference state.
    let mut seed=0x8123_4567_89ab_cdefu64;
    let t:Vec<u8>=(0..10001).map(|_|{seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;(seed&3) as u8}).collect();
    let mut q=t.clone();q.insert(3000,4);q.drain(7000..7010);
    for adaptive in [0,64] {
        for flags in [0,RIGHT_ALIGN,SCORE_ONLY,EXTENSION_ONLY|REV_CIGAR] {
            compare(&q,&t,&matrix(-1),(4,2,24,1),63,adaptive,0,flags,400,5,
                &format!("long exact register case adaptive={adaptive} flags={flags:#x}"));
        }
    }
}

/// Replay actual HG001 calls captured with exact kernel parameters. Kept ignored because
/// the external corpus is a local profiling artifact, not a repository test fixture.
#[test]
#[ignore]
fn benchmark_register64_real_calls() {
    if !available() || !is_x86_feature_detected!("avx512vbmi") { return; }
    let narrow = std::env::var_os("BONEFISH_REPLAY_NARROW").is_some();
    let narrow_vl = std::env::var_os("BONEFISH_REPLAY_NARROW_VL").is_some();
    let wide = std::env::var_os("BONEFISH_REPLAY_ARRAY_WIDE").is_some();
    assert!(!(narrow && wide));
    if narrow_vl { assert!(narrow && is_x86_feature_detected!("avx512vl")); }
    let path = std::env::var("BONEFISH_DP_REPLAY").expect("set BONEFISH_DP_REPLAY to the captured cases file");
    let data = std::fs::read(path).unwrap();
    assert_eq!(&data[..8], b"BFDPPR01");
    struct Case { q: Vec<u8>, t: Vec<u8>, mat: [i8;25], gaps: [i8;4], p: [i32;8] }
    let mut cases = Vec::new();
    let mut at = 8;
    while at < data.len() {
        let p: [i32;8] = std::array::from_fn(|i| i32::from_le_bytes(data[at+i*4..at+i*4+4].try_into().unwrap()));
        at += 32;
        let gaps = std::array::from_fn(|i| data[at+i] as i8); at += 4;
        let mat = std::array::from_fn(|i| data[at+i] as i8); at += 25;
        at += 8; // captured CPU nanoseconds
        let q = data[at..at+p[0] as usize].to_vec(); at += p[0] as usize;
        let t = data[at..at+p[1] as usize].to_vec(); at += p[1] as usize;
        let cells = super::super::dual_clean::clean_band_cells(q.len(), t.len(), p[2], p[3]);
        let eligible = if wide { cells > 64 } else if narrow {
            cells <= 32 && p[3] <= 0 && (0..=31).contains(&p[2]) && p[5] & APPROX_MAX != 0
        } else { (33..=64).contains(&cells) };
        if eligible { cases.push(Case{q,t,mat,gaps,p}); }
    }
    assert!(!cases.is_empty());
    fn run(c: &Case, registers: bool, narrow: bool, narrow_vl: bool, wide: bool) -> DpResult {
        let mut ez = DpResult::default();
        let (approx, cigar, right, has_n) = (c.p[5] & APPROX_MAX != 0, c.p[5] & SCORE_ONLY == 0,
            c.p[5] & RIGHT_ALIGN != 0, c.q.iter().chain(&c.t).any(|&x| x >= 4));
        unsafe {
            macro_rules! go {
                ($a:literal,$c:literal,$r:literal,$n:literal) => {
                    if registers && wide {
                        super::clean_impl512::<$a,$c,$r,$n>(&c.q,&c.t,5,&c.mat,
                            c.gaps[0],c.gaps[1],c.gaps[2],c.gaps[3],c.p[2],c.p[6],c.p[7],c.p[5],c.p[3],c.p[4],&mut ez)
                    } else if registers && narrow_vl {
                        super::super::dual_clean_reg32vl::reg32vl_impl::<$a,$c,$r,$n>(&c.q,&c.t,5,&c.mat,
                            c.gaps[0],c.gaps[1],c.gaps[2],c.gaps[3],c.p[2],c.p[6],c.p[7],c.p[5],c.p[3],c.p[4],&mut ez)
                    } else if registers {
                        super::super::dual_clean_reg64::reg64_impl::<$a,$c,$r,$n>(&c.q,&c.t,5,&c.mat,
                            c.gaps[0],c.gaps[1],c.gaps[2],c.gaps[3],c.p[2],c.p[6],c.p[7],c.p[5],c.p[3],c.p[4],&mut ez)
                    } else if narrow {
                        super::super::dual_clean::reg1_impl::<$c,$r,$n>(&c.q,&c.t,5,&c.mat,
                            c.gaps[0],c.gaps[1],c.gaps[2],c.gaps[3],c.p[2],c.p[6],c.p[7],c.p[5],c.p[4],&mut ez)
                    } else if wide {
                        super::super::dual_clean::extend_dual_affine_avx2_clean(&c.q,&c.t,5,&c.mat,
                            c.gaps[0],c.gaps[1],c.gaps[2],c.gaps[3],c.p[2],c.p[6],c.p[7],c.p[5],c.p[3],c.p[4],&mut ez)
                    } else {
                        super::clean_impl512::<$a,$c,$r,$n>(&c.q,&c.t,5,&c.mat,
                            c.gaps[0],c.gaps[1],c.gaps[2],c.gaps[3],c.p[2],c.p[6],c.p[7],c.p[5],c.p[3],c.p[4],&mut ez)
                    }
                };
            }
            match (approx,cigar,right,has_n) {
                (false,false,false,false)=>go!(false,false,false,false),
                (false,false,false,true)=>go!(false,false,false,true),
                (false,false,true,false)=>go!(false,false,true,false),
                (false,false,true,true)=>go!(false,false,true,true),
                (false,true,false,false)=>go!(false,true,false,false),
                (false,true,false,true)=>go!(false,true,false,true),
                (false,true,true,false)=>go!(false,true,true,false),
                (false,true,true,true)=>go!(false,true,true,true),
                (true,false,false,false)=>go!(true,false,false,false),
                (true,false,false,true)=>go!(true,false,false,true),
                (true,false,true,false)=>go!(true,false,true,false),
                (true,false,true,true)=>go!(true,false,true,true),
                (true,true,false,false)=>go!(true,true,false,false),
                (true,true,false,true)=>go!(true,true,false,true),
                (true,true,true,false)=>go!(true,true,true,false),
                (true,true,true,true)=>go!(true,true,true,true),
            }
        }
        ez
    }
    CHECK_SHIFTED_GAPS.with(|v|v.set(true));
    for (i,c) in cases.iter().enumerate() { assert_result(&run(c,false,narrow,narrow_vl,wide),&run(c,true,narrow,narrow_vl,wide),&format!("captured call {i}")); }
    CHECK_SHIFTED_GAPS.with(|v|v.set(false));
    eprintln!("ACTIVE_GAP_CHECK_ORACLE_PASSED cases={};CHECK=false for timing",cases.len());
    let mut csv = String::from("trial,mode,approximate_max,cases,repeats,wall_s,checksum\n");
    for approx in [false,true] {
        let selected: Vec<_> = cases.iter().filter(|c| (c.p[5] & APPROX_MAX != 0) == approx).collect();
        if selected.is_empty() { continue; }
        for trial in 1..=3 {
            let order = if trial % 2 == 1 { [false,true] } else { [true,false] };
            for registers in order {
                let begin = std::time::Instant::now(); let mut sum = 0i64;
                for _ in 0..8 {
                    for c in &selected {
                        let ez = std::hint::black_box(run(std::hint::black_box(c),registers,narrow,narrow_vl,wide));
                        sum = sum.wrapping_add(ez.score as i64 + ez.cigar.len() as i64);
                    }
                }
                let elapsed=begin.elapsed().as_secs_f64();
                let mode=if registers && wide {"array512"} else if registers && narrow_vl {"register32_vl"}
                    else if registers {"register64"} else if narrow {"register32_avx2"} else if wide {"array256"} else {"array512"};
                let line=format!("{trial},{mode},{approx},{},8,{elapsed:.9},{sum}\n",selected.len());
                print!("{line}");csv.push_str(&line);
            }
        }
    }
    if let Ok(path)=std::env::var("BONEFISH_DP_REPLAY_RESULTS") { std::fs::write(path,csv).unwrap(); }
}

struct TraceLimitGuard(Option<usize>);
impl TraceLimitGuard {
    fn set(limit:usize)->Self {Self(TEST_TRACE_STREAM_MIN.with(|v|v.replace(Some(limit))))}
}
impl Drop for TraceLimitGuard {fn drop(&mut self){TEST_TRACE_STREAM_MIN.with(|v|v.set(self.0));}}

fn run_array_trace(q:&[u8],t:&[u8],mat:&[i8],gaps:[i8;4],p:[i32;6],stream:bool)->DpResult {
    let _guard=TraceLimitGuard::set(if stream {1} else {usize::MAX});
    let mut out=DpResult::default();let [bw,adaptive,shift,flags,zdrop,bonus]=p;
    let has_n=q.iter().chain(t).any(|&x|x>=4);
    macro_rules! go {($a:literal,$c:literal,$r:literal,$n:literal)=>{unsafe {
        clean_impl512::<$a,$c,$r,$n>(q,t,5,mat,gaps[0],gaps[1],gaps[2],gaps[3],bw,zdrop,bonus,flags,adaptive,shift,&mut out)
    }}}
    match (flags&APPROX_MAX!=0,flags&SCORE_ONLY==0,flags&RIGHT_ALIGN!=0,has_n) {
        (false,false,false,false)=>go!(false,false,false,false),(false,false,false,true)=>go!(false,false,false,true),
        (false,false,true,false)=>go!(false,false,true,false),(false,false,true,true)=>go!(false,false,true,true),
        (false,true,false,false)=>go!(false,true,false,false),(false,true,false,true)=>go!(false,true,false,true),
        (false,true,true,false)=>go!(false,true,true,false),(false,true,true,true)=>go!(false,true,true,true),
        (true,false,false,false)=>go!(true,false,false,false),(true,false,false,true)=>go!(true,false,false,true),
        (true,false,true,false)=>go!(true,false,true,false),(true,false,true,true)=>go!(true,false,true,true),
        (true,true,false,false)=>go!(true,true,false,false),(true,true,false,true)=>go!(true,true,false,true),
        (true,true,true,false)=>go!(true,true,true,false),(true,true,true,true)=>go!(true,true,true,true),
    }
    out
}

#[test]
fn streamed_traceback_matches_cached_and_avx2_including_early_exits() {
    if !available(){return}
    let mut seed=0x5126_3417_abcd_098fu64;
    fn rng(x:&mut u64)->u64 {*x^=*x<<13;*x^=*x>>7;*x^=*x<<17;*x}
    let mut cigars=0;let mut dropped=0;
    for case in 0..240 {
        let n=100+(rng(&mut seed)%1100) as usize;
        let mut t:Vec<u8>=(0..n).map(|_|(rng(&mut seed)&3) as u8).collect();
        let mut q=t.clone();
        for _ in 0..n/17 {let j=rng(&mut seed) as usize%q.len();q[j]=(rng(&mut seed)&3) as u8;}
        if case%3==0 {q.insert(n/2,4);}if case%5==0 {t[n/3]=4;}
        if case%7==0 {for j in n/2..n {q[j]=t[j].wrapping_add(1)&3;}}
        let bw=[64,65,95,127,128,255,-1][case%7];
        let adaptive=if case%3==0 {[65,96,128][(case/3)%3]} else {0};
        let shift=if adaptive==0 {[-17,0,17][case%3]} else {0};
        let flags=(if case&1!=0 {APPROX_MAX|APPROX_DROP} else {0})
            |(if case&2!=0 {RIGHT_ALIGN} else {0})|(if case&4!=0 {SCORE_ONLY} else {0})
            |(if case&8!=0 {EXTENSION_ONLY|REV_CIGAR} else {0});
        let zdrop=if case%7==0 {10} else {[-1,10,40,400][case%4]};let bonus=[-1,0,5,20][case%4];
        let p=[bw,adaptive,shift,flags,zdrop,bonus];let gaps=[4,2,24,1];let mat=matrix(-1);
        let cached=run_array_trace(&q,&t,&mat,gaps,p,false);
        let streamed=run_array_trace(&q,&t,&mat,gaps,p,true);
        assert_result(&cached,&streamed,&format!("stream case={case} p={p:?}"));
        let mut avx2=DpResult::default();unsafe {super::super::dual_clean::extend_dual_affine_avx2_clean(
            &q,&t,5,&mat,gaps[0],gaps[1],gaps[2],gaps[3],bw,zdrop,bonus,flags,adaptive,shift,&mut avx2);}
        assert_result(&avx2,&streamed,&format!("AVX2 stream case={case}"));
        cigars+=usize::from(flags&SCORE_ONLY==0 && !streamed.cigar.is_empty());
        dropped+=usize::from(streamed.zdropped!=0);
    }
    assert!(cigars>20);assert!(dropped>10);
}

#[test]
#[ignore]
fn benchmark_trace_stores_real_calls() {
    if !available(){return}
    let bytes=std::fs::read(std::env::var("BONEFISH_DP_REPLAY").unwrap()).unwrap();assert_eq!(&bytes[..8],b"BFDPPR01");
    let min=std::env::var("BONEFISH_TRACE_REPLAY_MIN_KB").ok().and_then(|v|v.parse::<usize>().ok()).unwrap_or(1024)*1024;
    struct Case {q:Vec<u8>,t:Vec<u8>,mat:[i8;25],gaps:[i8;4],p:[i32;8]}
    let mut at=8;let mut cases=Vec::new();
    while at<bytes.len() {
        let p=std::array::from_fn(|i|i32::from_le_bytes(bytes[at+i*4..at+i*4+4].try_into().unwrap()));at+=32;
        let gaps=std::array::from_fn(|i|bytes[at+i] as i8);at+=4;
        let mat=std::array::from_fn(|i|bytes[at+i] as i8);at+=25+8;
        let q=bytes[at..at+p[0] as usize].to_vec();at+=p[0] as usize;
        let t=bytes[at..at+p[1] as usize].to_vec();at+=p[1] as usize;
        let cells=super::super::dual_clean::clean_band_cells(q.len(),t.len(),p[2],p[3]);
        let size=(q.len()+t.len()-1)*cells.div_ceil(64)*64;
        if cells>64 && p[5]&SCORE_ONLY==0 && size>=min {cases.push(Case{q,t,mat,gaps,p});}
    }
    assert!(!cases.is_empty());
    let run=|c:&Case,stream|run_array_trace(&c.q,&c.t,&c.mat,c.gaps,[c.p[2],c.p[3],c.p[4],c.p[5],c.p[6],c.p[7]],stream);
    for (i,c) in cases.iter().enumerate() {
        assert_result(&run(c,false),&run(c,true),&format!("captured trace {i}"));
    }
    let mut csv=String::from("trial,mode,approx,cases,repeats,wall_s,checksum\n");
    for approx in [false,true] {
        let selected:Vec<_>=cases.iter().filter(|c|(c.p[5]&APPROX_MAX!=0)==approx).collect();if selected.is_empty(){continue}
        for trial in 1..=3 {
            let order=if trial%2==1 {[false,true]} else {[true,false]};
            for stream in order {
                let start=std::time::Instant::now();let mut sum=0i64;
                for _ in 0..4 {for c in &selected {let out=std::hint::black_box(run(std::hint::black_box(c),stream));sum=sum.wrapping_add(out.score as i64+out.cigar.len() as i64);}}
                let mode=if stream {"stream512"} else {"cached512"};
                let line=format!("{trial},{mode},{approx},{},4,{:.9},{sum}\n",selected.len(),start.elapsed().as_secs_f64());print!("{line}");csv.push_str(&line);
            }
        }
    }
    if let Ok(p)=std::env::var("BONEFISH_DP_REPLAY_RESULTS"){std::fs::write(p,csv).unwrap();}
}

fn run_shifted_case<const PING:bool>(q:&[u8],t:&[u8],mat:&[i8],gaps:[i8;4],p:[i32;6],stream:bool)->DpResult {
    let _guard=TraceLimitGuard::set(if stream {1} else {usize::MAX});
    let mut out=DpResult::default();let [bw,adaptive,shift,flags,zdrop,bonus]=p;
    let has_n=q.iter().chain(t).any(|&x|x>=4);
    macro_rules! go {($a:literal,$c:literal,$r:literal,$n:literal)=>{unsafe {
        clean_impl512_state::<$a,$c,$r,$n,PING>(q,t,5,mat,gaps[0],gaps[1],gaps[2],gaps[3],bw,zdrop,bonus,flags,adaptive,shift,&mut out)
    }}}
    match (flags&APPROX_MAX!=0,flags&SCORE_ONLY==0,flags&RIGHT_ALIGN!=0,has_n) {
        (false,false,false,false)=>go!(false,false,false,false),(false,false,false,true)=>go!(false,false,false,true),
        (false,false,true,false)=>go!(false,false,true,false),(false,false,true,true)=>go!(false,false,true,true),
        (false,true,false,false)=>go!(false,true,false,false),(false,true,false,true)=>go!(false,true,false,true),
        (false,true,true,false)=>go!(false,true,true,false),(false,true,true,true)=>go!(false,true,true,true),
        (true,false,false,false)=>go!(true,false,false,false),(true,false,false,true)=>go!(true,false,false,true),
        (true,false,true,false)=>go!(true,false,true,false),(true,false,true,true)=>go!(true,false,true,true),
        (true,true,false,false)=>go!(true,true,false,false),(true,true,false,true)=>go!(true,true,false,true),
        (true,true,true,false)=>go!(true,true,true,false),(true,true,true,true)=>go!(true,true,true,true),
    }
    out
}

#[test]
fn shifted_previous_row_loads_match_in_place_and_avx2_all_fields() {
    if !available(){return}
    let mut seed=0x512d_1ec7_b0ff_3456u64;
    fn rng(x:&mut u64)->u64 {*x^=*x<<13;*x^=*x>>7;*x^=*x<<17;*x}
    let mut cigars=0;let mut drops=0;
    for case in 0..480 {
        let bounds=[1,31,32,33,63,64,65,127,128,129,255,256,257,511,512,513,1023,1024,1025];
        let n=if case%4==0 {bounds[(case/4)%bounds.len()]} else {100+(rng(&mut seed)%1500) as usize};
        let mut t:Vec<u8>=(0..n).map(|_|(rng(&mut seed)&3) as u8).collect();
        let mut q=t.clone();
        for _ in 0..n/19 {let j=rng(&mut seed) as usize%q.len();q[j]=(rng(&mut seed)&3) as u8;}
        if case%3==0 {q.insert(n/2,4);}if case%5==0 {t[n/3]=4;}
        if case%7==0 {for j in n/2..n {q[j]=t[j].wrapping_add(1)&3;}}
        if case%11==0 && n>50 {q.drain(n/3..n/3+17);}
        let bw=[31,32,63,64,65,95,127,128,129,255,752,-1][case%12];
        let adaptive=if case%3==0 {[33,65,96,128,255][(case/3)%5]} else {0};
        let shift=if adaptive==0 {[-33,-17,0,17,33][case%5]} else {0};
        let flags=(if case&1!=0 {APPROX_MAX|APPROX_DROP} else {0})
            |(if case&2!=0 {RIGHT_ALIGN} else {0})|(if case&4!=0 {SCORE_ONLY} else {0})
            |(if case&8!=0 {EXTENSION_ONLY|REV_CIGAR} else {0});
        let zdrop=if case%7==0 {10} else {[-1,10,40,400][case%4]};let bonus=[-1,0,5,20][case%4];
        let p=[bw,adaptive,shift,flags,zdrop,bonus];
        let gaps=[[4,2,24,1],[24,1,4,2],[4,1,8,1],[6,2,12,1]][case%4];
        let mat=matrix(if case%9==0 {0} else {-1});
        let inplace=run_shifted_case::<false>(&q,&t,&mat,gaps,p,case%2!=0);
        let direct=run_shifted_case::<true>(&q,&t,&mat,gaps,p,case%2!=0);
        assert_result(&inplace,&direct,&format!("shifted case={case} q={} t={} p={p:?}",q.len(),t.len()));
        let mut avx2=DpResult::default();unsafe {super::super::dual_clean::extend_dual_affine_avx2_clean(
            &q,&t,5,&mat,gaps[0],gaps[1],gaps[2],gaps[3],bw,zdrop,bonus,flags,adaptive,shift,&mut avx2);}
        assert_result(&avx2,&direct,&format!("AVX2 shifted case={case}"));
        cigars+=usize::from(flags&SCORE_ONLY==0 && !direct.cigar.is_empty());
        drops+=usize::from(direct.zdropped!=0);
    }
    assert!(cigars>20);assert!(drops>10);
}


#[test]
fn native_masked_continuation_bits_match_every_direction_and_mask_combination() {
    if !available(){return}
    #[target_feature(enable="avx2,avx512f,avx512bw")]
    unsafe fn run(d:&[u8;64],m:[u64;4])->[u8;64] {unsafe {
        let v=_mm512_loadu_si512(d.as_ptr().cast());
        let v=append_continuation_bits_native!(v,m[0],m[1],m[2],m[3],
            _mm512_set1_epi8(8),_mm512_set1_epi8(16),_mm512_set1_epi8(32),_mm512_set1_epi8(64));
        let mut out=[0u8;64];_mm512_storeu_si512(out.as_mut_ptr().cast(),v);out
    }}
    for direction in 0..8u8 {
        for bits in 0..16u8 {
            let d:[u8;64]=std::array::from_fn(|lane|(direction+lane as u8)%8);
            let mut masks=[0u64;4];
            let expected:[u8;64]=std::array::from_fn(|lane| {
                let b=bits^(lane as u8%16);
                for (i,m) in masks.iter_mut().enumerate(){if b&(1<<i)!=0 {*m|=1u64<<lane;}}
                d[lane]|(b<<3)
            });
            assert_eq!(unsafe {run(&d,masks)},expected,"direction{direction} mask{bits}");
        }
    }
}

fn run_plane_case<const PLANES:bool>(q:&[u8],t:&[u8],mat:&[i8],gaps:[i8;4],p:[i32;6],stream:bool)->DpResult {
    let _guard=TraceLimitGuard::set(if stream {1} else {usize::MAX});
    let mut out=DpResult::default();let [bw,adaptive,shift,flags,zdrop,bonus]=p;
    let has_n=q.iter().chain(t).any(|&x|x>=4);
    macro_rules! go {($a:literal,$c:literal,$r:literal,$n:literal)=>{unsafe {
        clean_impl512_planes::<$a,$c,$r,$n,true,PLANES>(q,t,5,mat,gaps[0],gaps[1],gaps[2],gaps[3],bw,zdrop,bonus,flags,adaptive,shift,&mut out)
    }}}
    match (flags&APPROX_MAX!=0,flags&SCORE_ONLY==0,flags&RIGHT_ALIGN!=0,has_n) {
        (false,false,false,false)=>go!(false,false,false,false),(false,false,false,true)=>go!(false,false,false,true),
        (false,false,true,false)=>go!(false,false,true,false),(false,false,true,true)=>go!(false,false,true,true),
        (false,true,false,false)=>go!(false,true,false,false),(false,true,false,true)=>go!(false,true,false,true),
        (false,true,true,false)=>go!(false,true,true,false),(false,true,true,true)=>go!(false,true,true,true),
        (true,false,false,false)=>go!(true,false,false,false),(true,false,false,true)=>go!(true,false,false,true),
        (true,false,true,false)=>go!(true,false,true,false),(true,false,true,true)=>go!(true,false,true,true),
        (true,true,false,false)=>go!(true,true,false,false),(true,true,false,true)=>go!(true,true,false,true),
        (true,true,true,false)=>go!(true,true,true,false),(true,true,true,true)=>go!(true,true,true,true),
    }
    out
}

#[test]
fn balanced_trace_planes_preserve_all_ties_masks_and_lanes_cached_and_streamed() {
    if !available(){return}
    #[target_feature(enable="avx2,avx512f,avx512bw")]
    unsafe fn run<const RIGHT:bool>(input:&[[i8;64];5],gaps:[u64;4],nt:bool)->([i8;64],[u8;64]) {unsafe {
        let v=[_mm512_loadu_si512(input[0].as_ptr().cast()),_mm512_loadu_si512(input[1].as_ptr().cast()),
            _mm512_loadu_si512(input[2].as_ptr().cast()),_mm512_loadu_si512(input[3].as_ptr().cast()),
            _mm512_loadu_si512(input[4].as_ptr().cast())];
        let (z,m0,m1,m2,m3)=trace_plane_argmax!(v[0],v[1],v[2],v[3],v[4],RIGHT);
        let mem=AlignedMemory::new(64,64);let p=mem.as_ptr();
        store_four_trace_planes!(p,nt,m0,m1,m2,m3);
        store_four_trace_planes!(p.add(32),nt,gaps[0],gaps[1],gaps[2],gaps[3]);
        if nt {_mm_sfence();}
        let mut scores=[0i8;64];_mm512_storeu_si512(scores.as_mut_ptr().cast(),z);
        let flags:[u8;64]=std::array::from_fn(|lane| {
            let mut d=trace_plane_direction64(p,lane) as u8;
            for bit in 0..4 {if trace_plane_bit64(p,lane,4+bit){d|=1<<(3+bit);}}
            d
        });
        (scores,flags)
    }}
    fn rng(s:&mut u64)->u64 {*s^=*s<<13;*s^=*s>>7;*s^=*s<<17;*s}
    let mut seed=0x5128_71ac_9123_abefu64;
    for case in 0..8192 {
        let mut input:[[i8;64];5]=std::array::from_fn(|_|std::array::from_fn(|_|rng(&mut seed) as i8));
        if case%4==0 {for lane in 0..64 {for state in 1..5 {input[state][lane]=input[0][lane];}}}
        if case%7==0 {input[2]=input[1];input[4]=input[3];}
        if case%11==0 {for state in &mut input {for v in state {*v=[i8::MIN,0,i8::MAX][rng(&mut seed) as usize%3];}}}
        let gaps=if case%13==0 {[0,u64::MAX,0xaaaa_aaaa_aaaa_aaaa,0x5555_5555_5555_5555]}
            else {std::array::from_fn(|_|rng(&mut seed))};
        for right in [false,true] {
            let mut scores=[0i8;64];let mut flags=[0u8;64];
            for lane in 0..64 {
                let mut winner=0usize;
                for state in 1..5 {
                    if input[state][lane]>input[winner][lane] || (right && input[state][lane]==input[winner][lane]) {winner=state;}
                }
                scores[lane]=input[winner][lane];flags[lane]=winner as u8;
                for bit in 0..4 {if (gaps[bit]>>lane)&1!=0 {flags[lane]|=1<<(3+bit);}}
            }
            for nt in [false,true] {
                let got=if right {unsafe {run::<true>(&input,gaps,nt)}} else {unsafe {run::<false>(&input,gaps,nt)}};
                assert_eq!(got,(scores,flags),"plane case{case} right{right} nt{nt}");
            }
        }
    }
}

#[test]
fn bitplane_traceback_matches_byte_trace_and_avx2_including_early_exits() {
    if !available(){return}
    let mut seed=0x512b_17a1_76ca_d203u64;
    fn rng(x:&mut u64)->u64 {*x^=*x<<13;*x^=*x>>7;*x^=*x<<17;*x}
    let mut cigars=0;let mut drops=0;
    for case in 0..320 {
        let bounds=[1,31,32,33,63,64,65,127,128,129,255,256,257,511,512,513,1023,1024,1025];
        let n=if case%4==0 {bounds[(case/4)%bounds.len()]} else {100+(rng(&mut seed)%1500) as usize};
        let mut t:Vec<u8>=(0..n).map(|_|(rng(&mut seed)&3) as u8).collect();
        let mut q=t.clone();
        for _ in 0..n/19 {let j=rng(&mut seed) as usize%q.len();q[j]=(rng(&mut seed)&3) as u8;}
        if case%3==0 {q.insert(n/2,4);}if case%5==0 {t[n/3]=4;}
        if case%7==0 {for j in n/2..n {q[j]=t[j].wrapping_add(1)&3;}}
        if case%11==0 && n>50 {q.drain(n/3..n/3+17);}
        let bw=[0,31,32,63,64,65,95,127,128,129,255,752,-1][case%13];
        let adaptive=if case%3==0 {[33,65,96,128,255][(case/3)%5]} else {0};
        let shift=if adaptive==0 {[-33,-17,0,17,33][case%5]} else {0};
        let flags=(if case&1!=0 {APPROX_MAX|APPROX_DROP} else {0})
            |(if case&2!=0 {RIGHT_ALIGN} else {0})|(if case&4!=0 {SCORE_ONLY} else {0})
            |(if case&8!=0 {EXTENSION_ONLY|REV_CIGAR} else {0});
        let zdrop=if case%7==0 {10} else {[-1,10,40,400][case%4]};let bonus=[-1,0,5,20][case%4];
        let p=[bw,adaptive,shift,flags,zdrop,bonus];
        let gaps=[[4,2,24,1],[24,1,4,2],[4,1,8,1],[6,2,12,1]][case%4];let mat=matrix(if case%9==0 {0} else {-1});
        let byte=run_plane_case::<false>(&q,&t,&mat,gaps,p,case%2!=0);
        let cached=run_plane_case::<true>(&q,&t,&mat,gaps,p,false);
        let streamed=run_plane_case::<true>(&q,&t,&mat,gaps,p,true);
        assert_result(&byte,&cached,&format!("plane cached case{case} p{p:?}"));
        assert_result(&byte,&streamed,&format!("plane NT case{case} p{p:?}"));
        let mut avx2=DpResult::default();unsafe {super::super::dual_clean::extend_dual_affine_avx2_clean(
            &q,&t,5,&mat,gaps[0],gaps[1],gaps[2],gaps[3],bw,zdrop,bonus,flags,adaptive,shift,&mut avx2);}
        assert_result(&avx2,&streamed,&format!("AVX2 plane case{case}"));
        cigars+=usize::from(flags&SCORE_ONLY==0 && !streamed.cigar.is_empty());drops+=usize::from(streamed.zdropped!=0);
    }
    assert!(cigars>20);assert!(drops>10);
}


#[test]
fn signed_gap_clip_matches_mask_selection_for_every_byte_and_alignment_side() {
    if !available(){return}
    #[target_feature(enable="avx2,avx512f,avx512bw")]
    unsafe fn run(bytes:&[i8;64],qe:i8,right:bool)->([i8;64],[i8;64]) {unsafe {
        let v=_mm512_loadu_si512(bytes.as_ptr().cast());let z=_mm512_setzero_si512();
        let m=if right {!_mm512_cmpgt_epi8_mask(z,v)} else {_mm512_cmpgt_epi8_mask(v,z)};
        let old=_mm512_sub_epi8(_mm512_maskz_mov_epi8(m,v),_mm512_set1_epi8(qe));
        let new=_mm512_sub_epi8(_mm512_max_epi8(v,z),_mm512_set1_epi8(qe));
        let mut a=[0i8;64];let mut b=[0i8;64];
        _mm512_storeu_si512(a.as_mut_ptr().cast(),old);_mm512_storeu_si512(b.as_mut_ptr().cast(),new);(a,b)
    }}
    for right in [false,true] {
        for qe in [-128i8,-127,-25,-6,-1,0,1,6,25,126,127] {
            for group in 0..4 {
                let bytes:[i8;64]=std::array::from_fn(|lane|(group*64+lane) as u8 as i8);
                let (old,new)=unsafe {run(&bytes,qe,right)};
                let expected:[i8;64]=std::array::from_fn(|lane|bytes[lane].max(0).wrapping_sub(qe));
                assert_eq!(old,expected);assert_eq!(new,expected);
            }
        }
    }
}


#[test]
fn unmasked_h16_blocks_preserve_four_lane_ties_and_unaligned_score_updates() {
    if !available(){return}
    #[target_feature(enable="avx2,avx512f,avx512bw")]
    unsafe fn run(h:&mut [i32],v:&[u8],pad:usize,n:usize)->[i32;16] {unsafe {
        let mut low=_mm512_setr_epi32(
            32767,24575,16383,8191,32766,24574,16382,8190,
            32765,24573,16381,8189,32764,24572,16380,8188);
        let four=_mm512_set1_epi32(4);let mut mk=_mm512_set1_epi32(i32::MIN);
        for at in (0..n).step_by(16) {
            mk=update_h16_packed(h.as_mut_ptr().add(pad+at),v.as_ptr().add(pad+at),low,mk);
            low=_mm512_sub_epi32(low,four);
        }
        let mut keys=[0i32;16];_mm512_storeu_si512(keys.as_mut_ptr().cast(),mk);keys
    }}
    fn rng(s:&mut u64)->u64 {*s^=*s<<13;*s^=*s>>7;*s^=*s<<17;*s}
    let mut seed=0x16ff_5120_8423_12afu64;
    for case in 0..512 {
        let n=[1,2,3,7,31,127,500][case%7]*16;let pad=[0,1,3,7,17][case%5];
        let mut h=vec![0x5a5a5a5ai32;n+pad+16];let mut v=vec![0u8;n+pad+16];
        for at in 0..n {
            h[pad+at]=if case%3==0 {17} else if case%5==0 {(rng(&mut seed)%5) as i32-2} else {(rng(&mut seed)%80001) as i32-40000};
            v[pad+at]=if case%3==0 || case%5==0 {0} else {rng(&mut seed) as u8};
        }
        if case%11==0 {
            for at in 0..n {h[pad+at]=-100;v[pad+at]=0;}
            h[pad+n-4]=100;h[pad+1]=100;
        } else if case%13==0 {
            for at in 0..n {h[pad+at]=-17;v[pad+at]=0;}
        }
        let mut expected=h.clone();let mut lane_h=[i32::MIN;4];let mut lane_t=[0usize;4];
        for at in 0..n {
            expected[pad+at]=expected[pad+at].wrapping_add(v[pad+at] as i8 as i32);
            let lane=at%4;
            if expected[pad+at]>lane_h[lane] {lane_h[lane]=expected[pad+at];lane_t[lane]=at;}
        }
        let mut winner=0;
        for lane in 1..4 {if lane_h[lane]>lane_h[winner] {winner=lane;}}
        let keys=unsafe {run(&mut h,&v,pad,n)};assert_eq!(h,expected,"all score bytes case{case}");
        let key=*keys.iter().max().unwrap();let rank=32767-(key&0x7fff);
        let got_position=((rank&0x1fff)*4+(rank>>13)) as usize;
        assert_eq!(key>>15,lane_h[winner],"score case{case}");
        assert_eq!(got_position,lane_t[winner],"tie position case{case}");
    }
}

#[test]
fn shifted_gap_clipping_and_ties_match_scalar_over_active_ranges() {
    if !available(){return}
    #[target_feature(enable="avx2,avx512f,avx512bw")]
    unsafe fn run(an:&[i8;64],qe:i8,right:bool)->([i8;64],u64) {unsafe {
        let old=_mm512_loadu_si512(an.as_ptr().cast());let penalty=_mm512_set1_epi8(qe);
        let shifted=_mm512_sub_epi8(old,penalty);let floor=_mm512_set1_epi8(-qe);
        let next=_mm512_max_epi8(shifted,floor);
        let mask=if right {!_mm512_cmpgt_epi8_mask(floor,shifted)} else {_mm512_cmpgt_epi8_mask(shifted,floor)};
        let mut out=[0i8;64];_mm512_storeu_si512(out.as_mut_ptr().cast(),next);(out,mask)
    }}
    for (low,high,qe) in [(-10i8,4i8,6i8),(-9,24,25)] {
        for start in low..=high {
            let values:[i8;64]=std::array::from_fn(|i|low+((((start-low) as usize+i)%((high-low+1) as usize)) as i8));
            for right in [false,true] {
                let expected:[i8;64]=std::array::from_fn(|i|values[i].max(0).wrapping_sub(qe));
                let expected_mask=values.iter().enumerate().fold(0u64,|m,(i,&v)|m|((if right {v>=0} else {v>0}) as u64)<<i);
                assert_eq!(unsafe {run(&values,qe,right)},(expected,expected_mask));
            }
        }
    }
}

fn checked_clip_case(q:&[u8],t:&[u8],mat:&[i8],g:[i8;4],p:[i32;6],stream:bool)->DpResult {
    let [bw,adaptive,shift,flags,zdrop,bonus]=p;assert_eq!(adaptive,0);
    let mut out=DpResult::default();let has_n=q.iter().chain(t).any(|&x|x>=4);
    macro_rules! go {($a:literal,$r:literal,$n:literal,$s:literal)=>{unsafe {
        clean_impl512_clip::<$a,true,$r,$n,true,true,$s,true,true>(q,t,5,mat,g[0],g[1],g[2],g[3],bw,zdrop,bonus,flags,0,shift,&mut out)
    }}}
    match(flags&APPROX_MAX!=0,flags&RIGHT_ALIGN!=0,has_n,stream) {
        (false,false,false,false)=>go!(false,false,false,false),(false,false,false,true)=>go!(false,false,false,true),
        (false,false,true,false)=>go!(false,false,true,false),(false,false,true,true)=>go!(false,false,true,true),
        (false,true,false,false)=>go!(false,true,false,false),(false,true,false,true)=>go!(false,true,false,true),
        (false,true,true,false)=>go!(false,true,true,false),(false,true,true,true)=>go!(false,true,true,true),
        (true,false,false,false)=>go!(true,false,false,false),(true,false,false,true)=>go!(true,false,false,true),
        (true,false,true,false)=>go!(true,false,true,false),(true,false,true,true)=>go!(true,false,true,true),
        (true,true,false,false)=>go!(true,true,false,false),(true,true,false,true)=>go!(true,true,false,true),
        (true,true,true,false)=>go!(true,true,true,false),(true,true,true,true)=>go!(true,true,true,true),
    }
    out
}

#[test]
fn checked_gap_shift_preserves_every_active_lane_and_full_avx2_result() {
    if !available(){return}
    fn rng(x:&mut u64)->u64 {*x^=*x<<13;*x^=*x>>7;*x^=*x<<17;*x}
    let mut seed=0x6a95_5120_b0ff_013fu64;
    for case in 0..800 {
        let n=if case%5==0 {[1,31,32,33,63,64,65,127,128,129,191,192,193,255,256,257,511,512,513,1023][(case/5)%20]} else {64+(rng(&mut seed)%2400) as usize};
        let mut t:Vec<u8>=(0..n).map(|_|(rng(&mut seed)&3) as u8).collect();let mut q=t.clone();
        if case%7==0 {q.fill(0);t.fill(0);}
        for _ in 0..n/17 {let i=rng(&mut seed) as usize%q.len();q[i]=(rng(&mut seed)&3) as u8;}
        if case%3==0 {q.insert(n/2,4);}if case%11==0 {t[n/3]=4;}
        if case%13==0 && n>200 {q.drain(n/3..n/3+129);}
        let bw=[0,1,31,32,33,63,64,65,127,128,129,191,192,193,255,256,752,-1][case%18];
        let shift=[-752,-33,-1,0,1,33,752][case%7];
        let flags=(if case&1!=0 {APPROX_MAX|APPROX_DROP} else {0})|(if case&2!=0 {RIGHT_ALIGN} else {0})
            |(if case&4!=0 {REV_CIGAR|EXTENSION_ONLY} else {0});
        let zdrop=[-1,10,40,400][case%4];let bonus=[-1,0,5,20][case%4];let params=[bw,0,shift,flags,zdrop,bonus];
        let g=if case%2==0 {[4,2,24,1]} else {[24,1,4,2]};let mat=matrix([0,-1,2][case%3]);
        assert!(shifted_clip_eligible(q.len(),t.len(),5,&mat,g[0],g[1],g[2],g[3],bw,0,shift));
        let mut expected=DpResult::default();unsafe {super::super::dual_clean::extend_dual_affine_avx2_clean(
            &q,&t,5,&mat,g[0],g[1],g[2],g[3],bw,zdrop,bonus,flags,0,shift,&mut expected);}
        for stream in [false,true] {
            let got=checked_clip_case(&q,&t,&mat,g,params,stream);
            assert_result(&expected,&got,&format!("checked gap shift case{case} stream{stream}"));
        }
    }
}
