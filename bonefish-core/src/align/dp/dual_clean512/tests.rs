//! The AVX-512 clean kernels (64-lane array kernel, `reg64_impl`, `reg32vl_impl`) against
//! `clean_impl_reference` (the AVX2 array kernel alone), plus unit tests of the asm and mask
//! helpers the array kernel is built from.

use super::super::dual_clean::{clean_band_cells, clean_impl_reference, extend_dual_affine_clean};
use super::super::dual_clean_reg32vl::reg32vl_impl;
use super::super::dual_clean_reg64::reg64_impl;
use super::*;
use crate::align::dp::simd_test_available;

fn available() -> bool {
    simd_test_available("avx2", is_x86_feature_detected!("avx2"))
        && simd_test_available("avx512bw", is_x86_feature_detected!("avx512bw"))
}

fn vbmi_available() -> bool {
    available() && simd_test_available("avx512vbmi", is_x86_feature_detected!("avx512vbmi"))
}

fn vbmi_vl_available() -> bool {
    vbmi_available() && simd_test_available("avx512vl", is_x86_feature_detected!("avx512vl"))
}

/// Gap parameter sets: the default (both orders, which use the shortened gap arithmetic) and
/// two that take the general path, one with equal extensions.
const GAPS: [(i8, i8, i8, i8); 4] = [(4, 2, 24, 1), (24, 1, 4, 2), (6, 2, 12, 1), (4, 1, 8, 1)];

type DpFn =
    unsafe fn(&[u8], &[u8], i8, &[i8], i8, i8, i8, i8, i32, i32, i32, i32, i32, i32, &mut DpResult);

/// Calls `$kernel::<APPROX, CIGAR, RIGHT, HAS_N, $extra...>$args` for the call's flags.
macro_rules! by_flags {
    ($kernel:ident, [$($extra:tt),*], $flags:expr, $has_n:expr, $args:tt) => {
        match (($flags) & APPROX_MAX != 0, ($flags) & SCORE_ONLY == 0, ($flags) & RIGHT_ALIGN != 0, $has_n) {
            (false, false, false, false) => $kernel::<false, false, false, false $(, $extra)*> $args,
            (false, false, false, true) => $kernel::<false, false, false, true $(, $extra)*> $args,
            (false, false, true, false) => $kernel::<false, false, true, false $(, $extra)*> $args,
            (false, false, true, true) => $kernel::<false, false, true, true $(, $extra)*> $args,
            (false, true, false, false) => $kernel::<false, true, false, false $(, $extra)*> $args,
            (false, true, false, true) => $kernel::<false, true, false, true $(, $extra)*> $args,
            (false, true, true, false) => $kernel::<false, true, true, false $(, $extra)*> $args,
            (false, true, true, true) => $kernel::<false, true, true, true $(, $extra)*> $args,
            (true, false, false, false) => $kernel::<true, false, false, false $(, $extra)*> $args,
            (true, false, false, true) => $kernel::<true, false, false, true $(, $extra)*> $args,
            (true, false, true, false) => $kernel::<true, false, true, false $(, $extra)*> $args,
            (true, false, true, true) => $kernel::<true, false, true, true $(, $extra)*> $args,
            (true, true, false, false) => $kernel::<true, true, false, false $(, $extra)*> $args,
            (true, true, false, true) => $kernel::<true, true, false, true $(, $extra)*> $args,
            (true, true, true, false) => $kernel::<true, true, true, false $(, $extra)*> $args,
            (true, true, true, true) => $kernel::<true, true, true, true $(, $extra)*> $args,
        }
    };
}

/// Defines a `DpFn` that runs a generic kernel directly, bypassing the band-width routing.
macro_rules! flag_dispatched {
    ($name:ident, $kernel:ident) => {
        #[allow(clippy::too_many_arguments)]
        unsafe fn $name(
            q: &[u8],
            t: &[u8],
            alphabet: i8,
            mat: &[i8],
            go: i8,
            ge: i8,
            go2: i8,
            ge2: i8,
            bw: i32,
            zdrop: i32,
            bonus: i32,
            flags: i32,
            adaptive: i32,
            shift: i32,
            out: &mut DpResult,
        ) {
            let has_n = q.iter().chain(t).any(|&c| c >= 4);
            unsafe {
                by_flags!(
                    $kernel,
                    [],
                    flags,
                    has_n,
                    (
                        q, t, alphabet, mat, go, ge, go2, ge2, bw, zdrop, bonus, flags, adaptive,
                        shift, out
                    )
                )
            }
        }
    };
}
flag_dispatched!(array_kernel, clean_impl512);
flag_dispatched!(reg64_kernel, reg64_impl);
flag_dispatched!(reg32vl_kernel, reg32vl_impl);

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn pick<T: Copy>(&mut self, xs: &[T]) -> T {
        xs[self.below(xs.len())]
    }
    fn bases(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| (self.next() & 3) as u8).collect()
    }
    /// `t` with about `per_32` edits per 32 bases (half substitutions, the rest split between
    /// deletions and insertions).
    fn mutate(&mut self, t: &[u8], per_32: usize) -> Vec<u8> {
        let mut q = Vec::with_capacity(t.len() + t.len() / 8);
        for &c in t {
            match self.below(32) {
                r if r < per_32 / 2 => q.push((self.next() & 3) as u8),
                r if r < per_32 * 3 / 4 => {}
                r if r < per_32 => {
                    q.push(c);
                    q.push((self.next() & 3) as u8);
                }
                _ => q.push(c),
            }
        }
        if q.is_empty() {
            q.push(0);
        }
        q
    }
}

/// Scoring matrix: match 2, mismatch -4, `n_score` for N.
fn matrix(n_score: i8) -> [i8; 25] {
    std::array::from_fn(|i| {
        let (a, b) = (i / 5, i % 5);
        if a == 4 || b == 4 {
            n_score
        } else if a == b {
            2
        } else {
            -4
        }
    })
}

/// Flags from the low four bits of `bits`.
fn flag_combo(bits: usize) -> i32 {
    (if bits & 1 != 0 {
        APPROX_MAX | APPROX_DROP
    } else {
        0
    }) | (if bits & 2 != 0 { RIGHT_ALIGN } else { 0 })
        | (if bits & 4 != 0 { SCORE_ONLY } else { 0 })
        | (if bits & 8 != 0 {
            EXTENSION_ONLY | REV_CIGAR
        } else {
            0
        })
}

/// Some N bases in some cases.
fn add_ns(rng: &mut Rng, case: usize, q: &mut [u8], t: &mut [u8]) {
    if case.is_multiple_of(7) {
        q[rng.below(q.len())] = 4;
    }
    if case.is_multiple_of(11) {
        t[rng.below(t.len())] = 4;
    }
}

struct Case {
    q: Vec<u8>,
    t: Vec<u8>,
    mat: [i8; 25],
    gaps: (i8, i8, i8, i8),
    bw: i32,
    adaptive: i32,
    shift: i32,
    flags: i32,
    zdrop: i32,
    bonus: i32,
}

impl Case {
    fn new(q: Vec<u8>, t: Vec<u8>, gaps: (i8, i8, i8, i8), bw: i32, flags: i32) -> Self {
        Case {
            q,
            t,
            mat: matrix(-1),
            gaps,
            bw,
            adaptive: 0,
            shift: 0,
            flags,
            zdrop: 400,
            bonus: 5,
        }
    }

    fn run(&self, kernel: DpFn) -> DpResult {
        let (go, ge, go2, ge2) = self.gaps;
        let mut out = DpResult::default();
        unsafe {
            kernel(
                &self.q,
                &self.t,
                5,
                &self.mat,
                go,
                ge,
                go2,
                ge2,
                self.bw,
                self.zdrop,
                self.bonus,
                self.flags,
                self.adaptive,
                self.shift,
                &mut out,
            );
        }
        out
    }

    fn cells(&self) -> usize {
        clean_band_cells(self.q.len(), self.t.len(), self.bw, self.adaptive)
    }

    fn has_n(&self) -> bool {
        self.q.iter().chain(&self.t).any(|&c| c >= 4)
    }

    fn fast_eligible(&self) -> bool {
        let (go, ge, go2, ge2) = self.gaps;
        self.flags & SCORE_ONLY == 0
            && shifted_clip_eligible(
                self.q.len(),
                self.t.len(),
                5,
                &self.mat,
                go,
                ge,
                go2,
                ge2,
                self.bw,
                self.adaptive,
                self.shift,
            )
    }

    fn context(&self) -> String {
        format!(
            "qlen {} tlen {} gaps {:?} bw {} adaptive {} shift {} flags {:#x} zdrop {} bonus {} n {}",
            self.q.len(),
            self.t.len(),
            self.gaps,
            self.bw,
            self.adaptive,
            self.shift,
            self.flags,
            self.zdrop,
            self.bonus,
            self.mat[24]
        )
    }
}

fn assert_result(want: &DpResult, got: &DpResult, context: &str) {
    assert_eq!(
        (
            got.max,
            got.max_score_query_pos,
            got.max_score_target_pos,
            got.max_query_end_score,
            got.max_query_end_target_pos,
            got.max_target_end_score,
            got.max_target_end_query_pos,
            got.score,
            got.reach_end,
            got.zdropped,
            got.cigar_len,
            &got.cigar
        ),
        (
            want.max,
            want.max_score_query_pos,
            want.max_score_target_pos,
            want.max_query_end_score,
            want.max_query_end_target_pos,
            want.max_target_end_score,
            want.max_target_end_query_pos,
            want.score,
            want.reach_end,
            want.zdropped,
            want.cigar_len,
            &want.cigar
        ),
        "{context}"
    );
}

/// The AVX-512 dispatch (`reg64_impl` for 33..=64 cells on VBMI CPUs), the array kernel alone,
/// and the top-level routing all reproduce the reference.
fn check_avx512_paths(c: &Case) {
    let want = c.run(clean_impl_reference);
    assert_result(
        &want,
        &c.run(extend_dual_affine_avx512_clean),
        &format!("dispatch: {}", c.context()),
    );
    assert_result(
        &want,
        &c.run(array_kernel),
        &format!("array kernel: {}", c.context()),
    );
    assert_result(
        &want,
        &c.run(extend_dual_affine_clean),
        &format!("routing: {}", c.context()),
    );
}

#[test]
fn clean_band_cells_counts_active_cells() {
    assert_eq!(clean_band_cells(1000, 1000, 31, 0), 32);
    assert_eq!(clean_band_cells(1000, 1000, 32, 0), 33);
    assert_eq!(clean_band_cells(32, 1000, 500, 0), 32);
    assert_eq!(clean_band_cells(1000, 31, -1, 0), 31);
    assert_eq!(clean_band_cells(1000, 1000, 500, 32), 32);
    assert_eq!(clean_band_cells(1000, 1000, 500, 33), 33);
    assert_eq!(clean_band_cells(1000, 1000, 16, 64), 16);
}

/// Lengths and band widths around every vector boundary, all flag combinations, with ties
/// (long runs of one base) in every cell.
#[test]
fn avx512_clean_matches_reference_boundaries_and_ties() {
    if !available() {
        return;
    }
    let lengths = [
        1, 15, 16, 31, 32, 33, 63, 64, 65, 95, 96, 127, 128, 129, 255, 256, 257,
    ];
    for &len in &lengths {
        for &bw in &[0, 15, 31, 32, 33, 63, 64, 65, 127, 128, -1] {
            for bits in 0..16 {
                let t = vec![0; len];
                let mut q = vec![0; len + len / 3];
                if len > 32 {
                    q[len / 2] = 1;
                }
                let mut c = Case::new(q, t, GAPS[0], bw, flag_combo(bits));
                c.zdrop = 40;
                check_avx512_paths(&c);
            }
        }
    }
}

#[test]
fn avx512_clean_matches_reference_random_static_shifted_adaptive() {
    if !available() {
        return;
    }
    let mut rng = Rng(0xBAD5_EED5_1234_6789);
    let widths = [
        0, 1, 15, 16, 30, 31, 32, 33, 62, 63, 64, 65, 95, 127, 128, 129, 255, -1,
    ];
    for case in 0..6000 {
        let len = 1 + rng.below(600);
        let mut t = rng.bases(len);
        let mut q = rng.mutate(&t, 7);
        add_ns(&mut rng, case, &mut q, &mut t);
        let mut flags = 0;
        for flag in [
            APPROX_MAX,
            APPROX_DROP,
            RIGHT_ALIGN,
            SCORE_ONLY,
            EXTENSION_ONLY,
            REV_CIGAR,
        ] {
            if rng.next() & 1 != 0 {
                flags |= flag;
            }
        }
        let mut c = Case::new(q, t, GAPS[case % 4], widths[case % widths.len()], flags);
        c.adaptive = if case % 3 == 0 {
            rng.pick(&[1, 16, 31, 32, 33, 63, 64, 65, 128])
        } else {
            0
        };
        c.shift = if c.adaptive == 0 && case % 4 == 0 {
            rng.below(31) as i32 - 15
        } else {
            0
        };
        c.zdrop = rng.pick(&[-1, 10, 40, 100, 400]);
        c.bonus = rng.pick(&[-1, 0, 5, 20]);
        c.mat = matrix(if case % 5 == 0 { 0 } else { -1 });
        check_avx512_paths(&c);
    }
}

/// Exact max tracking beyond the packed-key length guard (uncompressed H-score fallback).
#[test]
fn avx512_clean_long_exact_score_range() {
    if !available() {
        return;
    }
    let mut rng = Rng(0x8123_4567_89ab_cdef);
    let t = rng.bases(10001);
    let mut q = t.clone();
    q.insert(3000, 4);
    q.drain(7000..7010);
    for adaptive in [0, 64] {
        for flags in [0, RIGHT_ALIGN, SCORE_ONLY, EXTENSION_ONLY | REV_CIGAR] {
            let mut c = Case::new(q.clone(), t.clone(), GAPS[0], 63, flags);
            c.adaptive = adaptive;
            check_avx512_paths(&c);
        }
    }
}

/// `reg64_impl` on its whole domain (33..=64 active cells), static, shifted and adaptive bands.
#[test]
fn register64_matches_reference() {
    if !vbmi_available() {
        return;
    }
    let mut rng = Rng(0x6464_0f0f_1234_abcd);
    let mut cigars = 0;
    for case in 0..4000 {
        let len = 1 + rng.below(800);
        let mut t = rng.bases(len);
        let per_32 = 1 + rng.below(8);
        let mut q = rng.mutate(&t, per_32);
        add_ns(&mut rng, case, &mut q, &mut t);
        let bw = if case % 5 == 0 {
            -1
        } else {
            32 + rng.below(32) as i32
        };
        let mut c = Case::new(q, t, GAPS[case % 4], bw, flag_combo(case));
        if case % 3 == 0 {
            c.bw = rng.pick(&[64, 100, 500, -1]);
            c.adaptive = 33 + rng.below(32) as i32;
        } else if case % 4 == 0 {
            c.shift = rng.below(41) as i32 - 20;
        }
        c.zdrop = rng.pick(&[-1, 10, 40, 400]);
        c.bonus = rng.pick(&[-1, 0, 5, 20]);
        c.mat = matrix(if case % 9 == 0 { 0 } else { -1 });
        if !(33..=64).contains(&c.cells()) {
            continue;
        }
        let want = c.run(clean_impl_reference);
        assert_result(
            &want,
            &c.run(reg64_kernel),
            &format!("case {case}: {}", c.context()),
        );
        cigars += usize::from(!want.cigar.is_empty());
    }
    assert!(cigars > 500, "{cigars}");
}

/// `reg32vl_impl` on the calls `extend_dual_affine_avx2_clean` gives it: approximate max
/// tracking, static bands of bandwidth 0..=31.
#[test]
fn register32vl_matches_reference() {
    if !vbmi_vl_available() {
        return;
    }
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    for case in 0..20000 {
        let len = 1 + rng.below(400);
        let mut t = rng.bases(len);
        let per_32 = rng.below(10);
        let mut q = rng.mutate(&t, per_32);
        add_ns(&mut rng, case, &mut q, &mut t);
        let mut flags = APPROX_MAX;
        for flag in [
            APPROX_DROP,
            RIGHT_ALIGN,
            SCORE_ONLY,
            EXTENSION_ONLY | REV_CIGAR,
        ] {
            if rng.next() & 1 != 0 {
                flags |= flag;
            }
        }
        let mut c = Case::new(q, t, GAPS[case % 4], rng.below(32) as i32, flags);
        c.shift = if case % 5 == 0 {
            rng.below(9) as i32 - 4
        } else {
            0
        };
        c.zdrop = rng.pick(&[400, 100, 30, -1]);
        c.bonus = if case % 3 == 0 { 5 } else { -1 };
        c.mat = matrix(if case % 13 == 0 { 0 } else { -1 });
        assert_result(
            &c.run(clean_impl_reference),
            &c.run(reg32vl_kernel),
            &format!("case {case}: {}", c.context()),
        );
    }
}

/// The array kernel with a chosen traceback layout and store mode (the dispatch picks them
/// from the band width and matrix size); `CHECK` asserts the shortened gap arithmetic.
fn run_array<const PLANES: bool, const STREAM: bool, const CHECK: bool>(c: &Case) -> DpResult {
    let fast = c.fast_eligible();
    assert!(
        fast || !CHECK,
        "CHECK needs a call eligible for the shortened gap arithmetic"
    );
    let (go, ge, go2, ge2) = c.gaps;
    let mut out = DpResult::default();
    unsafe {
        if fast {
            by_flags!(
                clean_impl512_kernel,
                [PLANES, STREAM, true, CHECK],
                c.flags,
                c.has_n(),
                (
                    &c.q, &c.t, 5, &c.mat, go, ge, go2, ge2, c.bw, c.zdrop, c.bonus, c.flags,
                    c.adaptive, c.shift, &mut out
                )
            )
        } else {
            by_flags!(
                clean_impl512_kernel,
                [PLANES, STREAM, false, false],
                c.flags,
                c.has_n(),
                (
                    &c.q, &c.t, 5, &c.mat, go, ge, go2, ge2, c.bw, c.zdrop, c.bonus, c.flags,
                    c.adaptive, c.shift, &mut out
                )
            )
        }
    }
    out
}

/// Byte and bit-plane traceback, each with cached and non-temporal stores, on every band width;
/// sequences around vector boundaries, with N bases, divergent tails that end in z-drops, and
/// deletions.
#[test]
fn array_kernel_traceback_layouts_and_store_modes_match_reference() {
    if !available() {
        return;
    }
    let mut rng = Rng(0x512b_17a1_76ca_d203);
    let bounds = [
        1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257, 511, 512, 513, 1023, 1024, 1025,
    ];
    let (mut cigars, mut drops) = (0, 0);
    for case in 0..480 {
        let len = if case % 4 == 0 {
            bounds[(case / 4) % bounds.len()]
        } else {
            100 + rng.below(1500)
        };
        let mut t = rng.bases(len);
        let mut q = t.clone();
        for _ in 0..len / 19 {
            let j = rng.below(q.len());
            q[j] = (rng.next() & 3) as u8;
        }
        if case % 3 == 0 {
            q.insert(len / 2, 4);
        }
        if case % 5 == 0 {
            t[len / 3] = 4;
        }
        if case % 7 == 0 {
            for j in len / 2..len {
                q[j] = t[j].wrapping_add(1) & 3;
            }
        }
        if case % 11 == 0 && len > 50 {
            q.drain(len / 3..len / 3 + 17);
        }
        let bw = [0, 31, 32, 63, 64, 65, 95, 127, 128, 129, 255, 752, -1][case % 13];
        let mut c = Case::new(q, t, GAPS[case % 4], bw, flag_combo(case));
        c.adaptive = if case % 3 == 0 {
            [33, 65, 96, 128, 255][(case / 3) % 5]
        } else {
            0
        };
        c.shift = if c.adaptive == 0 {
            [-33, -17, 0, 17, 33][case % 5]
        } else {
            0
        };
        c.zdrop = if case % 7 == 0 {
            10
        } else {
            [-1, 10, 40, 400][case % 4]
        };
        c.bonus = [-1, 0, 5, 20][case % 4];
        c.mat = matrix(if case % 9 == 0 { 0 } else { -1 });
        let want = c.run(clean_impl_reference);
        let context = c.context();
        assert_result(
            &want,
            &run_array::<false, false, false>(&c),
            &format!("bytes: {context}"),
        );
        assert_result(
            &want,
            &run_array::<false, true, false>(&c),
            &format!("bytes, streamed: {context}"),
        );
        assert_result(
            &want,
            &run_array::<true, false, false>(&c),
            &format!("planes: {context}"),
        );
        assert_result(
            &want,
            &run_array::<true, true, false>(&c),
            &format!("planes, streamed: {context}"),
        );
        cigars += usize::from(!want.cigar.is_empty());
        drops += usize::from(want.zdropped != 0);
    }
    assert!(cigars > 20 && drops > 10, "{cigars} {drops}");
}

/// The shortened gap arithmetic keeps every active lane within its bounds and equal to the
/// general form (in-kernel assertions), and the full results equal the reference.
#[test]
fn array_kernel_checked_gap_arithmetic_matches_reference() {
    if !available() {
        return;
    }
    let mut rng = Rng(0x6a95_5120_b0ff_013f);
    let bounds = [
        1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256, 257, 511, 512, 513, 1023,
    ];
    for case in 0..800 {
        let len = if case % 5 == 0 {
            bounds[(case / 5) % bounds.len()]
        } else {
            64 + rng.below(2400)
        };
        let mut t = rng.bases(len);
        let mut q = t.clone();
        if case % 7 == 0 {
            q.fill(0);
            t.fill(0);
        }
        for _ in 0..len / 17 {
            let i = rng.below(q.len());
            q[i] = (rng.next() & 3) as u8;
        }
        if case % 3 == 0 {
            q.insert(len / 2, 4);
        }
        if case % 11 == 0 {
            t[len / 3] = 4;
        }
        if case % 13 == 0 && len > 200 {
            q.drain(len / 3..len / 3 + 129);
        }
        let bw = [
            0, 1, 31, 32, 33, 63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256, 752, -1,
        ][case % 18];
        let flags = (if case & 1 != 0 {
            APPROX_MAX | APPROX_DROP
        } else {
            0
        }) | (if case & 2 != 0 { RIGHT_ALIGN } else { 0 })
            | (if case & 4 != 0 {
                REV_CIGAR | EXTENSION_ONLY
            } else {
                0
            });
        let mut c = Case::new(q, t, GAPS[case % 2], bw, flags);
        c.shift = [-752, -33, -1, 0, 1, 33, 752][case % 7];
        c.zdrop = [-1, 10, 40, 400][case % 4];
        c.bonus = [-1, 0, 5, 20][case % 4];
        c.mat = matrix([0, -1, 2][case % 3]);
        let want = c.run(clean_impl_reference);
        assert_result(
            &want,
            &run_array::<true, false, true>(&c),
            &format!("checked: {}", c.context()),
        );
        assert_result(
            &want,
            &run_array::<true, true, true>(&c),
            &format!("checked, streamed: {}", c.context()),
        );
    }
}

#[test]
fn native_masked_continuation_bits_match_every_direction_and_mask_combination() {
    if !available() {
        return;
    }
    #[target_feature(enable = "avx2,avx512f,avx512bw")]
    unsafe fn run(d: &[u8; 64], m: [u64; 4]) -> [u8; 64] {
        unsafe {
            let v = _mm512_loadu_si512(d.as_ptr().cast());
            let v = append_continuation_bits_native!(
                v,
                m[0],
                m[1],
                m[2],
                m[3],
                _mm512_set1_epi8(8),
                _mm512_set1_epi8(16),
                _mm512_set1_epi8(32),
                _mm512_set1_epi8(64)
            );
            let mut out = [0u8; 64];
            _mm512_storeu_si512(out.as_mut_ptr().cast(), v);
            out
        }
    }
    for direction in 0..8u8 {
        for bits in 0..16u8 {
            let d: [u8; 64] = std::array::from_fn(|lane| (direction + lane as u8) % 8);
            let mut masks = [0u64; 4];
            let expected: [u8; 64] = std::array::from_fn(|lane| {
                let b = bits ^ (lane as u8 % 16);
                for (i, m) in masks.iter_mut().enumerate() {
                    if b & (1 << i) != 0 {
                        *m |= 1u64 << lane;
                    }
                }
                d[lane] | (b << 3)
            });
            assert_eq!(
                unsafe { run(&d, masks) },
                expected,
                "direction {direction} mask {bits}"
            );
        }
    }
}

/// `trace_plane_argmax` + `store_four_trace_planes` decode (with the traceback helpers) to the
/// byte layout's direction and continuation bits, including all ties, cached and streamed.
#[test]
#[allow(clippy::needless_range_loop)] // lanes index several parallel arrays
fn trace_planes_decode_to_byte_directions_with_ties() {
    if !available() {
        return;
    }
    #[target_feature(enable = "avx2,avx512f,avx512bw")]
    unsafe fn run<const RIGHT: bool>(
        input: &[[i8; 64]; 5],
        gaps: [u64; 4],
        nt: bool,
    ) -> ([i8; 64], [u8; 64]) {
        unsafe {
            let v = [
                _mm512_loadu_si512(input[0].as_ptr().cast()),
                _mm512_loadu_si512(input[1].as_ptr().cast()),
                _mm512_loadu_si512(input[2].as_ptr().cast()),
                _mm512_loadu_si512(input[3].as_ptr().cast()),
                _mm512_loadu_si512(input[4].as_ptr().cast()),
            ];
            let (z, m0, m1, m2, m3) = trace_plane_argmax!(v[0], v[1], v[2], v[3], v[4], RIGHT);
            let mem = AlignedMemory::new(64, 64);
            let p = mem.as_ptr();
            store_four_trace_planes!(p, nt, m0, m1, m2, m3);
            store_four_trace_planes!(p.add(32), nt, gaps[0], gaps[1], gaps[2], gaps[3]);
            if nt {
                _mm_sfence();
            }
            let mut scores = [0i8; 64];
            _mm512_storeu_si512(scores.as_mut_ptr().cast(), z);
            let flags: [u8; 64] = std::array::from_fn(|lane| {
                let mut d = trace_plane_direction64(p, lane) as u8;
                for bit in 0..4 {
                    if trace_plane_bit64(p, lane, 4 + bit) {
                        d |= 1 << (3 + bit);
                    }
                }
                d
            });
            (scores, flags)
        }
    }
    let mut rng = Rng(0x5128_71ac_9123_abef);
    for case in 0..8192 {
        let mut input: [[i8; 64]; 5] =
            std::array::from_fn(|_| std::array::from_fn(|_| rng.next() as i8));
        if case % 4 == 0 {
            for lane in 0..64 {
                for state in 1..5 {
                    input[state][lane] = input[0][lane];
                }
            }
        }
        if case % 7 == 0 {
            input[2] = input[1];
            input[4] = input[3];
        }
        if case % 11 == 0 {
            for state in &mut input {
                for v in state {
                    *v = rng.pick(&[i8::MIN, 0, i8::MAX]);
                }
            }
        }
        let gaps = if case % 13 == 0 {
            [0, u64::MAX, 0xaaaa_aaaa_aaaa_aaaa, 0x5555_5555_5555_5555]
        } else {
            std::array::from_fn(|_| rng.next())
        };
        for right in [false, true] {
            // Sequential argmax, as the byte-direction code computes it.
            let mut scores = [0i8; 64];
            let mut flags = [0u8; 64];
            for lane in 0..64 {
                let mut winner = 0usize;
                for state in 1..5 {
                    if input[state][lane] > input[winner][lane]
                        || (right && input[state][lane] == input[winner][lane])
                    {
                        winner = state;
                    }
                }
                scores[lane] = input[winner][lane];
                flags[lane] = winner as u8;
                for bit in 0..4 {
                    if (gaps[bit] >> lane) & 1 != 0 {
                        flags[lane] |= 1 << (3 + bit);
                    }
                }
            }
            for nt in [false, true] {
                let got = if right {
                    unsafe { run::<true>(&input, gaps, nt) }
                } else {
                    unsafe { run::<false>(&input, gaps, nt) }
                };
                assert_eq!(got, (scores, flags), "case {case} right {right} nt {nt}");
            }
        }
    }
}

/// `max(v, 0) - qe` (array kernel) equals the masked form `maskz(v > 0 or v >= 0, v) - qe`
/// (register kernels) for every byte, either alignment side.
#[test]
fn signed_gap_clip_matches_mask_selection_for_every_byte_and_alignment_side() {
    if !available() {
        return;
    }
    #[target_feature(enable = "avx2,avx512f,avx512bw")]
    unsafe fn run(bytes: &[i8; 64], qe: i8, right: bool) -> ([i8; 64], [i8; 64]) {
        unsafe {
            let v = _mm512_loadu_si512(bytes.as_ptr().cast());
            let z = _mm512_setzero_si512();
            let m = if right {
                !_mm512_cmpgt_epi8_mask(z, v)
            } else {
                _mm512_cmpgt_epi8_mask(v, z)
            };
            let masked = _mm512_sub_epi8(_mm512_maskz_mov_epi8(m, v), _mm512_set1_epi8(qe));
            let clipped = _mm512_sub_epi8(_mm512_max_epi8(v, z), _mm512_set1_epi8(qe));
            let mut a = [0i8; 64];
            let mut b = [0i8; 64];
            _mm512_storeu_si512(a.as_mut_ptr().cast(), masked);
            _mm512_storeu_si512(b.as_mut_ptr().cast(), clipped);
            (a, b)
        }
    }
    for right in [false, true] {
        for qe in [-128i8, -127, -25, -6, -1, 0, 1, 6, 25, 126, 127] {
            for group in 0..4 {
                let bytes: [i8; 64] = std::array::from_fn(|lane| (group * 64 + lane) as u8 as i8);
                let (masked, clipped) = unsafe { run(&bytes, qe, right) };
                let expected: [i8; 64] =
                    std::array::from_fn(|lane| bytes[lane].max(0).wrapping_sub(qe));
                assert_eq!(masked, expected);
                assert_eq!(clipped, expected);
            }
        }
    }
}

/// Over the operand ranges the in-kernel CHECK asserts (-10..=4 short, -9..=24 long), the
/// shortened form `max(v - qe, -qe)` and its threshold mask equal `max(v, 0) - qe` and `v > 0`
/// (`v >= 0` for right alignment).
#[test]
fn shifted_gap_clipping_and_ties_match_scalar_over_active_ranges() {
    if !available() {
        return;
    }
    #[target_feature(enable = "avx2,avx512f,avx512bw")]
    unsafe fn run(an: &[i8; 64], qe: i8, right: bool) -> ([i8; 64], u64) {
        unsafe {
            let old = _mm512_loadu_si512(an.as_ptr().cast());
            let shifted = _mm512_sub_epi8(old, _mm512_set1_epi8(qe));
            let floor = _mm512_set1_epi8(-qe);
            let next = _mm512_max_epi8(shifted, floor);
            let mask = if right {
                !_mm512_cmpgt_epi8_mask(floor, shifted)
            } else {
                _mm512_cmpgt_epi8_mask(shifted, floor)
            };
            let mut out = [0i8; 64];
            _mm512_storeu_si512(out.as_mut_ptr().cast(), next);
            (out, mask)
        }
    }
    for (low, high, qe) in [(-10i8, 4i8, 6i8), (-9, 24, 25)] {
        for start in low..=high {
            let span = (high - low + 1) as usize;
            let values: [i8; 64] =
                std::array::from_fn(|i| low + (((start - low) as usize + i) % span) as i8);
            for right in [false, true] {
                let expected: [i8; 64] = std::array::from_fn(|i| values[i].max(0).wrapping_sub(qe));
                let expected_mask = values.iter().enumerate().fold(0u64, |m, (i, &v)| {
                    m | ((if right { v >= 0 } else { v > 0 }) as u64) << i
                });
                assert_eq!(
                    unsafe { run(&values, qe, right) },
                    (expected, expected_mask)
                );
            }
        }
    }
}

#[test]
fn unmasked_h16_blocks_preserve_four_lane_ties_and_unaligned_score_updates() {
    if !available() {
        return;
    }
    #[target_feature(enable = "avx2,avx512f,avx512bw")]
    unsafe fn run(h: &mut [i32], v: &[u8], pad: usize, n: usize) -> [i32; 16] {
        unsafe {
            let mut low = _mm512_setr_epi32(
                32767, 24575, 16383, 8191, 32766, 24574, 16382, 8190, 32765, 24573, 16381, 8189,
                32764, 24572, 16380, 8188,
            );
            let four = _mm512_set1_epi32(4);
            let mut mk = _mm512_set1_epi32(i32::MIN);
            for at in (0..n).step_by(16) {
                mk = update_h16_packed(
                    h.as_mut_ptr().add(pad + at),
                    v.as_ptr().add(pad + at),
                    low,
                    mk,
                );
                low = _mm512_sub_epi32(low, four);
            }
            let mut keys = [0i32; 16];
            _mm512_storeu_si512(keys.as_mut_ptr().cast(), mk);
            keys
        }
    }
    let mut rng = Rng(0x16ff_5120_8423_12af);
    for case in 0..512 {
        let n = [1, 2, 3, 7, 31, 127, 500][case % 7] * 16;
        let pad = [0, 1, 3, 7, 17][case % 5];
        let mut h = vec![0x5a5a5a5ai32; n + pad + 16];
        let mut v = vec![0u8; n + pad + 16];
        for at in 0..n {
            h[pad + at] = if case % 3 == 0 {
                17
            } else if case % 5 == 0 {
                rng.below(5) as i32 - 2
            } else {
                rng.below(80001) as i32 - 40000
            };
            v[pad + at] = if case % 3 == 0 || case % 5 == 0 {
                0
            } else {
                rng.next() as u8
            };
        }
        if case % 11 == 0 {
            for at in 0..n {
                h[pad + at] = -100;
                v[pad + at] = 0;
            }
            h[pad + n - 4] = 100;
            h[pad + 1] = 100;
        } else if case % 13 == 0 {
            for at in 0..n {
                h[pad + at] = -17;
                v[pad + at] = 0;
            }
        }
        let mut expected = h.clone();
        let mut lane_h = [i32::MIN; 4];
        let mut lane_t = [0usize; 4];
        for at in 0..n {
            expected[pad + at] = expected[pad + at].wrapping_add(v[pad + at] as i8 as i32);
            let lane = at % 4;
            if expected[pad + at] > lane_h[lane] {
                lane_h[lane] = expected[pad + at];
                lane_t[lane] = at;
            }
        }
        let mut winner = 0;
        for lane in 1..4 {
            if lane_h[lane] > lane_h[winner] {
                winner = lane;
            }
        }
        let keys = unsafe { run(&mut h, &v, pad, n) };
        assert_eq!(h, expected, "score updates, case {case}");
        let key = *keys.iter().max().unwrap();
        let rank = 32767 - (key & 0x7fff);
        let position = ((rank & 0x1fff) * 4 + (rank >> 13)) as usize;
        assert_eq!(key >> 15, lane_h[winner], "score, case {case}");
        assert_eq!(position, lane_t[winner], "tie position, case {case}");
    }
}
