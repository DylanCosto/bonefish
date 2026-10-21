//! AVX-512 bitonic sorting networks for 16..=64-element radix-sort leaves.
//!
//! Leaves are sorted as 32-bit composites `(low key bits << 6) | input index`. The composites
//! are distinct, so the result is the stable order that `Minimizer::leaf_sort` defines.
use super::Minimizer;
use core::arch::x86_64::*;

/// Sorts `arr` stably by `x`. The caller ensures AVX-512F is available, that
/// `16 <= arr.len() <= 64` and `s <= 26`, and that all keys agree above bit `s`.
#[target_feature(enable = "avx512f")]
pub(super) unsafe fn sort_leaf(arr: &mut [Minimizer], s: u32) {
    unsafe {
        let n = arr.len();
        debug_assert!((16..=64).contains(&n) && s <= 26);
        let mask = (1u64 << s) - 1;
        let mut keys = [u32::MAX; 64];
        for (i, item) in arr.iter().enumerate() {
            keys[i] = (((item.x & mask) << 6) | i as u64) as u32;
        }
        // Padding (u32::MAX) sorts after every real key: a real composite can equal it only at
        // index 63, i.e. when n == 64 and there is no padding.
        match n {
            16 => sort_512_16(&mut keys),
            17..=32 => sort_512_32(&mut keys),
            _ => sort_512_64(&mut keys),
        }
        let mut tmp = [Minimizer { x: 0, y: 0 }; 64];
        for (i, &key) in keys[..n].iter().enumerate() {
            tmp[i] = arr[(key & 63) as usize];
        }
        arr.copy_from_slice(&tmp[..n]);
    }
}

// Bitonic networks over 16, 32 and 64 keys (generated). In each compare-exchange stage element
// i is paired with i ^ j and keeps the maximum iff (i & j == 0) != (i & k == 0); k doubles
// from stage group to stage group and j halves within each group.
#[target_feature(enable = "avx512f")]
unsafe fn sort_512_16(a: &mut [u32; 64]) {
    unsafe {
        let mut v0 = _mm512_loadu_si512(a.as_ptr().add(0) as *const _);
        let ix1 = _mm512_set_epi32(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
        let ix2 = _mm512_set_epi32(13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2);
        let ix4 = _mm512_set_epi32(11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4);
        let ix8 = _mm512_set_epi32(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 26214, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 15420, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 23130, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 4080, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 13260, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 21930, v0, t);
        let t = _mm512_permutexvar_epi32(ix8, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 65280, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 61680, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 52428, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 43690, v0, t);
        _mm512_storeu_si512(a.as_mut_ptr().add(0) as *mut _, v0);
    }
}
#[target_feature(enable = "avx512f")]
unsafe fn sort_512_32(a: &mut [u32; 64]) {
    unsafe {
        let mut v0 = _mm512_loadu_si512(a.as_ptr().add(0) as *const _);
        let mut v1 = _mm512_loadu_si512(a.as_ptr().add(16) as *const _);
        let ix1 = _mm512_set_epi32(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
        let ix2 = _mm512_set_epi32(13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2);
        let ix4 = _mm512_set_epi32(11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4);
        let ix8 = _mm512_set_epi32(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 26214, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 26214, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 15420, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 15420, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 23130, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 23130, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 4080, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 4080, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 13260, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 13260, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 21930, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 21930, v1, t);
        let t = _mm512_permutexvar_epi32(ix8, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 65280, v0, t);
        let t = _mm512_permutexvar_epi32(ix8, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 255, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 61680, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 3855, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 52428, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 13107, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 43690, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 21845, v1, t);
        let lo = _mm512_min_epu32(v0, v1);
        let hi = _mm512_max_epu32(v0, v1);
        v0 = lo;
        v1 = hi;
        let t = _mm512_permutexvar_epi32(ix8, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 65280, v0, t);
        let t = _mm512_permutexvar_epi32(ix8, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 65280, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 61680, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 61680, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 52428, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 52428, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 43690, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 43690, v1, t);
        _mm512_storeu_si512(a.as_mut_ptr().add(0) as *mut _, v0);
        _mm512_storeu_si512(a.as_mut_ptr().add(16) as *mut _, v1);
    }
}
#[target_feature(enable = "avx512f")]
unsafe fn sort_512_64(a: &mut [u32; 64]) {
    unsafe {
        let mut v0 = _mm512_loadu_si512(a.as_ptr().add(0) as *const _);
        let mut v1 = _mm512_loadu_si512(a.as_ptr().add(16) as *const _);
        let mut v2 = _mm512_loadu_si512(a.as_ptr().add(32) as *const _);
        let mut v3 = _mm512_loadu_si512(a.as_ptr().add(48) as *const _);
        let ix1 = _mm512_set_epi32(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1);
        let ix2 = _mm512_set_epi32(13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2);
        let ix4 = _mm512_set_epi32(11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4);
        let ix8 = _mm512_set_epi32(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 26214, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 26214, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 26214, v2, t);
        let t = _mm512_permutexvar_epi32(ix1, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 26214, v3, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 15420, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 15420, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 15420, v2, t);
        let t = _mm512_permutexvar_epi32(ix2, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 15420, v3, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 23130, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 23130, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 23130, v2, t);
        let t = _mm512_permutexvar_epi32(ix1, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 23130, v3, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 4080, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 4080, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 4080, v2, t);
        let t = _mm512_permutexvar_epi32(ix4, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 4080, v3, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 13260, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 13260, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 13260, v2, t);
        let t = _mm512_permutexvar_epi32(ix2, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 13260, v3, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 21930, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 21930, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 21930, v2, t);
        let t = _mm512_permutexvar_epi32(ix1, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 21930, v3, t);
        let t = _mm512_permutexvar_epi32(ix8, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 65280, v0, t);
        let t = _mm512_permutexvar_epi32(ix8, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 255, v1, t);
        let t = _mm512_permutexvar_epi32(ix8, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 65280, v2, t);
        let t = _mm512_permutexvar_epi32(ix8, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 255, v3, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 61680, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 3855, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 61680, v2, t);
        let t = _mm512_permutexvar_epi32(ix4, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 3855, v3, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 52428, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 13107, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 52428, v2, t);
        let t = _mm512_permutexvar_epi32(ix2, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 13107, v3, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 43690, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 21845, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 43690, v2, t);
        let t = _mm512_permutexvar_epi32(ix1, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 21845, v3, t);
        let lo = _mm512_min_epu32(v0, v1);
        let hi = _mm512_max_epu32(v0, v1);
        v0 = lo;
        v1 = hi;
        let lo = _mm512_min_epu32(v2, v3);
        let hi = _mm512_max_epu32(v2, v3);
        v3 = lo;
        v2 = hi;
        let t = _mm512_permutexvar_epi32(ix8, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 65280, v0, t);
        let t = _mm512_permutexvar_epi32(ix8, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 65280, v1, t);
        let t = _mm512_permutexvar_epi32(ix8, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 255, v2, t);
        let t = _mm512_permutexvar_epi32(ix8, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 255, v3, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 61680, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 61680, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 3855, v2, t);
        let t = _mm512_permutexvar_epi32(ix4, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 3855, v3, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 52428, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 52428, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 13107, v2, t);
        let t = _mm512_permutexvar_epi32(ix2, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 13107, v3, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 43690, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 43690, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 21845, v2, t);
        let t = _mm512_permutexvar_epi32(ix1, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 21845, v3, t);
        let lo = _mm512_min_epu32(v0, v2);
        let hi = _mm512_max_epu32(v0, v2);
        v0 = lo;
        v2 = hi;
        let lo = _mm512_min_epu32(v1, v3);
        let hi = _mm512_max_epu32(v1, v3);
        v1 = lo;
        v3 = hi;
        let lo = _mm512_min_epu32(v0, v1);
        let hi = _mm512_max_epu32(v0, v1);
        v0 = lo;
        v1 = hi;
        let lo = _mm512_min_epu32(v2, v3);
        let hi = _mm512_max_epu32(v2, v3);
        v2 = lo;
        v3 = hi;
        let t = _mm512_permutexvar_epi32(ix8, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 65280, v0, t);
        let t = _mm512_permutexvar_epi32(ix8, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 65280, v1, t);
        let t = _mm512_permutexvar_epi32(ix8, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 65280, v2, t);
        let t = _mm512_permutexvar_epi32(ix8, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 65280, v3, t);
        let t = _mm512_permutexvar_epi32(ix4, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 61680, v0, t);
        let t = _mm512_permutexvar_epi32(ix4, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 61680, v1, t);
        let t = _mm512_permutexvar_epi32(ix4, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 61680, v2, t);
        let t = _mm512_permutexvar_epi32(ix4, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 61680, v3, t);
        let t = _mm512_permutexvar_epi32(ix2, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 52428, v0, t);
        let t = _mm512_permutexvar_epi32(ix2, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 52428, v1, t);
        let t = _mm512_permutexvar_epi32(ix2, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 52428, v2, t);
        let t = _mm512_permutexvar_epi32(ix2, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 52428, v3, t);
        let t = _mm512_permutexvar_epi32(ix1, v0);
        v0 = _mm512_mask_max_epu32(_mm512_min_epu32(v0, t), 43690, v0, t);
        let t = _mm512_permutexvar_epi32(ix1, v1);
        v1 = _mm512_mask_max_epu32(_mm512_min_epu32(v1, t), 43690, v1, t);
        let t = _mm512_permutexvar_epi32(ix1, v2);
        v2 = _mm512_mask_max_epu32(_mm512_min_epu32(v2, t), 43690, v2, t);
        let t = _mm512_permutexvar_epi32(ix1, v3);
        v3 = _mm512_mask_max_epu32(_mm512_min_epu32(v3, t), 43690, v3, t);
        _mm512_storeu_si512(a.as_mut_ptr().add(0) as *mut _, v0);
        _mm512_storeu_si512(a.as_mut_ptr().add(16) as *mut _, v1);
        _mm512_storeu_si512(a.as_mut_ptr().add(32) as *mut _, v2);
        _mm512_storeu_si512(a.as_mut_ptr().add(48) as *mut _, v3);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rand(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }
    #[test]
    fn avx512_leaf_network_preserves_all_fields_ties_and_padding() {
        if !crate::align::dp::simd_test_available("avx512bw", is_x86_feature_detected!("avx512bw"))
        {
            return;
        }
        let mut seed = 0x5126_aaaa_1357_2468;
        for n in 16..=64 {
            for s in [0, 1, 8, 16, 24, 26] {
                for case in 0..256 {
                    let mask = (1u64 << s) - 1;
                    let upper = rand(&mut seed) & !mask;
                    let mut a: Vec<_> = (0..n)
                        .map(|i| {
                            let low = match case % 7 {
                                0 => 0,
                                1 => mask,
                                2 => i as u64 & mask,
                                3 => (n - i) as u64 & mask,
                                4 => rand(&mut seed) & 7 & mask,
                                _ => rand(&mut seed) & mask,
                            };
                            Minimizer {
                                x: upper | low,
                                y: rand(&mut seed),
                            }
                        })
                        .collect();
                    let mut expected = a.clone();
                    expected.sort_by_key(|m| m.x);
                    unsafe { sort_leaf(&mut a, s) };
                    assert_eq!(a, expected, "n={n} s={s} case={case}");
                }
            }
        }
    }
}
