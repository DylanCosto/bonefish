//! Minimizer sketching for DNA sequences.
//!
//! Computes (w,k)-minimizers: for each window of `w` consecutive k-mers, the
//! k-mer with the smallest hash is selected as a minimizer. Both forward and
//! reverse-complement k-mers are hashed; the canonical (smaller) strand is kept.
//!
//! Each minimizer is a 128-bit [`Minimizer`] with two packed fields:
//! - **x** (reference/hash): bits 0-7 kmer_span, bits 8+ minimizer hash.
//!   After seed collection, repacked as: bits 0-31 ref_pos, bit 32+ ref_id|strand.
//! - **y** (query): bits 0-31 query_pos (bit 0 = strand), bits 32-39 kmer_span,
//!   bits 48-55 segment ID, bit 43 self-mapping flag (ava mode).
//!
//! [`sketch_sequence`] clears and fills an output `Vec<Minimizer>`.


/// Anchor bit-field encoding constants for the `y` field of `Minimizer`.
/// Segment ID occupies bits 48-55 of y.
pub const SEED_SEG_SHIFT: u64 = 48;
pub const SEED_SEG_MASK: u64  = 0xFFu64 << SEED_SEG_SHIFT;
/// Flag for self-mapping anchors in overlap (ava) mode (bit 43 of y).
pub const SEED_SELF: u64      = 1u64 << 43;
/// Flag: anchor marked as long-join bridge between distant seeds (bit 40 of y).
pub const SEED_LONG_JOIN: u64 = 1u64 << 40;
/// Flag: anchor ignored during gap-fill alignment (bit 41 of y).
pub const SEED_IGNORE: u64    = 1u64 << 41;
/// Flag: anchor is part of a tandem repeat region (bit 42 of y).
pub const SEED_TANDEM: u64    = 1u64 << 42;

/// A minimizer anchor: packed 128-bit representation of a (reference, query) seed match.
///
/// **x field** (reference): `ref_id << 33 | ref_pos << 1 | strand` (upper) | position (lower 32)
/// **y field** (query): `seg_id << 48 | ... | kmer_span << 32 | query_pos << 1 | strand` (lower 32)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
  pub x: u64,
  pub y: u64,
}

pub type Mm128 = Minimizer;

impl Minimizer {
    /// Reference position (lower 32 bits of x, as signed integer).
    #[inline(always)]
    pub fn ref_pos(&self) -> i32 {
        self.x as u32 as i32
    }

    /// Set reference position (lower 32 bits of x), preserving upper bits (rid/strand).
    #[inline(always)]
    pub fn set_ref_pos(&mut self, pos: i32) {
        self.x = (self.x & 0xFFFFFFFF_00000000) | (pos as u32 as u64);
    }

    /// Reference ID and strand bits combined (upper 32 bits of x).
    /// Used for same-target comparisons without separating rid/strand.
    #[inline(always)]
    pub fn ref_id_strand(&self) -> u64 {
        self.x >> 32
    }

    /// Reference sequence ID (bits 33-62 of x, without strand bit).
    #[inline(always)]
    pub fn ref_id(&self) -> usize {
        ((self.x >> 32) & 0x7FFFFFFF) as usize
    }

    /// Query position (lower 32 bits of y, as signed integer).
    #[inline(always)]
    pub fn query_pos(&self) -> i32 {
        self.y as i32
    }

    /// K-mer span (bits 32-39 of y).
    #[inline(always)]
    pub fn query_span(&self) -> i32 {
        ((self.y >> 32) & 0xff) as i32
    }

    /// Segment ID for multi-segment reads (bits 48-55 of y).
    #[inline(always)]
    pub fn segment_id(&self) -> usize {
        ((self.y & SEED_SEG_MASK) >> SEED_SEG_SHIFT) as usize
    }

    /// Whether this is a self-mapping anchor in overlap (ava) mode.
    #[inline(always)]
    pub fn is_self_anchor(&self) -> bool {
        (self.y & SEED_SELF) != 0
    }

    /// Whether this anchor is ignored during gap-fill alignment.
    #[inline(always)]
    pub fn is_ignored(&self) -> bool {
        (self.y & SEED_IGNORE) != 0
    }

    /// Whether this anchor is part of a tandem repeat region.
    #[inline(always)]
    pub fn is_tandem(&self) -> bool {
        (self.y & SEED_TANDEM) != 0
    }

    /// Whether this anchor is a long-join bridge between distant seeds.
    #[inline(always)]
    pub fn is_long_join(&self) -> bool {
        (self.y & SEED_LONG_JOIN) != 0
    }
}

/// Trait for sequence sketching strategies.
///
/// A sketcher converts a DNA sequence into a set of [`Minimizer`] entries
/// that will be used for seed lookup in the index. Different strategies
/// select different subsets of k-mers (e.g., minimizers, syncmers).
///
/// The output `Minimizer` values must use the standard packing:
/// - `x`: `hash << 8 | kmer_span` (hash from [`kmer_hash`])
/// - `y`: `position << 1 | strand` (strand: 0=forward, 1=reverse)
pub trait Sketcher {
    /// Sketch a sequence, clearing and filling the output vector.
    fn sketch(&self, seq: &[u8], len: usize, rid: usize, out: &mut Vec<Minimizer>);
}

/// Lookup table for base encoding: A=0, C=1, G=2, T=3, other=4
/// Using a 256-byte table for branchless lookup (cache-friendly)
static BASE_ENCODING: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

#[inline(always)]
pub fn encode_base(b: u8) -> u8 {
    BASE_ENCODING[b as usize]
}

/// Invertible 64-bit hash for k-mer bits. Reduces nucleotide composition bias.
pub fn kmer_hash(mut bits: u64, mask: u64) -> u64 {
  bits = (!bits).wrapping_add(bits << 21) & mask; // bits = (bits << 21) - bits - 1;
  bits = bits ^ (bits >> 24);
  bits = (bits.wrapping_add(bits << 3)).wrapping_add(bits << 8) & mask; // bits * 265
  bits = bits ^ (bits >> 14);
  bits = (bits.wrapping_add(bits << 2)).wrapping_add(bits << 4) & mask; // bits * 21
  bits = bits ^ (bits >> 28);
  bits.wrapping_add(bits << 31) & mask
}

/// Standard (w,k)-minimizer sketcher. Wraps [`sketch_sequence`] to implement
/// the [`Sketcher`] trait.
pub struct MinimizerSketcher {
    pub w: usize,
    pub k: usize,
}

impl Sketcher for MinimizerSketcher {
    fn sketch(&self, seq: &[u8], len: usize, rid: usize, out: &mut Vec<Minimizer>) {
        sketch_sequence(seq, len, self.w, self.k, rid, out);
    }
}

/// Compute minimizers for a DNA sequence and write them to the output vector.
/// Clears `p` before writing.
pub fn sketch_sequence(
  seq: &[u8],
  len: usize,
  w: usize,
  k: usize,
  rid: usize,
  p: &mut Vec<Minimizer>,
) {
    p.clear();
    sketch_sequence_impl(seq, len, w, k, rid, p);
}

/// Core minimizer sketching algorithm. Appends minimizers to `p`.
/// sketch_sequence (clear + sketch) delegates here.
fn sketch_sequence_impl(
  seq: &[u8],
  len: usize,
  w: usize,
  k: usize,
  rid: usize,
  p: &mut Vec<Minimizer>,
) {
    #[cfg(target_arch = "x86_64")]
    if len > 0 && w > 0 && w < 256 && k > 0 && k <= 28 && is_x86_feature_detected!("avx2") && !*SKETCH_SCALAR {
        if len >= 1024 && k % 2 == 1 && super::dp::use_avx512() && !*SKETCH_AVX512_OFF {
            unsafe { sketch512::sketch(seq, len, w, k, rid, p) };
            return;
        }
        unsafe { sketch_chunked_avx2(seq, len, w, k, rid, p) };
        return;
    }
    sketch_sequence_scalar(seq, len, w, k, rid, p);
}

/// BONEFISH_SKETCH_SCALAR=1: always use the one-pass sketcher (comparisons).
static SKETCH_SCALAR: std::sync::LazyLock<bool> =
    std::sync::LazyLock::new(|| std::env::var("BONEFISH_SKETCH_SCALAR").is_ok_and(|v| v != "0"));

/// Sketching in chunks: (1) roll the k-mers and record, per position that the one-pass
/// loop would process (symmetric k-mers are skipped there too), the canonical k-mer, strand,
/// span and k-mer count `l`; (2) hash the chunk's k-mers in one tight loop (vectorised);
/// (3) run the one-pass loop's window logic on those values. Same minimizers, same order.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn sketch_chunked_avx2(seq: &[u8], len: usize, w: usize, k: usize, rid: usize, p: &mut Vec<Minimizer>) { unsafe {
    const CH: usize = 512;
    let shift1 = 2 * (k - 1);
    let mask = (1u64 << (2 * k)) - 1;
    let mut kmer: [u64; 2] = [0, 0];
    let mut l: usize = 0;
    let mut buf_pos: usize = 0;
    let mut min_pos: usize = 0;
    let mut buf: [Minimizer; 256] = [Minimizer { x: 0, y: 0 }; 256];
    let mut min: Minimizer = Minimizer { x: u64::MAX, y: u64::MAX };
    // Per event of the chunk: canonical k-mer (hashed in place), y, span, l; valid = l >= k.
    let mut ek = [0u64; CH];
    let mut ey = [0u64; CH];
    let mut esp = [0u8; CH];
    let mut el = [0u32; CH];
    let mut i = 0usize;
    while i < len {
        let end = (i + CH).min(len);
        let mut n = 0usize;
        // (1) serial roll
        while i < end {
            let c = encode_base(*seq.get_unchecked(i));
            if c < 4 {
                let span = if l + 1 < k { l + 1 } else { k };
                kmer[0] = (kmer[0] << 2 | c as u64) & mask;
                kmer[1] = (kmer[1] >> 2) | (3u64 ^ c as u64) << shift1;
                if kmer[0] != kmer[1] {
                    let z = (kmer[0] >= kmer[1]) as usize;
                    l += 1;
                    *ek.get_unchecked_mut(n) = kmer[z];
                    *ey.get_unchecked_mut(n) = (rid as u64) << 32 | (i as u64) << 1 | z as u64;
                    *esp.get_unchecked_mut(n) = span as u8;
                    *el.get_unchecked_mut(n) = l as u32;
                    n += 1;
                }
            } else {
                l = 0;
                *el.get_unchecked_mut(n) = 0;
                *esp.get_unchecked_mut(n) = 0;
                n += 1;
            }
            i += 1;
        }
        // (2) hashes
        for j in 0..n {
            let mut bits = *ek.get_unchecked(j);
            bits = (!bits).wrapping_add(bits << 21) & mask;
            bits = bits ^ (bits >> 24);
            bits = (bits.wrapping_add(bits << 3)).wrapping_add(bits << 8) & mask;
            bits = bits ^ (bits >> 14);
            bits = (bits.wrapping_add(bits << 2)).wrapping_add(bits << 4) & mask;
            bits = bits ^ (bits >> 28);
            *ek.get_unchecked_mut(j) = bits.wrapping_add(bits << 31) & mask;
        }
        // (3) window logic, as in sketch_sequence_scalar
        for j in 0..n {
            let lj = *el.get_unchecked(j) as usize;
            let info = if lj >= k {
                Minimizer { x: *ek.get_unchecked(j) << 8 | *esp.get_unchecked(j) as u64, y: *ey.get_unchecked(j) }
            } else {
                Minimizer { x: u64::MAX, y: u64::MAX }
            };
            *buf.get_unchecked_mut(buf_pos) = info;
            if lj == (w + k - 1) && min.x != u64::MAX {
                for item in &buf[(buf_pos + 1)..w] {
                    if min.x == item.x && item.y != min.y { p.push(*item); }
                }
                for item in &buf[..buf_pos] {
                    if min.x == item.x && item.y != min.y { p.push(*item); }
                }
            }
            if info.x <= min.x {
                if lj >= (w + k) && min.x != u64::MAX { p.push(min); }
                min = info;
                min_pos = buf_pos;
            } else if buf_pos == min_pos {
                if lj >= (w + k - 1) && min.x != u64::MAX { p.push(min); }
                min.x = u64::MAX;
                for (jj, item) in buf[(buf_pos + 1)..w].iter().enumerate() {
                    if min.x == u64::MAX || min.x >= item.x { min = *item; min_pos = buf_pos + 1 + jj; }
                }
                for (jj, item) in buf[..=buf_pos].iter().enumerate() {
                    if min.x == u64::MAX || min.x >= item.x { min = *item; min_pos = jj; }
                }
                if lj >= (w + k - 1) && min.x != u64::MAX {
                    for item in &buf[(buf_pos + 1)..w] {
                        if min.x == item.x && min.y != item.y { p.push(*item); }
                    }
                    for item in &buf[..=buf_pos] {
                        if min.x == item.x && min.y != item.y { p.push(*item); }
                    }
                }
            }
            if buf_pos + 1 == w { buf_pos = 0; } else { buf_pos += 1; }
        }
    }
    if min.x != u64::MAX {
        p.push(min);
    }
}}

/// The one-pass sketcher (for comparisons with the chunked path).
pub fn sketch_sequence_reference(seq: &[u8], len: usize, w: usize, k: usize, rid: usize, p: &mut Vec<Minimizer>) {
    p.clear();
    sketch_sequence_scalar(seq, len, w, k, rid, p);
}

/// One-pass sketcher (reference for the chunked path).
fn sketch_sequence_scalar(
  seq: &[u8],
  len: usize,
  w: usize,
  k: usize,
  rid: usize,
  p: &mut Vec<Minimizer>,
) {
    if len == 0 {
        return;
    }
    assert!(w > 0 && w < 256 && k > 0 && k <= 28);

    let shift1 = 2 * (k - 1);
    let mask = (1u64 << (2 * k)) - 1;
    let mut kmer: [u64; 2] = [0, 0];
    let mut i: usize = 0;
    let mut l: usize = 0;
    let mut buf_pos: usize = 0;
    let mut min_pos: usize = 0;
    let mut buf: Vec<Minimizer> = vec![Minimizer { x: 0, y: 0 }; 256];
    let mut min: Minimizer = Minimizer { x: u64::MAX, y: u64::MAX };

    while i < len {
        let c = encode_base(seq[i]);
        let mut info = Minimizer { x: u64::MAX, y: u64::MAX };
        if c < 4 { // not an ambiguous base
            let kmer_span = if l + 1 < k { l + 1 } else { k };
            kmer[0] = (kmer[0] << 2 | c as u64) & mask;           // forward k-mer
            kmer[1] = (kmer[1] >> 2) | (3u64 ^ c as u64) << shift1; // reverse k-mer
            if kmer[0] == kmer[1] { // skip symmetric k-mers
                i += 1;
                continue;
            }
            let z: usize = if kmer[0] < kmer[1] { 0 } else { 1 }; // strand
            l += 1;
            if l >= k && kmer_span < 256 {
                info.x = kmer_hash(kmer[z], mask) << 8 | kmer_span as u64;
                info.y = (rid as u64) << 32 | (i as u64) << 1 | z as u64;
            }
        } else {
            l = 0;
        }

        buf[buf_pos] = info;

        // Special case for the first window — push identical k-mers
        if l == (w + k - 1) && min.x != u64::MAX {
            for item in &buf[(buf_pos + 1)..w] {
                if min.x == item.x && item.y != min.y { p.push(*item); }
            }
            for item in &buf[..buf_pos] {
                if min.x == item.x && item.y != min.y { p.push(*item); }
            }
        }

        if info.x <= min.x { // a new minimum; write the old min
            if l >= (w + k) && min.x != u64::MAX { p.push(min); }
            min = info;
            min_pos = buf_pos;
        } else if buf_pos == min_pos { // old min has moved outside the window
            if l >= (w + k - 1) && min.x != u64::MAX { p.push(min); }
            min.x = u64::MAX;
            for (j, item) in buf[(buf_pos + 1)..w].iter().enumerate() {
                if min.x == u64::MAX || min.x >= item.x { min = *item; min_pos = buf_pos + 1 + j; }
            }
            for (j, item) in buf[..=buf_pos].iter().enumerate() {
                if min.x == u64::MAX || min.x >= item.x { min = *item; min_pos = j; }
            }
            if l >= (w + k - 1) && min.x != u64::MAX { // write identical k-mers
                for item in &buf[(buf_pos + 1)..w] {
                    if min.x == item.x && min.y != item.y { p.push(*item); }
                }
                for item in &buf[..=buf_pos] {
                    if min.x == item.x && min.y != item.y { p.push(*item); }
                }
            }
        }

        if buf_pos + 1 == w { buf_pos = 0; } else { buf_pos += 1; }
        i += 1;
    }
    if min.x != u64::MAX {
        p.push(min);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sketch_basic() {
        let seq = "GTTGATAATCACTCACTGAGTGACATCCAAATCATGGCGTCCCAAGGCACCAAACGGTCTTATGAACAGATGGAAACTGATGGAGATCGCCAGAATGCAACTGAGATTAGGGCATCCGTCGGAAAGATGATTGATGGAATTGGGAGATTCTACATCCAAATGTGCACTGAACTTAAACTCAGTGATCATGAAGGACGGTTGATCCAAAACAGCTTGACAATAGAGAAAATGGTGCTTTCTGCTTTTGATGAAAGAAGGAATAAATACCTGGAAGAACACCCCAGCGCGGGGAAAGATCCCAAGAAAACCGGGGGGCCCATATACAGGAGAGTCGATGGGAAATGGATGAGAGAACTCGTCCTTTATGACAAAGAAGAAATAAGGCGAATCTGGCGCCAAGCCAACAATGGTGAGGATGCTACATCTGGTCTAACCCACCTAATGATTTGGCATTCCAATTTGAATGATGCAACATACCAAAGGACAAGAGCTCTTGTTCGGACTGGAATGGACCCCAGAATGTGCTCTCTGATGCAGGGCTCGACTCTCCCTAGAAGGTCCGGAGCTGCCGGTGCTGCAGTCAAAGGAATCGGAACAATGGTGATGGAACTGATCAGAATGATCAAACGGGGGATCAACGATCGAAATTTTTGGAGAGGTGAGAATGGGCGGAAAACAAGAAGTGCTTATGAGAGAATGTGCAACATTCTCAAAGGAAAATTTCAAACAGCTGCACAAAAAGCAATGGTGGATCAAGTTAGAGAAAGCCGGAATCCAGGAAACGCTGAGATCGAAGATCTCATATTTTTAGCAAGATCTGCACTGATATTGAGAGGATCAGTTGCTCACAAATCTTGCCTACCTGCCTGTGCATATGGACCTGCAGTATCCAGTGGTTATGACTTTGAAAAAGAGGGATATTCCTTGGTGGGAATAGACCCTTTCAAACTACTTCAAAATAGCCAAATATACAGCTTAATCAGACCTAATGAGAATCCAGCACACAAGAGTCAGCTGGTGTGGATGGCATGTCATTCTGCTGCATTTGAAGATTTAAGATTGTTAAGCTTCATCAGAGGAACAAAAGTATCTCCTCGGGGGAAACTGTCAACTAGAGGAGTACAAATTGCTTCAAATGAGAACATGGATAATATGGGATCAAGCACTCTTGAACTGAGAAGCGGGTACTGGGCCATAAGGACCAGGAGTGGAGGAAACACTAATCAGCAGAGGGCCTCCGCAGGCCAAACCAGTGTGCAACCAACGTTTTCTGTACAAAGAAACCTCCCATTTGAAAAGTCAACCATCATGGCAGCATTCACTGGAAATACGGAAGGAAGAACTTCAGACATGAGGGCAGAAATTATAAGGATGATGGAAGGTGCAAAACCAGAAGAAGTGTCATTCCGGGGGAGGGGAGTTTTCGAGCTCTCTGACGAGAAGGCAGCGAACCCGATCGTGCCCTCTTTTGATATGAGTAACGAAGGATCTTATTTCTTCGGAGACAATGCAGAAGAATACGACAATTAAGAAAAAANNNN";
        let mut mins = Vec::new();
        sketch_sequence(seq.as_bytes(), seq.len(), 10, 15, 0, &mut mins);
        assert!(!mins.is_empty());
    }

    /// The chunked path must give exactly the one-pass sketcher's minimizers.
    #[test]
    fn chunked_sketch_matches_one_pass() {
        let mut x = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = move || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for case in 0..3000 {
            let len = (rnd() % 3000) as usize;
            let mut seq: Vec<u8> = (0..len).map(|_| b"ACGT"[(rnd() & 3) as usize]).collect();
            // N runs, lowercase, low-complexity stretches and palindromes (symmetric k-mers)
            if case % 3 == 0 && len > 10 { let s = (rnd() as usize) % len; for c in seq[s..(s + 1 + (rnd() % 40) as usize).min(len)].iter_mut() { *c = b'N'; } }
            if case % 5 == 0 && len > 10 { let s = (rnd() as usize) % len; for c in seq[s..(s + 60).min(len)].iter_mut() { *c = b'A'; } }
            if case % 7 == 0 && len > 40 { for j in 0..20 { seq[j] = b"ACGT"[j % 4]; seq[39 - j] = match seq[j] { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' }; } }
            if case % 11 == 0 { for c in seq.iter_mut().step_by(3) { *c = c.to_ascii_lowercase(); } }
            let k = [15usize, 17, 19, 21, 28, 5][(rnd() % 6) as usize];
            let w = [1usize, 2, 10, 19, 50, 255][(rnd() % 6) as usize];
            sketch_sequence(&seq, len, w, k, (case % 3) as usize, &mut a);
            sketch_sequence_reference(&seq, len, w, k, (case % 3) as usize, &mut b);
            assert_eq!(a, b, "case {case}: len {len} k {k} w {w}");
        }
    }

    #[test]
    fn test_sketch_empty() {
        let seq = "";
        let mut mins = Vec::new();
        sketch_sequence(seq.as_bytes(), seq.len(), 10, 15, 0, &mut mins);
        assert_eq!(mins.len(), 0);
    }

    /// The hash wraps by design. These k-mers make an intermediate exceed
    /// u64, which panics in overflow-checked builds unless every step wraps.
    #[test]
    fn test_kmer_hash_wraps_instead_of_overflowing() {
        assert_eq!(kmer_hash(5210431, (1 << 42) - 1), 401579441571); // k=21
        assert_eq!(kmer_hash(68568, (1 << 56) - 1), 636782588647065); // k=28
        // Every k-mer for the largest supported k must be hashable.
        for k in [15usize, 17, 19, 21, 25, 28] {
            let mask = (1u64 << (2 * k)) - 1;
            for b in 0..200_000u64 {
                let _ = kmer_hash(b, mask);
            }
        }
    }

}

#[cfg(target_arch = "x86_64")]
#[path = "sketch512.rs"]
mod sketch512;
#[cfg(target_arch = "x86_64")]
static SKETCH_AVX512_OFF:std::sync::LazyLock<bool>=std::sync::LazyLock::new(||
    std::env::var("BONEFISH_SKETCH_AVX512_OFF").is_ok_and(|v|v!="0"));
