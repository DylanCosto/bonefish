//! Per-bucket hash table index backend.
//!
//! Each of 2^B buckets has its own open-addressing hash table that maps
//! hash suffixes to (offset, count) ranges in a shared flat `positions` array.
//! During build, each bucket's temporary (hash, position) pairs are processed
//! and freed independently — the key to controlling peak build memory.
//!
//! All positions (including singletons) live in `positions[]` so that the
//! existing `get_range`/`get_by_range` API works without changes.

use serde::{Serialize, Deserialize, Serializer, Deserializer};
use serde::de::SeqAccess;
use super::index::SeedLookup;
use super::index_storage::Storage;

/// Multiplier that maps a hash suffix to its sort/probe key. Odd, so the map is a bijection
/// on u64 (the key is recoverable with [`MIX_INV`]); the product's high bits are well mixed.
const MIX: u64 = 0x9E37_79B9_7F4A_7C15;
/// Multiplicative inverse of [`MIX`] mod 2^64.
const MIX_INV: u64 = {
    let mut x = MIX; // Newton iteration: each step doubles the number of correct bits.
    let mut i = 0;
    while i < 6 { x = x.wrapping_mul(2u64.wrapping_sub(MIX.wrapping_mul(x))); i += 1; }
    x
};
/// Empty-slot marker; compares greater than every stored probe key.
const EMPTY: u64 = u64::MAX;

/// One hash-table slot: `m` = hash suffix × [`MIX`] (or [`EMPTY`]), `v` = (offset << 32) | count.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Slot { pub m: u64, pub v: u64 }

/// A bucket's table region within `slots`: home slots are `base .. base + cap`; entries that
/// overflow past the last home slot follow, then one [`EMPTY`] sentinel ends the region.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BucketHdr { pub base: u64, pub cap: u64 }

/// Home slot of probe key `m` in a region with `cap` home slots.
#[inline(always)]
fn home(m: u64, cap: u64) -> u64 {
    ((m as u128 * cap as u128) >> 64) as u64
}

/// Lay out one bucket's entries (`(probe key, value)`, consumed) as a region: sorted by probe
/// key and placed at `max(home, previous + 1)`. Every run of occupied slots is therefore in
/// ascending key order, so a lookup scans from its home slot until it meets a key >= its own
/// (a hit, another key's run, or an empty slot) — a miss stops as early as a hit.
fn layout_region(entries: Vec<(u64, u64)>) -> (u64, Vec<Slot>) {
    if entries.is_empty() {
        return (0, vec![Slot { m: EMPTY, v: 0 }]);
    }
    let n = entries.len() as u64;
    let cap = n + n / 3 + 1; // load factor <= 0.75
    // Ordered linear probing: insert each entry at the first slot of its probe sequence whose
    // key is larger, shifting that run one slot right. The final table doesn't depend on the
    // insertion order: it equals placing the entries in key order at max(home, previous + 1).
    let mut region = vec![Slot { m: EMPTY, v: 0 }; cap as usize + 64];
    for (m, v) in entries {
        assert!(m != EMPTY, "hash key collides with the empty-slot marker");
        let mut p = home(m, cap) as usize;
        while region[p].m < m { p += 1; }
        debug_assert!(region[p].m != m, "duplicate hash key");
        let mut q = p;
        while region[q].m != EMPTY { q += 1; }
        if q + 1 >= region.len() { region.resize(q + 64, Slot { m: EMPTY, v: 0 }); }
        region.copy_within(p..q, p + 1);
        region[p] = Slot { m, v };
    }
    // Trim to the last occupied slot (at least the home range), then add the sentinel.
    let used = region.iter().rposition(|s| s.m != EMPTY).map_or(0, |i| i + 1).max(cap as usize);
    region.truncate(used);
    region.push(Slot { m: EMPTY, v: 0 }); // sentinel
    (cap, region)
}

/// Appends bucket regions in order. On native targets the regions are laid out by two worker
/// threads, overlapping with whatever produces the entries (deserialization or index build).
#[derive(Default)]
struct TableBuilder {
    hdrs: Vec<BucketHdr>,
    slots: Vec<Slot>,
    #[cfg(not(target_arch = "wasm32"))]
    workers: Option<LayoutWorkers>,
}

#[cfg(not(target_arch = "wasm32"))]
struct LayoutWorkers {
    tx: Vec<std::sync::mpsc::SyncSender<Vec<(u64, u64)>>>,
    rx: Vec<std::sync::mpsc::Receiver<(u64, Vec<Slot>)>>,
    sent: usize,
    received: usize,
}

#[cfg(not(target_arch = "wasm32"))]
impl LayoutWorkers {
    const N: usize = 2;
    fn new() -> Self {
        let (mut tx, mut rx) = (Vec::new(), Vec::new());
        for _ in 0..Self::N {
            let (in_tx, in_rx) = std::sync::mpsc::sync_channel::<Vec<(u64, u64)>>(2);
            let (out_tx, out_rx) = std::sync::mpsc::sync_channel(2);
            std::thread::spawn(move || {
                for entries in in_rx {
                    if out_tx.send(layout_region(entries)).is_err() { break; }
                }
            });
            tx.push(in_tx);
            rx.push(out_rx);
        }
        Self { tx, rx, sent: 0, received: 0 }
    }
}

impl TableBuilder {
    /// Queue one bucket's `(hash suffix, value)` entries.
    fn push(&mut self, mut entries: Vec<(u64, u64)>) {
        for e in entries.iter_mut() { e.0 = e.0.wrapping_mul(MIX); }
        #[cfg(not(target_arch = "wasm32"))]
        {
            let w = self.workers.get_or_insert_with(LayoutWorkers::new);
            w.tx[w.sent % LayoutWorkers::N].send(entries).expect("layout worker exited");
            w.sent += 1;
            // Keep a bounded number of buckets in flight; collect results in bucket order.
            while w.sent - w.received > 2 * LayoutWorkers::N {
                let r = w.rx[w.received % LayoutWorkers::N].recv().expect("layout worker exited");
                w.received += 1;
                Self::append(&mut self.hdrs, &mut self.slots, r);
            }
        }
        #[cfg(target_arch = "wasm32")]
        Self::append(&mut self.hdrs, &mut self.slots, layout_region(entries));
    }

    fn append(hdrs: &mut Vec<BucketHdr>, slots: &mut Vec<Slot>, (cap, region): (u64, Vec<Slot>)) {
        hdrs.push(BucketHdr { base: slots.len() as u64, cap });
        slots.extend_from_slice(&region);
    }

    fn finish(mut self) -> (Vec<BucketHdr>, Vec<Slot>) {
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(mut w) = self.workers.take() {
            while w.received < w.sent {
                let r = w.rx[w.received % LayoutWorkers::N].recv().expect("layout worker exited");
                w.received += 1;
                Self::append(&mut self.hdrs, &mut self.slots, r);
            }
        }
        self.slots.shrink_to_fit();
        (self.hdrs, self.slots)
    }
}

/// Minimap2-style per-bucket hash index.
///
/// Each of the 2^`bucket_bits` buckets owns a static open-addressing table region in the flat
/// `slots` array (see [`layout_region`]); all positions live in `positions`. All three arrays
/// are plain `#[repr(C)]` data, so a raw index file can be mapped and used in place.
#[derive(Debug, Clone)]
pub struct BucketHashLookup {
    bucket_bits: u32,
    hdrs: Storage<BucketHdr>,
    slots: Storage<Slot>,
    /// Flat position array for ALL hashes (singletons and multi-occurrence).
    positions: Storage<u64>,
    /// Optional precomputed histogram of occurrence counts (`(count, n_hashes)`, ascending),
    /// stored in raw index files so `cal_mid_occ` need not scan every slot.
    count_hist: Option<Vec<(u32, u64)>>,
}

impl BucketHashLookup {
    /// Create an empty lookup (for placeholder/header-only indices).
    pub fn empty() -> Self {
        Self { bucket_bits: 0, hdrs: Storage::default(), slots: Storage::default(), positions: Storage::default(), count_hist: None }
    }

    fn from_builder(bucket_bits: u32, tb: TableBuilder, positions: Vec<u64>) -> Self {
        let (hdrs, slots) = tb.finish();
        Self { bucket_bits, hdrs: hdrs.into(), slots: slots.into(), positions: positions.into(), count_hist: None }
    }

    /// Load from the `.mmi` per-bucket hash table format.
    ///
    /// Each bucket contributes: `p[]` (multi-occurrence positions) + hash entries
    /// (key/value pairs where key bit 0 = singleton flag). We convert to our
    /// own bucket hash format with a unified positions array.
    pub fn from_minimap2_buckets(
        bucket_bits: u32,
        bucket_data: Vec<(Vec<u64>, Vec<(u64, u64)>)>, // per-bucket: (p_array, hash_entries)
    ) -> Self {
        let mut positions: Vec<u64> = Vec::new();
        let mut tb = TableBuilder::default();

        for (p, hash_entries) in &bucket_data {
            let mut entries = Vec::with_capacity(hash_entries.len());
            for &(key, value) in hash_entries {
                // Source-format key: (hash_suffix << 1) | singleton_flag.
                let our_key = key >> 1;

                if key & 1 != 0 {
                    // Singleton: value is the position directly
                    let offset = positions.len() as u64;
                    positions.push(value);
                    entries.push((our_key, (offset << 32) | 1));
                } else {
                    // Multi-occurrence: value = (offset_in_p << 32) | count
                    let count = (value & 0xFFFFFFFF) as usize;
                    let start = (value >> 32) as usize;
                    let offset = positions.len() as u64;
                    positions.extend_from_slice(&p[start..start + count]);
                    entries.push((our_key, (offset << 32) | count as u64));
                }
            }
            tb.push(entries);
        }

        Self::from_builder(bucket_bits, tb, positions)
    }

    /// Build from pre-distributed, pre-sorted bucket data (consumed bucket by bucket so the
    /// source entries are freed as the tables are laid out).
    pub fn build(bucket_bits: u32, sorted_buckets: &mut [Vec<(u64, u64)>], max_occ: usize) -> Self {
        let mut total_positions = 0usize;
        for b in sorted_buckets.iter() {
            let mut i = 0;
            while i < b.len() {
                let h = b[i].0;
                let start = i;
                while i < b.len() && b[i].0 == h { i += 1; }
                if max_occ < usize::MAX && (i - start) > max_occ { continue; }
                total_positions += i - start;
            }
        }

        let mut positions: Vec<u64> = Vec::with_capacity(total_positions);
        let mut tb = TableBuilder::default();
        for bucket_entries in sorted_buckets.iter_mut() {
            let b = std::mem::take(bucket_entries);
            let mut entries = Vec::new();
            let mut i = 0;
            while i < b.len() {
                let h = b[i].0;
                let start = i;
                while i < b.len() && b[i].0 == h { i += 1; }
                let count = i - start;
                if max_occ < usize::MAX && count > max_occ { continue; }
                let offset = positions.len() as u64;
                positions.extend(b[start..i].iter().map(|e| e.1));
                entries.push((h >> bucket_bits, (offset << 32) | count as u64));
            }
            tb.push(entries);
        }

        Self::from_builder(bucket_bits, tb, positions)
    }

    /// Raw arrays, for writing a raw (mappable) index file.
    pub fn raw_parts(&self) -> (u32, &[BucketHdr], &[Slot], &[u64]) {
        (self.bucket_bits, &self.hdrs, &self.slots, &self.positions)
    }

    /// Reassemble from raw arrays (e.g. mapped from a raw index file), checking the table
    /// invariants that lookups rely on.
    pub fn from_raw_parts(bucket_bits: u32, hdrs: Storage<BucketHdr>, slots: Storage<Slot>, positions: Storage<u64>,
                          count_hist: Option<Vec<(u32, u64)>>) -> Result<Self, String> {
        if bucket_bits >= 32 || hdrs.len() != 1usize << bucket_bits {
            return Err(format!("bad bucket table: {} headers for {} bucket bits", hdrs.len(), bucket_bits));
        }
        for (b, h) in hdrs.iter().enumerate() {
            let end = hdrs.get(b + 1).map_or(slots.len() as u64, |n| n.base);
            // The home slots lie inside the region and the region's last slot is an empty
            // sentinel, so every probe scan terminates inside the region.
            if h.base.checked_add(h.cap).is_none_or(|e| e >= end) || end > slots.len() as u64
                || slots[end as usize - 1].m != EMPTY {
                return Err(format!("corrupt hash table region for bucket {}", b));
            }
        }
        Ok(Self { bucket_bits, hdrs, slots, positions, count_hist })
    }

    /// Histogram of occurrence counts, `(count, n_hashes)` in ascending count order.
    pub fn count_histogram(&self) -> Vec<(u32, u64)> {
        if let Some(h) = &self.count_hist { return h.clone(); }
        let mut counts: Vec<u32> = self.slots.iter().filter(|s| s.m != EMPTY).map(|s| s.v as u32).collect();
        counts.sort_unstable();
        let mut hist: Vec<(u32, u64)> = Vec::new();
        for c in counts {
            match hist.last_mut() {
                Some(last) if last.0 == c => last.1 += 1,
                _ => hist.push((c, 1)),
            }
        }
        hist
    }

    #[inline(always)]
    fn find(&self, hash: u64) -> Option<u64> {
        let hdr = self.hdrs[(hash & ((1u64 << self.bucket_bits) - 1)) as usize];
        let m = (hash >> self.bucket_bits).wrapping_mul(MIX);
        let mut i = (hdr.base + home(m, hdr.cap)) as usize;
        loop {
            // SAFETY: the scan starts at a home slot of this bucket's region and stops at the
            // latest at the region's EMPTY sentinel (checked at construction / load).
            let s = unsafe { *self.slots.get_unchecked(i) };
            if s.m >= m {
                return if s.m == m { Some(s.v) } else { None };
            }
            i += 1;
        }
    }
}

// Serialized in the original bincode index layout: bucket bits, then per bucket a sequence of
// (hash suffix, value) pairs, then the positions array.
impl Serialize for BucketHashLookup {
    fn serialize<S: Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
        use serde::ser::{SerializeSeq, SerializeStruct};
        struct BucketSer<'a>(&'a [Slot]);
        impl Serialize for BucketSer<'_> {
            fn serialize<S: Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
                let n = self.0.iter().filter(|s| s.m != EMPTY).count();
                let mut seq = ser.serialize_seq(Some(n))?;
                for s in self.0.iter().filter(|s| s.m != EMPTY) {
                    seq.serialize_element(&(s.m.wrapping_mul(MIX_INV), s.v))?;
                }
                seq.end()
            }
        }
        struct BucketsSer<'a>(&'a BucketHashLookup);
        impl Serialize for BucketsSer<'_> {
            fn serialize<S: Serializer>(&self, ser: S) -> Result<S::Ok, S::Error> {
                let l = self.0;
                let mut seq = ser.serialize_seq(Some(l.hdrs.len()))?;
                for (b, h) in l.hdrs.iter().enumerate() {
                    let end = l.hdrs.get(b + 1).map_or(l.slots.len(), |n| n.base as usize);
                    seq.serialize_element(&BucketSer(&l.slots[h.base as usize..end]))?;
                }
                seq.end()
            }
        }
        let mut st = ser.serialize_struct("BucketHashLookup", 3)?;
        st.serialize_field("bucket_bits", &self.bucket_bits)?;
        st.serialize_field("buckets", &BucketsSer(self))?;
        st.serialize_field("positions", &self.positions)?;
        st.end()
    }
}

impl<'de> Deserialize<'de> for BucketHashLookup {
    fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
        /// Streams the per-bucket entry sequences straight into a table builder.
        struct Tables(TableBuilder);
        impl<'de> Deserialize<'de> for Tables {
            fn deserialize<D: Deserializer<'de>>(de: D) -> Result<Self, D::Error> {
                struct V;
                impl<'de> serde::de::Visitor<'de> for V {
                    type Value = Tables;
                    fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                        f.write_str("seq of buckets")
                    }
                    fn visit_seq<A: SeqAccess<'de>>(self, mut s: A) -> Result<Tables, A::Error> {
                        let mut tb = TableBuilder::default();
                        while let Some(entries) = s.next_element::<Vec<(u64, u64)>>()? {
                            tb.push(entries);
                        }
                        Ok(Tables(tb))
                    }
                }
                de.deserialize_seq(V)
            }
        }
        #[derive(Deserialize)]
        struct Raw { bucket_bits: u32, buckets: Tables, positions: Vec<u64> }
        let r = Raw::deserialize(de)?;
        Ok(BucketHashLookup::from_builder(r.bucket_bits, r.buckets.0, r.positions))
    }
}

impl BucketHashLookup {
    /// Prefetch the home slot for a hash value into L1 cache.
    /// Call this 5-10 iterations ahead of the actual lookup.
    #[inline]
    pub fn prefetch(&self, hash: u64) {
        if self.hdrs.is_empty() { return; }
        // aarch64's prefetch intrinsic (`stdarch_aarch64_prefetch`) is still
        // nightly-only as of Rust 1.96, so we skip the hint there to stay on
        // stable. It costs nothing on Apple M2 and 1-2% on a Raspberry Pi.
        #[cfg(target_arch = "x86_64")]
        unsafe {
            let hdr = self.hdrs[(hash & ((1u64 << self.bucket_bits) - 1)) as usize];
            let m = (hash >> self.bucket_bits).wrapping_mul(MIX);
            let ptr = self.slots.as_ptr().add((hdr.base + home(m, hdr.cap)) as usize) as *const i8;
            std::arch::x86_64::_mm_prefetch(ptr, std::arch::x86_64::_MM_HINT_T0);
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = hash;
    }

    /// Prefetch a positions range into L1 cache.
    #[inline]
    pub fn prefetch_positions(&self, range: (u32, u32)) {
        if (range.0 as usize) < self.positions.len() {
            unsafe {
                let ptr = self.positions.as_ptr().add(range.0 as usize) as *const u8;
                #[cfg(target_arch = "x86_64")]
                std::arch::x86_64::_mm_prefetch(ptr as *const i8, std::arch::x86_64::_MM_HINT_T0);
                #[cfg(not(target_arch = "x86_64"))]
                { let _ = ptr; }
            }
        }
    }
}

impl SeedLookup for BucketHashLookup {
    #[inline]
    fn get(&self, hash: u64) -> Option<&[u64]> {
        let encoded = self.find(hash)?;
        let offset = (encoded >> 32) as usize;
        let count = (encoded & 0xFFFFFFFF) as usize;
        Some(&self.positions[offset..offset + count])
    }

    #[inline]
    fn get_range(&self, hash: u64) -> Option<(u32, u32)> {
        let encoded = self.find(hash)?;
        let offset = (encoded >> 32) as u32;
        let count = (encoded & 0xFFFFFFFF) as u32;
        Some((offset, offset + count))
    }

    #[inline]
    fn get_by_range(&self, range: (u32, u32)) -> &[u64] {
        &self.positions[range.0 as usize..range.1 as usize]
    }

    fn occurrence_counts(&self) -> Box<dyn Iterator<Item = u32> + '_> {
        Box::new(self.slots.iter().filter(|s| s.m != EMPTY).map(|s| (s.v & 0xFFFFFFFF) as u32))
    }

    fn is_empty(&self) -> bool {
        self.slots.iter().all(|s| s.m == EMPTY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bucket_hash_basic() {
        let bucket_bits = 2u32; // 4 buckets
        let mask = (1u64 << bucket_bits) - 1;

        // Create test entries: 3 hashes, each in a different bucket
        let mut buckets = vec![Vec::new(); 4];
        let entries = vec![
            (100u64, 1000u64), // bucket 100 & 3 = 0
            (101, 2000),       // bucket 101 & 3 = 1
            (101, 2001),       // same hash, bucket 1
            (102, 3000),       // bucket 102 & 3 = 2
        ];
        for &(h, p) in &entries {
            buckets[(h & mask) as usize].push((h, p));
        }
        for b in &mut buckets { b.sort_unstable(); }

        let lut = BucketHashLookup::build(bucket_bits, &mut buckets, usize::MAX);

        // Singleton lookup
        assert_eq!(lut.get(100), Some(&[1000u64][..]));
        // Multi-occurrence lookup
        let r = lut.get(101).unwrap();
        assert_eq!(r.len(), 2);
        assert!(r.contains(&2000));
        assert!(r.contains(&2001));
        // Another singleton
        assert_eq!(lut.get(102), Some(&[3000u64][..]));
        // Missing hash
        assert_eq!(lut.get(999), None);

        // get_range / get_by_range
        let (s, e) = lut.get_range(101).unwrap();
        assert_eq!((e - s) as usize, 2);
        let slice = lut.get_by_range((s, e));
        assert_eq!(slice.len(), 2);
    }
}

#[cfg(test)]
mod bench_layout {
    #[test]
    #[ignore]
    fn bench_layout_region() {
        let mut x = 88172645463325252u64;
        let mut v: Vec<(u64, u64)> = (0..210_000).map(|i| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; ((x >> 40).wrapping_mul(super::MIX), i) }).collect();
        v.sort_unstable(); v.dedup_by_key(|e| e.0);
        let t = std::time::Instant::now();
        for _ in 0..20 { std::hint::black_box(super::layout_region(v.clone())); }
        eprintln!("layout: {:.1} ns/entry", t.elapsed().as_nanos() as f64 / (20 * v.len()) as f64);
    }
}
