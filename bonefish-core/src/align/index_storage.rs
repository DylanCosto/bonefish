//! Array storage for the large index arrays: an owned `Vec`, or (on native targets) a region
//! of a memory-mapped raw index file, so a pre-built index can be used straight from the page
//! cache instead of being copied into freshly allocated memory.
//!
//! Serializes exactly like `Vec<T>` (length, then elements), so the bincode index format is
//! unchanged.

use std::fmt;
use std::marker::PhantomData;
use std::ops::Deref;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

pub enum Storage<T: Copy + 'static> {
    Owned(Vec<T>),
    #[cfg(not(target_arch = "wasm32"))]
    Mapped {
        map: std::sync::Arc<memmap2::Mmap>,
        /// Byte offset of the first element (aligned for `T`).
        off: usize,
        len: usize,
        _t: PhantomData<T>,
    },
}

impl<T: Copy + 'static> Storage<T> {
    /// A `len`-element view of `map` starting at byte `off` (`off` must be aligned for `T`
    /// and the range must lie inside the map).
    #[cfg(not(target_arch = "wasm32"))]
    pub fn mapped(map: std::sync::Arc<memmap2::Mmap>, off: usize, len: usize) -> Self {
        assert!(off.is_multiple_of(std::mem::align_of::<T>()), "misaligned mapped array");
        assert!(off + len * std::mem::size_of::<T>() <= map.len(), "mapped array out of range");
        Storage::Mapped { map, off, len, _t: PhantomData }
    }
}

impl<T: Copy + 'static> Storage<T> {
    /// A copy in anonymous memory advised for transparent huge pages (fewer TLB misses on
    /// random lookups). Linux only; elsewhere a plain copy.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn to_huge_pages(&self) -> std::io::Result<Self> {
        let bytes = std::mem::size_of_val(&self[..]);
        if bytes == 0 { return Ok(Storage::Owned(Vec::new())); }
        let mut m = memmap2::MmapMut::map_anon(bytes)?;
        #[cfg(target_os = "linux")]
        let _ = m.advise(memmap2::Advice::HugePage);
        // SAFETY: plain-old-data copy into a fresh mapping of exactly `bytes` bytes.
        unsafe { std::ptr::copy_nonoverlapping(self.as_ptr() as *const u8, m.as_mut_ptr(), bytes); }
        let map = std::sync::Arc::new(m.make_read_only()?);
        Ok(Storage::Mapped { map, off: 0, len: self.len(), _t: PhantomData })
    }
}

impl<T: Copy + 'static> Deref for Storage<T> {
    type Target = [T];
    #[inline]
    fn deref(&self) -> &[T] {
        match self {
            Storage::Owned(v) => v,
            #[cfg(not(target_arch = "wasm32"))]
            Storage::Mapped { map, off, len, .. } => unsafe {
                std::slice::from_raw_parts(map.as_ptr().add(*off) as *const T, *len)
            },
        }
    }
}

impl<T: Copy + 'static> Default for Storage<T> {
    fn default() -> Self {
        Storage::Owned(Vec::new())
    }
}

impl<T: Copy + 'static> From<Vec<T>> for Storage<T> {
    fn from(v: Vec<T>) -> Self {
        Storage::Owned(v)
    }
}

impl<T: Copy + 'static> Clone for Storage<T> {
    fn clone(&self) -> Self {
        match self {
            Storage::Owned(v) => Storage::Owned(v.clone()),
            #[cfg(not(target_arch = "wasm32"))]
            Storage::Mapped { map, off, len, .. } => Storage::Mapped { map: map.clone(), off: *off, len: *len, _t: PhantomData },
        }
    }
}

impl<T: Copy + 'static> fmt::Debug for Storage<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Storage({} elements)", self.len())
    }
}

impl<T: Copy + Serialize + 'static> Serialize for Storage<T> {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.collect_seq(self.iter())
    }
}

impl<'de, T: Copy + Deserialize<'de> + 'static> Deserialize<'de> for Storage<T> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Vec::<T>::deserialize(d).map(Storage::Owned)
    }
}
