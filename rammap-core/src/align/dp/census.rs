//! DP census (Gate B): records every extension/gap-fill DP call for offline analysis.
//!
//! Compiled only with `--features census`; enabled at run time by `RAMMAP_CENSUS=<prefix>`.
//! Writes `<prefix>.calls.bin` (one fixed 68-byte record per DP call) and
//! `<prefix>.seqs.bin` (full inputs for every `RAMMAP_CENSUS_SEQ_EVERY`-th call, default
//! 200, so kernels can be replayed exactly in the Gate C bake-off).
//!
//! Records are buffered per thread and appended to a shared registry; `flush()` must be
//! called once at the end of the run (rayon pool threads never run TLS destructors).
//!
//! calls.bin record (little-endian, 68 bytes):
//!   u8 kind, u8 flags_hi (bit0 zdropped, bit1 reach_end, bit2 seq dumped, bit3 dual-affine,
//!   bit4 skipped by the max_sw_mat guard),
//!   u16 pad, i32 dp_flags, i32 bandwidth, i32 zdrop, i32 end_bonus,
//!   u32 qlen, u32 tlen, i32 score, i32 max, u32 q_used, u32 t_used,
//!   u32 n_mismatch, u32 n_ins_ops, u32 ins_bases, u32 n_del_ops, u32 del_bases, u32 time_ns
//!
//! seqs.bin record: u32 call_index_in_thread_stream, u64 global_seq_id, the same 68-byte
//!   call record, i8[25] score matrix, i8 q, i8 e, i8 q2, i8 e2, u32 qlen, u32 tlen,
//!   qseq (nt4 bytes), tseq (nt4 bytes).

use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;

use super::common::DpResult;

/// Which part of the alignment a DP call belongs to.
#[derive(Clone, Copy)]
#[repr(u8)]
pub enum Kind {
    LeftExt = 0,
    Gap = 1,
    /// Gap-fill re-run without APPROX_MAX after a z-drop test failed.
    GapRerun = 2,
    RightExt = 3,
}

const RECORD_BYTES: usize = 68;
const FLUSH_BYTES: usize = 4 << 20;

struct Config {
    seq_every: u64,
    calls: Mutex<File>,
    seqs: Mutex<File>,
}

static CONFIG: OnceLock<Option<Config>> = OnceLock::new();
static SEQ_COUNTER: AtomicU64 = AtomicU64::new(0);
static REGISTRY: Mutex<Vec<Arc<Mutex<Buffers>>>> = Mutex::new(Vec::new());

#[derive(Default)]
struct Buffers {
    calls: Vec<u8>,
    seqs: Vec<u8>,
    n: u32,
}

thread_local! {
    static LOCAL: Arc<Mutex<Buffers>> = {
        let b = Arc::new(Mutex::new(Buffers::default()));
        REGISTRY.lock().unwrap().push(b.clone());
        b
    };
}

fn config() -> Option<&'static Config> {
    CONFIG
        .get_or_init(|| {
            let prefix = std::env::var("RAMMAP_CENSUS").ok()?;
            let seq_every = std::env::var("RAMMAP_CENSUS_SEQ_EVERY").ok().and_then(|s| s.parse().ok()).unwrap_or(200);
            let open = |ext: &str| File::create(format!("{prefix}.{ext}")).unwrap_or_else(|e| panic!("census: {prefix}.{ext}: {e}"));
            log::info!("DP census enabled: {prefix}.calls.bin / .seqs.bin (seq every {seq_every})");
            Some(Config { seq_every, calls: Mutex::new(open("calls.bin")), seqs: Mutex::new(open("seqs.bin")) })
        })
        .as_ref()
}

/// Gap-penalty context needed to replay a call exactly.
pub struct Params<'a> {
    pub mat: &'a [i8],
    pub q: i8,
    pub e: i8,
    pub q2: i8,
    pub e2: i8,
    pub bandwidth: i32,
    pub zdrop: i32,
    pub end_bonus: i32,
    pub flags: i32,
}

/// Start timing a DP call. Returns None (and costs one atomic load) when disabled.
#[inline]
pub fn start() -> Option<Instant> {
    config().map(|_| Instant::now())
}

/// Record a finished DP call started with `start()`.
pub fn finish(t0: Option<Instant>, kind: Kind, qseq: &[u8], tseq: &[u8], p: &Params, ez: &DpResult) {
    let (Some(t0), Some(cfg)) = (t0, config()) else { return };
    let ns = t0.elapsed().as_nanos().min(u32::MAX as u128) as u32;

    let (q_used, t_used, mism, ni, bi, nd, bd) = cigar_stats(&ez.cigar, qseq, tseq);
    let dump = SEQ_COUNTER.fetch_add(1, Ordering::Relaxed) % cfg.seq_every == 0;
    let dual = p.q != p.q2 || p.e != p.e2;
    // skipped: the max_sw_mat guard fired and no DP ran.
    let skipped = ez.cigar.is_empty() && ez.zdropped != 0 && ez.max_score_query_pos == -1;
    let hi = (ez.zdropped != 0) as u8 | ((ez.reach_end != 0) as u8) << 1 | (dump as u8) << 2 | (dual as u8) << 3 | (skipped as u8) << 4;

    let mut rec = [0u8; RECORD_BYTES];
    let mut w = Cursor(&mut rec, 0);
    w.u8(kind as u8);
    w.u8(hi);
    w.u16(0);
    for v in [p.flags, p.bandwidth, p.zdrop, p.end_bonus] { w.i32(v); }
    w.u32(qseq.len() as u32);
    w.u32(tseq.len() as u32);
    w.i32(ez.score);
    w.i32(ez.max);
    for v in [q_used, t_used, mism, ni, bi, nd, bd, ns] { w.u32(v); }

    LOCAL.with(|b| {
        let mut b = b.lock().unwrap();
        let idx = b.n;
        b.n += 1;
        b.calls.extend_from_slice(&rec);
        if dump {
            let id = SEQ_COUNTER.load(Ordering::Relaxed);
            b.seqs.extend_from_slice(&idx.to_le_bytes());
            b.seqs.extend_from_slice(&id.to_le_bytes());
            b.seqs.extend_from_slice(&rec);
            let mut m = [0i8; 25];
            m[..p.mat.len().min(25)].copy_from_slice(&p.mat[..p.mat.len().min(25)]);
            b.seqs.extend(m.iter().map(|&x| x as u8));
            b.seqs.extend([p.q as u8, p.e as u8, p.q2 as u8, p.e2 as u8]);
            b.seqs.extend_from_slice(&(qseq.len() as u32).to_le_bytes());
            b.seqs.extend_from_slice(&(tseq.len() as u32).to_le_bytes());
            b.seqs.extend_from_slice(qseq);
            b.seqs.extend_from_slice(tseq);
        }
        if b.calls.len() >= FLUSH_BYTES || b.seqs.len() >= FLUSH_BYTES {
            drain(cfg, &mut b);
        }
    });
}

fn drain(cfg: &Config, b: &mut Buffers) {
    if !b.calls.is_empty() {
        cfg.calls.lock().unwrap().write_all(&b.calls).expect("census: write calls");
        b.calls.clear();
    }
    if !b.seqs.is_empty() {
        cfg.seqs.lock().unwrap().write_all(&b.seqs).expect("census: write seqs");
        b.seqs.clear();
    }
}

/// Write out every thread's remaining records. Call once after mapping finishes.
pub fn flush() {
    let Some(cfg) = config() else { return };
    for b in REGISTRY.lock().unwrap().iter() {
        drain(cfg, &mut b.lock().unwrap());
    }
    cfg.calls.lock().unwrap().flush().ok();
    cfg.seqs.lock().unwrap().flush().ok();
}

/// Walk a DP CIGAR ((len << 4) | op, op 0=M 1=I 2=D) against the inputs.
/// Returns (query bases used, target bases used, mismatches, ins ops, ins bases, del ops, del bases).
fn cigar_stats(cigar: &[u32], q: &[u8], t: &[u8]) -> (u32, u32, u32, u32, u32, u32, u32) {
    let (mut qi, mut ti) = (0usize, 0usize);
    let (mut mism, mut ni, mut bi, mut nd, mut bd) = (0u32, 0u32, 0u32, 0u32, 0u32);
    for &c in cigar {
        let len = (c >> 4) as usize;
        match c & 0xf {
            0 => {
                let end = (qi + len).min(q.len()).min(qi + t.len().saturating_sub(ti));
                mism += q[qi..end].iter().zip(&t[ti..ti + (end - qi)]).filter(|(a, b)| a != b).count() as u32;
                qi += len;
                ti += len;
            }
            1 => { ni += 1; bi += len as u32; qi += len; }
            2 => { nd += 1; bd += len as u32; ti += len; }
            _ => {}
        }
    }
    (qi as u32, ti as u32, mism, ni, bi, nd, bd)
}

struct Cursor<'a>(&'a mut [u8; RECORD_BYTES], usize);
impl Cursor<'_> {
    fn put(&mut self, b: &[u8]) { self.0[self.1..self.1 + b.len()].copy_from_slice(b); self.1 += b.len(); }
    fn u8(&mut self, v: u8) { self.put(&[v]); }
    fn u16(&mut self, v: u16) { self.put(&v.to_le_bytes()); }
    fn u32(&mut self, v: u32) { self.put(&v.to_le_bytes()); }
    fn i32(&mut self, v: i32) { self.put(&v.to_le_bytes()); }
}
