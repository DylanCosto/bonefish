//! GPU offload of large dual-affine extension / gap-fill DP calls (`--features gpu`).
//!
//! Enabled at run time with `RAMMAP_GPU=1`. Speculate-and-retry design, so aligner threads
//! never wait on the GPU:
//!   1. A read is aligned on the CPU as usual. When it reaches an eligible DP call, the result
//!      is looked up in a cache keyed by a 128-bit hash of the exact inputs.
//!   2. On the first miss the request is queued and the read switches to *speculative* mode:
//!      every later DP call (GPU-eligible or not) returns a placeholder at once, eligible misses
//!      are queued too, and z-drop tests are skipped. Gap/extension inputs depend only on the
//!      anchors, so one speculative pass predicts nearly all of the read's GPU requests.
//!      The read's output is discarded ([`was_speculative`]) and it is re-run next round.
//!   3. Background GPU workers drain the queue in large batches *while the round runs*.
//!   4. [`wait_pending`] waits for the queue to empty; pending reads are re-run exactly, now
//!      hitting the cache; repeat until none are pending. [`clear_cache`] frees the chunk.
//! Results are bit-identical to the AVX2 path (`cuda/gpu_dp.cu`), so output is unchanged.
//!
//! Offload only happens inside [`retry_scope`], which callers that implement the retry loop
//! wrap around each read; everywhere else (paired-end, split index) DP stays on the CPU.
//!
//! Eligible: dual-affine (q != q2 or e != e2), q+e <= q2+e2, no splice flags, scoring equal to
//! the first job's, max(qlen, tlen) >= `RAMMAP_GPU_MIN` (default 4096), and band diagonal
//! <= 16384 cells. `RAMMAP_GPU_WORKERS` (default 2) contexts/streams of `RAMMAP_GPU_MEM_GB`
//! (default 5) each share every batch.

use std::collections::{HashMap, HashSet};
use std::ffi::{c_int, c_void};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Condvar, Mutex, OnceLock};
use std::time::Duration;
use std::time::Instant;

use super::common::{DpResult, NEG_INF};

#[repr(C)]
#[derive(Clone, Copy)]
struct SegIn {
    q_off: i64,
    t_off: i64,
    state_off: i64,
    tb_off: i64,
    qlen: i32,
    tlen: i32,
    bw: i32,
    zdrop: i32,
    flags: i32,
    end_bonus: i32,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct SegOut {
    max: i32,
    max_q: i32,
    max_t: i32,
    zdropped: i32,
    score: i32,
    mqe: i32,
    mqe_t: i32,
    mte: i32,
    mte_q: i32,
    n_diag: i32,
    reach_end: i32,
    n_cigar: i32,
    cells: i64,
}

unsafe extern "C" {
    fn gpu_ctx_create(budget: i64, mat: *const i8, q: i8, e: i8, q2: i8, e2: i8) -> *mut c_void;
    fn gpu_ctx_run(ctx: *mut c_void, n: c_int, segs: *const SegIn, qbuf: *const u8, qbytes: i64, tbuf: *const u8,
                   tbytes: i64, with_tb: c_int, cig_off: *const i64, cig_total: i64, cigars: *mut u32,
                   outs: *mut SegOut, kernel_ms: *mut f64, n_batches: *mut c_int) -> c_int;
}

const MAX_DIAG: usize = 64 * 256;
const SPLICE_FLAGS: i32 = 0x100 | 0x200 | 0x400 | 0x800 | 0x1000;
const SHARDS: usize = 64;


thread_local! {
    static IN_SCOPE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static SPECULATING: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether the current read is in speculative mode (its DP results are placeholders).
#[inline]
pub fn speculating() -> bool {
    SPECULATING.with(|c| c.get())
}

/// After a read's attempt (normal return or panic), whether it went speculative and must be
/// re-run. Resets the flag.
pub fn was_speculative() -> bool {
    SPECULATING.with(|c| c.replace(false))
}

/// Placeholder result for speculative mode: an extension that extends nothing, or a gap filled
/// with matches plus one indel for the length difference. Never reaches output.
pub fn placeholder(ql: usize, tl: usize, flags: i32, ez: &mut DpResult) {
    ez.max_query_end_score = NEG_INF;
    ez.max_target_end_score = NEG_INF;
    if flags & 0x40 != 0 {   // EXTENSION_ONLY
        ez.max = 0;
        ez.max_score_query_pos = -1;
        ez.max_score_target_pos = -1;
        ez.score = NEG_INF;
        ez.cigar.clear();
    } else {
        ez.max = 0;
        ez.score = 0;
        ez.max_score_query_pos = ql as i32 - 1;
        ez.max_score_target_pos = tl as i32 - 1;
        ez.cigar.clear();
        let m = ql.min(tl) as u32;
        if m > 0 { ez.cigar.push(m << 4); }
        if ql > tl { ez.cigar.push((((ql - tl) as u32) << 4) | 1); }
        if tl > ql { ez.cigar.push((((tl - ql) as u32) << 4) | 2); }
    }
}

/// Run `f` (one read's alignment) with GPU offload allowed. Afterwards (also after a caught
/// panic) the caller must check [`was_speculative`] and, if set, discard the output and re-run
/// the read after [`wait_pending`].
pub fn retry_scope<R>(f: impl FnOnce() -> R) -> R {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) { IN_SCOPE.with(|c| c.set(false)); }
    }
    IN_SCOPE.with(|c| c.set(true));
    SPECULATING.with(|c| c.set(false));
    let _reset = Reset;
    f()
}

/// Whether the GPU path is enabled for this run (`RAMMAP_GPU` set and built with the feature).
pub fn enabled() -> bool {
    service().is_some()
}

struct Job {
    key: u128,
    q: Vec<u8>,
    t: Vec<u8>,
    bw: i32,
    zdrop: i32,
    end_bonus: i32,
    flags: i32,
}

struct Scoring {
    mat: [i8; 25],
    q: i8,
    e: i8,
    q2: i8,
    e2: i8,
}

#[derive(Default)]
struct Pending {
    jobs: Vec<Job>,
    keys: HashSet<u128>,   // queued or in flight, so a request is never run twice
    in_flight: usize,
    flush: bool,           // a round has ended: run whatever is queued now
}

struct Ctx(*mut c_void);
// SAFETY: a context is only ever used by one thread at a time (each sits behind its own Mutex).
unsafe impl Send for Ctx {}

struct Service {
    min_len: usize,
    mem_bytes: i64,
    n_workers: usize,
    scoring: OnceLock<Scoring>,
    cache: Vec<Mutex<HashMap<u128, DpResult>>>,
    pending: Mutex<Pending>,
    work_cv: Condvar,
    done_cv: Condvar,
    ctxs: Vec<Mutex<Ctx>>,
    jobs: AtomicU64,
    rounds: AtomicU64,
    hits: AtomicU64,
    cells: AtomicU64,
    gpu_us: AtomicU64,
}

static SERVICE: OnceLock<Option<&'static Service>> = OnceLock::new();

fn service() -> Option<&'static Service> {
    *SERVICE.get_or_init(|| {
        if std::env::var_os("RAMMAP_GPU").is_none() {
            return None;
        }
        let env = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
        let n_workers = (env("RAMMAP_GPU_WORKERS", 2.0) as usize).max(1);
        let svc: &'static Service = Box::leak(Box::new(Service {
            min_len: env("RAMMAP_GPU_MIN", 4096.0) as usize,
            mem_bytes: (env("RAMMAP_GPU_MEM_GB", 5.0) * (1u64 << 30) as f64) as i64,
            n_workers,
            scoring: OnceLock::new(),
            cache: (0..SHARDS).map(|_| Mutex::new(HashMap::new())).collect(),
            pending: Mutex::new(Pending::default()),
            work_cv: Condvar::new(),
            done_cv: Condvar::new(),
            ctxs: (0..n_workers).map(|_| Mutex::new(Ctx(std::ptr::null_mut()))).collect(),
            jobs: AtomicU64::new(0),
            rounds: AtomicU64::new(0),
            hits: AtomicU64::new(0),
            cells: AtomicU64::new(0),
            gpu_us: AtomicU64::new(0),
        }));
        for w in 0..n_workers {
            std::thread::Builder::new().name(format!("gpu-dp-{w}")).spawn(move || worker(svc, w)).expect("spawn GPU worker");
        }
        log::info!("GPU DP offload (speculate-and-retry): {n_workers} workers, {:.1} GB each, min segment {}",
                   svc.mem_bytes as f64 / (1u64 << 30) as f64, svc.min_len);
        Some(svc)
    })
}

fn key(qseq: &[u8], tseq: &[u8], bw: i32, zdrop: i32, end_bonus: i32, flags: i32) -> u128 {
    let mut lo = DefaultHasher::new();
    let mut hi = DefaultHasher::new();
    0x5eed_u64.hash(&mut hi);
    for h in [&mut lo, &mut hi] {
        (qseq, tseq, bw, zdrop, end_bonus, flags).hash(h);
    }
    ((hi.finish() as u128) << 64) | lo.finish() as u128
}

/// Use the GPU result for this DP call if eligible. Returns false (leaving `ez` untouched)
/// when the call is not eligible. On a cache miss, queues the request, switches the read to
/// speculative mode and fills `ez` with a placeholder.
#[allow(clippy::too_many_arguments)]
pub fn try_offload(qseq: &[u8], tseq: &[u8], mat: &[i8], q: i8, e: i8, q2: i8, e2: i8, bw: i32, zdrop: i32,
                   end_bonus: i32, flags: i32, ez: &mut DpResult) -> bool {
    if !IN_SCOPE.with(|c| c.get()) {
        return false;
    }
    let Some(svc) = service() else { return false };
    let (ql, tl) = (qseq.len(), tseq.len());
    if ql.max(tl) < svc.min_len || ql == 0 || tl == 0 || flags & SPLICE_FLAGS != 0 {
        return false;
    }
    if (q == q2 && e == e2) || (q2 as i32 + e2 as i32) < (q as i32 + e as i32) || mat.len() < 25 {
        return false;
    }
    let wl = if bw < 0 { ql.max(tl) } else { bw as usize };
    if ql.min(tl).min(wl + 1) > MAX_DIAG {
        return false;
    }
    let sc = svc.scoring.get_or_init(|| Scoring { mat: mat[..25].try_into().unwrap(), q, e, q2, e2 });
    if sc.mat[..] != mat[..25] || sc.q != q || sc.e != e || sc.q2 != q2 || sc.e2 != e2 {
        return false;
    }

    let k = key(qseq, tseq, bw, zdrop, end_bonus, flags);
    let hit = svc.cache[(k as usize) % SHARDS].lock().unwrap().get(&k).cloned();
    let Some(r) = hit else {
        let mut p = svc.pending.lock().unwrap();
        if p.keys.insert(k) {
            p.jobs.push(Job { key: k, q: qseq.to_vec(), t: tseq.to_vec(), bw, zdrop, end_bonus, flags });
            if p.jobs.len() >= BATCH_TARGET { svc.work_cv.notify_one(); }
        }
        drop(p);
        SPECULATING.with(|c| c.set(true));
        placeholder(ql, tl, flags, ez);
        return true;
    };
    svc.hits.fetch_add(1, Ordering::Relaxed);
    // Mirror the CPU kernels: fields the DP never writes keep the caller's values.
    ez.max = r.max;
    ez.max_score_query_pos = r.max_score_query_pos;
    ez.max_score_target_pos = r.max_score_target_pos;
    ez.max_query_end_score = r.max_query_end_score;
    if r.max_query_end_score != NEG_INF { ez.max_query_end_target_pos = r.max_query_end_target_pos; }
    ez.max_target_end_score = r.max_target_end_score;
    if r.max_target_end_score != NEG_INF { ez.max_target_end_query_pos = r.max_target_end_query_pos; }
    ez.score = r.score;
    if r.zdropped != 0 { ez.zdropped = 1; }
    if r.reach_end != 0 { ez.reach_end = 1; }
    ez.cigar = r.cigar;   // every call site passes a fresh DpResult, whose CIGAR is empty
    true
}

/// Jobs a worker waits to accumulate before launching mid-round (bigger batches fill the GPU).
const BATCH_TARGET: usize = 1024;
/// Most jobs per launch, so several workers share a large backlog.
const BATCH_MAX: usize = 4096;

fn worker(svc: &'static Service, w: usize) {
    loop {
        let mut jobs = {
            let mut p = svc.pending.lock().unwrap();
            loop {
                if p.jobs.len() >= BATCH_TARGET || (p.flush && !p.jobs.is_empty()) {
                    break;
                }
                p = svc.work_cv.wait_timeout(p, Duration::from_millis(50)).unwrap().0;
            }
            let take = p.jobs.len().min(BATCH_MAX);
            let batch: Vec<Job> = p.jobs.drain(..take).collect();
            p.in_flight += batch.len();
            batch
        };
        let t0 = Instant::now();
        // Largest first: memory-bounded sub-batches then hold similar sizes.
        jobs.sort_by_key(|j| std::cmp::Reverse(j.q.len().max(j.t.len())));
        let n = jobs.len();
        let keys: Vec<u128> = jobs.iter().map(|j| j.key).collect();
        run_share(svc, w, jobs);
        svc.jobs.fetch_add(n as u64, Ordering::Relaxed);
        svc.gpu_us.fetch_add(t0.elapsed().as_micros() as u64, Ordering::Relaxed);
        let mut p = svc.pending.lock().unwrap();
        p.in_flight -= n;
        for k in keys { p.keys.remove(&k); }
        drop(p);
        svc.done_cv.notify_all();
    }
}

/// End of a round: run everything still queued and wait until all requests are cached.
/// Returns the number of GPU jobs completed since the previous call.
pub fn wait_pending() -> usize {
    let Some(svc) = service() else { return 0 };
    static LAST: AtomicU64 = AtomicU64::new(0);
    let mut p = svc.pending.lock().unwrap();
    p.flush = true;
    svc.work_cv.notify_all();
    while !p.jobs.is_empty() || p.in_flight > 0 {
        p = svc.done_cv.wait(p).unwrap();
    }
    p.flush = false;
    drop(p);
    svc.rounds.fetch_add(1, Ordering::Relaxed);
    let now = svc.jobs.load(Ordering::Relaxed);
    (now - LAST.swap(now, Ordering::Relaxed)) as usize
}

fn run_share(svc: &Service, w: usize, jobs: Vec<Job>) {
    if jobs.is_empty() {
        return;
    }
    let mut ctx = svc.ctxs[w].lock().unwrap();
    if ctx.0.is_null() {
        let sc = svc.scoring.get().expect("scoring set before first job");
        ctx.0 = unsafe { gpu_ctx_create(svc.mem_bytes, sc.mat.as_ptr(), sc.q, sc.e, sc.q2, sc.e2) };
        assert!(!ctx.0.is_null(), "GPU context creation failed (lower RAMMAP_GPU_MEM_GB?)");
    }
    let (mut qbuf, mut tbuf) = (Vec::new(), Vec::new());
    let mut ins = Vec::with_capacity(jobs.len());
    let mut cig_off = Vec::with_capacity(jobs.len());
    let mut cig_total = 0i64;
    for j in &jobs {
        ins.push(SegIn { q_off: qbuf.len() as i64, t_off: tbuf.len() as i64, state_off: 0, tb_off: 0,
                         qlen: j.q.len() as i32, tlen: j.t.len() as i32, bw: j.bw, zdrop: j.zdrop,
                         flags: j.flags, end_bonus: j.end_bonus });
        qbuf.extend_from_slice(&j.q);
        tbuf.extend_from_slice(&j.t);
        cig_off.push(cig_total);
        cig_total += (j.q.len() + j.t.len() + 2) as i64;
    }
    let mut outs = vec![SegOut::default(); jobs.len()];
    let mut cigars = vec![0u32; cig_total as usize];
    let rc = unsafe {
        gpu_ctx_run(ctx.0, jobs.len() as c_int, ins.as_ptr(), qbuf.as_ptr(), qbuf.len() as i64, tbuf.as_ptr(),
                    tbuf.len() as i64, 1, cig_off.as_ptr(), cig_total, cigars.as_mut_ptr(), outs.as_mut_ptr(),
                    std::ptr::null_mut(), std::ptr::null_mut())
    };
    assert_eq!(rc, 0, "GPU DP batch failed");
    let mut cells = 0u64;
    for (k, (j, o)) in jobs.iter().zip(&outs).enumerate() {
        cells += o.cells as u64;
        let off = cig_off[k] as usize;
        let r = DpResult {
            max: o.max,
            max_score_query_pos: o.max_q,
            max_score_target_pos: o.max_t,
            max_query_end_score: o.mqe,
            max_query_end_target_pos: o.mqe_t,
            max_target_end_score: o.mte,
            max_target_end_query_pos: o.mte_q,
            score: o.score,
            reach_end: o.reach_end,
            zdropped: o.zdropped,
            cigar: cigars[off..off + o.n_cigar as usize].to_vec(),
            ..Default::default()
        };
        svc.cache[(j.key as usize) % SHARDS].lock().unwrap().insert(j.key, r);
    }
    svc.cells.fetch_add(cells, Ordering::Relaxed);
}

/// Drop all cached results (call when a chunk is fully aligned).
pub fn clear_cache() {
    if let Some(svc) = service() {
        for sh in &svc.cache {
            sh.lock().unwrap().clear();
        }
    }
}

/// One-line summary of GPU offload activity (None when the GPU path is off).
pub fn report() -> Option<String> {
    let svc = (*SERVICE.get()?)?;
    Some(format!("GPU DP offload ({} workers): {} segments over {} rounds, {} cache hits, {:.1} G cells, {:.1} s GPU worker time", svc.n_workers,
                 svc.jobs.load(Ordering::Relaxed), svc.rounds.load(Ordering::Relaxed), svc.hits.load(Ordering::Relaxed),
                 svc.cells.load(Ordering::Relaxed) as f64 / 1e9, svc.gpu_us.load(Ordering::Relaxed) as f64 / 1e6))
}
