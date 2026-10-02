//! GPU offload of large dual-affine extension / gap-fill DP calls (`--features gpu`).
//!
//! Enabled at run time with `RAMMAP_GPU=1`. Aligner threads call [`try_offload`] from the DP
//! dispatch; eligible calls are queued and the calling thread blocks until a GPU worker has
//! run them. Workers (`RAMMAP_GPU_WORKERS`, default 2) each own a CUDA stream and device pool
//! (`RAMMAP_GPU_MEM_GB` per worker, default 5), gather everything queued into one batch and run
//! `cuda/gpu_dp.cu` (forward pass + traceback), whose results are bit-identical to the AVX2 path.
//!
//! Because callers block, throughput depends on oversubscribing aligner threads (`-t` well
//! above the core count) so the CPU stays busy while many jobs wait on the GPU.
//!
//! Eligible: dual-affine (q != q2 or e != e2), q+e <= q2+e2, no splice flags, scoring matrix
//! equal to the first job's, max(qlen, tlen) >= `RAMMAP_GPU_MIN` (default 4096), and band
//! diagonal <= 16384 cells (the kernel's per-thread register bound).

use std::ffi::{c_int, c_void};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock};
use std::time::{Duration, Instant};

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
/// Most jobs one worker takes per batch, so several workers share a large queue.
const MAX_BATCH: usize = 512;

struct Job {
    q: Vec<u8>,
    t: Vec<u8>,
    bw: i32,
    zdrop: i32,
    end_bonus: i32,
    flags: i32,
    reply: Arc<(Mutex<Option<DpResult>>, Condvar)>,
}

struct Scoring {
    mat: [i8; 25],
    q: i8,
    e: i8,
    q2: i8,
    e2: i8,
}

struct Service {
    min_len: usize,
    mem_bytes: i64,
    scoring: OnceLock<Scoring>,
    queue: Mutex<Vec<Job>>,
    cv: Condvar,
    jobs: AtomicU64,
    batches: AtomicU64,
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
        let svc: &'static Service = Box::leak(Box::new(Service {
            min_len: env("RAMMAP_GPU_MIN", 4096.0) as usize,
            mem_bytes: (env("RAMMAP_GPU_MEM_GB", 5.0) * (1u64 << 30) as f64) as i64,
            scoring: OnceLock::new(),
            queue: Mutex::new(Vec::new()),
            cv: Condvar::new(),
            jobs: AtomicU64::new(0),
            batches: AtomicU64::new(0),
            cells: AtomicU64::new(0),
            gpu_us: AtomicU64::new(0),
        }));
        let workers = env("RAMMAP_GPU_WORKERS", 2.0) as usize;
        for w in 0..workers {
            std::thread::Builder::new()
                .name(format!("gpu-dp-{w}"))
                .spawn(move || worker(svc))
                .expect("spawn GPU worker");
        }
        log::info!("GPU DP offload: {workers} workers, {:.1} GB each, min segment {}",
                   svc.mem_bytes as f64 / (1u64 << 30) as f64, svc.min_len);
        Some(svc)
    })
}

/// Run this DP call on the GPU if eligible; returns false (leaving `ez` untouched) otherwise.
#[allow(clippy::too_many_arguments)]
pub fn try_offload(qseq: &[u8], tseq: &[u8], mat: &[i8], q: i8, e: i8, q2: i8, e2: i8, bw: i32, zdrop: i32,
                   end_bonus: i32, flags: i32, ez: &mut DpResult) -> bool {
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

    let reply = Arc::new((Mutex::new(None), Condvar::new()));
    {
        let mut qu = svc.queue.lock().unwrap();
        qu.push(Job { q: qseq.to_vec(), t: tseq.to_vec(), bw, zdrop, end_bonus, flags, reply: reply.clone() });
    }
    svc.cv.notify_one();
    let (lock, cv) = &*reply;
    let mut slot = lock.lock().unwrap();
    while slot.is_none() {
        slot = cv.wait(slot).unwrap();
    }
    let r = slot.take().unwrap();
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

fn worker(svc: &'static Service) {
    let mut ctx: *mut c_void = std::ptr::null_mut();
    loop {
        let jobs: Vec<Job> = {
            let mut qu = svc.queue.lock().unwrap();
            while qu.is_empty() {
                qu = svc.cv.wait(qu).unwrap();
            }
            // Give concurrent submitters a moment to pile on, so batches are not tiny.
            if qu.len() < 64 {
                qu = svc.cv.wait_timeout(qu, Duration::from_millis(2)).unwrap().0;
            }
            let take = qu.len().min(MAX_BATCH);
            qu.drain(..take).collect()
        };
        if jobs.is_empty() {
            continue;
        }
        if ctx.is_null() {
            let sc = svc.scoring.get().expect("scoring set before first job");
            ctx = unsafe { gpu_ctx_create(svc.mem_bytes, sc.mat.as_ptr(), sc.q, sc.e, sc.q2, sc.e2) };
            assert!(!ctx.is_null(), "GPU context creation failed (lower RAMMAP_GPU_MEM_GB?)");
        }
        run_batch(svc, ctx, jobs);
    }
}

fn run_batch(svc: &Service, ctx: *mut c_void, jobs: Vec<Job>) {
    let t0 = Instant::now();
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
        gpu_ctx_run(ctx, jobs.len() as c_int, ins.as_ptr(), qbuf.as_ptr(), qbuf.len() as i64, tbuf.as_ptr(),
                    tbuf.len() as i64, 1, cig_off.as_ptr(), cig_total, cigars.as_mut_ptr(), outs.as_mut_ptr(),
                    std::ptr::null_mut(), std::ptr::null_mut())
    };
    assert_eq!(rc, 0, "GPU DP batch failed");
    let mut cells = 0u64;
    for (k, (j, o)) in jobs.into_iter().zip(&outs).enumerate() {
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
        let (lock, cv) = &*j.reply;
        *lock.lock().unwrap() = Some(r);
        cv.notify_one();
    }
    svc.jobs.fetch_add(outs.len() as u64, Ordering::Relaxed);
    svc.batches.fetch_add(1, Ordering::Relaxed);
    svc.cells.fetch_add(cells, Ordering::Relaxed);
    svc.gpu_us.fetch_add(t0.elapsed().as_micros() as u64, Ordering::Relaxed);
}

/// One-line summary of GPU offload activity (empty when the GPU path is off).
pub fn report() -> Option<String> {
    let svc = (*SERVICE.get()?)?;
    let (jobs, batches) = (svc.jobs.load(Ordering::Relaxed), svc.batches.load(Ordering::Relaxed));
    Some(format!("GPU DP offload: {jobs} segments in {batches} batches (mean {:.1}/batch), {:.1} G cells, {:.1} s worker time",
                 jobs as f64 / batches.max(1) as f64, svc.cells.load(Ordering::Relaxed) as f64 / 1e9,
                 svc.gpu_us.load(Ordering::Relaxed) as f64 / 1e6))
}
