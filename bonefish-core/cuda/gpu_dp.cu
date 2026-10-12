// GPU port of bonefish's dual-affine banded extension DP forward pass (ksw2 extd2 semantics).
//
// One thread block per segment. Cells on an anti-diagonal depend only on the previous
// anti-diagonal, so a block computes each diagonal in parallel and synchronises between
// diagonals. Bit-exact with bonefish's AVX2 path (`extend_dual_affine_avx2_fn` in
// bonefish-core/src/align/dp/dual.rs), which reproduces minimap2's SSE ksw2 quirks:
//   - each diagonal computes every cell of the 16-aligned band [st0/16*16, ceil16(en0+1)-1];
//   - a persistent score array is refreshed only in 16-byte blocks from st0, so cells outside
//     that range use stale scores; reads past either sequence end see zero padding (base A);
//   - the boundary carry-in, new-entry init and last_st/last_en use the rounded band;
// and otherwise follows the same recurrence as `extend_dual_affine_scalar`:
//   - difference-encoded int8 state (u, v, x, y, x2, y2) with wrapping arithmetic;
//   - x, x2 and v are read at t-1 from the previous diagonal, so they are double-buffered;
//   - exact max tracking reproduces the CPU's 4-lane reduction tie-break order;
//   - approximate max (APPROX_MAX) follows the same single-path walk;
//   - z-drop, query-end and target-end tracking and the final-corner score as on the CPU;
//   - optional traceback bytes in the CPU's layout (traceback itself still runs on the CPU).
#include <climits>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <cuda_runtime.h>

#define NEG_INF (-0x40000000)
#define APPROX_MAX 0x08
#define APPROX_DROP 0x10
#define GENERIC_SCORING 0x04
#define RIGHT_ALIGN 0x02
#define SCORE_ONLY 0x01
#define EXTENSION_ONLY 0x40
#define REV_CIGAR 0x80

struct SegIn {
  int64_t q_off, t_off;           // byte offsets into the packed query/target buffers
  int64_t state_off, tb_off;      // byte offsets into the device pool (set by the host batcher)
  int32_t qlen, tlen, bw, zdrop, flags, end_bonus;
};

struct SegOut {
  int32_t max, max_q, max_t, zdropped, score;
  int32_t mqe, mqe_t, mte, mte_q;
  int32_t n_diag;                 // diagonals computed (for cell accounting)
  int32_t reach_end, n_cigar;     // set by the traceback kernel
  int64_t cells;
};

struct Scoring {
  int8_t mat[25];
  int8_t q, e, q2, e2;
};

__device__ __forceinline__ int8_t w8(int v) { return (int8_t)v; }  // wrap to int8

constexpr int BLOCK = 256;

__global__ void __launch_bounds__(BLOCK, 4)
dp_fwd(const SegIn* __restrict__ segs, SegOut* __restrict__ outs,
       const uint8_t* __restrict__ qbuf, const uint8_t* __restrict__ tbuf,
       uint8_t* __restrict__ tb_pool, Scoring sc, int with_tb) {
  const SegIn s = segs[blockIdx.x];
  const int tid = threadIdx.x;
  const int qlen = s.qlen, tlen = s.tlen;
  const uint8_t* q = qbuf + s.q_off;
  const uint8_t* tq = tbuf + s.t_off;

  const int ge = sc.e, ge2 = sc.e2, go = sc.q, go2 = sc.q2;
  const int qe = go + ge, qe2 = go2 + ge2;
  const int8_t neg_qe = w8(-qe), neg_qe2 = w8(-qe2);
  const bool approx = (s.flags & APPROX_MAX) != 0;
  const bool generic = (s.flags & GENERIC_SCORING) != 0;
  const bool right = (s.flags & RIGHT_ALIGN) != 0;
  const bool tb = with_tb && !(s.flags & SCORE_ONLY);

  int long_thres = (ge != ge2) ? (go2 - go) / (ge - ge2) - 1 : 0;
  if (go2 + ge2 + long_thres * ge2 > go + ge + long_thres * ge) long_thres += 1;
  const int8_t long_diff = w8(long_thres * (ge - ge2) - (go2 - go) - ge2);
  const int wl = s.bw < 0 ? max(qlen, tlen) : s.bw;
  const int8_t sc_mch = sc.mat[0], sc_mis = sc.mat[1];
  const int8_t sc_n = sc.mat[24] == 0 ? w8(-ge2) : sc.mat[24];
  const int n_col = (min(min(qlen, tlen), wl + 1) + 31) / 32 + 1;   // AVX2 traceback layout
  const int stride = n_col * 32;

  // State slab: 10 int8 arrays u, y, y2, x[2], x2[2], v[2], s (each padded) + int32 h.
  // Padding covers the rounded band (<= ceil16(tlen)-1) and score writes up to en0+15.
  const int pad = ((tlen + 15) & ~15) + 16;
  int8_t* base = (int8_t*)(tb_pool + s.state_off);
  int8_t* u = base;
  int8_t* y = u + pad;
  int8_t* y2 = y + pad;
  int8_t* xb[2] = {y2 + pad, y2 + 2 * pad};
  int8_t* x2b[2] = {y2 + 3 * pad, y2 + 4 * pad};
  int8_t* vb[2] = {y2 + 5 * pad, y2 + 6 * pad};
  int8_t* sarr = base + 9 * pad;
  int32_t* h = (int32_t*)(base + 10 * pad);
  uint8_t* p = tb ? tb_pool + s.tb_off : nullptr;
  int32_t* band_off = tb ? (int32_t*)(p + (int64_t)(qlen + tlen - 1) * stride) : nullptr;
  int32_t* band_end = tb ? band_off + (qlen + tlen - 1) : nullptr;

  for (int t = tid; t < pad; t += BLOCK) {
    u[t] = neg_qe; y[t] = neg_qe; y2[t] = neg_qe2;
    xb[0][t] = xb[1][t] = neg_qe; x2b[0][t] = x2b[1][t] = neg_qe2; vb[0][t] = vb[1][t] = neg_qe;
    sarr[t] = 0;
    if (!approx) h[t] = 0;
  }

  __shared__ long long red[BLOCK / 32];
  __shared__ int sh_h_en0, sh_h_st0, sh_stop;
  // Result state (meaningful in thread 0; broadcast what others need).
  int r_max = 0, r_max_q = -1, r_max_t = -1, r_zd = 0, r_score = NEG_INF;
  int r_mqe = NEG_INF, r_mqe_t = -1, r_mte = NEG_INF, r_mte_q = -1;
  int last_st = -1, last_en = -1, h0 = 0, last_h0_t = 0, n_diag = 0;
  long long cells = 0;
  if (tid == 0) sh_stop = 0;
  __syncthreads();

  const int valid = qlen + tlen - 1;
  for (int r = 0; r < valid; ++r) {
    int st0 = max(0, r - qlen + 1), en0 = min(tlen - 1, r);
    st0 = max(st0, (r - wl + 1) >> 1);
    en0 = min(en0, (r + wl) >> 1);
    if (st0 > en0) { r_zd = 1; break; }
    const int st = st0 / 16 * 16, en = (en0 + 16) / 16 * 16 - 1;
    const int cur = r & 1, prv = cur ^ 1;
    int8_t* xc = xb[cur]; int8_t* x2c = x2b[cur]; int8_t* vc = vb[cur];
    const int8_t* xp = xb[prv]; const int8_t* x2p = x2b[prv]; const int8_t* vp = vb[prv];

    // Score refresh, as the SIMD code does: 16-byte blocks from st0 while the block start <= en0.
    {
      const int s_end = generic ? en0 : st0 + 16 * ((en0 - st0) / 16 + 1) - 1;
      for (int t = st0 + tid; t <= s_end; t += BLOCK) {
        const uint8_t sq = t < tlen ? tq[t] : 0;
        const uint8_t sv = t <= r ? q[r - t] : 0;   // reversed-query padding is zero
        sarr[t] = generic ? sc.mat[sq * 5 + sv] : ((sq == 4 || sv == 4) ? sc_n : (sq == sv ? sc_mch : sc_mis));
      }
    }
    __syncthreads();

    int8_t x1, x21, v1;
    if (st > 0) {
      if (st > last_st && st - 1 <= last_en) { x1 = xp[st - 1]; x21 = x2p[st - 1]; v1 = vp[st - 1]; }
      else { x1 = neg_qe; x21 = neg_qe2; v1 = neg_qe; }
    } else {
      x1 = neg_qe; x21 = neg_qe2;
      v1 = r == 0 ? neg_qe : r < long_thres ? w8(-ge) : r == long_thres ? long_diff : w8(-ge2);
    }
    if (tb && tid == 0) { band_off[r] = st; band_end[r] = en; }

    for (int t = st + tid; t <= en; t += BLOCK) {
      if (t == r) {  // new diagonal entry (when the rounded en >= r)
        y[t] = neg_qe; y2[t] = neg_qe2;
        u[t] = r == 0 ? neg_qe : r < long_thres ? w8(-ge) : r == long_thres ? long_diff : w8(-ge2);
      }
      const int8_t zs = sarr[t];
      const int8_t xt1 = t == st ? x1 : xp[t - 1];
      const int8_t x2t1 = t == st ? x21 : x2p[t - 1];
      const int8_t vt1 = t == st ? v1 : vp[t - 1];
      const int8_t ut = u[t];
      const int8_t a = w8(xt1 + vt1), b = w8(y[t] + ut), a2 = w8(x2t1 + vt1), b2 = w8(y2[t] + ut);
      int8_t z = zs; uint8_t d = 0;
      if (!right) {
        if (a > z) { d = 1; z = a; }
        if (b > z) { d = 2; z = b; }
        if (a2 > z) { d = 3; z = a2; }
        if (b2 > z) { d = 4; z = b2; }
      } else {
        if (z > a) d = 0; else { d = 1; z = a; }
        if (z <= b) { d = 2; z = b; }
        if (z <= a2) { d = 3; z = a2; }
        if (z <= b2) { d = 4; z = b2; }
      }
      if (z > sc_mch) z = sc_mch;
      u[t] = w8(z - vt1);
      vc[t] = w8(z - ut);
      const int8_t a_n = w8(a - w8(z - go)), b_n = w8(b - w8(z - go));
      const int8_t a2_n = w8(a2 - w8(z - go2)), b2_n = w8(b2 - w8(z - go2));
      const bool ka = right ? a_n >= 0 : a_n > 0, kb = right ? b_n >= 0 : b_n > 0;
      const bool ka2 = right ? a2_n >= 0 : a2_n > 0, kb2 = right ? b2_n >= 0 : b2_n > 0;
      xc[t] = w8((ka ? a_n : 0) - qe);
      y[t] = w8((kb ? b_n : 0) - qe);
      x2c[t] = w8((ka2 ? a2_n : 0) - qe2);
      y2[t] = w8((kb2 ? b2_n : 0) - qe2);
      if (tb) {
        d |= (ka ? 0x08 : 0) | (kb ? 0x10 : 0) | (ka2 ? 0x20 : 0) | (kb2 ? 0x40 : 0);
        p[(int64_t)r * stride + (t - st)] = d;
      }
    }
    // No carry-over needed: st0 and en0 never decrease, so every t-1 read by diagonal r+1
    // (and its st0-1 boundary read) lies in diagonal r's band and is in the buffer just written.
    cells += (tid == 0) ? (en - st + 1) : 0;
    ++n_diag;
    __syncthreads();

    if (!approx) {
      // h update: read old values, sync, write; reduce max with the CPU's tie-break order.
      const int en1 = st0 + (en0 - st0) / 4 * 4;
      long long best = LLONG_MIN;
      int nh[64]; int cnt = 0;   // per-thread new h values; host routes diagonals > 64*BLOCK cells to CPU
      for (int t = st0 + tid; t <= en0; t += BLOCK) {
        int hn;
        if (r == 0) hn = (int)vc[0] - qe;
        else if (t == en0) hn = en0 > 0 ? h[en0 - 1] + (int)u[en0] : h[en0] + (int)vc[en0];
        else hn = h[t] + (int)vc[t];
        unsigned prio;
        if (t == en0) prio = 0;
        else if (t < en1) prio = 1u + ((unsigned)((t - st0) & 3) << 20) + (unsigned)((t - st0) >> 2);
        else prio = 1u + (4u << 20) + (unsigned)(t - en1);
        const long long key = ((long long)hn << 32) | (long long)(0xFFFFFFFFu - prio);
        if (key > best) best = key;
        nh[cnt++] = hn;
      }
      __syncthreads();
      cnt = 0;
      for (int t = st0 + tid; t <= en0; t += BLOCK) {
        h[t] = nh[cnt++];
        if (t == en0) sh_h_en0 = h[t];
        if (t == st0) sh_h_st0 = h[t];
      }
      for (int o = 16; o > 0; o >>= 1) { long long v2 = __shfl_down_sync(0xffffffff, best, o); if (v2 > best) best = v2; }
      if ((tid & 31) == 0) red[tid >> 5] = best;
      __syncthreads();
      if (tid == 0) {
        for (int i = 0; i < BLOCK / 32; ++i) if (red[i] > best) best = red[i];
        const int max_h = (int)(best >> 32);
        const unsigned prio = 0xFFFFFFFFu - (unsigned)(best & 0xFFFFFFFF);
        int max_t;
        if (prio == 0) max_t = en0;
        else if (prio < 1u + (4u << 20)) { unsigned k = prio - 1u; max_t = st0 + (int)((k & 0xFFFFF) << 2) + (int)(k >> 20); }
        else max_t = en1 + (int)(prio - 1u - (4u << 20));
        if (en0 == tlen - 1 && sh_h_en0 > r_mte) { r_mte = sh_h_en0; r_mte_q = r - en0; }
        if (r - st0 == qlen - 1 && sh_h_st0 > r_mqe) { r_mqe = sh_h_st0; r_mqe_t = st0; }
        if (max_h > r_max) { r_max = max_h; r_max_t = max_t; r_max_q = r - max_t; }
        else if (max_t >= r_max_t && r - max_t >= r_max_q) {
          const int tl = max_t - r_max_t, ql = (r - max_t) - r_max_q;
          const int l = tl > ql ? tl - ql : ql - tl;
          if (s.zdrop >= 0 && r_max - max_h > s.zdrop + l * ge2) { r_zd = 1; sh_stop = 1; }
        }
        if (r == qlen + tlen - 2 && en0 == tlen - 1) r_score = sh_h_en0;
      }
    } else if (tid == 0) {
      if (r > 0) {
        if (last_h0_t >= st0 && last_h0_t <= en0 && last_h0_t + 1 >= st0 && last_h0_t < en0) {
          const int d0 = vc[last_h0_t], d1 = u[last_h0_t + 1];
          if (d0 > d1) h0 += d0; else { h0 += d1; ++last_h0_t; }
        } else if (last_h0_t >= st0 && last_h0_t <= en0) h0 += vc[last_h0_t];
        else { ++last_h0_t; h0 += u[last_h0_t]; }
      } else { h0 = (int)vc[0] - qe; last_h0_t = 0; }
      if (s.flags & APPROX_DROP) {
        if (h0 > r_max) { r_max = h0; r_max_t = last_h0_t; r_max_q = r - last_h0_t; }
        else if (last_h0_t >= r_max_t && r - last_h0_t >= r_max_q) {
          const int tl = last_h0_t - r_max_t, ql = (r - last_h0_t) - r_max_q;
          const int l = tl > ql ? tl - ql : ql - tl;
          if (s.zdrop >= 0 && r_max - h0 > s.zdrop + l * ge2) { r_zd = 1; sh_stop = 1; }
        }
      }
      if (r == qlen + tlen - 2 && en0 == tlen - 1) r_score = h0;
    }
    last_st = st; last_en = en;
    __syncthreads();
    if (sh_stop) break;
  }
  if (tid == 0) {
    if (approx && r_score == NEG_INF) r_score = r_max;
    SegOut o; o.max = r_max; o.max_q = r_max_q; o.max_t = r_max_t; o.zdropped = r_zd; o.score = r_score;
    o.mqe = r_mqe; o.mqe_t = r_mqe_t; o.mte = r_mte; o.mte_q = r_mte_q; o.n_diag = n_diag; o.cells = cells;
    o.reach_end = 0; o.n_cigar = 0;
    outs[blockIdx.x] = o;
  }
}

// Traceback, one thread per segment: port of `traceback_start_position` + `traceback_dual_affine`
// (bonefish-core/src/align/dp/common.rs) over the AVX2-layout traceback bytes written by dp_fwd.
// CIGAR ops are (len << 4 | op), op 0=M 1=I 2=D, merged like `push_cigar`.
__device__ __forceinline__ void push_cig(uint32_t* c, int& n, uint32_t op, uint32_t len) {
  if (n > 0 && (c[n - 1] & 0xf) == op) c[n - 1] += len << 4;
  else c[n++] = (len << 4) | op;
}

__global__ void dp_tb(const SegIn* __restrict__ segs, SegOut* __restrict__ outs, const uint8_t* __restrict__ pool,
                      uint32_t* __restrict__ cig, const int64_t* __restrict__ cig_off, int n) {
  const int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n) return;
  const SegIn s = segs[k];
  SegOut o = outs[k];
  const int qlen = s.qlen, tlen = s.tlen;
  const int wl = s.bw < 0 ? max(qlen, tlen) : s.bw;
  const int stride = ((min(min(qlen, tlen), wl + 1) + 31) / 32 + 1) * 32;
  const int vr = qlen + tlen - 1;
  const uint8_t* p = pool + s.tb_off;
  const int32_t* band_off = (const int32_t*)(p + (int64_t)vr * stride);
  const int32_t* band_end = band_off + vr;
  uint32_t* c = cig + cig_off[k];
  int nc = 0, i, j;
  o.reach_end = 0;
  if (!o.zdropped && !(s.flags & EXTENSION_ONLY)) { i = tlen - 1; j = qlen - 1; }
  else if (!o.zdropped && (s.flags & EXTENSION_ONLY) && o.mqe + s.end_bonus > o.max) { o.reach_end = 1; i = o.mqe_t; j = qlen - 1; }
  else if (o.max_t >= 0 && o.max_q >= 0) { i = o.max_t; j = o.max_q; }
  else { o.n_cigar = 0; outs[k] = o; return; }
  int state = 0;
  while (i >= 0 && j >= 0) {
    const int r = i + j;
    const int off_r = band_off[r], off_end_r = band_end[r];
    int force = -1;
    if (i < off_r) force = 2;
    if (i > off_end_r) force = 1;
    const uint8_t tmp = force < 0 ? p[(int64_t)r * stride + (i - off_r)] : 0;
    if (state == 0) state = tmp & 7;
    else if (((tmp >> (state + 2)) & 1) == 0) state = 0;
    if (state == 0) state = tmp & 7;
    if (force >= 0) state = force;
    switch (state) {
      case 1: case 3: push_cig(c, nc, 2, 1); --i; break;           // D1, D2
      case 2: case 4: push_cig(c, nc, 1, 1); --j; break;           // I1, I2
      default: push_cig(c, nc, 0, 1); --i; --j; break;             // M
    }
  }
  if (i >= 0) push_cig(c, nc, 2, i + 1);
  if (j >= 0) push_cig(c, nc, 1, j + 1);
  if (!(s.flags & REV_CIGAR))
    for (int a = 0, b = nc - 1; a < b; ++a, --b) { uint32_t t = c[a]; c[a] = c[b]; c[b] = t; }
  o.n_cigar = nc;
  outs[k] = o;
}

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return -1; } } while (0)

// ---------------------------------------------------------------------------------------------
// Host API. A context owns a CUDA stream, a device pool (DP state + traceback bytes) and growable
// I/O buffers, so steady-state batches do no cudaMalloc (which synchronises the whole device).
// One context per host worker thread; contexts on different streams overlap on the GPU.
// ---------------------------------------------------------------------------------------------
struct GpuCtx {
  cudaStream_t st;
  uint8_t* pool; int64_t budget;
  uint8_t *dq, *dt; int64_t qcap, tcap;
  SegIn* dsegs; int64_t ncap;
  SegOut* douts; int64_t ocap;
  uint32_t* dcig; int64_t cigcap;
  int64_t* dcigoff; int64_t offcap;
  Scoring sc;
  std::vector<SegIn> batch;
};

template <class T> static int grow(T** p, int64_t* cap, int64_t need) {
  if (need <= *cap) return 0;
  if (*p) CK(cudaFree(*p));
  const int64_t c = need + need / 2;
  CK(cudaMalloc(p, sizeof(T) * c));
  *cap = c;
  return 0;
}

// Per-segment device bytes: DP state slab + (optional) traceback matrix and band offsets.
static void seg_bytes(const SegIn& s, int with_tb, int64_t* st_bytes, int64_t* tb_bytes) {
  const int pad = ((s.tlen + 15) & ~15) + 16;
  *st_bytes = (int64_t)14 * pad;   // 10 int8 arrays + int32 h
  *tb_bytes = 0;
  if (with_tb) {
    const int wl = s.bw < 0 ? (s.qlen > s.tlen ? s.qlen : s.tlen) : s.bw;
    int mn = s.qlen < s.tlen ? s.qlen : s.tlen; if (mn > wl + 1) mn = wl + 1;
    const int64_t stride = ((mn + 31) / 32 + 1) * 32;
    const int64_t vr = s.qlen + s.tlen - 1;
    *tb_bytes = (vr * stride + 8 * vr + 15) & ~15LL;
  }
}

extern "C" int64_t gpu_seg_bytes(int qlen, int tlen, int bw, int with_tb) {
  SegIn s{}; s.qlen = qlen; s.tlen = tlen; s.bw = bw;
  int64_t a, b; seg_bytes(s, with_tb, &a, &b); return a + b;
}

extern "C" void* gpu_ctx_create(int64_t budget, const int8_t* mat, int8_t q, int8_t e, int8_t q2, int8_t e2) {
  GpuCtx* c = new GpuCtx();
  for (int i = 0; i < 25; ++i) c->sc.mat[i] = mat[i];
  c->sc.q = q; c->sc.e = e; c->sc.q2 = q2; c->sc.e2 = e2;
  c->budget = budget;
  c->dq = c->dt = nullptr; c->qcap = c->tcap = 0; c->dsegs = nullptr; c->ncap = 0;
  c->douts = nullptr; c->ocap = 0; c->dcig = nullptr; c->cigcap = 0; c->dcigoff = nullptr; c->offcap = 0;
  if (cudaStreamCreateWithFlags(&c->st, cudaStreamNonBlocking) != cudaSuccess ||
      cudaMalloc(&c->pool, budget) != cudaSuccess) {
    fprintf(stderr, "gpu_ctx_create: CUDA allocation of %lld bytes failed\n", (long long)budget);
    delete c; return nullptr;
  }
  return c;
}

extern "C" void gpu_ctx_destroy(void* h) {
  GpuCtx* c = (GpuCtx*)h; if (!c) return;
  cudaFree(c->pool); cudaFree(c->dq); cudaFree(c->dt); cudaFree(c->dsegs); cudaFree(c->douts);
  cudaFree(c->dcig); cudaFree(c->dcigoff); cudaStreamDestroy(c->st);
  delete c;
}

// Run n segments (sub-batched by the context's memory budget). With with_tb, CIGARs land in
// `cigars` at cig_off[i] (capacity qlen+tlen+2 each). kernel_ms may be null.
extern "C" int gpu_ctx_run(void* h, int n, const SegIn* segs_in, const uint8_t* qbuf, int64_t qbytes,
                           const uint8_t* tbuf, int64_t tbytes, int with_tb,
                           const int64_t* cig_off, int64_t cig_total, uint32_t* cigars,
                           SegOut* outs, double* kernel_ms, int* n_batches) {
  GpuCtx* c = (GpuCtx*)h;
  if (n == 0) return 0;
  if (grow(&c->dq, &c->qcap, qbytes) || grow(&c->dt, &c->tcap, tbytes) ||
      grow(&c->dsegs, &c->ncap, (int64_t)n) || grow(&c->douts, &c->ocap, (int64_t)n)) return -1;
  if (with_tb) {
    if (grow(&c->dcig, &c->cigcap, cig_total) || grow(&c->dcigoff, &c->offcap, (int64_t)n)) return -1;
    CK(cudaMemcpyAsync(c->dcigoff, cig_off, sizeof(int64_t) * n, cudaMemcpyHostToDevice, c->st));
  }
  CK(cudaMemcpyAsync(c->dq, qbuf, qbytes, cudaMemcpyHostToDevice, c->st));
  CK(cudaMemcpyAsync(c->dt, tbuf, tbytes, cudaMemcpyHostToDevice, c->st));
  cudaEvent_t k0, k1; CK(cudaEventCreate(&k0)); CK(cudaEventCreate(&k1));
  CK(cudaEventRecord(k0, c->st));
  int nb = 0, i = 0;
  // Each sub-batch's SegIn array is written to its own slice of dsegs, so the async copies of
  // successive sub-batches never overwrite one another while earlier kernels still read them.
  while (i < n) {
    c->batch.clear(); int64_t used = 0; const int start = i;
    while (i < n) {
      SegIn s = segs_in[i];
      int64_t stb, tbb; seg_bytes(s, with_tb, &stb, &tbb);
      if (used + stb + tbb > c->budget && !c->batch.empty()) break;
      if (stb + tbb > c->budget) { fprintf(stderr, "segment needs %lld bytes > budget\n", (long long)(stb + tbb)); return -1; }
      s.state_off = used; s.tb_off = used + stb; used += stb + tbb;
      c->batch.push_back(s); ++i;
    }
    const int nbt = (int)c->batch.size();
    CK(cudaMemcpyAsync(c->dsegs + start, c->batch.data(), sizeof(SegIn) * nbt, cudaMemcpyHostToDevice, c->st));
    CK(cudaStreamSynchronize(c->st));   // batch vector is reused next iteration
    dp_fwd<<<nbt, BLOCK, 0, c->st>>>(c->dsegs + start, c->douts + start, c->dq, c->dt, c->pool, c->sc, with_tb);
    CK(cudaGetLastError());
    if (with_tb) {
      dp_tb<<<(nbt + 63) / 64, 64, 0, c->st>>>(c->dsegs + start, c->douts + start, c->pool, c->dcig, c->dcigoff + start, nbt);
      CK(cudaGetLastError());
    }
    ++nb;
  }
  CK(cudaEventRecord(k1, c->st));
  CK(cudaMemcpyAsync(outs, c->douts, sizeof(SegOut) * n, cudaMemcpyDeviceToHost, c->st));
  if (with_tb) CK(cudaMemcpyAsync(cigars, c->dcig, sizeof(uint32_t) * cig_total, cudaMemcpyDeviceToHost, c->st));
  CK(cudaStreamSynchronize(c->st));
  if (kernel_ms) { float ms; CK(cudaEventElapsedTime(&ms, k0, k1)); *kernel_ms = ms; }
  cudaEventDestroy(k0); cudaEventDestroy(k1);
  if (n_batches) *n_batches = nb;
  return 0;
}

// One-shot wrapper used by the kernels bake-off: create a context, run, destroy.
extern "C" int gpu_dp_run(int n, const SegIn* segs_in, const uint8_t* qbuf, int64_t qbytes,
                          const uint8_t* tbuf, int64_t tbytes, const int8_t* mat, int8_t q, int8_t e,
                          int8_t q2, int8_t e2, int with_tb, int64_t mem_budget,
                          const int64_t* cig_off, int64_t cig_total, uint32_t* cigars,
                          SegOut* outs, double* kernel_ms, double* total_ms, int* n_batches) {
  cudaEvent_t t0, t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
  CK(cudaEventRecord(t0));
  void* c = gpu_ctx_create(mem_budget, mat, q, e, q2, e2);
  if (!c) return -1;
  const int rc = gpu_ctx_run(c, n, segs_in, qbuf, qbytes, tbuf, tbytes, with_tb, cig_off, cig_total, cigars,
                             outs, kernel_ms, n_batches);
  gpu_ctx_destroy(c);
  CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
  float ms; CK(cudaEventElapsedTime(&ms, t0, t1)); *total_ms = ms;
  return rc;
}
