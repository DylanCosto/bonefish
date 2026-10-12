//! Gap-affine 2-piece wavefront alignment (WFA, Marco-Sola et al. 2021) for global gap fills.
//!
//! Scores: match `a`, mismatch `-b`, a gap of length `l` costs `min(q + e*l, q2 + e2*l)` (the
//! ksw2 extd2 model). Global alignment maximising that score is the same as minimising the
//! penalty `x*X + sum over gaps of min(o1 + e1*l, o2 + e2*l)` with `x = 2(a+b)`, `o = 2q`,
//! `e = 2e + a` (Eizenga & Paten 2022), and `score = (a*(n+m) - penalty) / 2`.
//!
//! Work grows with the penalty, not with the sequence lengths, so near-identical gaps (most of
//! a high-accuracy read) align far faster than with banded DP. `align` gives up (returns None)
//! once the penalty passes `max_pen`; callers then use the DP. The sequences are aligned
//! reversed so that forward match extension, which is greedy, pushes gaps to the left, like
//! ksw2's left-aligned traceback.

const NONE: i32 = i32::MIN / 2;

#[derive(Clone, Copy)]
pub struct Penalties { pub x: i32, pub o1: i32, pub e1: i32, pub o2: i32, pub e2: i32 }

impl Penalties {
    /// From ksw2 dual-affine scores (match `a`, mismatch `b` > 0, gap pieces (q, e), (q2, e2)).
    pub fn from_scores(a: i32, b: i32, q: i32, e: i32, q2: i32, e2: i32) -> Self {
        Penalties { x: 2 * (a + b), o1: 2 * q, e1: 2 * e + a, o2: 2 * q2, e2: 2 * e2 + a }
    }
}

/// One wavefront: furthest-reaching offsets (text position h) per diagonal k = h - v in lo..=hi.
#[derive(Clone, Default)]
struct Wf { lo: i32, hi: i32, off: Vec<i32> }

impl Wf {
    #[inline(always)]
    fn get(&self, k: i32) -> i32 {
        if k < self.lo || k > self.hi { NONE } else { self.off[(k - self.lo) as usize] }
    }
}

#[derive(Default)]
struct Wfs { m: Vec<Option<Wf>>, i1: Vec<Option<Wf>>, d1: Vec<Option<Wf>>, i2: Vec<Option<Wf>>, d2: Vec<Option<Wf>> }

#[inline(always)]
fn at(v: &[Option<Wf>], s: i32) -> Option<&Wf> {
    if s < 0 { None } else { v.get(s as usize).and_then(|w| w.as_ref()) }
}

#[inline(always)]
fn get(v: &[Option<Wf>], s: i32, k: i32) -> i32 {
    at(v, s).map_or(NONE, |w| w.get(k))
}

/// Length of the common prefix of `a` and `b` (nt4 codes).
#[inline(always)]
fn lcp(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    while i + 8 <= n {
        let x = u64::from_le_bytes(a[i..i + 8].try_into().unwrap()) ^ u64::from_le_bytes(b[i..i + 8].try_into().unwrap());
        if x != 0 { return i + (x.trailing_zeros() / 8) as usize; }
        i += 8;
    }
    while i < n && a[i] == b[i] { i += 1; }
    i
}

/// Global alignment of `q` (pattern, read) against `t` (text, reference), both nt4 with no
/// ambiguous bases. Returns (penalty, CIGAR as ksw2 ops `len << 4 | op`, 0 M / 1 I / 2 D) in
/// forward orientation, or None if the penalty would exceed `max_pen`.
pub fn align(q: &[u8], t: &[u8], p: &Penalties, max_pen: i32) -> Option<(i32, Vec<u32>)> {
    // Align the reversed sequences (see the module doc).
    let qr: Vec<u8> = q.iter().rev().copied().collect();
    let tr: Vec<u8> = t.iter().rev().copied().collect();
    let (n, m) = (qr.len() as i32, tr.len() as i32); // v in 0..=n, h in 0..=m
    let k_end = m - n;
    let mut w = Wfs::default();
    // s = 0: M on diagonal 0 at offset 0, extended.
    let e0 = lcp(&qr, &tr) as i32;
    w.m.push(Some(Wf { lo: 0, hi: 0, off: vec![e0] }));
    for v in [&mut w.i1, &mut w.d1, &mut w.i2, &mut w.d2] { v.push(None); }
    if k_end == 0 && e0 >= m {
        return Some((0, backtrace(&w, p, 0, k_end, n, m, &qr, &tr)));
    }
    let mut s = 0;
    loop {
        s += 1;
        if s > max_pen { return None; }
        // Diagonal range from the source wavefronts.
        let srcs = [(&w.m, s - p.x), (&w.m, s - p.o1 - p.e1), (&w.i1, s - p.e1), (&w.d1, s - p.e1),
                    (&w.m, s - p.o2 - p.e2), (&w.i2, s - p.e2), (&w.d2, s - p.e2)];
        let (mut lo, mut hi, mut any) = (i32::MAX, i32::MIN, false);
        for (v, ss) in srcs {
            if let Some(wf) = at(v, ss) { lo = lo.min(wf.lo); hi = hi.max(wf.hi); any = true; }
        }
        if !any {
            for v in [&mut w.m, &mut w.i1, &mut w.d1, &mut w.i2, &mut w.d2] { v.push(None); }
            continue;
        }
        let (lo, hi) = (lo - 1, hi + 1);
        let len = (hi - lo + 1) as usize;
        let mut mo = vec![NONE; len];
        let mut i1o = vec![NONE; len];
        let mut d1o = vec![NONE; len];
        let mut i2o = vec![NONE; len];
        let mut d2o = vec![NONE; len];
        for k in lo..=hi {
            let idx = (k - lo) as usize;
            // Insertion (consumes the pattern/read: v + 1, h same): from diagonal k + 1.
            let a1 = get(&w.m, s - p.o1 - p.e1, k + 1).max(get(&w.i1, s - p.e1, k + 1));
            let a2 = get(&w.m, s - p.o2 - p.e2, k + 1).max(get(&w.i2, s - p.e2, k + 1));
            // Deletion (consumes the text/reference: h + 1): from diagonal k - 1.
            let b1 = get(&w.m, s - p.o1 - p.e1, k - 1).max(get(&w.d1, s - p.e1, k - 1));
            let b2 = get(&w.m, s - p.o2 - p.e2, k - 1).max(get(&w.d2, s - p.e2, k - 1));
            let ins1 = a1;
            let ins2 = a2;
            let del1 = if b1 > NONE { b1 + 1 } else { NONE };
            let del2 = if b2 > NONE { b2 + 1 } else { NONE };
            let mis = { let x = get(&w.m, s - p.x, k); if x > NONE { x + 1 } else { NONE } };
            // Clip to the matrix: v = h - k in 0..=n, h in 0..=m.
            let ok = |h: i32| h > NONE && h <= m && h - k >= 0 && h - k <= n;
            let ins1 = if ok(ins1) { ins1 } else { NONE };
            let ins2 = if ok(ins2) { ins2 } else { NONE };
            let del1 = if ok(del1) { del1 } else { NONE };
            let del2 = if ok(del2) { del2 } else { NONE };
            let mis = if ok(mis) { mis } else { NONE };
            i1o[idx] = ins1; i2o[idx] = ins2; d1o[idx] = del1; d2o[idx] = del2;
            let mut h = mis.max(ins1).max(ins2).max(del1).max(del2);
            if h > NONE {
                let v = h - k;
                h += lcp(&qr[v as usize..], &tr[h as usize..]) as i32;
            }
            mo[idx] = h;
        }
        let done = k_end >= lo && k_end <= hi && mo[(k_end - lo) as usize] >= m;
        w.m.push(Some(Wf { lo, hi, off: mo }));
        w.i1.push(Some(Wf { lo, hi, off: i1o }));
        w.d1.push(Some(Wf { lo, hi, off: d1o }));
        w.i2.push(Some(Wf { lo, hi, off: i2o }));
        w.d2.push(Some(Wf { lo, hi, off: d2o }));
        if done {
            return Some((s, backtrace(&w, p, s, k_end, n, m, &qr, &tr)));
        }
    }
}

/// Walk back from (s, k_end) at the end of M, building the CIGAR of the reversed alignment; it
/// is emitted forward-oriented for the original sequences.
#[allow(clippy::too_many_arguments)]
fn backtrace(w: &Wfs, p: &Penalties, s_end: i32, k_end: i32, n: i32, m: i32, qr: &[u8], tr: &[u8]) -> Vec<u32> {
    let _ = (qr, tr);
    // ops in backtrace order (end of reversed alignment first) = start of the forward alignment
    // first, since the reversal flips the order back.
    let mut ops: Vec<u32> = Vec::new();
    let push = |ops: &mut Vec<u32>, op: u32, len: u32| {
        if len == 0 { return; }
        if let Some(last) = ops.last_mut() && (*last & 0xf) == op { *last += len << 4; return; }
        ops.push(len << 4 | op);
    };
    #[derive(Clone, Copy, PartialEq)]
    enum C { M, I1, D1, I2, D2 }
    let (mut s, mut k, mut c) = (s_end, k_end, C::M);
    let mut h = m;
    loop {
        match c {
            C::M => {
                // Matches back to where this M value was set.
                let from_mis = { let x = get(&w.m, s - p.x, k); if x > NONE && x + 1 <= m && x + 1 - k <= n { x + 1 } else { NONE } };
                let i1 = get(&w.i1, s, k); let i2 = get(&w.i2, s, k);
                let d1 = get(&w.d1, s, k); let d2 = get(&w.d2, s, k);
                let src = from_mis.max(i1).max(i2).max(d1).max(d2);
                if s == 0 || src <= NONE {
                    // s = 0 start: matches from the origin.
                    push(&mut ops, 0, h as u32);
                    break;
                }
                push(&mut ops, 0, (h - src) as u32);
                h = src;
                // Prefer a gap closing here over a mismatch (left-aligned gaps), then piece 1.
                if src == d1 { c = C::D1; }
                else if src == d2 { c = C::D2; }
                else if src == i1 { c = C::I1; }
                else if src == i2 { c = C::I2; }
                else { push(&mut ops, 0, 1); h -= 1; s -= p.x; }
            }
            C::I1 | C::I2 => {
                let (oe, e, wv) = if c == C::I1 { (p.o1 + p.e1, p.e1, &w.i1) } else { (p.o2 + p.e2, p.e2, &w.i2) };
                push(&mut ops, 1, 1);
                // came from diagonal k + 1 at the same h
                let ext = get(wv, s - e, k + 1);
                let open = get(&w.m, s - oe, k + 1);
                k += 1;
                if ext >= open && ext > NONE && ext == h { s -= e; } else { s -= oe; c = C::M; }
            }
            C::D1 | C::D2 => {
                let (oe, e, wv) = if c == C::D1 { (p.o1 + p.e1, p.e1, &w.d1) } else { (p.o2 + p.e2, p.e2, &w.d2) };
                push(&mut ops, 2, 1);
                let ext = get(wv, s - e, k - 1);
                let open = get(&w.m, s - oe, k - 1);
                k -= 1;
                h -= 1;
                if ext >= open && ext > NONE && ext == h { s -= e; } else { s -= oe; c = C::M; }
            }
        }
        if h <= 0 && k == 0 && c == C::M && s == 0 { break; }
    }
    // Backtrace order on the reversed sequences = forward order on the originals.
    ops
}


// ---- Arena implementation (the one used for speed) -------------------------------------------

const C_M: usize = 0;
const C_I1: usize = 1;
const C_D1: usize = 2;
const C_I2: usize = 3;
const C_D2: usize = 4;

/// Per score: shared diagonal range and, per component, its start in `data` (u32::MAX = none).
#[derive(Clone, Copy)]
struct Meta { lo: i32, hi: i32, base: [u32; 5] }

const NO_META: Meta = Meta { lo: 1, hi: 0, base: [u32::MAX; 5] };

/// Reusable buffers for `align_fast`.
#[derive(Default)]
pub struct WfaBuf { meta: Vec<Meta>, data: Vec<i32>, qr: Vec<u8>, tr: Vec<u8> }

#[inline(always)]
fn aget(meta: &[Meta], data: &[i32], s: i32, c: usize, k: i32) -> i32 {
    if s < 0 { return NONE; }
    let m = &meta[s as usize];
    let b = m.base[c];
    if b == u32::MAX || k < m.lo || k > m.hi { return NONE; }
    unsafe { *data.get_unchecked(b as usize + (k - m.lo) as usize) }
}

fn gcd(a: i32, b: i32) -> i32 { if b == 0 { a } else { gcd(b, a % b) } }

/// As `align`, with flat reusable storage and penalties divided by their gcd.
pub fn align_fast(q: &[u8], t: &[u8], p0: &Penalties, max_pen: i32, buf: &mut WfaBuf) -> Option<(i32, Vec<u32>)> {
    let g = [p0.x, p0.o1, p0.e1, p0.o2, p0.e2].into_iter().fold(0, gcd).max(1);
    let p = Penalties { x: p0.x / g, o1: p0.o1 / g, e1: p0.e1 / g, o2: p0.o2 / g, e2: p0.e2 / g };
    let max_s = max_pen / g;
    let (oe1, oe2) = (p.o1 + p.e1, p.o2 + p.e2);
    buf.qr.clear(); buf.qr.extend(q.iter().rev());
    buf.tr.clear(); buf.tr.extend(t.iter().rev());
    let (qr, tr) = (&buf.qr[..], &buf.tr[..]);
    let (n, m) = (qr.len() as i32, tr.len() as i32);
    let k_end = m - n;
    buf.meta.clear();
    buf.data.clear();
    let e0 = lcp(qr, tr) as i32;
    buf.meta.push(Meta { lo: 0, hi: 0, base: [0, u32::MAX, u32::MAX, u32::MAX, u32::MAX] });
    buf.data.push(e0);
    let mut s_done = if k_end == 0 && e0 >= m { Some(0) } else { None };
    let mut s = 0;
    while s_done.is_none() {
        s += 1;
        if s > max_s { return None; }
        let (meta, data) = (&buf.meta, &buf.data);
        let has = |ss: i32, c: usize| ss >= 0 && meta[ss as usize].base[c] != u32::MAX;
        let (sx, so1, se1, so2, se2) = (s - p.x, s - oe1, s - p.e1, s - oe2, s - p.e2);
        let need_i1 = has(so1, C_M) || has(se1, C_I1);
        let need_d1 = has(so1, C_M) || has(se1, C_D1);
        let need_i2 = has(so2, C_M) || has(se2, C_I2);
        let need_d2 = has(so2, C_M) || has(se2, C_D2);
        let need_x = has(sx, C_M);
        if !(need_i1 || need_d1 || need_i2 || need_d2 || need_x) {
            buf.meta.push(NO_META);
            continue;
        }
        let (mut lo, mut hi) = (i32::MAX, i32::MIN);
        for (ss, c) in [(sx, C_M), (so1, C_M), (se1, C_I1), (se1, C_D1), (so2, C_M), (se2, C_I2), (se2, C_D2)] {
            if has(ss, c) { let mm = &meta[ss as usize]; lo = lo.min(mm.lo); hi = hi.max(mm.hi); }
        }
        let (lo, hi) = ((lo - 1).max(-n), (hi + 1).min(m));
        let len = (hi - lo + 1).max(0) as usize;
        let start = data.len() as u32;
        let mut base = [u32::MAX; 5];
        let mut nb = 0u32;
        for (c, need) in [(C_M, true), (C_I1, need_i1), (C_D1, need_d1), (C_I2, need_i2), (C_D2, need_d2)] {
            if need { base[c] = start + nb * len as u32; nb += 1; }
        }
        let total = start as usize + nb as usize * len;
        buf.data.resize(total, NONE);
        let (meta, data) = (&buf.meta[..], &mut buf.data[..]);
        let mut done = false;
        for k in lo..=hi {
            let idx = (k - lo) as usize;
            let v_ok = |h: i32| h > NONE && h <= m && h - k >= 0 && h - k <= n;
            let mut best = NONE;
            if need_i1 {
                let h = aget(meta, data, so1, C_M, k + 1).max(aget(meta, data, se1, C_I1, k + 1));
                let h = if v_ok(h) { h } else { NONE };
                data[base[C_I1] as usize + idx] = h; best = best.max(h);
            }
            if need_i2 {
                let h = aget(meta, data, so2, C_M, k + 1).max(aget(meta, data, se2, C_I2, k + 1));
                let h = if v_ok(h) { h } else { NONE };
                data[base[C_I2] as usize + idx] = h; best = best.max(h);
            }
            if need_d1 {
                let h0 = aget(meta, data, so1, C_M, k - 1).max(aget(meta, data, se1, C_D1, k - 1));
                let h = if h0 > NONE && v_ok(h0 + 1) { h0 + 1 } else { NONE };
                data[base[C_D1] as usize + idx] = h; best = best.max(h);
            }
            if need_d2 {
                let h0 = aget(meta, data, so2, C_M, k - 1).max(aget(meta, data, se2, C_D2, k - 1));
                let h = if h0 > NONE && v_ok(h0 + 1) { h0 + 1 } else { NONE };
                data[base[C_D2] as usize + idx] = h; best = best.max(h);
            }
            if need_x {
                let h0 = aget(meta, data, sx, C_M, k);
                if h0 > NONE && v_ok(h0 + 1) { best = best.max(h0 + 1); }
            }
            if best > NONE {
                let v = best - k;
                best += lcp(&qr[v as usize..], &tr[best as usize..]) as i32;
            }
            data[base[C_M] as usize + idx] = best;
            if k == k_end && best >= m { done = true; }
        }
        buf.meta.push(Meta { lo, hi, base });
        if done { s_done = Some(s); }
    }
    let s_end = s_done.unwrap();
    // Backtrace (as in `backtrace`, on the arena).
    let (meta, data) = (&buf.meta[..], &buf.data[..]);
    let get = |s: i32, c: usize, k: i32| aget(meta, data, s, c, k);
    let mut ops: Vec<u32> = Vec::new();
    let push = |ops: &mut Vec<u32>, op: u32, len: u32| {
        if len == 0 { return; }
        if let Some(last) = ops.last_mut() && (*last & 0xf) == op { *last += len << 4; return; }
        ops.push(len << 4 | op);
    };
    let (mut s, mut k, mut c, mut h) = (s_end, k_end, C_M, m);
    loop {
        if c == C_M {
            let fx = { let x = get(s - p.x, C_M, k); if x > NONE && x + 1 <= m && x + 1 - k <= n { x + 1 } else { NONE } };
            let (i1, i2, d1, d2) = (get(s, C_I1, k), get(s, C_I2, k), get(s, C_D1, k), get(s, C_D2, k));
            let src = fx.max(i1).max(i2).max(d1).max(d2);
            if s == 0 || src <= NONE { push(&mut ops, 0, h as u32); break; }
            push(&mut ops, 0, (h - src) as u32);
            h = src;
            if src == d1 { c = C_D1; } else if src == d2 { c = C_D2; }
            else if src == i1 { c = C_I1; } else if src == i2 { c = C_I2; }
            else { push(&mut ops, 0, 1); h -= 1; s -= p.x; }
        } else if c == C_I1 || c == C_I2 {
            let (oe, e) = if c == C_I1 { (oe1, p.e1) } else { (oe2, p.e2) };
            push(&mut ops, 1, 1);
            let ext = get(s - e, c, k + 1);
            k += 1;
            if ext > NONE && ext == h { s -= e; } else { s -= oe; c = C_M; }
        } else {
            let (oe, e) = if c == C_D1 { (oe1, p.e1) } else { (oe2, p.e2) };
            push(&mut ops, 2, 1);
            let ext = get(s - e, c, k - 1);
            k -= 1;
            h -= 1;
            if ext > NONE && ext == h { s -= e; } else { s -= oe; c = C_M; }
        }
    }
    Some((s_end * g, ops))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dp_score(q: &[u8], t: &[u8], a: i32, b: i32, go: i32, ge: i32, go2: i32, ge2: i32) -> i32 {
        // Plain O(nm) dual-affine global DP (max score) as the reference.
        let (n, m) = (q.len(), t.len());
        let ninf = i32::MIN / 4;
        let mut h = vec![vec![ninf; m + 1]; n + 1];
        let mut e1 = h.clone(); let mut f1 = h.clone(); let mut e2 = h.clone(); let mut f2 = h.clone();
        h[0][0] = 0;
        for i in 0..=n {
            for j in 0..=m {
                if i == 0 && j == 0 { continue; }
                if j > 0 {
                    e1[i][j] = (h[i][j - 1] - go - ge).max(e1[i][j - 1] - ge);
                    e2[i][j] = (h[i][j - 1] - go2 - ge2).max(e2[i][j - 1] - ge2);
                }
                if i > 0 {
                    f1[i][j] = (h[i - 1][j] - go - ge).max(f1[i - 1][j] - ge);
                    f2[i][j] = (h[i - 1][j] - go2 - ge2).max(f2[i - 1][j] - ge2);
                }
                let mut best = e1[i][j].max(e2[i][j]).max(f1[i][j]).max(f2[i][j]);
                if i > 0 && j > 0 { best = best.max(h[i - 1][j - 1] + if q[i - 1] == t[j - 1] { a } else { -b }); }
                h[i][j] = best;
            }
        }
        h[n][m]
    }

    fn cigar_score(c: &[u32], q: &[u8], t: &[u8], a: i32, b: i32, go: i32, ge: i32, go2: i32, ge2: i32) -> i32 {
        let (mut i, mut j, mut sc) = (0usize, 0usize, 0i32);
        for &op in c {
            let l = (op >> 4) as usize;
            match op & 0xf {
                0 => { for _ in 0..l { sc += if q[i] == t[j] { a } else { -b }; i += 1; j += 1; } }
                1 => { sc -= (go + ge * l as i32).min(go2 + ge2 * l as i32); i += l; }
                _ => { sc -= (go + ge * l as i32).min(go2 + ge2 * l as i32); j += l; }
            }
        }
        assert_eq!((i, j), (q.len(), t.len()), "cigar does not span the sequences");
        sc
    }

    #[test]
    fn wfa_matches_dp_score() {
        let (a, b, go, ge, go2, ge2) = (2, 4, 4, 2, 24, 1);
        let p = Penalties::from_scores(a, b, go, ge, go2, ge2);
        let mut x = 0x1234_5678_9abc_def1u64;
        let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        for case in 0..3000 {
            let len = 1 + (rnd() % 300) as usize;
            let t: Vec<u8> = (0..len).map(|_| (rnd() & 3) as u8).collect();
            let rate = [1u64, 3, 6, 12, 25][case % 5];
            let mut q = Vec::new();
            let mut i = 0;
            while i < t.len() {
                let r = rnd() % 100;
                if r < rate { q.push((rnd() & 3) as u8); i += 1; }
                else if r < rate + rate / 2 { i += 1 + (rnd() % 3) as usize; }
                else if r < 2 * rate { q.push((rnd() & 3) as u8); }
                else if r == 99 && case % 7 == 0 { i += 40; } // occasional long deletion
                else { q.push(t[i]); i += 1; }
            }
            if q.is_empty() { continue; }
            let want = dp_score(&q, &t, a, b, go, ge, go2, ge2);
            let (pen, cig) = align(&q, &t, &p, i32::MAX / 4).expect("no limit");
            let got = (a * (q.len() + t.len()) as i32 - pen) / 2;
            assert_eq!(got, want, "case {case}: penalty-derived score");
            assert_eq!(cigar_score(&cig, &q, &t, a, b, go, ge, go2, ge2), want, "case {case}: cigar score");
            let mut buf = WfaBuf::default();
            let (pen2, cig2) = align_fast(&q, &t, &p, i32::MAX / 4, &mut buf).expect("no limit");
            assert_eq!(pen2, pen, "case {case}: arena penalty");
            assert_eq!(cigar_score(&cig2, &q, &t, a, b, go, ge, go2, ge2), want, "case {case}: arena cigar score");
        }
    }
}

#[cfg(test)]
mod vs_ksw2 {
    use super::*;
    use crate::align::dp::{self, DpResult};

    /// Compare against the gap-fill DP (clean kernel, band |D| + 32 widened to cover the path):
    /// score agreement, CIGAR agreement and speed. `cargo test --release -- --ignored --nocapture vs_ksw2`
    #[test]
    #[ignore]
    fn wfa_vs_gap_dp() {
        let (a, b, go, ge, go2, ge2) = (2i8, 4i8, 4i8, 2i8, 24i8, 1i8);
        let p = Penalties::from_scores(a as i32, b as i32, go as i32, ge as i32, go2 as i32, ge2 as i32);
        let mat = { let mut m = [0i8; 25]; for i in 0..5 { for j in 0..5 { m[i * 5 + j] = if i == 4 || j == 4 { -1 } else if i == j { a } else { -b }; } } m };
        dp::set_clean_dual(true);
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let mut cases = Vec::new();
        for c in 0..4000 {
            let len = 100 + (rnd() % 700) as usize;
            let t: Vec<u8> = (0..len).map(|_| (rnd() & 3) as u8).collect();
            let rate = [1u64, 2, 4, 6, 10][c % 5];
            let mut q = Vec::new();
            let mut i = 0;
            while i < t.len() {
                let r = rnd() % 200;
                if r < rate { q.push((rnd() & 3) as u8); i += 1; }
                else if r < rate * 2 { i += 1; }
                else if r < rate * 3 { q.push((rnd() & 3) as u8); }
                else { q.push(t[i]); i += 1; }
            }
            cases.push((q, t));
        }
        let (mut same_score, mut same_cigar, mut n, mut wfa_none) = (0, 0, 0, 0);
        let t0 = std::time::Instant::now();
        let mut dp_res = Vec::new();
        for (q, t) in &cases {
            let nb = (t.len() as i32 - q.len() as i32).abs() + 32;
            let mut ez = DpResult::default();
            dp::extend_dual_affine_dispatch(q, t, 5, &mat, go, ge, go2, ge2, nb, 400, -1, dp::APPROX_MAX, &mut ez);
            dp_res.push(ez);
        }
        let t_dp = t0.elapsed();
        let t1 = std::time::Instant::now();
        let mut wfa_res = Vec::new();
        let mut buf = WfaBuf::default();
        for (q, t) in &cases { wfa_res.push(align_fast(q, t, &p, 1000, &mut buf)); }
        let t_wfa = t1.elapsed();
        for ((q, t), (ez, w)) in cases.iter().zip(dp_res.iter().zip(wfa_res.iter())) {
            n += 1;
            let Some((pen, cig)) = w else { wfa_none += 1; continue; };
            let s = (2 * (q.len() + t.len()) as i32 - pen) / 2;
            if s == ez.score { same_score += 1; }
            if *cig == ez.cigar { same_cigar += 1; }
        }
        eprintln!("{n} gaps: dp {:.1} us/gap, wfa {:.1} us/gap; same score {same_score}, same cigar {same_cigar}, wfa gave up {wfa_none}",
                  t_dp.as_secs_f64() * 1e6 / n as f64, t_wfa.as_secs_f64() * 1e6 / n as f64);
        // Per divergence class (rate index = case % 5), WFA capped at penalty 100 vs the DP.
        for cls in 0..5 {
            let sel: Vec<usize> = (0..cases.len()).filter(|i| i % 5 == cls).collect();
            let t2 = std::time::Instant::now();
            let mut solved = 0;
            for &i in &sel { if align_fast(&cases[i].0, &cases[i].1, &p, 100, &mut buf).is_some() { solved += 1; } }
            let tw = t2.elapsed().as_secs_f64() * 1e6 / sel.len() as f64;
            let t3 = std::time::Instant::now();
            for &i in &sel {
                let (q, t) = &cases[i];
                let nb = (t.len() as i32 - q.len() as i32).abs() + 32;
                let mut ez = DpResult::default();
                dp::extend_dual_affine_dispatch(q, t, 5, &mat, go, ge, go2, ge2, nb, 400, -1, dp::APPROX_MAX, &mut ez);
            }
            let td = t3.elapsed().as_secs_f64() * 1e6 / sel.len() as f64;
            eprintln!("  class {cls}: wfa(cap 100) {tw:.1} us/gap, solved {solved}/{}, dp {td:.1} us/gap", sel.len());
        }
    }
}
