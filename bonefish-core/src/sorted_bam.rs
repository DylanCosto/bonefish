//! Coordinate-sorted, indexed BAM output (`bonefish --sort -o out.bam`), in place of
//! `bonefish -a | samtools sort --write-index`.
//!
//! The mapper's SAM text is parsed into BAM records as it is written (one writer thread; parsing
//! is cheap next to compression). Records collect in a memory batch. A full batch is sorted and
//! handed to a background thread that writes it as a temporary run, in chunks compressed with
//! libdeflate level 1 (a private format, much cheaper than the BGZF level 6 that samtools uses for
//! its temporary files), while the next batch fills. After the last read, the runs and the last
//! batch are merged and written as BGZF blocks compressed in parallel, and the BAI index is built
//! from the block offsets as they are written.
//!
//! Order is samtools sort's coordinate order: reference id (unmapped records last), position, then
//! reverse strand after forward, ties in input order. Records are encoded as htslib's `sam_parse1`
//! does (smallest integer type for `i` tags, CIGARs over 65535 operations in a `CG` tag), so
//! `samtools view` of the result equals that of `bonefish -a | samtools sort` record for record.

use std::collections::{BTreeMap, BinaryHeap, HashMap};
use std::cmp::Reverse;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread::JoinHandle;

/// Uncompressed size limit of a BGZF block (htslib's BGZF_BLOCK_SIZE).
const BLOCK: usize = 0xff00;
/// Raw size of a compressed chunk in a temporary run.
const RUN_CHUNK: usize = 4 << 20;
/// BGZF blocks per compression job in the final merge.
const BLOCKS_PER_JOB: usize = 64;

pub struct SortConfig {
    /// Memory for record batches (two are alive at a time: one filling, one being spilled).
    pub mem: usize,
    /// Temporary run files are `<tmp_prefix>.NNNN.run`.
    pub tmp_prefix: String,
    /// Threads for compression.
    pub threads: usize,
    /// BGZF compression level of the output (htslib's default is 6).
    pub level: i32,
}

/// One batch of encoded records: `data` holds them back to back (each with its 4-byte
/// block_size), `ents` the sort key and offset of each, in input order.
#[derive(Default)]
struct Batch {
    data: Vec<u8>,
    ents: Vec<(u64, u64)>,
}

impl Batch {
    /// Stable sort by key (ties keep input order).
    fn sort(&mut self) {
        self.ents.sort_by_key(|e| e.0);
    }
    fn record(&self, off: u64) -> &[u8] {
        let o = off as usize;
        let n = u32::from_le_bytes(self.data[o..o + 4].try_into().unwrap()) as usize;
        &self.data[o..o + 4 + n]
    }
}

pub struct SortedBamWriter {
    cfg: SortConfig,
    out_path: String,
    partial: Vec<u8>,
    header: String,
    name2tid: HashMap<Vec<u8>, i32>,
    refs: Vec<(String, u32)>,
    header_done: bool,
    batch: Batch,
    spill_tx: Option<SyncSender<(Batch, PathBuf)>>,
    /// Batches back from the spill thread (to reuse their buffers); `false`: spill aborted.
    back_rx: Receiver<(Batch, bool)>,
    /// Set by `finish`: abandon a spill in progress, its batch is merged from memory instead.
    abort: std::sync::Arc<std::sync::atomic::AtomicBool>,
    spill_thread: Option<JoinHandle<io::Result<()>>>,
    runs: Vec<PathBuf>,
    scratch: Vec<u32>,
}

impl SortedBamWriter {
    pub fn new(out_path: &str, cfg: SortConfig) -> io::Result<Self> {
        // Fail early if the output can't be created.
        File::create(out_path)?;
        // Rendezvous: at most one batch filling and one being spilled.
        let (tx, rx) = mpsc::sync_channel::<(Batch, PathBuf)>(0);
        let (back_tx, back_rx) = mpsc::channel::<(Batch, bool)>();
        let abort = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let spill_threads = (cfg.threads / 4).max(1);
        let ab = abort.clone();
        let spill_thread = std::thread::Builder::new().name("sort-spill".into()).spawn(move || -> io::Result<()> {
            while let Ok((batch, path)) = rx.recv() {
                let done = write_run(&batch, &path, spill_threads, &ab)?;
                if !done { let _ = std::fs::remove_file(&path); }
                if back_tx.send((batch, done)).is_err() { break; }
            }
            Ok(())
        })?;
        Ok(SortedBamWriter {
            cfg, out_path: out_path.to_string(), partial: Vec::new(), header: String::new(),
            name2tid: HashMap::new(), refs: Vec::new(), header_done: false, batch: Batch::default(),
            spill_tx: Some(tx), back_rx, abort, spill_thread: Some(spill_thread), runs: Vec::new(), scratch: Vec::new(),
        })
    }

    fn line(&mut self, line: &[u8]) -> io::Result<()> {
        if line.is_empty() { return Ok(()); }
        if line[0] == b'@' && !self.header_done {
            let text = std::str::from_utf8(line).map_err(|_| bad("non-UTF-8 header line"))?;
            if text.starts_with("@HD") {
                // As samtools sort: SO:coordinate; the query grouping no longer holds.
                let mut f: Vec<String> = text.split('\t').filter(|t| !t.starts_with("SO:") && !t.starts_with("GO:")).map(String::from).collect();
                f.insert(2.min(f.len()), "SO:coordinate".into());
                self.header.push_str(&f.join("\t"));
            } else {
                if text.starts_with("@SQ") {
                    let (mut sn, mut ln) = (None, None);
                    for t in text.split('\t') {
                        if let Some(v) = t.strip_prefix("SN:") { sn = Some(v.to_string()); }
                        if let Some(v) = t.strip_prefix("LN:") { ln = v.parse::<u32>().ok(); }
                    }
                    let (sn, ln) = (sn.ok_or_else(|| bad("@SQ without SN"))?, ln.ok_or_else(|| bad("@SQ without LN"))?);
                    self.name2tid.insert(sn.as_bytes().to_vec(), self.refs.len() as i32);
                    self.refs.push((sn, ln));
                }
                self.header.push_str(text);
            }
            self.header.push('\n');
            return Ok(());
        }
        self.header_done = true;
        let off = self.batch.data.len() as u64;
        let key = encode_record(line, &self.name2tid, &mut self.batch.data, &mut self.scratch)?;
        self.batch.ents.push((key, off));
        if self.batch.data.len() >= self.cfg.mem / 2 {
            self.spill()?;
        }
        Ok(())
    }

    fn spill(&mut self) -> io::Result<()> {
        // Next batch: reuse a spilled batch's buffers if one is back (no 2 GB regrowth).
        let next = match self.back_rx.try_recv() {
            Ok((mut b, _)) => { b.data.clear(); b.ents.clear(); b }
            Err(_) => Batch { data: Vec::with_capacity(self.batch.data.capacity()), ents: Vec::with_capacity(self.batch.ents.capacity()) },
        };
        let mut b = std::mem::replace(&mut self.batch, next);
        b.sort();
        let path = PathBuf::from(format!("{}.{:04}.run", self.cfg.tmp_prefix, self.runs.len()));
        self.runs.push(path.clone());
        if self.spill_tx.as_ref().unwrap().send((b, path)).is_err() {
            // The spill thread stopped: report its error.
            return Err(self.join_spill().err().unwrap_or_else(|| bad("spill thread stopped")));
        }
        Ok(())
    }

    fn join_spill(&mut self) -> io::Result<()> {
        self.spill_tx = None;
        match self.spill_thread.take() {
            Some(h) => h.join().map_err(|_| bad("spill thread panicked"))?,
            None => Ok(()),
        }
    }

    /// Merge everything into the sorted BAM and its index. Returns the number of records.
    pub fn finish(mut self) -> io::Result<u64> {
        if !self.partial.is_empty() {
            let p = std::mem::take(&mut self.partial);
            self.line(&p)?;
        }
        let t0 = std::time::Instant::now();
        // A spill still being written is abandoned: its batch is still in memory.
        self.abort.store(true, std::sync::atomic::Ordering::Relaxed);
        self.join_spill()?;
        let mut mem_batches: Vec<Batch> = Vec::new();
        while let Ok((b, done)) = self.back_rx.try_recv() {
            if !done {
                self.runs.pop(); // the abandoned run was the last one started
                mem_batches.push(b);
            }
        }
        let t1 = std::time::Instant::now();
        let mut last = std::mem::take(&mut self.batch);
        last.sort();
        mem_batches.push(last);
        let res = merge_and_write(&self.out_path, &self.header, &self.refs, &self.runs, mem_batches, &self.cfg);
        if std::env::var_os("BONEFISH_SORT_TIMING").is_some() {
            eprintln!("[sort] wait for last spill {:.2}s, merge+write {:.2}s, {} runs", (t1 - t0).as_secs_f64(), t1.elapsed().as_secs_f64(), self.runs.len());
        }
        for r in &self.runs { let _ = std::fs::remove_file(r); }
        res
    }
}

impl Drop for SortedBamWriter {
    fn drop(&mut self) {
        self.abort.store(true, std::sync::atomic::Ordering::Relaxed);
        let _ = self.join_spill();
        for r in &self.runs { let _ = std::fs::remove_file(r); }
    }
}

impl Write for SortedBamWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut start = 0;
        while let Some(nl) = memchr::memchr(b'\n', &buf[start..]) {
            let end = start + nl;
            if self.partial.is_empty() {
                self.line(&buf[start..end])?;
            } else {
                let mut p = std::mem::take(&mut self.partial);
                p.extend_from_slice(&buf[start..end]);
                self.line(&p)?;
            }
            start = end + 1;
        }
        self.partial.extend_from_slice(&buf[start..]);
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> { Ok(()) }
}

/// `Write` handle that forwards to a shared writer, so the caller can `finish` it afterwards.
pub struct SharedSink(pub std::sync::Arc<std::sync::Mutex<Option<SortedBamWriter>>>);

impl Write for SharedSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().as_mut().ok_or_else(|| bad("sorted output already finished"))?.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> { Ok(()) }
}

struct Guard<F: FnMut()>(F);
impl<F: FnMut()> Drop for Guard<F> { fn drop(&mut self) { (self.0)() } }
fn scopeguard<F: FnMut()>(f: F) -> Guard<F> { Guard(f) }

fn bad(msg: &str) -> io::Error { io::Error::new(io::ErrorKind::InvalidData, msg.to_string()) }

// ───────────────────────── SAM text → BAM record ─────────────────────────

const SEQ_NT16: [u8; 256] = {
    let mut t = [15u8; 256];
    let s = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;
    while i < 16 {
        t[s[i] as usize] = i as u8;
        t[s[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    t
};

fn cigar_code(c: u8) -> Option<u32> {
    Some(match c { b'M' => 0, b'I' => 1, b'D' => 2, b'N' => 3, b'S' => 4, b'H' => 5, b'P' => 6, b'=' => 7, b'X' => 8, _ => return None })
}

fn parse_i64(s: &[u8]) -> io::Result<i64> {
    std::str::from_utf8(s).ok().and_then(|t| t.parse::<i64>().ok()).ok_or_else(|| bad("bad integer field"))
}

/// htslib's hts_reg2bin(beg, end, 14, 5) for BAI.
pub fn reg2bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    if beg >> 14 == end >> 14 { return (((1 << 15) - 1) / 7 + (beg >> 14)) as u32; }
    if beg >> 17 == end >> 17 { return (((1 << 12) - 1) / 7 + (beg >> 17)) as u32; }
    if beg >> 20 == end >> 20 { return (((1 << 9) - 1) / 7 + (beg >> 20)) as u32; }
    if beg >> 23 == end >> 23 { return (((1 << 6) - 1) / 7 + (beg >> 23)) as u32; }
    if beg >> 26 == end >> 26 { return (((1 << 3) - 1) / 7 + (beg >> 26)) as u32; }
    0
}

/// Sort key in samtools' coordinate order (tid as unsigned, so unmapped -1 sorts last).
fn sort_key(tid: i32, pos: i32, flag: u16) -> u64 {
    ((tid as u32 as u64) << 32) | (((pos as i64 + 1) as u64 & 0x7fff_ffff) << 1) | ((flag >> 4) & 1) as u64
}

/// Append one SAM line as a BAM record (with block_size) to `out`; returns its sort key.
fn encode_record(line: &[u8], name2tid: &HashMap<Vec<u8>, i32>, out: &mut Vec<u8>, cig: &mut Vec<u32>) -> io::Result<u64> {
    let mut f = line.split(|&b| b == b'\t');
    let mut next = || f.next().ok_or_else(|| bad("SAM line with fewer than 11 fields"));
    let qname = next()?;
    let flag = parse_i64(next()?)? as u16;
    let rname = next()?;
    let pos = parse_i64(next()?)? - 1;
    let mapq = parse_i64(next()?)? as u8;
    let cigar = next()?;
    let rnext = next()?;
    let pnext = parse_i64(next()?)? - 1;
    let tlen = parse_i64(next()?)? as i32;
    let seq = next()?;
    let qual = next()?;
    let tid = |n: &[u8]| -> io::Result<i32> {
        if n == b"*" { Ok(-1) } else { name2tid.get(n).copied().ok_or_else(|| bad("reference not in header")) }
    };
    let ref_id = tid(rname)?;
    let mate_id = if rnext == b"=" { ref_id } else { tid(rnext)? };

    cig.clear();
    let mut rlen: i64 = 0;
    if cigar != b"*" {
        let mut n: u32 = 0;
        for &c in cigar {
            if c.is_ascii_digit() { n = n * 10 + (c - b'0') as u32; continue; }
            let op = cigar_code(c).ok_or_else(|| bad("bad CIGAR operation"))?;
            if matches!(op, 0 | 2 | 3 | 7 | 8) { rlen += n as i64; }
            cig.push(n << 4 | op);
            n = 0;
        }
    }
    let l_seq = if seq == b"*" { 0 } else { seq.len() };
    // htslib: bin from the alignment end (pos + 1 when unmapped or without reference span).
    let end = if flag & 4 == 0 && rlen > 0 { pos + rlen } else { pos + 1 };
    let bin = reg2bin(pos, end);
    let long_cigar = cig.len() > 0xffff;

    let start = out.len();
    out.extend_from_slice(&[0u8; 4]); // block_size, patched below
    out.extend_from_slice(&ref_id.to_le_bytes());
    out.extend_from_slice(&(pos as i32).to_le_bytes());
    out.push((qname.len() + 1) as u8);
    out.push(mapq);
    out.extend_from_slice(&(bin as u16).to_le_bytes());
    out.extend_from_slice(&(if long_cigar { 2 } else { cig.len() as u16 }).to_le_bytes());
    out.extend_from_slice(&flag.to_le_bytes());
    out.extend_from_slice(&(l_seq as i32).to_le_bytes());
    out.extend_from_slice(&mate_id.to_le_bytes());
    out.extend_from_slice(&(pnext as i32).to_le_bytes());
    out.extend_from_slice(&tlen.to_le_bytes());
    out.extend_from_slice(qname);
    out.push(0);
    if long_cigar {
        // As htslib's bam_write1: placeholder <l_seq>S<rlen>N, real CIGAR in a CG:B:I tag.
        out.extend_from_slice(&((l_seq as u32) << 4 | 4).to_le_bytes());
        out.extend_from_slice(&((rlen as u32) << 4 | 3).to_le_bytes());
    } else {
        for &c in cig.iter() { out.extend_from_slice(&c.to_le_bytes()); }
    }
    if l_seq > 0 {
        for pair in seq.chunks(2) {
            let hi = SEQ_NT16[pair[0] as usize];
            let lo = if pair.len() > 1 { SEQ_NT16[pair[1] as usize] } else { 0 };
            out.push(hi << 4 | lo);
        }
        if qual == b"*" {
            out.extend(std::iter::repeat_n(0xffu8, l_seq));
        } else {
            if qual.len() != l_seq { return Err(bad("SEQ and QUAL lengths differ")); }
            out.extend(qual.iter().map(|&q| q.wrapping_sub(33)));
        }
    }
    for tag in f {
        encode_aux(tag, out)?;
    }
    if long_cigar {
        out.extend_from_slice(b"CGBI");
        out.extend_from_slice(&(cig.len() as u32).to_le_bytes());
        for &c in cig.iter() { out.extend_from_slice(&c.to_le_bytes()); }
    }
    let block_size = (out.len() - start - 4) as u32;
    out[start..start + 4].copy_from_slice(&block_size.to_le_bytes());
    Ok(sort_key(ref_id, pos as i32, flag))
}

fn push_int(v: i64, out: &mut Vec<u8>) {
    // htslib sam_parse1: smallest type that holds the value.
    if v < 0 {
        if v >= i8::MIN as i64 { out.push(b'c'); out.push(v as i8 as u8); }
        else if v >= i16::MIN as i64 { out.push(b's'); out.extend_from_slice(&(v as i16).to_le_bytes()); }
        else { out.push(b'i'); out.extend_from_slice(&(v as i32).to_le_bytes()); }
    } else if v <= u8::MAX as i64 { out.push(b'C'); out.push(v as u8); }
    else if v <= u16::MAX as i64 { out.push(b'S'); out.extend_from_slice(&(v as u16).to_le_bytes()); }
    else { out.push(b'I'); out.extend_from_slice(&(v as u32).to_le_bytes()); }
}

fn parse_f32(s: &[u8]) -> io::Result<f32> {
    // htslib: (float)strtod(...), i.e. through double.
    std::str::from_utf8(s).ok().and_then(|t| t.parse::<f64>().ok()).map(|x| x as f32).ok_or_else(|| bad("bad float"))
}

fn encode_aux(tag: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    if tag.len() < 5 || tag[2] != b':' || tag[4] != b':' { return Err(bad("bad optional field")); }
    out.extend_from_slice(&tag[..2]);
    let val = &tag[5..];
    match tag[3] {
        b'A' => { out.push(b'A'); out.push(*val.first().ok_or_else(|| bad("empty A field"))?); }
        b'i' => push_int(parse_i64(val)?, out),
        b'f' => { out.push(b'f'); out.extend_from_slice(&parse_f32(val)?.to_le_bytes()); }
        b'Z' | b'H' => { out.push(tag[3]); out.extend_from_slice(val); out.push(0); }
        b'B' => {
            let sub = *val.first().ok_or_else(|| bad("empty B field"))?;
            out.push(b'B');
            out.push(sub);
            let cnt_at = out.len();
            out.extend_from_slice(&[0u8; 4]);
            let mut n = 0u32;
            if val.len() > 1 {
                for v in val[2..].split(|&b| b == b',') {
                    match sub {
                        b'c' => out.push(parse_i64(v)? as i8 as u8),
                        b'C' => out.push(parse_i64(v)? as u8),
                        b's' => out.extend_from_slice(&(parse_i64(v)? as i16).to_le_bytes()),
                        b'S' => out.extend_from_slice(&(parse_i64(v)? as u16).to_le_bytes()),
                        b'i' => out.extend_from_slice(&(parse_i64(v)? as i32).to_le_bytes()),
                        b'I' => out.extend_from_slice(&(parse_i64(v)? as u32).to_le_bytes()),
                        b'f' => out.extend_from_slice(&parse_f32(v)?.to_le_bytes()),
                        _ => return Err(bad("bad B array type")),
                    }
                    n += 1;
                }
            }
            out[cnt_at..cnt_at + 4].copy_from_slice(&n.to_le_bytes());
        }
        _ => return Err(bad("unknown optional field type")),
    }
    Ok(())
}

// ───────────────────────── temporary runs ─────────────────────────

/// Write a sorted batch as a run: chunks of whole records, each `[raw_len u32][comp_len u32][deflate]`.
/// Returns false (file incomplete) if `abort` was set meanwhile.
fn write_run(b: &Batch, path: &PathBuf, threads: usize, abort: &std::sync::atomic::AtomicBool) -> io::Result<bool> {
    let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let mut group: Vec<Vec<u8>> = Vec::new();
    let mut cur = Vec::with_capacity(RUN_CHUNK + 65536);
    let flush = |group: &mut Vec<Vec<u8>>, w: &mut BufWriter<File>| -> io::Result<bool> {
        if abort.load(std::sync::atomic::Ordering::Relaxed) { return Ok(false); }
        let comp = compress_all(group, 1, threads);
        for (raw, c) in group.iter().zip(&comp) {
            w.write_all(&(raw.len() as u32).to_le_bytes())?;
            w.write_all(&(c.len() as u32).to_le_bytes())?;
            w.write_all(c)?;
        }
        group.clear();
        Ok(true)
    };
    for &(_, off) in &b.ents {
        cur.extend_from_slice(b.record(off));
        if cur.len() >= RUN_CHUNK {
            group.push(std::mem::replace(&mut cur, Vec::with_capacity(RUN_CHUNK + 65536)));
            if group.len() >= 4 * threads && !flush(&mut group, &mut w)? { return Ok(false); }
        }
    }
    if !cur.is_empty() { group.push(cur); }
    if !group.is_empty() && !flush(&mut group, &mut w)? { return Ok(false); }
    w.flush()?;
    Ok(true)
}

/// Deflate every buffer (level `lvl`) on up to `threads` threads, in order.
fn compress_all(bufs: &[Vec<u8>], lvl: i32, threads: usize) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new(); bufs.len()];
    let per = bufs.len().div_ceil(threads.max(1)).max(1);
    std::thread::scope(|s| {
        for (bs, os) in bufs.chunks(per).zip(out.chunks_mut(per)) {
            s.spawn(move || {
                let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(lvl).unwrap());
                for (b, o) in bs.iter().zip(os.iter_mut()) {
                    let mut v = vec![0u8; c.deflate_compress_bound(b.len())];
                    let n = c.deflate_compress(b, &mut v).expect("deflate");
                    v.truncate(n);
                    *o = v;
                }
            });
        }
    });
    out
}

static RUN_WAIT_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static RUN_READ_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static RUN_DEC_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Reads a run back, decompressing chunks on a prefetch thread.
struct RunReader {
    rx: Receiver<Receiver<io::Result<Vec<u8>>>>,
    buf: Vec<u8>,
    /// Current record: buf[cur.0..cur.1].
    cur: (usize, usize),
}

impl RunReader {
    /// `dec` threads decompress chunks (in order) ahead of the merge.
    fn open(path: &PathBuf, dec: usize) -> io::Result<RunReader> {
        let mut r = BufReader::with_capacity(1 << 20, File::open(path)?);
        let dec = dec.max(1);
        // In-order futures: the feeder sends each chunk's result receiver, decoders fill them.
        let (tx, rx) = mpsc::sync_channel::<Receiver<io::Result<Vec<u8>>>>(2 * dec);
        let (job_tx, job_rx) = mpsc::sync_channel::<(Vec<u8>, usize, SyncSender<io::Result<Vec<u8>>>)>(dec);
        let job_rx = std::sync::Arc::new(std::sync::Mutex::new(job_rx));
        for _ in 0..dec {
            let jr = job_rx.clone();
            std::thread::Builder::new().name("sort-dec".into()).spawn(move || {
                let mut d = libdeflater::Decompressor::new();
                loop {
                    let job = jr.lock().unwrap().recv();
                    let Ok((c, raw, done)) = job else { break };
                    let td = std::time::Instant::now();
                    let mut out = vec![0u8; raw];
                    let res = d.deflate_decompress(&c, &mut out).map(|_| out).map_err(|_| bad("corrupt temporary run"));
                    RUN_DEC_NS.fetch_add(td.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
                    let _ = done.send(res);
                }
            })?;
        }
        std::thread::Builder::new().name("sort-run".into()).spawn(move || {
            let mut hdr = [0u8; 8];
            loop {
                let tr = std::time::Instant::now();
                match r.read_exact(&mut hdr) {
                    Ok(()) => {}
                    Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                    Err(e) => { let (t, f) = mpsc::sync_channel(1); let _ = t.send(Err(e)); let _ = tx.send(f); break; }
                }
                let raw = u32::from_le_bytes(hdr[..4].try_into().unwrap()) as usize;
                let cl = u32::from_le_bytes(hdr[4..].try_into().unwrap()) as usize;
                let mut c = vec![0u8; cl];
                if let Err(e) = r.read_exact(&mut c) { let (t, f) = mpsc::sync_channel(1); let _ = t.send(Err(e)); let _ = tx.send(f); break; }
                RUN_READ_NS.fetch_add(tr.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
                let (t, f) = mpsc::sync_channel(1);
                if tx.send(f).is_err() || job_tx.send((c, raw, t)).is_err() { break; }
            }
        })?;
        Ok(RunReader { rx, buf: Vec::new(), cur: (0, 0) })
    }

    /// Move to the next record; returns its sort key, or None at the end.
    fn advance(&mut self) -> io::Result<Option<u64>> {
        if self.cur.1 >= self.buf.len() {
            let tw = std::time::Instant::now();
            let got = self.rx.recv();
            let chunk = match got {
                Ok(fut) => fut.recv().map_err(|_| bad("run decoder stopped"))?,
                Err(_) => return Ok(None),
            };
            RUN_WAIT_NS.fetch_add(tw.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            self.buf = chunk?;
            self.cur = (0, 0);
        }
        let o = self.cur.1;
        let n = u32::from_le_bytes(self.buf[o..o + 4].try_into().unwrap()) as usize;
        self.cur = (o, o + 4 + n);
        Ok(Some(record_key(self.current())))
    }

    /// The current record (with block_size).
    fn current(&self) -> &[u8] {
        &self.buf[self.cur.0..self.cur.1]
    }
}

/// Sort key of an encoded record.
fn record_key(r: &[u8]) -> u64 {
    let tid = i32::from_le_bytes(r[4..8].try_into().unwrap());
    let pos = i32::from_le_bytes(r[8..12].try_into().unwrap());
    let flag = u16::from_le_bytes(r[18..20].try_into().unwrap());
    sort_key(tid, pos, flag)
}

// ───────────────────────── merge, BGZF output, BAI ─────────────────────────

/// Index fields of one record, with its position in the uncompressed stream:
/// (block number, offset in block) at its start and end.
struct RecMeta { tid: i32, beg: i64, end: i64, bin: u32, mapped: bool, sb: u64, so: u32, eb: u64, eo: u32 }

fn rec_meta(r: &[u8], sb: u64, so: u32, eb: u64, eo: u32) -> RecMeta {
    let tid = i32::from_le_bytes(r[4..8].try_into().unwrap());
    let pos = i32::from_le_bytes(r[8..12].try_into().unwrap()) as i64;
    let l_name = r[12] as usize;
    let n_cig = u16::from_le_bytes(r[16..18].try_into().unwrap()) as usize;
    let flag = u16::from_le_bytes(r[18..20].try_into().unwrap());
    let bin = u16::from_le_bytes(r[14..16].try_into().unwrap()) as u32;
    let mapped = flag & 4 == 0;
    let mut rlen: i64 = 0;
    let cig_at = 36 + l_name;
    let ops: Vec<u32> = (0..n_cig).map(|i| u32::from_le_bytes(r[cig_at + 4 * i..cig_at + 4 * i + 4].try_into().unwrap())).collect();
    // A long CIGAR is stored as <l_seq>S<rlen>N with the real one in CG: the N carries rlen.
    for &c in &ops {
        if matches!(c & 15, 0 | 2 | 3 | 7 | 8) { rlen += (c >> 4) as i64; }
    }
    let end = if mapped && rlen > 0 { pos + rlen } else { pos + 1 };
    RecMeta { tid, beg: pos, end, bin, mapped, sb, so, eb, eo }
}

/// One BGZF block from raw data.
fn bgzf_block(raw: &[u8], deflated: &[u8]) -> Vec<u8> {
    let total = 18 + deflated.len() + 8;
    let mut b = Vec::with_capacity(total);
    b.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0]);
    b.extend_from_slice(&((total - 1) as u16).to_le_bytes());
    b.extend_from_slice(deflated);
    let mut crc = libdeflater::Crc::new();
    crc.update(raw);
    b.extend_from_slice(&crc.sum().to_le_bytes());
    b.extend_from_slice(&(raw.len() as u32).to_le_bytes());
    b
}

const BGZF_EOF: [u8; 28] = [0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0];

/// Per-reference BAI state, filled as records are written.
#[derive(Default)]
struct RefIndex {
    bins: HashMap<u32, Vec<(u64, u64)>>,
    linear: Vec<u64>,
    off_beg: u64,
    off_end: u64,
    n_mapped: u64,
    n_unmapped: u64,
    seen: bool,
}

struct Indexer { refs: Vec<RefIndex>, n_no_coor: u64 }

impl Indexer {
    fn push(&mut self, m: &RecMeta, vbeg: u64, vend: u64) {
        if m.tid < 0 { self.n_no_coor += 1; return; }
        let r = &mut self.refs[m.tid as usize];
        if !r.seen { r.seen = true; r.off_beg = vbeg; }
        r.off_end = vend;
        if m.mapped { r.n_mapped += 1; } else { r.n_unmapped += 1; }
        let chunks = r.bins.entry(m.bin).or_default();
        match chunks.last_mut() {
            Some(c) if c.1 == vbeg => c.1 = vend,
            _ => chunks.push((vbeg, vend)),
        }
        if m.mapped {
            let (b, e) = (m.beg.max(0), m.end.max(1));
            let (w0, w1) = ((b >> 14) as usize, ((e - 1) >> 14) as usize);
            if r.linear.len() <= w1 { r.linear.resize(w1 + 1, u64::MAX); }
            for w in w0..=w1 {
                if r.linear[w] == u64::MAX { r.linear[w] = vbeg; }
            }
        }
    }

    fn write(&mut self, path: &str) -> io::Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(b"BAI\x01")?;
        w.write_all(&(self.refs.len() as i32).to_le_bytes())?;
        for r in &mut self.refs {
            if !r.seen {
                w.write_all(&0i32.to_le_bytes())?; // n_bin
                w.write_all(&0i32.to_le_bytes())?; // n_intv
                continue;
            }
            let mut bins: Vec<(&u32, &Vec<(u64, u64)>)> = r.bins.iter().collect();
            bins.sort_by_key(|b| *b.0);
            w.write_all(&((bins.len() + 1) as i32).to_le_bytes())?;
            for (bin, chunks) in bins {
                w.write_all(&bin.to_le_bytes())?;
                w.write_all(&(chunks.len() as i32).to_le_bytes())?;
                for &(a, b) in chunks { w.write_all(&a.to_le_bytes())?; w.write_all(&b.to_le_bytes())?; }
            }
            // Pseudo-bin 37450: offsets of the reference's records and its mapped/unmapped counts.
            w.write_all(&37450u32.to_le_bytes())?;
            w.write_all(&2i32.to_le_bytes())?;
            for v in [r.off_beg, r.off_end, r.n_mapped, r.n_unmapped] { w.write_all(&v.to_le_bytes())?; }
            // Empty windows get the previous window's offset (the first: the reference's first record).
            let mut prev = r.off_beg;
            for v in r.linear.iter_mut() { if *v == u64::MAX { *v = prev; } else { prev = *v; } }
            w.write_all(&(r.linear.len() as i32).to_le_bytes())?;
            for v in &r.linear { w.write_all(&v.to_le_bytes())?; }
        }
        w.write_all(&self.n_no_coor.to_le_bytes())?;
        w.flush()
    }
}

/// A compression job: raw BGZF blocks (numbered from `first_block`) and the records they hold.
struct Job { seq: u64, blocks: Vec<Vec<u8>>, metas: Vec<RecMeta> }

fn merge_and_write(out_path: &str, header: &str, refs: &[(String, u32)], runs: &[PathBuf], mem: Vec<Batch>, cfg: &SortConfig) -> io::Result<u64> {
    let threads = cfg.threads.max(1);
    let level = cfg.level;
    let (job_tx, job_rx) = mpsc::sync_channel::<Job>(threads * 2);
    let job_rx = std::sync::Arc::new(std::sync::Mutex::new(job_rx));
    let (done_tx, done_rx) = mpsc::sync_channel::<(u64, Vec<Vec<u8>>, Vec<RecMeta>)>(threads * 2);

    // Compressors.
    let mut workers = Vec::new();
    for _ in 0..threads {
        let rx = job_rx.clone();
        let tx = done_tx.clone();
        workers.push(std::thread::spawn(move || {
            let mut c = libdeflater::Compressor::new(libdeflater::CompressionLvl::new(level).unwrap());
            let mut tmp = Vec::new();
            loop {
                let job = match rx.lock().unwrap().recv() { Ok(j) => j, Err(_) => break };
                let out: Vec<Vec<u8>> = job.blocks.iter().map(|raw| {
                    tmp.resize(c.deflate_compress_bound(raw.len()), 0);
                    let n = c.deflate_compress(raw, &mut tmp).expect("deflate");
                    bgzf_block(raw, &tmp[..n])
                }).collect();
                if tx.send((job.seq, out, job.metas)).is_err() { break; }
            }
        }));
    }
    drop(done_tx);

    // Writer: blocks in order, block file offsets → virtual offsets → index.
    let out_path_w = out_path.to_string();
    let n_refs = refs.len();
    let writer = std::thread::spawn(move || -> io::Result<()> {
        let mut w = BufWriter::with_capacity(4 << 20, File::create(&out_path_w)?);
        let mut idx = Indexer { refs: (0..n_refs).map(|_| RefIndex::default()).collect(), n_no_coor: 0 };
        let mut pending: BTreeMap<u64, (Vec<Vec<u8>>, Vec<RecMeta>)> = BTreeMap::new();
        let mut next_seq = 0u64;
        let mut addr: Vec<u64> = Vec::new(); // file offset of each block written so far
        let mut file_off = 0u64;
        // Metas can point one block past the job's last (a record ending exactly at a block end).
        let mut held: Vec<RecMeta> = Vec::new();
        let voff = |addr: &Vec<u64>, file_off: u64, b: u64, o: u32| -> u64 {
            let a = if (b as usize) < addr.len() { addr[b as usize] } else { file_off };
            a << 16 | o as u64
        };
        while let Ok((seq, blocks, metas)) = done_rx.recv() {
            pending.insert(seq, (blocks, metas));
            while let Some((blocks, metas)) = pending.remove(&next_seq) {
                for b in &blocks { addr.push(file_off); w.write_all(b)?; file_off += b.len() as u64; }
                held.extend(metas);
                let known = addr.len() as u64;
                let mut keep = Vec::new();
                for m in held.drain(..) {
                    if m.eb < known || (m.eb == known && m.eo == 0) {
                        let (vb, ve) = (voff(&addr, file_off, m.sb, m.so), voff(&addr, file_off, m.eb, m.eo));
                        idx.push(&m, vb, ve);
                    } else { keep.push(m); }
                }
                held = keep;
                next_seq += 1;
            }
        }
        for m in held.drain(..) {
            let (vb, ve) = (voff(&addr, file_off, m.sb, m.so), voff(&addr, file_off, m.eb, m.eo));
            idx.push(&m, vb, ve);
        }
        w.write_all(&BGZF_EOF)?;
        w.flush()?;
        idx.write(&format!("{}.bai", out_path_w))
    });

    // Merge (this thread): header block(s), then records in order.
    let mut seq = 0u64;
    let mut blocks: Vec<Vec<u8>> = Vec::new();
    let mut metas: Vec<RecMeta> = Vec::new();
    let mut cur: Vec<u8> = Vec::with_capacity(BLOCK);
    let mut block_no: u64 = 0;
    let t_merge_start = std::time::Instant::now();
    let send_blocked = std::cell::Cell::new(std::time::Duration::ZERO);
    let send = |seq: &mut u64, blocks: &mut Vec<Vec<u8>>, metas: &mut Vec<RecMeta>| -> io::Result<()> {
        let ts = std::time::Instant::now();
        let _g = scopeguard(|| send_blocked.set(send_blocked.get() + ts.elapsed()));
        job_tx.send(Job { seq: *seq, blocks: std::mem::take(blocks), metas: std::mem::take(metas) }).map_err(|_| bad("compressor stopped"))?;
        *seq += 1;
        Ok(())
    };
    // Header: htslib writes it in its own block(s).
    let mut h = Vec::new();
    h.extend_from_slice(b"BAM\x01");
    h.extend_from_slice(&(header.len() as i32).to_le_bytes());
    h.extend_from_slice(header.as_bytes());
    h.extend_from_slice(&(refs.len() as i32).to_le_bytes());
    for (name, len) in refs {
        h.extend_from_slice(&((name.len() + 1) as i32).to_le_bytes());
        h.extend_from_slice(name.as_bytes());
        h.push(0);
        h.extend_from_slice(&len.to_le_bytes());
    }
    for c in h.chunks(BLOCK) { blocks.push(c.to_vec()); block_no += 1; }

    // Sources in input order: runs (creation order), then the in-memory batches; ties go to the
    // earliest source, so equal keys keep input order.
    let dec = (threads / runs.len().max(1)).clamp(1, 4);
    let mut readers: Vec<RunReader> = runs.iter().map(|p| RunReader::open(p, dec)).collect::<io::Result<_>>()?;
    let n_runs = readers.len();
    let mut mem_i = vec![0usize; mem.len()];
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    for (s, r) in readers.iter_mut().enumerate() {
        if let Some(k) = r.advance()? { heap.push(Reverse((k, s))); }
    }
    for (m, b) in mem.iter().enumerate() {
        if let Some(&(k, _)) = b.ents.first() { heap.push(Reverse((k, n_runs + m))); }
    }
    let mut n_rec = 0u64;
    while let Some(Reverse((_, s))) = heap.pop() {
        let rec: &[u8] = if s >= n_runs {
            let m = s - n_runs;
            mem_i[m] += 1;
            mem[m].record(mem[m].ents[mem_i[m] - 1].1)
        } else {
            readers[s].current()
        };
        // As htslib's bam_write1: start a new block if the record does not fit in this one.
        if !cur.is_empty() && cur.len() + rec.len() > BLOCK {
            blocks.push(std::mem::replace(&mut cur, Vec::with_capacity(BLOCK)));
            block_no += 1;
            if blocks.len() >= BLOCKS_PER_JOB { send(&mut seq, &mut blocks, &mut metas)?; }
        }
        let (sb, so) = (block_no, cur.len() as u32);
        let mut rest = rec;
        while !rest.is_empty() {
            let take = (BLOCK - cur.len()).min(rest.len());
            cur.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if cur.len() == BLOCK {
                blocks.push(std::mem::replace(&mut cur, Vec::with_capacity(BLOCK)));
                block_no += 1;
            }
        }
        metas.push(rec_meta(rec, sb, so, block_no, cur.len() as u32));
        // Refill source s.
        if s >= n_runs {
            let m = s - n_runs;
            if let Some(&(k, _)) = mem[m].ents.get(mem_i[m]) { heap.push(Reverse((k, s))); }
        } else if let Some(k) = readers[s].advance()? {
            heap.push(Reverse((k, s)));
        }
        if blocks.len() >= BLOCKS_PER_JOB { send(&mut seq, &mut blocks, &mut metas)?; }
        n_rec += 1;
    }
    if !cur.is_empty() { blocks.push(cur); }
    if !blocks.is_empty() || !metas.is_empty() { send(&mut seq, &mut blocks, &mut metas)?; }
    let t_merge_end = std::time::Instant::now();
    let merge_busy = (t_merge_end - t_merge_start).saturating_sub(send_blocked.get());
    drop(job_tx);
    for h in workers { h.join().map_err(|_| bad("compressor panicked"))?; }
    writer.join().map_err(|_| bad("writer panicked"))??;
    if std::env::var_os("BONEFISH_SORT_TIMING").is_some() {
        eprintln!("[sort] merge loop done {:.2}s before compressors+writer finished; merge thread busy {:.2}s (of which waiting for runs {:.2}s); run readers: file reads {:.2}s, decompression {:.2}s", t_merge_end.elapsed().as_secs_f64(), merge_busy.as_secs_f64(), RUN_WAIT_NS.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, RUN_READ_NS.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9, RUN_DEC_NS.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9);
    }
    Ok(n_rec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reg2bin_matches_htslib() {
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(reg2bin(0, 1), 4681);
        assert_eq!(reg2bin(0, 16385), 585);
        assert_eq!(reg2bin(0, 1 << 29), 0);
    }

    #[test]
    fn int_tag_types() {
        let mut v = Vec::new();
        for (x, t) in [(0i64, b'C'), (255, b'C'), (256, b'S'), (65536, b'I'), (-1, b'c'), (-129, b's'), (-40000, b'i')] {
            v.clear();
            push_int(x, &mut v);
            assert_eq!(v[0], t, "{x}");
        }
    }
}
