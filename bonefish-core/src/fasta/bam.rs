//! BAM input (e.g. unaligned ONT basecaller output) as a FASTQ text stream.
//!
//! `BamFastq` decodes BGZF blocks with libdeflate and turns each primary or unmapped record
//! (flags 0x100 and 0x800 skipped, like `samtools fastq -F 0x900`) into a FASTQ record. Reverse-
//! strand records are reverse-complemented back to the read's original orientation. With `-y`
//! the record's tags are written as the FASTQ comment in SAM text form, except alignment tags
//! that the new alignment replaces (see [`keep_tag`]), so dorado's per-read tags (`MM`/`ML`/`MN`,
//! which refer to the original orientation, `RG`, `qs`, `ch`, ...) carry into the output as with
//! `samtools fastq -T '*' | bonefish -y`. [`header_lines`] returns the header's `@RG` lines.

use std::io::{self, Read};

/// Tags written by an aligner for a particular alignment: stale after re-alignment, or
/// duplicated by bonefish's own. Not copied from the input.
const ALIGNMENT_TAGS: [[u8; 2]; 27] = [
    *b"NM", *b"MD", *b"AS", *b"SA", *b"XS", *b"XA", *b"MC", *b"MQ", *b"HI", *b"NH", *b"H0", *b"H1",
    *b"H2", *b"OA", *b"OC", *b"cs", *b"cg", *b"ds", *b"tp", *b"cm", *b"s1", *b"s2", *b"de", *b"dv",
    *b"rl", *b"zd", *b"nn",
];

static COPY: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);
static DROP_RG: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Whether to write the tags at all (they are only used with `-y`; ML arrays are large).
pub fn set_copy_tags(on: bool) {
    COPY.store(on, std::sync::atomic::Ordering::Relaxed);
}

/// Drop the input's `RG` tags (`-R` sets the read group instead).
pub fn set_drop_rg(on: bool) {
    DROP_RG.store(on, std::sync::atomic::Ordering::Relaxed);
}

/// Is a read tag copied to the output with `-y`? Everything except alignment tags (and `RG`
/// when `-R` replaces it). `ms` (minimap2's DP score) is an alignment tag too.
pub fn keep_tag(tag: [u8; 2]) -> bool {
    !(ALIGNMENT_TAGS.contains(&tag) || tag == *b"ms"
        || (tag == *b"RG" && DROP_RG.load(std::sync::atomic::Ordering::Relaxed)))
}

/// The SAM-style tags (`XX:T:value`) of a FASTQ/FASTA comment that `-y` copies, tab-separated;
/// other words (free text, `key=value` pairs) are dropped, so the SAM stays valid. None if
/// nothing is left.
pub fn filter_comment_tags(comment: &str) -> Option<String> {
    let is_tag = |t: &str| {
        let b = t.as_bytes();
        b.len() >= 5 && b[0].is_ascii_alphabetic() && b[1].is_ascii_alphanumeric() && b[2] == b':'
            && matches!(b[3], b'A' | b'i' | b'f' | b'Z' | b'H' | b'B') && b[4] == b':'
            && keep_tag([b[0], b[1]])
    };
    let mut out = String::with_capacity(comment.len());
    for field in comment.split('\t') {
        let mut push = |t: &str| {
            if is_tag(t) {
                if !out.is_empty() { out.push('\t'); }
                out.push_str(t);
            }
        };
        if is_tag(field) { push(field) } else { field.split(' ').for_each(&mut push) }
    }
    (!out.is_empty()).then_some(out)
}

/// The `@RG` lines of a BAM header (read from the first BGZF blocks of `input`), or an empty
/// list if the input is not BAM.
pub fn header_lines<R: Read>(input: R) -> io::Result<Vec<String>> {
    let mut f = io::BufReader::new(input);
    let mut dec = libdeflater::Decompressor::new();
    let mut block = Vec::with_capacity(1 << 16);
    let mut data = Vec::new();
    let mut need = 8usize;
    while data.len() < need {
        if !next_block(&mut f, &mut dec, &mut block, &mut data)? {
            return Ok(Vec::new());
        }
        if data.len() >= 8 {
            if &data[..4] != b"BAM\x01" {
                return Ok(Vec::new());
            }
            need = 8 + i32::from_le_bytes(data[4..8].try_into().unwrap()).max(0) as usize;
        }
    }
    let text = String::from_utf8_lossy(&data[8..need]);
    Ok(text.lines().filter(|l| l.starts_with("@RG\t")).map(|l| l.trim_end_matches('\0').to_string()).collect())
}

/// Decimal text of 0..=255 with a leading comma, for `B:C` arrays (ML).
static DEC_C: std::sync::LazyLock<Vec<Vec<u8>>> =
    std::sync::LazyLock::new(|| (0..256).map(|v| format!(",{v}").into_bytes()).collect());

/// Is this the start of a BGZF file holding BAM data? `head` must hold the first BGZF block.
pub fn looks_like_bam(head: &[u8]) -> bool {
    if head.len() < 18 || head[0] != 0x1f || head[1] != 0x8b || head[2] != 8 || head[3] & 4 == 0 {
        return false;
    }
    let Some(bsize) = bgzf_block_size(head) else { return false };
    if head.len() < bsize {
        return false;
    }
    let isize = u32::from_le_bytes(head[bsize - 4..bsize].try_into().unwrap()) as usize;
    if isize > MAX_ISIZE {
        return false;
    }
    let mut out = vec![0u8; isize];
    let mut dec = libdeflater::Decompressor::new();
    matches!(dec.gzip_decompress(&head[..bsize], &mut out), Ok(n) if n >= 4 && &out[..4] == b"BAM\x01")
}

/// Smallest BGZF block: 18-byte header and 8-byte footer (CRC32, ISIZE).
const MIN_BLOCK: usize = 26;
/// Largest uncompressed BGZF block (ISIZE), as in htslib.
const MAX_ISIZE: usize = 1 << 16;

/// Total size of the BGZF block starting at `h[0]` (from the `BC` extra subfield); None if
/// there is no such field or the size cannot hold the header and footer.
fn bgzf_block_size(h: &[u8]) -> Option<usize> {
    if h.len() < 12 {
        return None;
    }
    let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
    let mut i = 12;
    while i + 4 <= 12 + xlen && i + 4 <= h.len() {
        let slen = u16::from_le_bytes([h[i + 2], h[i + 3]]) as usize;
        if h[i] == b'B' && h[i + 1] == b'C' && slen == 2 && i + 6 <= h.len() {
            let size = u16::from_le_bytes([h[i + 4], h[i + 5]]) as usize + 1;
            return (size >= MIN_BLOCK).then_some(size);
        }
        i += 4 + slen;
    }
    None
}

pub struct BamFastq {
    /// Decompressed data from the BGZF worker thread, in file order (empty Vec = end).
    rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    /// Decompressed BAM bytes not yet parsed: `raw[pos..]`.
    raw: Vec<u8>,
    pos: usize,
    /// FASTQ text not yet returned: `out[out_pos..]`.
    out: Vec<u8>,
    out_pos: usize,
    header_done: bool,
    eof: bool,
    input_done: bool,
}

impl BamFastq {
    /// Decompression runs on its own thread (it is over half the work: ONT BAMs carry large
    /// base-modification tags), record conversion on the caller's.
    pub fn new<R: Read + Send + 'static>(inner: R) -> Self {
        let (tx, rx) = std::sync::mpsc::sync_channel::<io::Result<Vec<u8>>>(8);
        std::thread::spawn(move || {
            let mut inner = inner;
            let mut dec = libdeflater::Decompressor::new();
            let mut block = Vec::with_capacity(1 << 16);
            loop {
                // ~1 MB of decompressed data per message
                let mut chunk: Vec<u8> = Vec::with_capacity(1 << 20);
                let mut end = false;
                while chunk.len() < (1 << 20) - (1 << 16) {
                    match next_block(&mut inner, &mut dec, &mut block, &mut chunk) {
                        Ok(true) => {}
                        Ok(false) => { end = true; break; }
                        Err(e) => { let _ = tx.send(Err(e)); return; }
                    }
                }
                if !chunk.is_empty() && tx.send(Ok(chunk)).is_err() {
                    return;
                }
                if end {
                    let _ = tx.send(Ok(Vec::new()));
                    return;
                }
            }
        });
        BamFastq {
            rx,
            raw: Vec::with_capacity(1 << 21),
            pos: 0,
            out: Vec::with_capacity(1 << 20),
            out_pos: 0,
            header_done: false,
            eof: false,
            input_done: false,
        }
    }

    /// Append the next chunk of decompressed data to `raw`; false at end of input.
    fn next_chunk(&mut self) -> io::Result<bool> {
        if self.input_done {
            return Ok(false);
        }
        let chunk = self.rx.recv().map_err(|_| io::Error::other("BGZF reader thread stopped"))??;
        if chunk.is_empty() {
            self.input_done = true;
            return Ok(false);
        }
        if self.pos > 0 && self.pos >= self.raw.len() / 2 {
            self.raw.drain(..self.pos);
            self.pos = 0;
        }
        self.raw.extend_from_slice(&chunk);
        Ok(true)
    }

    /// Make at least `n` unparsed bytes available; false if the input ends first.
    fn want(&mut self, n: usize) -> io::Result<bool> {
        while self.raw.len() - self.pos < n {
            if !self.next_chunk()? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn i32_at(&self, off: usize) -> i32 {
        i32::from_le_bytes(self.raw[self.pos + off..self.pos + off + 4].try_into().unwrap())
    }

    /// A length field (i32) at `off`; None if negative.
    fn len_at(&self, off: usize) -> Option<usize> {
        usize::try_from(self.i32_at(off)).ok()
    }

    fn skip_header(&mut self) -> io::Result<()> {
        let bad = || io::Error::new(io::ErrorKind::InvalidData, "truncated or invalid BAM header");
        if !self.want(12)? || &self.raw[self.pos..self.pos + 4] != b"BAM\x01" {
            return Err(bad());
        }
        let l_text = self.len_at(4).ok_or_else(bad)?;
        if !self.want(8 + l_text + 4)? {
            return Err(bad());
        }
        self.pos += 8 + l_text;
        let n_ref = self.len_at(0).ok_or_else(bad)?;
        self.pos += 4;
        for _ in 0..n_ref {
            if !self.want(4)? {
                return Err(bad());
            }
            let l_name = self.len_at(0).ok_or_else(bad)?;
            if !self.want(4 + l_name + 4)? {
                return Err(bad());
            }
            self.pos += 4 + l_name + 4;
        }
        self.header_done = true;
        Ok(())
    }

    /// Convert records into `out` until it holds at least `target` bytes or the input ends.
    fn fill(&mut self, target: usize) -> io::Result<()> {
        if !self.header_done {
            self.skip_header()?;
        }
        while self.out.len() - self.out_pos < target {
            if !self.want(4)? {
                if self.raw.len() > self.pos {
                    return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BAM record"));
                }
                self.eof = true;
                return Ok(());
            }
            let block_size = self.len_at(0)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "malformed BAM record"))?;
            if !self.want(4 + block_size)? {
                return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BAM record"));
            }
            let rec = &self.raw[self.pos + 4..self.pos + 4 + block_size];
            if self.out_pos > 0 && self.out_pos == self.out.len() {
                self.out.clear();
                self.out_pos = 0;
            }
            append_fastq(rec, &mut self.out)?;
            self.pos += 4 + block_size;
        }
        Ok(())
    }
}

impl Read for BamFastq {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.out_pos == self.out.len() {
            if self.eof {
                return Ok(0);
            }
            self.out.clear();
            self.out_pos = 0;
            self.fill(1 << 20)?;
        }
        let n = buf.len().min(self.out.len() - self.out_pos);
        buf[..n].copy_from_slice(&self.out[self.out_pos..self.out_pos + n]);
        self.out_pos += n;
        Ok(n)
    }
}

/// Decompress the next BGZF block from `inner`, appending to `out`; false at end of input.
fn next_block<R: Read>(inner: &mut R, dec: &mut libdeflater::Decompressor, block: &mut Vec<u8>, out: &mut Vec<u8>) -> io::Result<bool> {
    let mut head = [0u8; 18];
    match read_full(inner, &mut head)? {
        0 => return Ok(false),
        18 => {}
        _ => return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BGZF block header")),
    }
    let bsize = bgzf_block_size(&head)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "not a BGZF block (missing BC field)"))?;
    block.clear();
    block.extend_from_slice(&head);
    block.resize(bsize, 0);
    if read_full(inner, &mut block[18..])? != bsize - 18 {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BGZF block"));
    }
    let isize = u32::from_le_bytes(block[bsize - 4..].try_into().unwrap()) as usize;
    if isize > MAX_ISIZE {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "BGZF block larger than 64 KiB"));
    }
    let start = out.len();
    out.reserve(isize);
    // Decompress straight into the spare capacity (no zero-fill).
    let spare = unsafe { std::slice::from_raw_parts_mut(out.as_mut_ptr().add(start), isize) };
    let n = dec.gzip_decompress(block, spare)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("BGZF decompression failed: {e:?}")))?;
    unsafe { out.set_len(start + n) };
    Ok(true)
}

fn read_full<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut got = 0;
    while got < buf.len() {
        match r.read(&mut buf[got..]) {
            Ok(0) => break,
            Ok(n) => got += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(got)
}

/// Append one BAM record (without its block_size prefix) as FASTQ, unless it is secondary or
/// supplementary.
fn append_fastq(rec: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    let bad = || io::Error::new(io::ErrorKind::InvalidData, "malformed BAM record");
    if rec.len() < 32 {
        return Err(bad());
    }
    let l_read_name = rec[8] as usize;
    let n_cigar = u16::from_le_bytes([rec[12], rec[13]]) as usize;
    let flag = u16::from_le_bytes([rec[14], rec[15]]);
    let l_seq = i32::from_le_bytes(rec[16..20].try_into().unwrap());
    if flag & 0x900 != 0 {
        return Ok(());
    }
    let l_seq = usize::try_from(l_seq).map_err(|_| bad())?;
    let name_end = 32 + l_read_name;
    let seq_off = name_end + 4 * n_cigar;
    let qual_off = seq_off + l_seq.div_ceil(2);
    let tag_off = qual_off + l_seq;
    if tag_off > rec.len() || l_read_name == 0 {
        return Err(bad());
    }
    let name = &rec[32..name_end - 1]; // drop the NUL
    let rev = flag & 0x10 != 0;
    const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

    out.push(b'@');
    out.extend_from_slice(name);
    if COPY.load(std::sync::atomic::Ordering::Relaxed) {
        append_tags(&rec[tag_off..], out).map_err(|_| bad())?;
    }
    out.push(b'\n');

    let packed = &rec[seq_off..qual_off];
    let base = |i: usize| -> u8 {
        let b = packed[i >> 1];
        NT16[(if i & 1 == 0 { b >> 4 } else { b & 15 }) as usize]
    };
    if !rev {
        out.extend((0..l_seq).map(base));
    } else {
        out.extend((0..l_seq).rev().map(|i| complement(base(i))));
    }
    out.extend_from_slice(b"\n+\n");
    let qual = &rec[qual_off..tag_off];
    if l_seq > 0 && qual[0] == 0xff {
        // no qualities stored: `*`, as in SAM (the FASTQ reader maps it to no qualities)
        out.push(b'*');
    } else if !rev {
        out.extend(qual.iter().map(|&q| q.saturating_add(33)));
    } else {
        out.extend(qual.iter().rev().map(|&q| q.saturating_add(33)));
    }
    out.push(b'\n');
    Ok(())
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
        b'M' => b'K', b'K' => b'M', b'R' => b'Y', b'Y' => b'R',
        b'S' => b'S', b'W' => b'W', b'V' => b'B', b'B' => b'V',
        b'H' => b'D', b'D' => b'H', b'=' => b'=', _ => b'N',
    }
}

/// Write the tags `keep_tag` accepts among a record's aux fields as tab-separated SAM text.
fn append_tags(mut aux: &[u8], out: &mut Vec<u8>) -> Result<(), ()> {
    use std::io::Write;
    while aux.len() >= 3 {
        let tag = [aux[0], aux[1]];
        let ty = aux[2];
        let body = &aux[3..];
        let want = keep_tag(tag);
        let int_size = |t: u8| match t { b'c' | b'C' => Some(1), b's' | b'S' => Some(2), b'i' | b'I' | b'f' => Some(4), _ => None };
        let read_int = |t: u8, b: &[u8]| -> i64 {
            match t {
                b'c' => b[0] as i8 as i64, b'C' => b[0] as i64,
                b's' => i16::from_le_bytes([b[0], b[1]]) as i64, b'S' => u16::from_le_bytes([b[0], b[1]]) as i64,
                b'i' => i32::from_le_bytes(b[..4].try_into().unwrap()) as i64,
                _ => u32::from_le_bytes(b[..4].try_into().unwrap()) as i64,
            }
        };
        let used = match ty {
            b'Z' | b'H' => {
                let end = body.iter().position(|&c| c == 0).ok_or(())?;
                if want {
                    out.push(b'\t');
                    out.extend_from_slice(&tag);
                    out.push(b':');
                    out.push(ty);
                    out.push(b':');
                    out.extend_from_slice(&body[..end]);
                }
                end + 1
            }
            b'B' => {
                if body.len() < 5 { return Err(()); }
                let sub = body[0];
                let n = u32::from_le_bytes(body[1..5].try_into().unwrap()) as usize;
                let sz = int_size(sub).ok_or(())?;
                if body.len() < 5 + n * sz { return Err(()); }
                if want {
                    write!(out, "\t{}{}:B:{}", tag[0] as char, tag[1] as char, sub as char).map_err(|_| ())?;
                    if sub == b'C' {
                        let dec = &*DEC_C;
                        for &v in &body[5..5 + n] {
                            out.extend_from_slice(&dec[v as usize]);
                        }
                    } else {
                    for k in 0..n {
                        let e = &body[5 + k * sz..5 + (k + 1) * sz];
                        if sub == b'f' {
                            write!(out, ",{}", f32::from_le_bytes(e.try_into().unwrap())).map_err(|_| ())?;
                        } else {
                            write!(out, ",{}", read_int(sub, e)).map_err(|_| ())?;
                        }
                    }
                    }
                }
                5 + n * sz
            }
            b'A' => {
                if body.is_empty() { return Err(()); }
                if want {
                    write!(out, "\t{}{}:A:{}", tag[0] as char, tag[1] as char, body[0] as char).map_err(|_| ())?;
                }
                1
            }
            b'f' => {
                if body.len() < 4 { return Err(()); }
                if want {
                    write!(out, "\t{}{}:f:{}", tag[0] as char, tag[1] as char, f32::from_le_bytes(body[..4].try_into().unwrap())).map_err(|_| ())?;
                }
                4
            }
            t => {
                let sz = int_size(t).ok_or(())?;
                if body.len() < sz { return Err(()); }
                if want {
                    write!(out, "\t{}{}:i:{}", tag[0] as char, tag[1] as char, read_int(t, body)).map_err(|_| ())?;
                }
                sz
            }
        };
        aux = &body[used..];
    }
    Ok(())
}

/// BGZF decompression on several threads, as a `BufRead` of the decompressed stream (for
/// bgzipped FASTQ/FASTA). A splitter thread cuts the file into batches of whole blocks (~1 MiB
/// compressed), `workers` threads decompress them with libdeflate, and the reader puts them
/// back in file order.
pub struct BgzfParallel {
    rx: std::sync::mpsc::Receiver<(u64, io::Result<Vec<u8>>)>,
    pending: std::collections::BTreeMap<u64, io::Result<Vec<u8>>>,
    next: u64,
    cur: Vec<u8>,
    pos: usize,
    eof: bool,
}

/// Whether `head` starts with a BGZF block (gzip with the `BC` extra subfield).
pub fn looks_like_bgzf(head: &[u8]) -> bool {
    head.len() >= 18 && head[0] == 0x1f && head[1] == 0x8b && head[2] == 8 && head[3] & 4 != 0
        && bgzf_block_size(head).is_some()
}

impl BgzfParallel {
    pub fn new<R: Read + Send + 'static>(mut inner: R, workers: usize) -> Self {
        const BATCH: usize = 1 << 20;
        let workers = workers.max(1);
        let (job_tx, job_rx) = std::sync::mpsc::sync_channel::<(u64, io::Result<Vec<u8>>)>(2 * workers);
        let job_rx = std::sync::Arc::new(std::sync::Mutex::new(job_rx));
        let (res_tx, res_rx) = std::sync::mpsc::sync_channel::<(u64, io::Result<Vec<u8>>)>(2 * workers);
        // Splitter: whole blocks into batches.
        std::thread::spawn(move || {
            let mut seq = 0u64;
            let mut batch: Vec<u8> = Vec::with_capacity(BATCH + 65536);
            loop {
                let mut head = [0u8; 18];
                let got = match read_full(&mut inner, &mut head) { Ok(g) => g, Err(e) => { let _ = job_tx.send((seq, Err(e))); return; } };
                if got == 0 { break; }
                let bsize = if got == 18 { bgzf_block_size(&head) } else { None };
                let Some(bsize) = bsize else {
                    let _ = job_tx.send((seq, Err(io::Error::new(io::ErrorKind::InvalidData, "not a BGZF block (missing BC field)"))));
                    return;
                };
                let start = batch.len();
                batch.extend_from_slice(&head);
                batch.resize(start + bsize, 0);
                match read_full(&mut inner, &mut batch[start + 18..]) {
                    Ok(n) if n == bsize - 18 => {}
                    Ok(_) => { let _ = job_tx.send((seq, Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BGZF block")))); return; }
                    Err(e) => { let _ = job_tx.send((seq, Err(e))); return; }
                }
                if batch.len() >= BATCH {
                    if job_tx.send((seq, Ok(std::mem::take(&mut batch)))).is_err() { return; }
                    batch.reserve(BATCH + 65536);
                    seq += 1;
                }
            }
            if !batch.is_empty() { let _ = job_tx.send((seq, Ok(batch))); }
        });
        for _ in 0..workers {
            let job_rx = std::sync::Arc::clone(&job_rx);
            let res_tx = res_tx.clone();
            std::thread::spawn(move || {
                let mut dec = libdeflater::Decompressor::new();
                loop {
                    let job = { let rx = job_rx.lock().unwrap(); rx.recv() };
                    let Ok((seq, job)) = job else { return; };
                    let out = job.and_then(|b| {
                        let mut out: Vec<u8> = Vec::with_capacity(4 * b.len());
                        let mut i = 0;
                        while i < b.len() {
                            let bsize = bgzf_block_size(&b[i..]).unwrap_or(b.len() - i);
                            let blk = &b[i..i + bsize];
                            let isize = u32::from_le_bytes(blk[bsize - 4..].try_into().unwrap()) as usize;
                            if isize > MAX_ISIZE {
                                return Err(io::Error::new(io::ErrorKind::InvalidData, "BGZF block larger than 64 KiB"));
                            }
                            let start = out.len();
                            out.reserve(isize);
                            // SAFETY: decompress into reserved spare capacity; length set from the result.
                            let spare = unsafe { std::slice::from_raw_parts_mut(out.as_mut_ptr().add(start), isize) };
                            let n = dec.gzip_decompress(blk, spare)
                                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("BGZF decompression failed: {e:?}")))?;
                            unsafe { out.set_len(start + n) };
                            i += bsize;
                        }
                        Ok(out)
                    });
                    if res_tx.send((seq, out)).is_err() { return; }
                }
            });
        }
        BgzfParallel { rx: res_rx, pending: Default::default(), next: 0, cur: Vec::new(), pos: 0, eof: false }
    }
}

impl io::Read for BgzfParallel {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let avail = io::BufRead::fill_buf(self)?;
        let n = avail.len().min(out.len());
        out[..n].copy_from_slice(&avail[..n]);
        io::BufRead::consume(self, n);
        Ok(n)
    }
}

impl io::BufRead for BgzfParallel {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.pos >= self.cur.len() && !self.eof {
            if let Some(r) = self.pending.remove(&self.next) {
                self.next += 1;
                self.cur = r?;
                self.pos = 0;
                continue;
            }
            match self.rx.recv() {
                Ok((seq, r)) => { self.pending.insert(seq, r); }
                Err(_) => { self.eof = true; }
            }
        }
        Ok(&self.cur[self.pos.min(self.cur.len())..])
    }
    fn consume(&mut self, n: usize) {
        self.pos += n;
    }
}
