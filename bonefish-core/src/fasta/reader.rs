use std::io::{self, BufRead};
use std::path::Path;
use std::fs::File;
use flate2::read::MultiGzDecoder;

/// Convert RNA bases U->T and u->t in place
/// ASCII 'U' = 0x55, 'T' = 0x54, 'u' = 0x75, 't' = 0x74
/// so just decrement
#[inline]
fn rna_to_dna_inplace(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if *b == b'U' || *b == b'u' { *b -= 1; }
    }
}

use thiserror::Error;
use crate::fasta::record::{Record, RefRecord};

#[derive(Error, Debug)]
pub enum FastxError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("not FASTA, FASTQ or BAM (the first record must start with '>' or '@')")]
    UnknownFormat,
    #[error("this looks like SAM text; give the reads as FASTQ, FASTA or BAM")]
    SamText,
    #[error("line {line}: malformed FASTQ: {msg}")]
    MalformedFastq { line: usize, msg: String },
    #[error("line {line}: file ends inside FASTQ record '{name}'")]
    Truncated { line: usize, name: String },
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum RecordType {
    Fasta,
    Fastq,
}

pub struct Reader<R> {
    reader: R,
    record_type: Option<RecordType>, // Detected on first read
    buf: Vec<u8>,
    line_number: usize,
    done: bool,
    keep_description: bool,
}

impl<R: BufRead> Reader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            record_type: None,
            buf: Vec::with_capacity(64 * 1024),
            line_number: 1,
            done: false,
            keep_description: true,
        }
    }

    /// Whether FASTQ records keep their header comment (default true). Without it the comment
    /// bytes (e.g. basecaller MM/ML tags, often most of the file) are skipped without copying.
    pub fn set_keep_description(&mut self, keep: bool) {
        self.keep_description = keep;
    }

    /// Read a FASTQ header line whose '@' is the next byte, returning the name only and skipping
    /// the rest of the line in place. Returns None (consuming nothing) if the name holds a
    /// non-ASCII byte, so the caller can take the general path (whose whitespace split is
    /// Unicode-aware); Some(None) at EOF.
    fn read_fastq_name_only(&mut self) -> Result<Option<Option<String>>, FastxError> {
        // Peek: the whole name must be in the buffer to decide without consuming.
        let (name, consumed, at_eol) = {
            let b = self.reader.fill_buf()?;
            if b.is_empty() { return Ok(Some(None)); }
            debug_assert!(b[0] == b'@');
            let mut i = 1;
            while i < b.len() && !matches!(b[i], b'\t' | b'\n' | 0x0b | 0x0c | b'\r' | b' ') {
                if b[i] >= 0x80 { return Ok(None); }
                i += 1;
            }
            if i == b.len() { return Ok(None); } // name runs past the buffer: general path
            // SAFETY: bytes 1..i are ASCII.
            let name = unsafe { String::from_utf8_unchecked(b[1..i].to_vec()) };
            (name, i, b[i] == b'\n')
        };
        self.reader.consume(consumed);
        if at_eol {
            self.reader.consume(1);
        } else {
            // Skip the rest of the line (the comment) without copying it.
            loop {
                let b = self.reader.fill_buf()?;
                if b.is_empty() { break; }
                match b.iter().position(|&c| c == b'\n') {
                    Some(p) => { self.reader.consume(p + 1); break; }
                    None => { let n = b.len(); self.reader.consume(n); }
                }
            }
        }
        Ok(Some(Some(name)))
    }

    /// `Ok(None)` = empty input (EOF before any byte): zero records
    /// `Ok(Some(fmt))` = recognized
    /// `Err(UnknownFormat)` only for non-empty, unrecognized input
    fn detect_format(&mut self) -> Result<Option<RecordType>, FastxError> {
        // Leading blank lines are skipped.
        loop {
            let buffer = self.reader.fill_buf()?;
            if buffer.is_empty() {
                return Ok(None);
            }
            match buffer.iter().position(|&c| !matches!(c, b'\n' | b'\r' | b' ' | b'\t')) {
                Some(0) => break,
                Some(i) => { self.reader.consume(i); }
                None => { let n = buffer.len(); self.reader.consume(n); }
            }
        }
        let buffer = self.reader.fill_buf()?;
        if buffer[0] == b'>' {
            Ok(Some(RecordType::Fasta))
        } else if buffer[0] == b'@' {
            let sam = [b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"];
            if sam.iter().any(|t| buffer.starts_with(*t)) {
                return Err(FastxError::SamText);
            }
            Ok(Some(RecordType::Fastq))
        } else {
            Err(FastxError::UnknownFormat)
        }
    }

    /// Iterate over records using a callback (Zero-Copy)
    pub fn for_each_record<F>(&mut self, callback: F) -> Result<(), FastxError>
    where
        F: FnMut(RefRecord),
    {
        if self.record_type.is_none() {
            match self.detect_format()? {
                Some(fmt) => self.record_type = Some(fmt),
                None => return Ok(()), // empty input: nothing to do
            }
        }

        match self.record_type.unwrap() {
            RecordType::Fastq => self.parse_fastq(callback),
            RecordType::Fasta => self.parse_fasta(callback),
        }
    }
    
    /// Returns an iterator over owned Records.
    /// Note: This copies data. Use `for_each_record` for maximum performance.
    pub fn records(self) -> RecordIter<R> {
        RecordIter { reader: self }
    }

    fn parse_fastq<F>(&mut self, mut callback: F) -> Result<(), FastxError>
    where
        F: FnMut(RefRecord),
    {
        while let Some(r) = self.read_next_fastq()? {
            let head = match r.description() {
                Some(d) => format!("{} {}", r.name(), d),
                None => r.name().to_string(),
            };
            callback(RefRecord { head: head.as_bytes(), seq: r.sequence(), qual: r.quality() });
        }
        Ok(())
    }

    fn parse_fasta<F>(&mut self, mut callback: F) -> Result<(), FastxError>
    where
        F: FnMut(RefRecord),
    {
        let mut seq_buf = Vec::new(); // Buffer for accumulating sequence lines
        let mut header_buf = Vec::new(); 
        
        self.buf.clear();
        if self.reader.read_until(b'\n', &mut self.buf)? == 0 { return Ok(()); }
        self.line_number += 1;
        
        if self.buf[0] != b'>' { return Err(FastxError::UnknownFormat); }
        
        loop {
            header_buf.clear();
            header_buf.extend_from_slice(&self.buf); 
            
            seq_buf.clear();
            
            loop {
                self.buf.clear();
                let n = self.reader.read_until(b'\n', &mut self.buf)?;
                if n == 0 {
                    break;
                }
                self.line_number += 1;
                
                if self.buf[0] == b'>' {
                    break;
                }
                
                let l_end = trim_newline(&self.buf);
                seq_buf.extend_from_slice(&self.buf[..l_end]);
            }
            
            let l1_end = trim_newline(&header_buf);
            
            let r = RefRecord {
                head: &header_buf[1..l1_end],
                seq: &seq_buf,
                qual: None,
            };
            callback(r);
            
            if self.buf.is_empty() {
                break;
            }
        }
        Ok(())
    }
}

fn trim_newline(buf: &[u8]) -> usize {
    if buf.ends_with(b"\r\n") {
        buf.len().saturating_sub(2)
    } else if buf.ends_with(b"\n") {
        buf.len().saturating_sub(1)
    } else {
        buf.len()
    }
}

pub struct RecordIter<R> {
    reader: Reader<R>,
}

impl<R: BufRead> Iterator for RecordIter<R> {
    type Item = Result<Record, FastxError>;

    fn next(&mut self) -> Option<Self::Item> {
        // This is inefficient because we cannot yield from inside the `for_each_record` loop 
        // without coroutines.
        // We have to reimplement `parse_one` logic or use a state machine.
        // Or adapt `for_each_record` to be re-entrant? No.
        
        // For convenience in high-level code, we just need `next()`.
        // Let's implement a `read_next` method on Reader that returns `Option<Record>`.
        match self.reader.read_next() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

impl<R: BufRead> Reader<R> {
     pub fn read_next(&mut self) -> Result<Option<Record>, FastxError> {
        if self.done {
            return Ok(None);
        }
        if self.record_type.is_none() {
            match self.detect_format() {
                Ok(Some(fmt)) => self.record_type = Some(fmt),
                Ok(None) => { self.done = true; return Ok(None); } // empty input: 0 records
                Err(e) => { self.done = true; return Err(e); }     // report once, then done
            }
        }

        let r = match self.record_type.unwrap() {
             RecordType::Fastq => self.read_next_fastq(),
             RecordType::Fasta => self.read_next_fasta(),
        };
        if r.is_err() {
            self.done = true;
        }
        r
    }
    
    fn read_next_fastq(&mut self) -> Result<Option<Record>, FastxError> {
        // Skip blank lines between records.
        loop {
            let peek = self.reader.fill_buf()?;
            if peek.is_empty() { return Ok(None); }
            match peek[0] {
                b'\n' => { self.reader.consume(1); self.line_number += 1; }
                b'\r' => { self.reader.consume(1); }
                _ => break,
            }
        }
        let peek = self.reader.fill_buf()?;
        if peek[0] != b'@' {
            return Err(FastxError::MalformedFastq { line: self.line_number, msg: "expected a header line starting with '@'".into() });
        }
        // Fast path: name only, comment skipped in place.
        if !self.keep_description && let Some(name) = self.read_fastq_name_only()? {
            let Some(name) = name else { return Ok(None); };
            self.line_number += 1;
            return self.read_fastq_body(name, None);
        }
        self.buf.clear();
        self.reader.read_until(b'\n', &mut self.buf)?;
        self.line_number += 1;
        let h_end = trim_newline(&self.buf);
        let name_line = String::from_utf8_lossy(&self.buf[1..h_end]).to_string(); // Excluding @
        let mut parts = name_line.splitn(2, char::is_whitespace);
        let name = parts.next().unwrap_or("").to_string();
        let desc = if self.keep_description { parts.next().map(|s| s.trim_start().to_string()).filter(|s| !s.is_empty()) } else { None };
        self.read_fastq_body(name, desc)
    }

    /// Sequence, '+' and quality lines of a FASTQ record whose header has been read. Sequence
    /// and quality may span several lines (the quality is read until it is as long as the
    /// sequence). A quality of exactly `*` means none (as SAM; written by the BAM reader).
    fn read_fastq_body(&mut self, name: String, desc: Option<String>) -> Result<Option<Record>, FastxError> {
        if name.is_empty() || name.len() > 254 {
            return Err(FastxError::MalformedFastq { line: self.line_number - 1, msg: format!("read name empty or longer than 254 characters (the SAM limit): '{name}'") });
        }
        let truncated = |s: &Self, name: &str| FastxError::Truncated { line: s.line_number, name: name.to_string() };
        // Sequence lines, up to the '+' line.
        let mut seq = Vec::new();
        loop {
            self.buf.clear();
            if self.reader.read_until(b'\n', &mut self.buf)? == 0 { return Err(truncated(self, &name)); }
            self.line_number += 1;
            if self.buf[0] == b'+' { break; }
            if self.buf[0] == b'@' && seq.is_empty() && trim_newline(&self.buf) > 1 {
                return Err(FastxError::MalformedFastq { line: self.line_number, msg: format!("record '{name}' has no '+' line") });
            }
            let e = trim_newline(&self.buf);
            seq.extend_from_slice(&self.buf[..e]);
        }
        rna_to_dna_inplace(&mut seq);
        // Quality lines until the quality is as long as the sequence.
        let mut qual = Vec::with_capacity(seq.len());
        while qual.len() < seq.len() {
            self.buf.clear();
            if self.reader.read_until(b'\n', &mut self.buf)? == 0 { return Err(truncated(self, &name)); }
            self.line_number += 1;
            let e = trim_newline(&self.buf);
            if e == 0 { break; }
            qual.extend_from_slice(&self.buf[..e]);
            if qual.as_slice() == b"*" && seq.len() > 1 { break; }
        }
        if qual.as_slice() == b"*" && seq.len() > 1 {
            return Ok(Some(Record::new(name, desc, seq, None)));
        }
        if qual.len() != seq.len() {
            return Err(FastxError::MalformedFastq {
                line: self.line_number,
                msg: format!("record '{name}': {} bases but {} quality values", seq.len(), qual.len()),
            });
        }
        Ok(Some(Record::new(name, desc, seq, Some(qual))))
    }

    fn read_next_fasta(&mut self) -> Result<Option<Record>, FastxError> {
        // We need a persistent buffer for the "next header" line in FASTA
        // because we read until '>' to finish current record.
        // Let's assume `self.buf` holds the header line entering this function
        // IF we are in a sequence.
        // But `self.buf` is dynamic.
        
        // Strategy: 
        // 1. If `self.buf` is empty, read header.
        // 2. If `self.buf` is not empty (contains last >Line), use it.
        
        if self.buf.is_empty() {
             if self.reader.read_until(b'\n', &mut self.buf)? == 0 { return Ok(None); }
             self.line_number += 1;
        }
        
        // Now self.buf should have header
        if self.buf[0] != b'>' {
            // maybe we finished? or bad format?
             // try reading one more?
             return Ok(None); 
        }
        
        let h_end = trim_newline(&self.buf);
        let name_line = String::from_utf8_lossy(&self.buf[1..h_end]).to_string();
        let mut parts = name_line.splitn(2, char::is_whitespace);
        let name = parts.next().unwrap_or("").to_string();
        let desc = parts.next().map(|s| s.to_string());
        if name.is_empty() || name.len() > 254 {
            return Err(FastxError::MalformedFastq { line: self.line_number - 1, msg: format!("sequence name empty or longer than 254 characters (the SAM limit): '{name}'") });
        }

        let mut seq = Vec::new();
        loop {
            self.buf.clear();
             let n = self.reader.read_until(b'\n', &mut self.buf)?;
             if n == 0 { 
                 // EOF
                 self.buf.clear(); // Clear to signal EOF for next call
                 break; 
             }
             self.line_number += 1;
             
             if self.buf[0] == b'>' {
                 // Next record found. `self.buf` has >Header.
                 // Break loop, return current record.
                 // Keep `self.buf` populated for next call.
                 break;
             }
             
             let l_end = trim_newline(&self.buf);
             seq.extend_from_slice(&self.buf[..l_end]);
        }
        rna_to_dna_inplace(&mut seq);

        Ok(Some(Record::new(name, desc, seq, None)))
    }

    /// Read sequences until cumulative bases exceed batch_size.
    /// Returns (seqs, is_eof). Caller can call again for the next batch.
    /// Check before the next read,
    /// meaning the last sequence that pushes past the limit IS included.
    #[allow(clippy::type_complexity)]
    pub fn read_batch(&mut self, batch_size: u64) -> Result<(Vec<(String, Vec<u8>)>, bool), FastxError> {
        let mut seqs = Vec::new();
        let mut sum_len: u64 = 0;
        loop {
            if sum_len > batch_size && !seqs.is_empty() {
                return Ok((seqs, false));
            }
            match self.read_next()? {
                Some(r) => {
                    sum_len += r.sequence().len() as u64;
                    seqs.push((r.name().to_string(), r.sequence().to_vec()));
                }
                None => return Ok((seqs, true)),
            }
        }
    }
}

/// Inputs whose first bytes were already consumed by [`input_read_groups`] (stdin, pipes):
/// the bytes read, replayed ahead of the rest of the stream by the next [`open`] of that path.
#[allow(clippy::type_complexity)]
static PREREAD: std::sync::Mutex<Vec<(std::path::PathBuf, Box<dyn io::Read + Send>)>> = std::sync::Mutex::new(Vec::new());

fn open_raw(path: &Path) -> io::Result<Box<dyn io::Read + Send>> {
    {
        let mut pre = PREREAD.lock().unwrap();
        if let Some(i) = pre.iter().position(|(p, _)| p == path) {
            return Ok(pre.swap_remove(i).1);
        }
    }
    Ok(if path.as_os_str() == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| io::Error::new(e.kind(), format!("cannot open '{}': {}", path.display(), e)))?)
    })
}

/// Copies everything read from `inner` into `seen`.
#[cfg(not(target_arch = "wasm32"))]
struct Recording<'a, R> { inner: &'a mut R, seen: Vec<u8> }
#[cfg(not(target_arch = "wasm32"))]
impl<R: io::Read> io::Read for Recording<'_, R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.seen.extend_from_slice(&buf[..n]);
        Ok(n)
    }
}

/// The `@RG` header lines of a BAM input (empty for FASTA/FASTQ). A regular file is read
/// again by [`open`]; stdin or a pipe can be read only once, so the bytes consumed here are
/// kept and handed back to [`open`].
#[cfg(not(target_arch = "wasm32"))]
pub fn input_read_groups<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
    let path = path.as_ref();
    if path.as_os_str() != "-" && std::fs::metadata(path)?.is_file() {
        return super::bam::header_lines(File::open(path)?);
    }
    let mut raw = open_raw(path)?;
    let mut rec = Recording { inner: &mut raw, seen: Vec::new() };
    let lines = super::bam::header_lines(&mut rec);
    let seen = std::mem::take(&mut rec.seen);
    PREREAD.lock().unwrap().push((path.to_path_buf(), Box::new(io::Read::chain(io::Cursor::new(seen), raw))));
    lines
}

pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Reader<Box<dyn BufRead + Send>>> {
    let file = open_raw(path.as_ref())?;

    // Detect by content (gzip magic 0x1f 0x8b) rather than trusting the extension.
    // `fill_buf` peeks without consuming, so the magic bytes remain available to
    // whichever decoder we wire up below. An empty file falls through to plain text.
    // Large enough to peek at a whole BGZF block (<= 64 KiB) when checking for BAM.
    let mut buf_reader = io::BufReader::with_capacity(1 << 17, file);
    let is_gz_content = {
        let peek = buf_reader.fill_buf()?;
        peek.len() >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
    };
    #[cfg(not(target_arch = "wasm32"))]
    if is_gz_content && super::bam::looks_like_bam(buf_reader.fill_buf()?) {
        return Ok(Reader::new(Box::new(io::BufReader::with_capacity(1 << 20, super::bam::BamFastq::new(buf_reader)))));
    }
    let ext_says_gz = path.as_ref().to_string_lossy().ends_with(".gz");
    let stdin = path.as_ref().as_os_str() == "-";
    if stdin {
    } else if ext_says_gz && !is_gz_content {
        eprintln!("Warning: '{}' has a .gz extension but is not gzipped; reading as plain text.", path.as_ref().display());
    } else if !ext_says_gz && is_gz_content {
        eprintln!("Warning: '{}' is gzipped but lacks a .gz extension; decompressing.", path.as_ref().display());
    }

    let reader: Box<dyn BufRead + Send> = if is_gz_content {
        // bgzip output (independent blocks) decompresses on several threads; plain gzip on one.
        #[cfg(not(target_arch = "wasm32"))]
        if super::bam::looks_like_bgzf(buf_reader.fill_buf()?) {
            Box::new(super::bam::BgzfParallel::new(buf_reader, 3))
        } else {
            Box::new(ThreadedRead::new(MultiGzDecoder::new(buf_reader)))
        }
        #[cfg(target_arch = "wasm32")]
        { Box::new(io::BufReader::new(MultiGzDecoder::new(buf_reader))) }
    } else {
        Box::new(buf_reader)
    };
    Ok(Reader::new(reader))
}

/// Reads `inner` (e.g. a gzip decoder) on a worker thread in 1 MiB blocks, so decompression
/// overlaps with parsing on the consuming thread.
#[cfg(not(target_arch = "wasm32"))]
struct ThreadedRead {
    rx: std::sync::mpsc::Receiver<io::Result<Vec<u8>>>,
    cur: Vec<u8>,
    pos: usize,
    eof: bool,
}

#[cfg(not(target_arch = "wasm32"))]
impl ThreadedRead {
    fn new<R: io::Read + Send + 'static>(mut inner: R) -> Self {
        const BLOCK: usize = 1 << 20;
        let (tx, rx) = std::sync::mpsc::sync_channel::<io::Result<Vec<u8>>>(8);
        std::thread::spawn(move || loop {
            let mut b = vec![0u8; BLOCK];
            let mut n = 0;
            while n < BLOCK {
                match inner.read(&mut b[n..]) {
                    Ok(0) => break,
                    Ok(k) => n += k,
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                    Err(e) => { let _ = tx.send(Err(e)); return; }
                }
            }
            b.truncate(n);
            let last = n == 0;
            if tx.send(Ok(b)).is_err() || last { return; }
        });
        ThreadedRead { rx, cur: Vec::new(), pos: 0, eof: false }
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl io::Read for ThreadedRead {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let avail = self.fill_buf()?;
        let n = avail.len().min(out.len());
        out[..n].copy_from_slice(&avail[..n]);
        self.consume(n);
        Ok(n)
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl BufRead for ThreadedRead {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.pos >= self.cur.len() && !self.eof {
            match self.rx.recv() {
                Ok(Ok(b)) => {
                    if b.is_empty() { self.eof = true; }
                    self.cur = b;
                    self.pos = 0;
                }
                Ok(Err(e)) => { self.eof = true; return Err(e); }
                // The worker ends with an empty block; a closed channel without it means the
                // thread died, which must not look like the end of the input.
                Err(_) => { self.eof = true; return Err(io::Error::other("decompression thread stopped before the end of the input")); }
            }
        }
        Ok(&self.cur[self.pos.min(self.cur.len())..])
    }
    fn consume(&mut self, n: usize) {
        self.pos += n;
    }
}

// ------ Memory-Mapped FASTA Reader for Index Building -----
/// Read all sequences from a FASTA file into memory.
/// Returns Vec of (name, sequence) pairs.
/// Used for index building where all sequences are needed at once.
pub fn read_fasta<P: AsRef<Path>>(path: P) -> io::Result<Vec<(String, Vec<u8>)>> {
    let data = std::fs::read(path)?;
    parse_fasta_bytes(&data)
}

/// Parse FASTA from a byte slice (works with mmap or regular buffers)
pub fn parse_fasta_bytes(data: &[u8]) -> io::Result<Vec<(String, Vec<u8>)>> {
    let mut sequences = Vec::new();

    if data.is_empty() {
        return Ok(sequences);
    }

    // The streaming reader (Reader::detect_format) requires the first byte to be
    // '>' before accepting input as FASTA. This mmap/batch path used to skip
    // forward to the first '>' *anywhere* in the buffer instead of checking the
    // first byte, which meant any non-FASTA binary file containing a stray 0x3E
    // byte (near-guaranteed for any file of a few hundred bytes or more — e.g. a
    // wrong path, a corrupted download, or an index file with an unrecognized
    // extension that fell through the index/sequence-file dispatch) was silently
    // "parsed" as one bogus sequence with binary garbage as its name, with no
    // error at all. Enforce the same first-byte check here for consistency.
    if data[0] != b'>' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "not a valid FASTA file: expected '>' as the first byte, found {:?}. \
                 Check that the correct file was provided — this is not a FASTA/FASTQ \
                 file, and if it's meant to be a pre-built index it may have an \
                 unrecognized extension (expected .mmi/.idx/.rmmi).",
                data[0] as char
            ),
        ));
    }

    // Use memchr for fast newline scanning
    let mut pos = 0;
    let len = data.len();

    while pos < len {
        // Skip to next '>'
        while pos < len && data[pos] != b'>' {
            pos += 1;
        }
        if pos >= len {
            break;
        }

        // Parse header line
        let header_start = pos + 1; // Skip '>'

        // Find end of header line using memchr
        let header_end = memchr::memchr(b'\n', &data[header_start..])
            .map(|i| header_start + i)
            .unwrap_or(len);

        // Extract name (first whitespace-delimited token)
        let header_bytes = &data[header_start..header_end];
        let name_end = header_bytes.iter()
            .position(|&b| b == b' ' || b == b'\t' || b == b'\r')
            .unwrap_or(header_bytes.len());

        let name = String::from_utf8_lossy(&header_bytes[..name_end]).to_string();

        pos = header_end + 1;

        // Collect sequence lines until next '>' or EOF
        // Pre-estimate capacity based on typical sequence density
        let mut seq = Vec::with_capacity(1024);

        while pos < len && data[pos] != b'>' {
            // Find end of line
            let line_end = memchr::memchr(b'\n', &data[pos..])
                .map(|i| pos + i)
                .unwrap_or(len);

            // Copy sequence bytes (skip newlines and handle \r\n)
            let line = &data[pos..line_end];
            let line = if line.ends_with(b"\r") {
                &line[..line.len() - 1]
            } else {
                line
            };

            // Normalize to uppercase while copying
            seq.extend(line.iter().map(|&b| b.to_ascii_uppercase()));

            pos = line_end + 1;
        }

        if name.is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "reference sequence with an empty name"));
        }
        if seq.is_empty() {
            eprintln!("Warning: skipping reference sequence '{name}' of length 0");
            continue;
        }
        sequences.push((name, seq));
    }

    if sequences.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no sequences in the reference"));
    }
    let mut seen = std::collections::HashSet::with_capacity(sequences.len());
    for (name, _) in &sequences {
        if !seen.insert(name.as_str()) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("duplicate reference sequence name '{name}'")));
        }
    }
    Ok(sequences)
}

// ------ Indexed Reader (Keep or minimal refactor) -----
// Needs logic for Indexed Reading.
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_detect_format_fasta() {
        let data = b">seq1\nACGT";
        let cursor = Cursor::new(data);
        let mut reader = Reader::new(cursor);
        let format = reader.detect_format().unwrap();
        assert_eq!(format, Some(RecordType::Fasta));
    }

    #[test]
    fn test_detect_format_fastq() {
        let data = b"@seq1\nACGT\n+\nIIII";
        let cursor = Cursor::new(data);
        let mut reader = Reader::new(cursor);
        let format = reader.detect_format().unwrap();
        assert_eq!(format, Some(RecordType::Fastq));
    }

    #[test]
    fn test_read_fasta_owned() {
        let data = b">seq1\nACGT\n>seq2 desc\nAAAA\nTTTT";
        let cursor = Cursor::new(data);
        let reader = Reader::new(cursor);
        let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), "seq1");
        assert_eq!(records[0].sequence(), b"ACGT");
        
        assert_eq!(records[1].name(), "seq2");
        assert_eq!(records[1].description(), Some("desc"));
        assert_eq!(records[1].sequence(), b"AAAATTTT");
    }

    #[test]
    fn test_read_fastq_owned() {
        let data = b"@seq1\nACGT\n+\nIIII\n@seq2\nAA\n+\nKK";
        let cursor = Cursor::new(data);
        let reader = Reader::new(cursor);
        let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name(), "seq1");
        assert_eq!(records[0].sequence(), b"ACGT");
        assert_eq!(records[0].quality().unwrap(), b"IIII");
        
        assert_eq!(records[1].name(), "seq2");
        assert_eq!(records[1].sequence(), b"AA");
        assert_eq!(records[1].quality().unwrap(), b"KK");
    }

    #[test]
    fn test_for_each_record_fasta() {
        let data = b">seq1\nACGT";
        let cursor = Cursor::new(data);
        let mut reader = Reader::new(cursor);
        let mut count = 0;
        reader.for_each_record(|r| {
            assert_eq!(r.head, b"seq1");
            assert_eq!(r.seq, b"ACGT");
            count += 1;
        }).unwrap();
        assert_eq!(count, 1);
    }

     #[test]
    fn test_for_each_record_fastq() {
        let data = b"@s1\nA\n+\nI";
        let cursor = Cursor::new(data);
        let mut reader = Reader::new(cursor);
        let mut count = 0;
        let _ = reader.for_each_record(|r| {
            assert_eq!(r.head, b"s1");
            assert_eq!(r.seq, b"A");
            assert_eq!(r.qual.unwrap(), b"I");
            count += 1;
        });
        assert_eq!(count, 1);
    }

    /// Empty input (e.g. an empty stdin query when only one positional is
    /// given, as in `-x ava-ont reads.fq`) must yield ZERO records and
    /// terminate cleanly — no error
    #[test]
    fn test_empty_input_is_clean_eof() {
        let mut it = Reader::new(Cursor::new(b"" as &[u8])).records();
        assert!(it.next().is_none(), "empty input must yield no records");
        // And it must STAY terminated (fused), not resurrect an error.
        assert!(it.next().is_none());
        assert!(it.next().is_none());
    }

    /// Non-empty but unrecognized input reports the error exactly ONCE
    #[test]
    fn test_bad_input_errors_once_then_fused() {
        let mut it = Reader::new(Cursor::new(b"not a fastx file\nsecond line\n" as &[u8])).records();
        assert!(matches!(it.next(), Some(Err(_))), "first read should surface the error");
        for _ in 0..1000 {
            assert!(it.next().is_none(), "iterator must not re-emit the error (no spin/spam)");
        }
    }

    /// Simulate a log-and-continue consumer (like the CLI read loops) over bad
    /// input: it must terminate with a bounded number of warnings, not forever
    #[test]
    fn test_consumer_loop_terminates_on_bad_input() {
        let it = Reader::new(Cursor::new(b"garbage not fastx\n" as &[u8])).records();
        let mut warnings = 0usize;
        let mut records = 0usize;
        for item in it {
            match item {
                Ok(_) => records += 1,
                Err(_) => { warnings += 1; continue; } // mirrors the CLI's warn+continue
            }
        }
        assert_eq!(records, 0);
        assert!(warnings <= 1, "bad input should warn at most once, got {warnings}");
    }
}
