//! Chunk-fed FASTA / FASTQ parsers.
//!
//! Bytes are pushed in arbitrarily-sized chunks (e.g. from a `ReadableStream`
//! reader in JS, or a chunked file read) and complete records pop out the
//! other side as `(name, sequence_bytes)`. Designed for WASM streaming so the
//! caller never has to hold a multi-GB buffer.
//!
//! No allocations for the sequence in the steady state beyond the final
//! `Vec<u8>` that gets handed off — partial sequences accumulate into the
//! pending record's vec, which becomes the yielded record at end of input.

use std::collections::VecDeque;

/// Streaming FASTA parser. Push bytes via [`push`](Self::push); pop completed
/// `(name, seq)` pairs via [`next_record`](Self::next_record); flush the
/// in-flight record (if any) at end-of-input via [`finalize`](Self::finalize).
pub struct FastaStreamer {
    line_buf: Vec<u8>,
    cur_name: Option<String>,
    cur_seq: Vec<u8>,
    completed: VecDeque<(String, Vec<u8>)>,
    rna_to_dna: bool,
}

impl Default for FastaStreamer {
    fn default() -> Self { Self::new() }
}

impl FastaStreamer {
    pub fn new() -> Self {
        Self {
            line_buf: Vec::new(),
            cur_name: None,
            cur_seq: Vec::new(),
            completed: VecDeque::new(),
            rna_to_dna: true,
        }
    }

    /// If `false`, leave 'U'/'u' bases alone (default: `true`, U→T mapping).
    pub fn with_rna_to_dna(mut self, enabled: bool) -> Self {
        self.rna_to_dna = enabled;
        self
    }

    /// Feed a chunk of bytes. Records that complete inside this chunk are
    /// queued for [`next_record`](Self::next_record).
    pub fn push(&mut self, chunk: &[u8]) {
        let mut start = 0;
        while let Some(off) = memchr(b'\n', &chunk[start..]) {
            let line_end = start + off;
            if self.line_buf.is_empty() {
                self.process_line(&chunk[start..line_end]);
            } else {
                self.line_buf.extend_from_slice(&chunk[start..line_end]);
                let buf = std::mem::take(&mut self.line_buf);
                self.process_line(&buf);
            }
            start = line_end + 1;
        }
        if start < chunk.len() {
            self.line_buf.extend_from_slice(&chunk[start..]);
        }
    }

    fn process_line(&mut self, line: &[u8]) {
        let line = strip_cr(line);
        if line.is_empty() { return; }
        if line[0] == b'>' {
            // Finish previous record.
            if let Some(name) = self.cur_name.take() {
                let seq = std::mem::take(&mut self.cur_seq);
                self.completed.push_back((name, seq));
            }
            // New name = first whitespace-separated token after '>'.
            let after_gt = &line[1..];
            let name_end = after_gt.iter().position(u8::is_ascii_whitespace).unwrap_or(after_gt.len());
            let name = String::from_utf8_lossy(&after_gt[..name_end]).into_owned();
            self.cur_name = Some(name);
        } else if self.cur_name.is_some() {
            append_seq(&mut self.cur_seq, line, self.rna_to_dna);
        }
        // else: pre-header noise (e.g. blank lines), ignored
    }

    pub fn next_record(&mut self) -> Option<(String, Vec<u8>)> {
        self.completed.pop_front()
    }

    /// Flush the trailing partial line + the in-flight record. Call once after
    /// the last `push`. The flushed record (if any) is appended to the queue
    /// — drain with `next_record`.
    pub fn finalize(&mut self) {
        if !self.line_buf.is_empty() {
            let buf = std::mem::take(&mut self.line_buf);
            self.process_line(&buf);
        }
        if let Some(name) = self.cur_name.take() {
            let seq = std::mem::take(&mut self.cur_seq);
            self.completed.push_back((name, seq));
        }
    }
}

// ─── FASTQ ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FqLine { Header, Sequence, Plus, Quality }

/// Streaming FASTQ parser. 4-line records: `@name`, sequence, `+`, quality.
/// Quality scores are consumed and discarded; only `(name, seq)` is yielded,
/// to keep the WASM API identical between FASTA and FASTQ inputs.
pub struct FastqStreamer {
    line_buf: Vec<u8>,
    expect: FqLine,
    cur_name: Option<String>,
    cur_seq: Vec<u8>,
    qual_remaining: usize,
    completed: VecDeque<(String, Vec<u8>)>,
    rna_to_dna: bool,
}

impl Default for FastqStreamer {
    fn default() -> Self { Self::new() }
}

impl FastqStreamer {
    pub fn new() -> Self {
        Self {
            line_buf: Vec::new(),
            expect: FqLine::Header,
            cur_name: None,
            cur_seq: Vec::new(),
            qual_remaining: 0,
            completed: VecDeque::new(),
            rna_to_dna: true,
        }
    }

    pub fn with_rna_to_dna(mut self, enabled: bool) -> Self {
        self.rna_to_dna = enabled;
        self
    }

    pub fn push(&mut self, chunk: &[u8]) {
        let mut start = 0;
        while let Some(off) = memchr(b'\n', &chunk[start..]) {
            let line_end = start + off;
            if self.line_buf.is_empty() {
                self.process_line(&chunk[start..line_end]);
            } else {
                self.line_buf.extend_from_slice(&chunk[start..line_end]);
                let buf = std::mem::take(&mut self.line_buf);
                self.process_line(&buf);
            }
            start = line_end + 1;
        }
        if start < chunk.len() {
            self.line_buf.extend_from_slice(&chunk[start..]);
        }
    }

    fn process_line(&mut self, line: &[u8]) {
        let line = strip_cr(line);
        match self.expect {
            FqLine::Header => {
                if line.is_empty() || line[0] != b'@' { return; }  // tolerate blanks pre-record
                let after_at = &line[1..];
                let name_end = after_at.iter().position(u8::is_ascii_whitespace).unwrap_or(after_at.len());
                self.cur_name = Some(String::from_utf8_lossy(&after_at[..name_end]).into_owned());
                self.cur_seq.clear();
                self.expect = FqLine::Sequence;
            }
            FqLine::Sequence => {
                append_seq(&mut self.cur_seq, line, self.rna_to_dna);
                self.expect = FqLine::Plus;
            }
            FqLine::Plus => {
                let _ = line;
                self.qual_remaining = self.cur_seq.len();
                self.expect = FqLine::Quality;
                if self.qual_remaining == 0 {
                    // empty sequence: skip the quality line and emit immediately
                    self.emit_record();
                    self.expect = FqLine::Header;
                }
            }
            FqLine::Quality => {
                // Consume quality chars equal to sequence length, but accept any.
                if line.len() >= self.qual_remaining {
                    self.qual_remaining = 0;
                } else {
                    self.qual_remaining -= line.len();
                }
                if self.qual_remaining == 0 {
                    self.emit_record();
                    self.expect = FqLine::Header;
                }
            }
        }
    }

    fn emit_record(&mut self) {
        if let Some(name) = self.cur_name.take() {
            let seq = std::mem::take(&mut self.cur_seq);
            self.completed.push_back((name, seq));
        }
    }

    pub fn next_record(&mut self) -> Option<(String, Vec<u8>)> {
        self.completed.pop_front()
    }

    pub fn finalize(&mut self) {
        if !self.line_buf.is_empty() {
            let buf = std::mem::take(&mut self.line_buf);
            self.process_line(&buf);
        }
        // If we have a name + seq but didn't see the quality line, emit anyway.
        if matches!(self.expect, FqLine::Quality | FqLine::Plus) {
            self.emit_record();
            self.expect = FqLine::Header;
        }
    }
}

// ─── shared helpers ────────────────────────────────────────────────────────

#[inline]
fn strip_cr(line: &[u8]) -> &[u8] {
    if line.last() == Some(&b'\r') { &line[..line.len() - 1] } else { line }
}

#[inline]
fn memchr(needle: u8, haystack: &[u8]) -> Option<usize> {
    haystack.iter().position(|&b| b == needle)
}

#[inline]
fn append_seq(dst: &mut Vec<u8>, line: &[u8], rna_to_dna: bool) {
    dst.reserve(line.len());
    for &b in line {
        if b.is_ascii_whitespace() { continue; }
        let c = if rna_to_dna && (b == b'U' || b == b'u') { b - 1 } else { b };
        dst.push(c);
    }
}

// ─── tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn collect_fasta<I: IntoIterator<Item = &'static [u8]>>(chunks: I) -> Vec<(String, Vec<u8>)> {
        let mut s = FastaStreamer::new();
        for c in chunks { s.push(c); }
        s.finalize();
        let mut out = Vec::new();
        while let Some(r) = s.next_record() { out.push(r); }
        out
    }

    fn collect_fastq<I: IntoIterator<Item = &'static [u8]>>(chunks: I) -> Vec<(String, Vec<u8>)> {
        let mut s = FastqStreamer::new();
        for c in chunks { s.push(c); }
        s.finalize();
        let mut out = Vec::new();
        while let Some(r) = s.next_record() { out.push(r); }
        out
    }

    #[test]
    fn fasta_single_record() {
        let r = collect_fasta([b">chr1\nACGT\n".as_slice()]);
        assert_eq!(r, vec![("chr1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fasta_no_trailing_newline() {
        let r = collect_fasta([b">chr1\nACGT".as_slice()]);
        assert_eq!(r, vec![("chr1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fasta_two_records() {
        let r = collect_fasta([b">chr1\nACGT\n>chr2\nGGGG\n".as_slice()]);
        assert_eq!(r, vec![
            ("chr1".to_string(), b"ACGT".to_vec()),
            ("chr2".to_string(), b"GGGG".to_vec()),
        ]);
    }

    #[test]
    fn fasta_multiline_sequence() {
        let r = collect_fasta([b">chr1\nACG\nTAC\nGT\n".as_slice()]);
        assert_eq!(r, vec![("chr1".to_string(), b"ACGTACGT".to_vec())]);
    }

    #[test]
    fn fasta_chunk_splits_header() {
        let r = collect_fasta([b">chr".as_slice(), b"1\nACGT\n".as_slice()]);
        assert_eq!(r, vec![("chr1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fasta_chunk_splits_sequence() {
        let r = collect_fasta([b">chr1\nAC".as_slice(), b"GT\n>chr2\nTTTT\n".as_slice()]);
        assert_eq!(r, vec![
            ("chr1".to_string(), b"ACGT".to_vec()),
            ("chr2".to_string(), b"TTTT".to_vec()),
        ]);
    }

    #[test]
    fn fasta_byte_by_byte() {
        let input = b">chr1\nACGT\n>chr2\nGGGG\n";
        let mut s = FastaStreamer::new();
        for &b in input { s.push(&[b]); }
        s.finalize();
        let mut out = Vec::new();
        while let Some(r) = s.next_record() { out.push(r); }
        assert_eq!(out, vec![
            ("chr1".to_string(), b"ACGT".to_vec()),
            ("chr2".to_string(), b"GGGG".to_vec()),
        ]);
    }

    #[test]
    fn fasta_crlf_line_endings() {
        let r = collect_fasta([b">chr1\r\nACGT\r\n".as_slice()]);
        assert_eq!(r, vec![("chr1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fasta_name_strips_at_whitespace() {
        let r = collect_fasta([b">chr1 description here\nACGT\n".as_slice()]);
        assert_eq!(r[0].0, "chr1");
    }

    #[test]
    fn fasta_rna_to_dna() {
        let r = collect_fasta([b">rna\nACGU\n".as_slice()]);
        assert_eq!(r[0].1, b"ACGT");
    }

    #[test]
    fn fastq_single_record() {
        let r = collect_fastq([b"@r1\nACGT\n+\n!!!!\n".as_slice()]);
        assert_eq!(r, vec![("r1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fastq_two_records() {
        let r = collect_fastq([b"@r1\nACGT\n+\n!!!!\n@r2\nGGGG\n+\nIIII\n".as_slice()]);
        assert_eq!(r, vec![
            ("r1".to_string(), b"ACGT".to_vec()),
            ("r2".to_string(), b"GGGG".to_vec()),
        ]);
    }

    #[test]
    fn fastq_chunk_splits_record() {
        let r = collect_fastq([b"@r1\nAC".as_slice(), b"GT\n+\n!!!!\n".as_slice()]);
        assert_eq!(r, vec![("r1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fastq_byte_by_byte() {
        let input = b"@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\nIIII\n";
        let mut s = FastqStreamer::new();
        for &b in input { s.push(&[b]); }
        s.finalize();
        let mut out = Vec::new();
        while let Some(r) = s.next_record() { out.push(r); }
        assert_eq!(out, vec![
            ("r1".to_string(), b"ACGT".to_vec()),
            ("r2".to_string(), b"TTTT".to_vec()),
        ]);
    }

    #[test]
    fn fastq_no_trailing_newline() {
        let r = collect_fastq([b"@r1\nACGT\n+\n!!!!".as_slice()]);
        assert_eq!(r, vec![("r1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fastq_crlf_line_endings() {
        let r = collect_fastq([b"@r1\r\nACGT\r\n+\r\n!!!!\r\n".as_slice()]);
        assert_eq!(r, vec![("r1".to_string(), b"ACGT".to_vec())]);
    }

    #[test]
    fn fastq_name_strips_at_whitespace() {
        let r = collect_fastq([b"@r1 some comment\nACGT\n+\n!!!!\n".as_slice()]);
        assert_eq!(r[0].0, "r1");
    }
}
