//! Robustness stress test for the read/reference parsers. Random and mutated inputs (FASTA and
//! FASTQ text, BGZF containers, BAM records, gzip headers with garbage) go through every reader:
//! `Reader` (several buffer sizes, `read_next`, `for_each_record`, `read_batch`),
//! `parse_fasta_bytes`, the chunk-fed streamers, `BamFastq`, `BgzfParallel`,
//! `looks_like_bam`/`looks_like_bgzf` and, for a sample of cases, `open` and `header_lines` on a
//! file. Each must return Ok or Err: no panic on any thread and no hang (a watchdog exits the
//! process if a case makes no progress for 60 s).
//!
//! `parser_stress_fast` runs with the normal test suite. The long run is ignored by default:
//!
//!   BONEFISH_STRESS_ITERS=3000000 cargo test --release -p bonefish-core --test parser_stress -- --ignored
//!
//! `BONEFISH_STRESS_SEED` picks another case stream; `BONEFISH_STRESS_CASE=<case seed>` replays one
//! case. Failing inputs are written to `$CARGO_TARGET_TMPDIR/parser-stress-<case seed>.bin`.
//! Uses a SplitMix64 generator, so cases are deterministic and no dependency is needed.

use std::io::{BufRead, BufReader, Cursor};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, Once};

use bonefish::fasta::bam::{self, BamFastq, BgzfParallel};
use bonefish::fasta::{self, FastaStreamer, FastqStreamer, Reader};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize { (self.next() % n.max(1) as u64) as usize }
    fn chance(&mut self, num: usize, den: usize) -> bool { self.below(den) < num }
    fn pick<'a, T>(&mut self, xs: &'a [T]) -> &'a T { &xs[self.below(xs.len())] }
    fn bytes(&mut self, n: usize) -> Vec<u8> { (0..n).map(|_| self.next() as u8).collect() }
}

// ---------------------------------------------------------------------------------------------
// Generators

fn name(rng: &mut Rng) -> Vec<u8> {
    let len = match rng.below(100) {
        0..=3 => 0,
        4..=7 => 254,
        8..=10 => 255,
        11 => 1000 + rng.below(100_000),
        _ => 1 + rng.below(30),
    };
    (0..len).map(|_| match rng.below(40) {
        0 => rng.next() as u8,
        1 => b' ',
        2 => b'\t',
        _ => b'a' + rng.below(26) as u8,
    }).collect()
}

fn eol(rng: &mut Rng) -> &'static [u8] {
    if rng.chance(1, 4) { b"\r\n" } else { b"\n" }
}

fn bases(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| if rng.chance(1, 20) { *rng.pick(b"NnUuacgtRY") } else { *rng.pick(b"ACGT") }).collect()
}

/// `s` as lines of one random width (an empty `s` gives one empty line).
fn wrapped(rng: &mut Rng, o: &mut Vec<u8>, s: &[u8]) {
    if s.is_empty() {
        o.extend_from_slice(eol(rng));
        return;
    }
    let w = if rng.chance(3, 4) { s.len() } else { 1 + rng.below(80) };
    for line in s.chunks(w) {
        o.extend_from_slice(line);
        o.extend_from_slice(eol(rng));
    }
}

const COMMENTS: [&[u8]; 6] = [
    b"MM:Z:C+m?,1,2;\tML:B:C,200,10",
    b"runid=abc read=12 ch=5 start_time=2024-01-01T00:00:00Z",
    b"RG:Z:grp\tqs:i:12\tNM:i:3",
    b"",
    b"x",
    b"  \t ",
];

fn gen_fastq(rng: &mut Rng) -> Vec<u8> {
    let mut o = Vec::new();
    for _ in 0..rng.below(6) {
        if rng.chance(1, 8) { o.extend_from_slice(eol(rng)); }
        o.push(b'@');
        o.extend(name(rng));
        if rng.chance(1, 3) {
            o.push(*rng.pick(b" \t"));
            o.extend_from_slice(rng.pick(&COMMENTS));
        }
        o.extend_from_slice(eol(rng));
        let n = if rng.chance(1, 10) { 0 } else { rng.below(300) };
        let s = bases(rng, n);
        wrapped(rng, &mut o, &s);
        o.push(b'+');
        if rng.chance(1, 4) { o.extend(name(rng)); }
        o.extend_from_slice(eol(rng));
        let qn = match rng.below(10) { 0 => usize::MAX, 1 => rng.below(300), _ => n };
        let q: Vec<u8> = if qn == usize::MAX { b"*".to_vec() } else { (0..qn).map(|_| b'!' + rng.below(60) as u8).collect() };
        wrapped(rng, &mut o, &q);
    }
    o
}

fn gen_fasta(rng: &mut Rng) -> Vec<u8> {
    let mut o = Vec::new();
    for _ in 0..rng.below(6) {
        if rng.chance(1, 8) { o.extend_from_slice(eol(rng)); }
        o.push(b'>');
        o.extend(name(rng));
        if rng.chance(1, 3) {
            o.push(b' ');
            o.extend_from_slice(rng.pick(&COMMENTS));
        }
        o.extend_from_slice(eol(rng));
        let n = if rng.chance(1, 10) { 0 } else { rng.below(500) };
        let s = bases(rng, n);
        wrapped(rng, &mut o, &s);
    }
    o
}

fn put_i32(o: &mut Vec<u8>, v: i32) { o.extend_from_slice(&v.to_le_bytes()); }

/// One BAM aux field of a random type.
fn gen_tag(rng: &mut Rng, r: &mut Vec<u8>) {
    r.extend_from_slice(rng.pick(&[*b"MM", *b"ML", *b"RG", *b"qs", *b"NM", *b"ch", *b"XY", *b"ms"]));
    match rng.below(10) {
        0 => { r.push(b'A'); r.push(b'!' + rng.below(90) as u8); }
        1 | 2 => {
            r.push(*rng.pick(b"ZH"));
            let n = rng.below(40);
            r.extend((0..n).map(|_| b'0' + rng.below(40) as u8));
            r.push(0);
        }
        3 => { r.push(b'f'); r.extend(rng.bytes(4)); }
        4..=6 => {
            let t = *rng.pick(b"cCsSiI");
            r.push(t);
            let sz = match t { b'c' | b'C' => 1, b's' | b'S' => 2, _ => 4 };
            r.extend(rng.bytes(sz));
        }
        _ => {
            let sub = *rng.pick(b"cCsSiIf");
            let sz = match sub { b'c' | b'C' => 1, b's' | b'S' => 2, _ => 4 };
            let n = rng.below(20);
            r.push(b'B');
            r.push(sub);
            r.extend_from_slice(&(n as u32).to_le_bytes());
            r.extend(rng.bytes(n * sz));
        }
    }
}

/// An uncompressed BAM stream: header, references and a few unaligned/aligned records.
fn gen_bam(rng: &mut Rng) -> Vec<u8> {
    let mut o = b"BAM\x01".to_vec();
    let text: &[u8] = rng.pick(&[
        b"@HD\tVN:1.6\tSO:unknown\n@RG\tID:rg1\tSM:s\n@PG\tID:basecaller\n" as &[u8],
        b"",
        b"@RG\tID:a\n\0\0",
    ]);
    put_i32(&mut o, text.len() as i32);
    o.extend_from_slice(text);
    let n_ref = rng.below(3);
    put_i32(&mut o, n_ref as i32);
    for k in 0..n_ref {
        let nm = format!("chr{k}\0");
        put_i32(&mut o, nm.len() as i32);
        o.extend_from_slice(nm.as_bytes());
        put_i32(&mut o, 1000);
    }
    for _ in 0..rng.below(5) {
        let mut r = Vec::new();
        put_i32(&mut r, if rng.chance(1, 2) { -1 } else { 0 });
        put_i32(&mut r, rng.below(1000) as i32);
        let nm: Vec<u8> = (0..1 + rng.below(40)).map(|_| b'a' + rng.below(26) as u8).collect();
        r.push(nm.len() as u8 + 1);
        r.push(rng.below(61) as u8);
        r.extend_from_slice(&4680u16.to_le_bytes());
        let n_cigar = if rng.chance(1, 2) { 0 } else { 1 + rng.below(3) };
        r.extend_from_slice(&(n_cigar as u16).to_le_bytes());
        let flag: u16 = *rng.pick(&[0, 4, 16, 20, 0x100, 0x800, 0x900]);
        r.extend_from_slice(&flag.to_le_bytes());
        let l_seq = rng.below(200);
        put_i32(&mut r, l_seq as i32);
        put_i32(&mut r, -1);
        put_i32(&mut r, -1);
        put_i32(&mut r, 0);
        r.extend_from_slice(&nm);
        r.push(0);
        for _ in 0..n_cigar { r.extend_from_slice(&((l_seq as u32) << 4).to_le_bytes()); }
        r.extend(rng.bytes(l_seq.div_ceil(2)));
        if rng.chance(1, 5) {
            r.extend(std::iter::repeat_n(0xff, l_seq));
        } else {
            r.extend((0..l_seq).map(|_| rng.below(60) as u8));
        }
        for _ in 0..rng.below(6) { gen_tag(rng, &mut r); }
        put_i32(&mut o, r.len() as i32);
        o.extend_from_slice(&r);
    }
    o
}

/// One BGZF block holding `data`.
fn bgzf_block(c: &mut libdeflater::Compressor, data: &[u8], o: &mut Vec<u8>) {
    let mut def = vec![0u8; c.deflate_compress_bound(data.len())];
    let n = c.deflate_compress(data, &mut def).unwrap();
    let bsize = 18 + n + 8;
    o.extend_from_slice(&[0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0]);
    o.extend_from_slice(&((bsize - 1) as u16).to_le_bytes());
    o.extend_from_slice(&def[..n]);
    o.extend_from_slice(&libdeflater::crc32(data).to_le_bytes());
    o.extend_from_slice(&(data.len() as u32).to_le_bytes());
}

fn bgzf(rng: &mut Rng, raw: &[u8]) -> Vec<u8> {
    let mut c = libdeflater::Compressor::default();
    let max_block = if rng.chance(1, 2) { 65280 } else { 1 + rng.below(512) };
    let mut o = Vec::new();
    for chunk in raw.chunks(max_block) { bgzf_block(&mut c, chunk, &mut o); }
    if rng.chance(3, 4) { bgzf_block(&mut c, b"", &mut o); } // EOF marker block
    o
}

/// Plain (non-BGZF) gzip, which `open` hands to the threaded MultiGzDecoder.
fn gzip(raw: &[u8]) -> Vec<u8> {
    let mut c = libdeflater::Compressor::default();
    let mut o = vec![0u8; c.gzip_compress_bound(raw.len())];
    let n = c.gzip_compress(raw, &mut o).unwrap();
    o.truncate(n);
    o
}

const BYTES: [u8; 14] = [0, 0xff, 0x7f, 0x80, b'\n', b'\r', b'>', b'@', b'+', b'*', b' ', b'\t', 0x1f, 0x8b];
const INSERTS: [&[u8]; 10] = [b"\n", b"\r", b"\r\n", b"\n\n", b"@", b">", b"+", b"\n+\n", b"\n@\n", b"*"];
const INTS: [i32; 12] = [-1, 0, 1, -100, 255, 256, 0xffff, 0x10000, 1 << 24, i32::MIN, i32::MAX, i32::MAX - 3];

/// 0-4 random edits: truncation, bit flips, special bytes, inserted line breaks, deleted or
/// duplicated ranges, random bytes, extreme 16/32-bit fields (lengths, sizes), donor splices.
fn mutate(rng: &mut Rng, d: &mut Vec<u8>, donor: &[u8]) {
    for _ in 0..rng.below(5) {
        let len = d.len();
        match rng.below(11) {
            0 if len > 0 => d.truncate(rng.below(len)),
            1 if len > 0 => { let i = rng.below(len); d[i] ^= 1 << rng.below(8); }
            2 if len > 0 => { let i = rng.below(len); d[i] = *rng.pick(&BYTES); }
            3 => { let i = rng.below(len + 1); let s = *rng.pick(&INSERTS); d.splice(i..i, s.iter().copied()); }
            4 if len > 0 => { let i = rng.below(len); let j = (i + 1 + rng.below(64)).min(len); d.drain(i..j); }
            5 if len > 0 => {
                let i = rng.below(len);
                let j = (i + 1 + rng.below(64)).min(len);
                let c = d[i..j].to_vec();
                d.splice(j..j, c);
            }
            6 => { let i = rng.below(len + 1); let n = 1 + rng.below(16); let b = rng.bytes(n); d.splice(i..i, b); }
            7 if len >= 4 => { let i = rng.below(len - 3); d[i..i + 4].copy_from_slice(&rng.pick(&INTS).to_le_bytes()); }
            8 if !donor.is_empty() => {
                let a = rng.below(donor.len());
                let b = (a + rng.below(256)).min(donor.len());
                let i = rng.below(len + 1);
                d.splice(i..i, donor[a..b].iter().copied());
            }
            9 if len >= 2 => { let i = rng.below(len - 1); d[i..i + 2].copy_from_slice(&(rng.next() as u16).to_le_bytes()); }
            10 if len > 0 => { let i = rng.below(len); d[i] = rng.next() as u8; }
            _ => {}
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Kind { Fastq, Fasta, Random, BamPayload, BamContainer, TextBgzf, GzipGarbage, PlainGzip }
const N_KINDS: usize = 8;

fn make_case(rng: &mut Rng) -> (Kind, Vec<u8>) {
    let donor = if rng.chance(1, 2) { gen_fastq(rng) } else { gen_bam(rng) };
    let text = |rng: &mut Rng| if rng.chance(1, 2) { gen_fastq(rng) } else { gen_fasta(rng) };
    match rng.below(12) {
        0..=2 => { let mut d = gen_fastq(rng); mutate(rng, &mut d, &donor); (Kind::Fastq, d) }
        3 => { let mut d = gen_fasta(rng); mutate(rng, &mut d, &donor); (Kind::Fasta, d) }
        4 => {
            let prefix: &[u8] = rng.pick(&[
                b"" as &[u8], b"@", b">", b"\x1f\x8b", b"\x1f\x8b\x08\x04", b"BAM\x01",
                b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0",
            ]);
            let mut d = prefix.to_vec();
            let n = rng.below(512);
            d.extend(rng.bytes(n));
            (Kind::Random, d)
        }
        5 | 6 => {
            let mut raw = gen_bam(rng);
            mutate(rng, &mut raw, &donor);
            (Kind::BamPayload, bgzf(rng, &raw))
        }
        7 | 8 => {
            let raw = gen_bam(rng);
            let mut d = bgzf(rng, &raw);
            mutate(rng, &mut d, &donor);
            (Kind::BamContainer, d)
        }
        9 => {
            let mut raw = text(rng);
            mutate(rng, &mut raw, &donor);
            let mut d = bgzf(rng, &raw);
            if rng.chance(1, 4) { mutate(rng, &mut d, &donor); }
            (Kind::TextBgzf, d)
        }
        10 => {
            // gzip magic, optionally with a BC subfield of random BSIZE, then garbage
            let mut d = vec![0x1f, 0x8b, 8, *rng.pick(&[0u8, 4, 0xff])];
            d.extend(rng.bytes(6));
            if rng.chance(3, 4) {
                d.extend_from_slice(&[6, 0, b'B', b'C', 2, 0]);
                d.extend_from_slice(&(rng.next() as u16 % *rng.pick(&[8u16, 32, 64, 0xffff])).to_le_bytes());
            }
            let n = rng.below(256);
            d.extend(rng.bytes(n));
            (Kind::GzipGarbage, d)
        }
        _ => {
            let mut raw = if rng.chance(1, 2) { text(rng) } else { gen_bam(rng) };
            mutate(rng, &mut raw, &donor);
            let mut d = gzip(&raw);
            if rng.chance(1, 2) { mutate(rng, &mut d, &donor); }
            (Kind::PlainGzip, d)
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Consumers

/// Read every record; after an error the reader must stay finished.
fn drain<R: BufRead>(mut r: Reader<R>) {
    loop {
        match r.read_next() {
            Ok(Some(rec)) => { std::hint::black_box(rec); }
            Ok(None) => break,
            Err(_) => {
                assert!(matches!(r.read_next(), Ok(None)), "reader not finished after an error");
                break;
            }
        }
    }
}

fn feed_text(rng: &mut Rng, d: &[u8]) {
    let cap = *rng.pick(&[1usize, 2, 3, 7, 64, 1 << 16]);
    let mut r = Reader::new(BufReader::with_capacity(cap, Cursor::new(d)));
    r.set_keep_description(rng.chance(1, 2));
    drain(r);

    let mut r = Reader::new(Cursor::new(d));
    let _ = r.for_each_record(|rec| { std::hint::black_box(rec); });

    let mut r = Reader::new(Cursor::new(d));
    let batch = rng.below(500) as u64;
    for _ in 0..d.len() + 2 {
        match r.read_batch(batch) {
            Ok((_, false)) => {}
            Ok((_, true)) | Err(_) => break,
        }
    }

    let _ = std::hint::black_box(fasta::parse_fasta_bytes(d));

    let mut fa = FastaStreamer::new().with_rna_to_dna(rng.chance(1, 2));
    let mut fq = FastqStreamer::new().with_rna_to_dna(rng.chance(1, 2));
    let step = 1 + rng.below(64);
    for chunk in d.chunks(step) {
        fa.push(chunk);
        fq.push(chunk);
        while fa.next_record().is_some() {}
        while fq.next_record().is_some() {}
    }
    fa.finalize();
    fq.finalize();
    while fa.next_record().is_some() {}
    while fq.next_record().is_some() {}

    if d.len() < 4096 {
        let _ = std::hint::black_box(bam::filter_comment_tags(&String::from_utf8_lossy(d)));
    }
}

fn feed_binary(rng: &mut Rng, d: &[u8]) {
    std::hint::black_box(bam::looks_like_bam(d));
    std::hint::black_box(bam::looks_like_bgzf(d));
    let cap = *rng.pick(&[1usize, 5, 1 << 16]);
    let mut r = Reader::new(BufReader::with_capacity(cap, BamFastq::new(Cursor::new(d.to_vec()))));
    r.set_keep_description(rng.chance(1, 2));
    drain(r);
    drain(Reader::new(BgzfParallel::new(Cursor::new(d.to_vec()), 1 + rng.below(3))));
}

fn feed_file(rng: &mut Rng, d: &[u8], tmp: &Path) {
    let gz = d.starts_with(&[0x1f, 0x8b]);
    let tid = format!("{:?}", std::thread::current().id()).replace(|c: char| !c.is_ascii_alphanumeric(), "");
    let path = tmp.join(format!("parser-stress-{}-{tid}.{}", std::process::id(), if gz { "fq.gz" } else { "fq" }));
    std::fs::write(&path, d).unwrap();
    if let Ok(mut r) = fasta::open(&path) {
        r.set_keep_description(rng.chance(1, 2));
        drain(r);
    }
    let _ = std::hint::black_box(bam::header_lines(&path));
    let _ = std::fs::remove_file(&path);
}

fn feed(rng: &mut Rng, kind: Kind, d: &[u8], file: bool, tmp: &Path) {
    feed_text(rng, d);
    let binary = !matches!(kind, Kind::Fastq | Kind::Fasta);
    if binary || rng.chance(1, 16) { feed_binary(rng, d); }
    if file { feed_file(rng, d, tmp); }
}

// ---------------------------------------------------------------------------------------------
// Driver

static PANICS: AtomicUsize = AtomicUsize::new(0);
static FIRST_PANIC: Mutex<Option<(String, String)>> = Mutex::new(None); // (location, report)

/// Count panics on every thread (the BAM/BGZF readers decompress on their own threads, where a
/// panic would otherwise only show up as an early end or a generic error).
fn install_panic_counter() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let default = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            default(info);
            PANICS.fetch_add(1, Ordering::SeqCst);
            let loc = info.location().map(|l| l.to_string()).unwrap_or_default();
            let msg = format!("{info} (thread {:?})", std::thread::current().name());
            let mut first = FIRST_PANIC.lock().unwrap_or_else(|e| e.into_inner());
            first.get_or_insert((loc, msg));
        }));
    });
}

fn case_seed(seed: u64, i: u64) -> u64 {
    Rng(seed ^ i.wrapping_mul(0xD1B5_4A32_D192_ED03)).next()
}

fn tmp_dir() -> PathBuf { PathBuf::from(env!("CARGO_TARGET_TMPDIR")) }

/// Run one case; Err((panic location, report)) if anything panicked.
fn run_case(cs: u64, file: bool, counts: &mut [u64; N_KINDS]) -> Result<(), (String, String)> {
    let mut rng = Rng(cs);
    let (kind, d) = make_case(&mut rng);
    counts[kind as usize] += 1;
    let before = PANICS.load(Ordering::SeqCst);
    let res = catch_unwind(AssertUnwindSafe(|| feed(&mut rng, kind, &d, file, &tmp_dir())));
    if res.is_ok() && PANICS.load(Ordering::SeqCst) == before {
        return Ok(());
    }
    let path = tmp_dir().join(format!("parser-stress-{cs:#018x}.bin"));
    let _ = std::fs::write(&path, &d);
    let (loc, msg) = FIRST_PANIC.lock().unwrap_or_else(|e| e.into_inner()).take().unwrap_or_default();
    Err((loc, format!("case seed {cs:#018x} ({kind:?}, {} bytes, saved to {}): {msg}", d.len(), path.display())))
}

fn stress(seed: u64, iters: u64, file_every: u64) {
    install_panic_counter();
    let progress = Arc::new(AtomicU64::new(0));
    let current = Arc::new(AtomicU64::new(0));
    let done = Arc::new(AtomicBool::new(false));
    {
        let (progress, current, done) = (progress.clone(), current.clone(), done.clone());
        std::thread::spawn(move || {
            let (mut last, mut still) = (u64::MAX, 0);
            while !done.load(Ordering::Relaxed) {
                std::thread::sleep(std::time::Duration::from_secs(1));
                let p = progress.load(Ordering::Relaxed);
                if p != last { (last, still) = (p, 0); continue; }
                still += 1;
                if still >= 60 && !done.load(Ordering::Relaxed) {
                    let cs = current.load(Ordering::Relaxed);
                    println!("parser_stress: no progress for 60 s at case seed {cs:#018x} (hang)");
                    eprintln!("parser_stress: no progress for 60 s at case seed {cs:#018x} (hang)");
                    std::process::exit(2);
                }
            }
        });
    }

    let t0 = std::time::Instant::now();
    let mut counts = [0u64; N_KINDS];
    // Panic sites with their first case and number of cases.
    let mut failures: std::collections::BTreeMap<String, (String, u64)> = Default::default();
    let mut n_failed = 0u64;
    let replay = env_u64("BONEFISH_STRESS_CASE");
    let iters = if replay.is_some() { 1 } else { iters };
    for i in 0..iters {
        let cs = replay.unwrap_or_else(|| case_seed(seed, i));
        current.store(cs, Ordering::Relaxed);
        progress.store(i + 1, Ordering::Relaxed);
        if let Err((loc, e)) = run_case(cs, replay.is_some() || i % file_every == 0, &mut counts) {
            n_failed += 1;
            let site = failures.entry(loc).or_insert_with(|| {
                println!("parser_stress: PANIC {e}");
                (e, 0)
            });
            site.1 += 1;
            if failures.len() >= 20 || n_failed >= 10_000 { break; }
        }
    }
    done.store(true, Ordering::Relaxed);
    println!(
        "parser_stress: seed {seed:#x}: {} cases in {:.1} s ({}), {} with a panic",
        counts.iter().sum::<u64>(),
        t0.elapsed().as_secs_f64(),
        ["fastq", "fasta", "random", "bam-payload", "bam-container", "text-bgzf", "gzip-garbage", "plain-gzip"]
            .iter().zip(counts).map(|(k, n)| format!("{k} {n}")).collect::<Vec<_>>().join(", "),
        n_failed,
    );
    let report: Vec<String> = failures.values().map(|(e, n)| format!("[{n} cases] {e}")).collect();
    assert!(failures.is_empty(), "parser panics at {} sites:\n{}", failures.len(), report.join("\n"));
}

fn env_u64(key: &str) -> Option<u64> {
    std::env::var(key).ok().and_then(|s| {
        let s = s.replace('_', "");
        match s.strip_prefix("0x") { Some(h) => u64::from_str_radix(h, 16).ok(), None => s.parse().ok() }
    })
}

#[test]
fn parser_stress_fast() {
    stress(env_u64("BONEFISH_STRESS_SEED").unwrap_or(0x5eed), 1500, 16);
}

#[test]
#[ignore = "long run: BONEFISH_STRESS_ITERS (default 2M) cases, release build recommended"]
fn parser_stress_long() {
    let iters = env_u64("BONEFISH_STRESS_ITERS").unwrap_or(2_000_000);
    stress(env_u64("BONEFISH_STRESS_SEED").unwrap_or(0xb0_7e_f1_54), iters, 64);
}

/// Read all records: Ok(number of records) or Err(message of the first error).
fn read_all<R: BufRead>(mut r: Reader<R>) -> Result<usize, String> {
    let mut n = 0;
    while let Some(_rec) = r.read_next().map_err(|e| e.to_string())? { n += 1; }
    Ok(n)
}

/// A BGZF stream (with EOF block) holding a BAM header without references, then `records`:
/// (record without its block_size, block_size if not the record's length).
fn bam_bgzf(l_text: i32, records: &[(&[u8], Option<i32>)]) -> Vec<u8> {
    let mut raw = b"BAM\x01".to_vec();
    put_i32(&mut raw, l_text);
    put_i32(&mut raw, 0);
    for (rec, size) in records {
        put_i32(&mut raw, size.unwrap_or(rec.len() as i32));
        raw.extend_from_slice(rec);
    }
    let mut c = libdeflater::Compressor::default();
    let mut o = Vec::new();
    bgzf_block(&mut c, &raw, &mut o);
    bgzf_block(&mut c, b"", &mut o);
    o
}

/// An unmapped BAM record named `r1` with `l_seq` bases and the given aux bytes.
fn bam_record(l_seq: i32, aux: &[u8]) -> Vec<u8> {
    let mut r = Vec::new();
    put_i32(&mut r, -1);
    put_i32(&mut r, -1);
    r.extend_from_slice(&[3, 0, 0x48, 0x12, 0, 0, 4, 0]); // l_read_name, mapq, bin, n_cigar, flag
    put_i32(&mut r, l_seq);
    put_i32(&mut r, -1);
    put_i32(&mut r, -1);
    put_i32(&mut r, 0);
    r.extend_from_slice(b"r1\0");
    let n = l_seq.max(0) as usize;
    r.extend(std::iter::repeat_n(0x11, n.div_ceil(2)));
    r.extend(std::iter::repeat_n(20, n));
    r.extend_from_slice(aux);
    r
}

/// Inputs that panicked before the BAM/BGZF readers were hardened (found by the stress run).
#[test]
fn parser_regressions() {
    install_panic_counter();
    let before = PANICS.load(Ordering::SeqCst);
    let bam = |d: &[u8]| read_all(Reader::new(BufReader::new(BamFastq::new(Cursor::new(d.to_vec())))));
    let bgzf_text = |d: &[u8]| read_all(Reader::new(BgzfParallel::new(Cursor::new(d.to_vec()), 2)));

    // Sanity: a valid record reads.
    assert_eq!(bam(&bam_bgzf(0, &[(&bam_record(4, b"qsC\x05"), None)])), Ok(1));

    // BGZF header whose BC subfield gives BSIZE 0 (block size 1): looks_like_bam sliced
    // head[bsize - 4..]; the BAM reader thread and the bgzip splitter sliced block[18..]. The
    // splitter's panic ended the stream early with no error (silent truncation).
    let mut tiny = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0, 0];
    tiny.extend_from_slice(&[0; 32]);
    assert!(!bam::looks_like_bam(&tiny));
    assert!(!bam::looks_like_bgzf(&tiny));
    assert!(bam(&tiny).is_err());
    assert!(bgzf_text(&tiny).is_err());

    // ISIZE far above the 64 KiB BGZF limit (allocated before decompressing).
    let mut big = bam_bgzf(0, &[]);
    let bsize = u16::from_le_bytes([big[16], big[17]]) as usize + 1;
    big[bsize - 4..bsize].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(!bam::looks_like_bam(&big));
    assert!(bam(&big).is_err());
    assert!(bgzf_text(&big).is_err());

    // Negative i32 lengths: l_text, block_size, l_seq (overflowed, then sliced out of range).
    assert!(bam(&bam_bgzf(-1, &[])).is_err());
    assert!(bam(&bam_bgzf(0, &[(&bam_record(4, b""), Some(-1))])).is_err());
    assert!(bam(&bam_bgzf(0, &[(&bam_record(4, b""), Some(i32::MIN))])).is_err());
    assert!(bam(&bam_bgzf(0, &[(&bam_record(-2, b""), None)])).is_err());

    // B array of subtype A: was read as a 4-byte integer from a 1-byte element.
    assert!(bam(&bam_bgzf(0, &[(&bam_record(4, b"XYBA\x01\0\0\0z"), None)])).is_err());

    assert_eq!(PANICS.load(Ordering::SeqCst), before, "a parser panicked: {:?}", FIRST_PANIC.lock().unwrap());
}
