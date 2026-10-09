//! End-to-end regression tests on a small synthetic ONT-like dataset (`tests/data`, made by
//! `scripts/make_regression_data.py`).
//!
//! * `lr:hq` (and its output-option variants) must stay byte-identical to minimap2 2.31: the
//!   expected hashes in `tests/data/expected.tsv` marked `minimap2-2.31` were computed from
//!   minimap2's own output, and `lrhq.minimap2.sam.gz` holds that SAM with SEQ/QUAL masked so a
//!   failure can name the first differing record.
//! * The output-changing presets (`lr:hq-fast`, `lr:hq-balanced`) are pinned to golden hashes
//!   (marked `rammap`). After an intended change, `RAMMAP_BLESS=1 cargo test --release --test
//!   regression` prints the new lines for `expected.tsv`.
//! * A raw (memory-mapped) index and BAM input must give the same output as FASTA/FASTQ.
//!
//! Hashes are 64-bit FNV-1a over the output with `@PG` lines removed.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn tmp_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("rammap-regression-{}-{}", std::process::id(), name));
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run rammap and return stdout without `@PG` lines.
fn run(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_rammap")).args(args).output().expect("failed to run rammap");
    assert!(out.status.success(), "rammap {:?} failed:\n{}", args, String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).expect("non-UTF-8 output");
    text.lines().filter(|l| !l.starts_with("@PG")).map(|l| format!("{}\n", l)).collect()
}

fn fnv1a(s: &str) -> String {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in s.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{:016x}", h)
}

/// expected.tsv rows: (case name, source, arguments, hash).
fn expected() -> Vec<(String, String, String, String)> {
    std::fs::read_to_string(data("expected.tsv")).unwrap().lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].into(), f[1].into(), f[2].into(), f[3].into())
        })
        .collect()
}

fn map_args(case_args: &str) -> Vec<String> {
    let mut a: Vec<String> = case_args.split_whitespace().map(String::from).collect();
    a.extend(["-t", "2"].map(String::from));
    a.push(data("ref.fa.gz").to_str().unwrap().into());
    a.push(data("reads.fq.gz").to_str().unwrap().into());
    a
}

fn mask_seq_qual(sam: &str) -> String {
    sam.lines().map(|l| {
        if l.starts_with('@') { return format!("{}\n", l); }
        let mut f: Vec<&str> = l.split('\t').collect();
        if f.len() > 10 { f[9] = "*"; f[10] = "*"; }
        format!("{}\n", f.join("\t"))
    }).collect()
}

#[test]
fn outputs_match_expected() {
    let bless = std::env::var_os("RAMMAP_BLESS").is_some();
    let mut failures = Vec::new();
    for (name, source, args, hash) in expected() {
        let a = map_args(&args);
        let out = run(&a.iter().map(String::as_str).collect::<Vec<_>>());
        let got = fnv1a(&out);
        if bless && source == "rammap" {
            println!("{}\t{}\t{}\t{}", name, source, args, got);
            continue;
        }
        if got != hash {
            let path = tmp_dir("out").join(format!("{}.out", name));
            std::fs::write(&path, &out).unwrap();
            let mut msg = format!("{} ({}; expected from {}): hash {} != {}; output saved to {}", name, args, source, got, hash, path.display());
            if name == "lrhq_sam" {
                let mut exp = String::new();
                flate2::read::GzDecoder::new(std::fs::File::open(data("lrhq.minimap2.sam.gz")).unwrap())
                    .read_to_string(&mut exp).unwrap();
                let masked = mask_seq_qual(&out);
                if let Some((i, (e, g))) = exp.lines().zip(masked.lines()).enumerate().find(|(_, (e, g))| e != g) {
                    msg += &format!("\n  first difference at line {}:\n  minimap2: {}\n  rammap:   {}", i + 1, e, g);
                } else if exp.lines().count() != masked.lines().count() {
                    msg += &format!("\n  line counts differ: minimap2 {} rammap {}", exp.lines().count(), masked.lines().count());
                } else {
                    msg += "\n  (records match with SEQ/QUAL masked: SEQ or QUAL differs)";
                }
            }
            failures.push(msg);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn raw_index_gives_same_output() {
    let dir = tmp_dir("raw");
    let raw = dir.join("ref.raw");
    let raw_s = raw.to_str().unwrap();
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    let (reference, reads) = (reference.to_str().unwrap(), reads.to_str().unwrap());
    run(&["-x", "lr:hq", "-d", raw_s, "--idx-raw", reference]);
    for preset in ["lr:hq", "lr:hq-fast"] {
        let direct = run(&["-ax", preset, "-t", "2", reference, reads]);
        let mapped = run(&["-ax", preset, "-t", "2", raw_s, reads]);
        assert!(direct == mapped, "{}: raw index output differs from FASTA-built index", preset);
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// Write the reads as an unaligned BAM (as a basecaller would), to check BAM input.
fn write_ubam(fq_gz: &Path, out: &Path) {
    let mut fq = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(fq_gz).unwrap()).read_to_string(&mut fq).unwrap();
    let mut raw: Vec<u8> = Vec::new();
    let header = "@HD\tVN:1.6\tSO:unknown\n";
    raw.extend_from_slice(b"BAM\x01");
    raw.extend_from_slice(&(header.len() as i32).to_le_bytes());
    raw.extend_from_slice(header.as_bytes());
    raw.extend_from_slice(&0i32.to_le_bytes()); // no reference sequences
    let lines: Vec<&str> = fq.lines().collect();
    for rec in lines.chunks(4) {
        let (name, seq, qual) = (&rec[0][1..], rec[1].as_bytes(), rec[3].as_bytes());
        let mut r: Vec<u8> = Vec::new();
        r.extend_from_slice(&(-1i32).to_le_bytes()); // refID
        r.extend_from_slice(&(-1i32).to_le_bytes()); // pos
        r.push((name.len() + 1) as u8); // l_read_name
        r.push(0); // mapq
        r.extend_from_slice(&4680u16.to_le_bytes()); // bin
        r.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op
        r.extend_from_slice(&4u16.to_le_bytes()); // flag: unmapped
        r.extend_from_slice(&(seq.len() as i32).to_le_bytes());
        r.extend_from_slice(&(-1i32).to_le_bytes()); // next refID
        r.extend_from_slice(&(-1i32).to_le_bytes()); // next pos
        r.extend_from_slice(&0i32.to_le_bytes()); // tlen
        r.extend_from_slice(name.as_bytes());
        r.push(0);
        let code = |c: u8| -> u8 { match c { b'A' => 1, b'C' => 2, b'G' => 4, b'T' => 8, _ => 15 } };
        for pair in seq.chunks(2) {
            r.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |&c| code(c)));
        }
        r.extend(qual.iter().map(|&q| q - 33));
        raw.extend_from_slice(&(r.len() as i32).to_le_bytes());
        raw.extend_from_slice(&r);
    }
    write_bgzf(&raw, out);
}

/// BGZF: deflate blocks of < 64 KiB, each with the BC extra field, then the EOF block.
fn write_bgzf(raw: &[u8], out: &Path) {
    let mut f = std::fs::File::create(out).unwrap();
    for chunk in raw.chunks(60_000) {
        let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
        enc.write_all(chunk).unwrap();
        let cdata = enc.finish().unwrap();
        let mut crc = flate2::Crc::new();
        crc.update(chunk);
        let bsize = (cdata.len() + 25) as u16;
        f.write_all(&[31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, b'B', b'C', 2, 0]).unwrap();
        f.write_all(&bsize.to_le_bytes()).unwrap();
        f.write_all(&cdata).unwrap();
        f.write_all(&crc.sum().to_le_bytes()).unwrap();
        f.write_all(&(chunk.len() as u32).to_le_bytes()).unwrap();
    }
    f.write_all(&[31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]).unwrap();
}

#[test]
fn bam_input_gives_same_output() {
    let dir = tmp_dir("bam");
    let bam = dir.join("reads.bam");
    write_ubam(&data("reads.fq.gz"), &bam);
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    for preset in ["lr:hq", "lr:hq-fast"] {
        let from_fq = run(&["-ax", preset, "-t", "2", reference.to_str().unwrap(), reads.to_str().unwrap()]);
        let from_bam = run(&["-ax", preset, "-t", "2", reference.to_str().unwrap(), bam.to_str().unwrap()]);
        assert!(from_fq == from_bam, "{}: BAM input output differs from FASTQ input", preset);
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn bgzipped_fastq_gives_same_output() {
    let dir = tmp_dir("bgzf");
    let mut fq = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(data("reads.fq.gz")).unwrap()).read_to_string(&mut fq).unwrap();
    let bgz = dir.join("reads.fq.bgz.gz");
    write_bgzf(fq.as_bytes(), &bgz);
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    for preset in ["lr:hq", "lr:hq-fast"] {
        let a = run(&["-ax", preset, "-t", "2", reference.to_str().unwrap(), reads.to_str().unwrap()]);
        let b = run(&["-ax", preset, "-t", "2", reference.to_str().unwrap(), bgz.to_str().unwrap()]);
        assert!(a == b, "{}: bgzipped FASTQ output differs from gzipped", preset);
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// `--sort -o x.bam`: same records as the SAM output, in samtools sort's coordinate order
/// (reference id with unmapped last, position, reverse strand after forward, ties in input
/// order), plus a BAI index.
#[test]
fn sorted_bam_output() {
    let dir = tmp_dir("sort");
    let bam = dir.join("out.bam");
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    let (reference, reads) = (reference.to_str().unwrap(), reads.to_str().unwrap());
    for preset in ["lr:hq", "lr:hq-fast"] {
        let sam = run(&["-ax", preset, "-t", "2", reference, reads]);
        let names: Vec<&str> = sam.lines().filter(|l| l.starts_with("@SQ")).map(|l| &l.split('\t').nth(1).unwrap()[3..]).collect();
        let tid = |n: &str| if n == "*" { u32::MAX } else { names.iter().position(|&x| x == n).unwrap() as u32 };
        let mut expect: Vec<(u64, String, u16)> = sam.lines().filter(|l| !l.starts_with('@')).map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let flag: u16 = f[1].parse().unwrap();
            let pos: i64 = f[3].parse::<i64>().unwrap() - 1;
            ((tid(f[2]) as u64) << 32 | ((pos + 1) as u64) << 1 | ((flag >> 4) & 1) as u64, f[0].to_string(), flag)
        }).collect();
        expect.sort_by_key(|e| e.0);

        run(&["-ax", preset, "-t", "2", "--sort", "--sort-mem", "100K", "-o", bam.to_str().unwrap(), reference, reads]);
        let mut raw = Vec::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(&bam).unwrap()).read_to_end(&mut raw).unwrap();
        assert_eq!(&raw[..4], b"BAM\x01");
        let l_text = i32::from_le_bytes(raw[4..8].try_into().unwrap()) as usize;
        let mut p = 8 + l_text;
        let n_ref = i32::from_le_bytes(raw[p..p + 4].try_into().unwrap()) as usize;
        p += 4;
        for _ in 0..n_ref { p += 4 + i32::from_le_bytes(raw[p..p + 4].try_into().unwrap()) as usize + 4; }
        let mut got: Vec<(u64, String, u16)> = Vec::new();
        while p < raw.len() {
            let n = u32::from_le_bytes(raw[p..p + 4].try_into().unwrap()) as usize;
            let r = &raw[p + 4..p + 4 + n];
            let tid = i32::from_le_bytes(r[0..4].try_into().unwrap()) as u32;
            let pos = i32::from_le_bytes(r[4..8].try_into().unwrap()) as i64;
            let flag = u16::from_le_bytes(r[14..16].try_into().unwrap());
            let name = String::from_utf8(r[32..32 + r[8] as usize - 1].to_vec()).unwrap();
            got.push(((tid as u64) << 32 | ((pos + 1) as u64) << 1 | ((flag >> 4) & 1) as u64, name, flag));
            p += 4 + n;
        }
        assert_eq!(got, expect, "{}: sorted BAM records differ from the sorted SAM output", preset);
        let bai = std::fs::read(dir.join("out.bam.bai")).unwrap();
        assert_eq!(&bai[..4], b"BAI\x01");
    }
    std::fs::remove_dir_all(&dir).ok();
}
