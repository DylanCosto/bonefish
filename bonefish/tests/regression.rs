//! End-to-end regression tests on a small synthetic ONT-like dataset (`tests/data`, made by
//! `scripts/make_regression_data.py`).
//!
//! * The output (SAM, PAF, and the --cs / --MD --eqx variants) is pinned to golden hashes in
//!   `tests/data/expected.tsv`, so speed-only changes can be checked for identical output.
//!   Accuracy-changing work is judged on the GIAB benchmarks instead (see the project README);
//!   after an intended output change, `BONEFISH_BLESS=1 cargo test --release --test regression`
//!   prints the new lines for `expected.tsv`.
//! * A raw (memory-mapped) index and BAM input must give the same output as FASTA/FASTQ.
//! * Input handling: dorado BAM tags and read groups with -y, wrapped FASTQ, malformed input
//!   failing with a non-zero exit, `-o x.bam` as --sort, options after the file arguments.
//!
//! Hashes are 64-bit FNV-1a over the output with `@PG` lines removed.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn tmp_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("bonefish-regression-{}-{}", std::process::id(), name));
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Run bonefish and return stdout without `@PG` lines.
fn run(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_bonefish")).args(args).output().expect("failed to run bonefish");
    assert!(out.status.success(), "bonefish {:?} failed:\n{}", args, String::from_utf8_lossy(&out.stderr));
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

#[test]
fn outputs_match_expected() {
    let bless = std::env::var_os("BONEFISH_BLESS").is_some();
    let mut failures = Vec::new();
    for (name, source, args, hash) in expected() {
        let a = map_args(&args);
        let out = run(&a.iter().map(String::as_str).collect::<Vec<_>>());
        let got = fnv1a(&out);
        if bless && source == "bonefish" {
            println!("{}\t{}\t{}\t{}", name, source, args, got);
            continue;
        }
        if got != hash {
            let path = tmp_dir("out").join(format!("{}.out", name));
            std::fs::write(&path, &out).unwrap();
            let msg = format!("{} ({}; expected from {}): hash {} != {}; output saved to {}", name, args, source, got, hash, path.display());
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
    run(&["-d", raw_s, "--idx-raw", reference]);
    let direct = run(&["-a", "-t", "2", reference, reads]);
    let mapped = run(&["-a", "-t", "2", raw_s, reads]);
    assert!(direct == mapped, "raw index output differs from FASTA-built index");
    std::fs::remove_dir_all(&dir).ok();
}

/// Write the reads as an unaligned BAM (as a basecaller would), to check BAM input.
fn write_ubam(fq_gz: &Path, out: &Path) {
    write_ubam_with(fq_gz, out, "", |_| Vec::new());
}

/// `write_ubam` with extra header lines and per-record binary aux data (record index -> bytes).
fn write_ubam_with(fq_gz: &Path, out: &Path, header_extra: &str, aux: impl Fn(usize) -> Vec<u8>) {
    let mut fq = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(fq_gz).unwrap()).read_to_string(&mut fq).unwrap();
    let mut raw: Vec<u8> = Vec::new();
    let header = format!("@HD\tVN:1.6\tSO:unknown\n{header_extra}");
    raw.extend_from_slice(b"BAM\x01");
    raw.extend_from_slice(&(header.len() as i32).to_le_bytes());
    raw.extend_from_slice(header.as_bytes());
    raw.extend_from_slice(&0i32.to_le_bytes()); // no reference sequences
    let lines: Vec<&str> = fq.lines().collect();
    for (i, rec) in lines.chunks(4).enumerate() {
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
        r.extend_from_slice(&aux(i));
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
    let from_fq = run(&["-a", "-t", "2", reference.to_str().unwrap(), reads.to_str().unwrap()]);
    let from_bam = run(&["-a", "-t", "2", reference.to_str().unwrap(), bam.to_str().unwrap()]);
    assert!(from_fq == from_bam, "BAM input output differs from FASTQ input");
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
    let a = run(&["-a", "-t", "2", reference.to_str().unwrap(), reads.to_str().unwrap()]);
    let b = run(&["-a", "-t", "2", reference.to_str().unwrap(), bgz.to_str().unwrap()]);
    assert!(a == b, "bgzipped FASTQ output differs from gzipped");
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
    {
        let sam = run(&["-a", "-t", "2", reference, reads]);
        let names: Vec<&str> = sam.lines().filter(|l| l.starts_with("@SQ")).map(|l| &l.split('\t').nth(1).unwrap()[3..]).collect();
        let tid = |n: &str| if n == "*" { u32::MAX } else { names.iter().position(|&x| x == n).unwrap() as u32 };
        let mut expect: Vec<(u64, String, u16)> = sam.lines().filter(|l| !l.starts_with('@')).map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let flag: u16 = f[1].parse().unwrap();
            let pos: i64 = f[3].parse::<i64>().unwrap() - 1;
            ((tid(f[2]) as u64) << 32 | ((pos + 1) as u64) << 1 | ((flag >> 4) & 1) as u64, f[0].to_string(), flag)
        }).collect();
        expect.sort_by_key(|e| e.0);

        run(&["-a", "-t", "2", "--sort", "--sort-mem", "100K", "-o", bam.to_str().unwrap(), reference, reads]);
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
        assert_eq!(got, expect, "sorted BAM records differ from the sorted SAM output");
        let bai = std::fs::read(dir.join("out.bam.bai")).unwrap();
        assert_eq!(&bai[..4], b"BAI\x01");
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// Run bonefish expecting failure; returns stderr.
fn run_fails(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_bonefish")).args(args).output().expect("failed to run bonefish");
    assert!(!out.status.success(), "bonefish {:?} should have failed", args);
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn reads_text() -> String {
    let mut fq = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(data("reads.fq.gz")).unwrap()).read_to_string(&mut fq).unwrap();
    fq
}

/// dorado BAM input with -y: the read tags (here RG, qs, MM/ML) and the @RG header line reach
/// the output; an old alignment tag (NM) is replaced, not duplicated; -R replaces the read group.
#[test]
fn bam_tags_and_read_groups_with_y() {
    let dir = tmp_dir("tags");
    let bam = dir.join("reads.bam");
    let rg = "@RG\tID:run1_dna_r10.4.1_e8.2_400bps_sup@v5.0.0\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0\n";
    write_ubam_with(&data("reads.fq.gz"), &bam, rg, |i| {
        let mut a = Vec::new();
        a.extend_from_slice(b"RGZrun1_dna_r10.4.1_e8.2_400bps_sup@v5.0.0\0");
        a.extend_from_slice(b"qsi"); a.extend_from_slice(&(10 + i as i32 % 20).to_le_bytes());
        a.extend_from_slice(b"MMZC+m?,0,1;\0");
        a.extend_from_slice(b"MLBC"); a.extend_from_slice(&2u32.to_le_bytes()); a.extend_from_slice(&[200, 13]);
        a.extend_from_slice(b"NMi"); a.extend_from_slice(&999i32.to_le_bytes());
        a
    });
    let reference = data("ref.fa.gz");
    let (reference, bam_s) = (reference.to_str().unwrap(), bam.to_str().unwrap());
    let sam = run(&["-a", "-y", "-t", "2", reference, bam_s]);
    assert!(sam.lines().any(|l| l == rg.trim_end()), "input @RG line missing from the header");
    let recs: Vec<&str> = sam.lines().filter(|l| !l.starts_with('@')).collect();
    assert!(!recs.is_empty());
    for l in &recs {
        let tags: Vec<&str> = l.split('\t').skip(11).collect();
        for t in ["RG:Z:run1_dna_r10.4.1_e8.2_400bps_sup@v5.0.0", "MM:Z:C+m?,0,1;", "ML:B:C,200,13"] {
            assert!(tags.contains(&t), "tag {t} missing: {l}");
        }
        assert!(tags.iter().any(|t| t.starts_with("qs:i:")), "qs missing: {l}");
        assert!(!tags.contains(&"NM:i:999"), "stale NM copied: {l}");
        let mut names: Vec<&str> = tags.iter().map(|t| &t[..2]).collect();
        let n = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), n, "duplicate tags: {l}");
    }
    let sam_r = run(&["-a", "-y", "-t", "2", "-R", "@RG\\tID:mine\\tSM:s1", reference, bam_s]);
    assert!(!sam_r.contains("run1_dna"), "-R must replace the input read group");
    assert!(sam_r.lines().filter(|l| !l.starts_with('@')).all(|l| l.contains("\tRG:Z:mine")));
    std::fs::remove_dir_all(&dir).ok();
}

/// Reads given on stdin or through a pipe (`<(zcat reads.fq.gz)`) give the same output as the
/// file, and a BAM streamed this way keeps its `@RG` header lines (they are read before mapping,
/// so the bytes consumed for that must be replayed).
#[cfg(unix)]
#[test]
fn stdin_and_pipe_input() {
    use std::process::Stdio;
    let dir = tmp_dir("pipe");
    let bam = dir.join("reads.bam");
    let rg = "@RG\tID:run1\tPL:ONT\n";
    write_ubam_with(&data("reads.fq.gz"), &bam, rg, |_| b"RGZrun1\0".to_vec());
    let reference = data("ref.fa.gz");
    let reference = reference.to_str().unwrap();
    let reads = data("reads.fq.gz");
    for (input, extra) in [(&reads, &[][..]), (&bam, &["-y"][..])] {
        assert!(extra.is_empty() || run(&["-a", "-y", "-t", "2", reference, bam.to_str().unwrap()]).contains(rg.trim_end()));
        let mut args = vec!["-a", "-t", "2"];
        args.extend_from_slice(extra);
        let expected = run(&[&args[..], &[reference, input.to_str().unwrap()]].concat());
        // stdin
        let out = Command::new(env!("CARGO_BIN_EXE_bonefish")).args(&args).args([reference, "-"])
            .stdin(Stdio::from(std::fs::File::open(input).unwrap())).output().unwrap();
        assert!(out.status.success(), "stdin run failed: {}", String::from_utf8_lossy(&out.stderr));
        let text: String = String::from_utf8(out.stdout).unwrap().lines()
            .filter(|l| !l.starts_with("@PG")).map(|l| format!("{}\n", l)).collect();
        assert!(text == expected, "stdin output differs for {}", input.display());
        // named pipe
        let fifo = dir.join("fifo");
        std::fs::remove_file(&fifo).ok();
        assert!(Command::new("mkfifo").arg(&fifo).status().unwrap().success());
        let (src, dst) = (input.clone(), fifo.clone());
        let writer = std::thread::spawn(move || {
            std::io::copy(&mut std::fs::File::open(src).unwrap(), &mut std::fs::File::create(dst).unwrap()).ok();
        });
        let out = run(&[&args[..], &[reference, fifo.to_str().unwrap()]].concat());
        writer.join().unwrap();
        assert!(out == expected, "pipe output differs for {}", input.display());
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// Wrapped (multi-line) FASTQ gives the same alignments as the four-line form.
#[test]
fn multiline_fastq_gives_same_output() {
    let dir = tmp_dir("wrap");
    let fq = reads_text();
    let mut wrapped = String::new();
    for rec in fq.lines().collect::<Vec<_>>().chunks(4) {
        wrapped.push_str(rec[0]); wrapped.push('\n');
        for c in rec[1].as_bytes().chunks(60) { wrapped.push_str(std::str::from_utf8(c).unwrap()); wrapped.push('\n'); }
        wrapped.push_str("+\n");
        for c in rec[3].as_bytes().chunks(60) { wrapped.push_str(std::str::from_utf8(c).unwrap()); wrapped.push('\n'); }
    }
    let p = dir.join("wrapped.fq");
    std::fs::write(&p, wrapped).unwrap();
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    let a = run(&["-a", "-t", "2", reference.to_str().unwrap(), reads.to_str().unwrap()]);
    let b = run(&["-a", "-t", "2", reference.to_str().unwrap(), p.to_str().unwrap()]);
    assert!(a == b, "wrapped FASTQ output differs");
    std::fs::remove_dir_all(&dir).ok();
}

/// Malformed input fails (non-zero exit) instead of producing silently incomplete output.
#[test]
fn malformed_input_fails() {
    let dir = tmp_dir("bad");
    let reference = data("ref.fa.gz");
    let reference = reference.to_str().unwrap();
    let cases: [(&str, &str); 5] = [
        ("qual_mismatch.fq", "@r1\nACGTACGTAC\n+\nIIIII\n@r2\nACGT\n+\nIIII\n"),
        ("truncated.fq", "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nACGTACGT\n"),
        ("garbage.txt", "this is not a sequence file\n"),
        ("sam.txt", "@HD\tVN:1.6\n@SQ\tSN:x\tLN:10\n"),
        ("longname.fq", &format!("@{}\nACGT\n+\nIIII\n", "n".repeat(255))),
    ];
    for (name, text) in cases {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        let err = run_fails(&["-a", "-t", "2", reference, p.to_str().unwrap()]);
        assert!(err.contains(name), "error should name the file: {err}");
    }
    let err = run_fails(&["-a", "-k", "40", reference, data("reads.fq.gz").to_str().unwrap()]);
    assert!(err.contains("-k"), "{err}");
    std::fs::remove_dir_all(&dir).ok();
}

/// `-o x.bam` alone writes the sorted BAM (as --sort), options may follow the file arguments,
/// and a failing run (missing reference) leaves an existing output file untouched.
#[test]
fn bam_output_by_extension_and_option_order() {
    let dir = tmp_dir("ext");
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    let (reference, reads) = (reference.to_str().unwrap(), reads.to_str().unwrap());
    let a = dir.join("a.bam");
    let b = dir.join("b.bam");
    run(&["-t", "2", "-o", a.to_str().unwrap(), reference, reads]);
    run(&["-a", "--sort", "-t", "2", "-o", b.to_str().unwrap(), reference, reads]);
    let records = |p: &Path| {
        let mut raw = Vec::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(p).unwrap()).read_to_end(&mut raw).unwrap();
        let l_text = i32::from_le_bytes(raw[4..8].try_into().unwrap()) as usize;
        raw.split_off(8 + l_text) // references and records; the header text holds @PG CL
    };
    assert!(records(&a) == records(&b), "-o x.bam differs from --sort");
    assert!(dir.join("a.bam.bai").exists());
    let before = run(&["-a", "-t", "2", reference, reads]);
    let after = run(&[reference, reads, "-a", "-t", "2"]);
    assert!(before == after, "options after the file arguments changed the output");
    let keep = dir.join("keep.bam");
    std::fs::write(&keep, b"precious").unwrap();
    run_fails(&["-o", keep.to_str().unwrap(), "missing_reference.fa", reads]);
    assert_eq!(std::fs::read(&keep).unwrap(), b"precious");
    std::fs::remove_dir_all(&dir).ok();
}

/// Every SIMD path gives the same output: the kernels picked at run time (AVX-512 where the CPU
/// has it, e.g. on CI runners) against forced AVX2, SSE4.1 and scalar code, for SAM and PAF with
/// --cs. On arm64 only NEON exists, so this is x86_64-only.
#[cfg(target_arch = "x86_64")]
#[test]
fn simd_paths_give_same_output() {
    let reference = data("ref.fa.gz");
    let reads = data("reads.fq.gz");
    let (reference, reads) = (reference.to_str().unwrap(), reads.to_str().unwrap());
    let has_avx512 = std::arch::is_x86_feature_detected!("avx512bw");
    eprintln!("simd_paths_give_same_output: CPU has AVX-512BW: {has_avx512}");
    for args in [vec!["-a"], vec!["-c", "--cs"]] {
        let run_env = |env: Option<&str>| {
            let mut a: Vec<&str> = args.clone();
            a.extend(["-t", "2", reference, reads]);
            let mut c = Command::new(env!("CARGO_BIN_EXE_bonefish"));
            c.args(&a);
            if let Some(e) = env { c.env(e, "1"); }
            let out = c.output().expect("failed to run bonefish");
            assert!(out.status.success(), "{env:?} failed: {}", String::from_utf8_lossy(&out.stderr));
            String::from_utf8(out.stdout).unwrap().lines().filter(|l| !l.starts_with("@PG")).map(|l| format!("{l}\n")).collect::<String>()
        };
        let auto = run_env(None);
        for env in ["BONEFISH_FORCE_AVX2", "BONEFISH_FORCE_SSE", "BONEFISH_FORCE_SCALAR"] {
            assert!(run_env(Some(env)) == auto, "{env} output differs from the run-time kernels ({args:?}; AVX-512: {has_avx512})");
        }
    }
}
