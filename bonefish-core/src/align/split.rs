//! Split-index mode: multi-part reference processing and merge.
//!
//! When the reference is too large for a single index, it is split into parts.
//! Each part is aligned independently, results are written to temporary files,
//! and `split_merge` combines them with re-filtering.

use std::io::{self, BufReader, BufWriter, Read, Write};
use std::fs::File;
use crate::align::index::{Index, TargetSequence};
use crate::align::pipeline::{AlnResult, OutputConfig, ReadInfo, refilter_merged_results, format_output};
use crate::align::map::MapOptions;
use crate::align::stats::AlignmentStats;

/// Create a temp file for one index part, write the header.
/// Header format:
///   k: u32, n_seq: u32, then per-seq: name_len(u32), name([u8]), seq_len(u32)
pub fn split_init(prefix: &str, part_index: usize, mi: &Index) -> io::Result<BufWriter<File>> {
    let path = format!("{}.{:04}.tmp", prefix, part_index);
    let f = File::create(&path).map_err(|e| {
        io::Error::new(e.kind(), format!("failed to write to temporary file '{}': {}", path, e))
    })?;
    let mut w = BufWriter::new(f);
    w.write_all(&(mi.kmer_size as u32).to_le_bytes())?;
    w.write_all(&(mi.seqs.len() as u32).to_le_bytes())?;
    for s in &mi.seqs {
        let name_bytes = s.name.as_bytes();
        w.write_all(&(name_bytes.len() as u32).to_le_bytes())?;
        w.write_all(name_bytes)?;
        w.write_all(&(s.len as u32).to_le_bytes())?;
    }
    Ok(w)
}

/// Write one query's results to a temp file.
/// Per-query format: n_results(i32), rep_len(i32),
/// then for each result: bincode-serialized AlnResult.
pub fn split_write_query(
    w: &mut impl Write,
    results: &[AlnResult],
    rep_len: i32,
) -> io::Result<()> {
    w.write_all(&(results.len() as i32).to_le_bytes())?;
    w.write_all(&rep_len.to_le_bytes())?;
    for r in results {
        let encoded = bincode::serialize(r).map_err(|e| {
            io::Error::other(format!("bincode serialize error: {}", e))
        })?;
        w.write_all(&(encoded.len() as u32).to_le_bytes())?;
        w.write_all(&encoded)?;
    }
    Ok(())
}

/// Read headers from all part files, build global seqs list + rid_shifts.
/// Returns (k, all_seqs, rid_shifts, readers).
#[allow(clippy::type_complexity)]
pub fn split_merge_prep(
    prefix: &str,
    n_parts: usize,
) -> io::Result<(usize, Vec<TargetSequence>, Vec<usize>, Vec<BufReader<File>>)> {
    let mut readers = Vec::with_capacity(n_parts);
    let mut n_seq_parts = Vec::with_capacity(n_parts);
    let mut k_val = 0usize;

    for i in 0..n_parts {
        let path = format!("{}.{:04}.tmp", prefix, i);
        let f = File::open(&path).map_err(|e| {
            io::Error::new(e.kind(), format!("failed to open temporary file '{}': {}", path, e))
        })?;
        let mut r = BufReader::new(f);

        let mut buf4 = [0u8; 4];
        r.read_exact(&mut buf4)?;
        let k = u32::from_le_bytes(buf4) as usize;
        if i == 0 { k_val = k; }

        r.read_exact(&mut buf4)?;
        let n_seq = u32::from_le_bytes(buf4) as usize;
        n_seq_parts.push(n_seq);
        readers.push(r);
    }

    // Compute rid_shifts
    let mut rid_shifts = vec![0usize; n_parts];
    for i in 1..n_parts {
        rid_shifts[i] = rid_shifts[i - 1] + n_seq_parts[i - 1];
    }

    // Read all sequence headers
    let total_seqs: usize = n_seq_parts.iter().sum();
    let mut all_seqs = Vec::with_capacity(total_seqs);
    for (i, reader) in readers.iter_mut().enumerate() {
        for _ in 0..n_seq_parts[i] {
            let mut buf4 = [0u8; 4];
            reader.read_exact(&mut buf4)?;
            let name_len = u32::from_le_bytes(buf4) as usize;
            let mut name_buf = vec![0u8; name_len];
            reader.read_exact(&mut name_buf)?;
            let name = String::from_utf8(name_buf).map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("invalid UTF-8 in seq name: {}", e))
            })?;
            reader.read_exact(&mut buf4)?;
            let seq_len = u32::from_le_bytes(buf4) as usize;
            all_seqs.push(TargetSequence {
                name,
                len: seq_len,
                offset: 0,
                is_alt: false,
            });
        }
    }

    Ok((k_val, all_seqs, rid_shifts, readers))
}

/// Read one query's results from a part file. Returns None on EOF.
pub fn split_read_query(r: &mut impl Read) -> io::Result<Option<(Vec<AlnResult>, i32)>> {
    let mut buf4 = [0u8; 4];
    match r.read_exact(&mut buf4) {
        Ok(()) => {},
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e),
    }
    let n_results = i32::from_le_bytes(buf4);
    r.read_exact(&mut buf4)?;
    let rep_len = i32::from_le_bytes(buf4);

    let mut results = Vec::with_capacity(n_results as usize);
    for _ in 0..n_results {
        r.read_exact(&mut buf4)?;
        let enc_len = u32::from_le_bytes(buf4) as usize;
        let mut enc_buf = vec![0u8; enc_len];
        r.read_exact(&mut enc_buf)?;
        let aln: AlnResult = bincode::deserialize(&enc_buf).map_err(|e| {
            io::Error::other(format!("bincode deserialize error: {}", e))
        })?;
        results.push(aln);
    }

    Ok(Some((results, rep_len)))
}

/// Delete temp files.
pub fn split_rm_tmp(prefix: &str, n_parts: usize) {
    for i in 0..n_parts {
        let path = format!("{}.{:04}.tmp", prefix, i);
        let _ = std::fs::remove_file(&path);
    }
}

/// Configuration for split-index merge.
pub struct MergeConfig<'a> {
    pub prefix: &'a str,
    pub n_parts: usize,
    pub k: usize,
    pub w: usize,
    pub query_path: &'a str,
    pub sam_rg: Option<&'a str>,
    pub copy_comment: bool,
}

/// Merge results from all split parts and output.
/// Merge alignment results from multiple index parts and output.
pub fn merge_split_results(
    config: &MergeConfig,
    opt: &MapOptions,
    out_cfg: &OutputConfig,
    handle: &mut dyn Write,
) -> anyhow::Result<AlignmentStats> {
    let prefix = config.prefix;
    let n_parts = config.n_parts;
    let k = config.k;
    let w = config.w;
    let query_path = config.query_path;
    let sam_rg = config.sam_rg;
    let copy_comment = config.copy_comment;
    // 1. Read all temp file headers, build global seqs + rid_shifts
    let (k_from_file, all_seqs, rid_shifts, mut readers) = split_merge_prep(prefix, n_parts)
        .map_err(|e| anyhow::anyhow!("Split merge prep failed: {}", e))?;

    let _ = k_from_file; // k already known from caller

    // 2. Build header-only index for output formatting
    let global_mi = Index::header_only(k, w, all_seqs);

    // 3. Write SAM header if SAM output
    if out_cfg.output_sam {
        writeln!(handle, "@HD\tVN:1.6\tSO:unsorted\tGO:query")?;
        for s in &global_mi.seqs {
            writeln!(handle, "@SQ\tSN:{}\tLN:{}", s.name, s.len)?;
        }
        let version = env!("CARGO_PKG_VERSION");
        let cmd_line = std::env::args().collect::<Vec<String>>().join(" ");
        writeln!(handle, "@PG\tID:bonefish\tPN:bonefish\tVN:{}\tCL:{}", version, cmd_line)?;
        if let Some(rg_line) = sam_rg {
            let expanded = rg_line.replace("\\t", "\t");
            writeln!(handle, "{}", expanded)?;
        }
    }

    // 4. Re-open query files for output
    let query_reader = crate::fasta::open(query_path)
        .map_err(|e| anyhow::anyhow!("Error re-opening query for merge: {}", e))?;

    let mut total_stats = AlignmentStats::default();
    let mut record_iter = query_reader.records();
    let mut output_buffer = String::new();

    {
        // Single-end merge
        loop {
            let rec = match record_iter.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => { eprintln!("Warning: Error reading record: {}", e); continue; },
                None => break,
            };

            let qname = rec.name().to_string();
            let qseq = rec.sequence().to_vec();
            let qual = rec.quality().map(|qs| String::from_utf8_lossy(qs).to_string());
            let comment = if copy_comment { rec.description().map(|s| s.to_string()) } else { None };
            let qlen = qseq.len();

            // Read results from all parts
            let mut all_results: Vec<AlnResult> = Vec::new();
            let mut max_rep_len = 0i32;

            for (j, reader) in readers.iter_mut().enumerate() {
                if let Some((mut results, rep_len)) = split_read_query(reader)
                    .map_err(|e| anyhow::anyhow!("Error reading split part {}: {}", j, e))? {
                    for r in results.iter_mut() {
                        r.ref_id += rid_shifts[j];
                    }
                    all_results.extend(results);
                    if rep_len > max_rep_len { max_rep_len = rep_len; }
                }
            }

            // Refilter (uses max_rep_len for MAPQ, matching mm_set_mapq2 in merge_hits)
            let mut pq = refilter_merged_results(all_results, opt, k, qlen, max_rep_len, out_cfg.do_cigar);
            // merge doesn't update rep_len, so rl:i: tag is 0
            pq.rep_len = 0;

            // Format output
            output_buffer.clear();
            let ri = ReadInfo { qname: &qname, qseq: &qseq, qual: qual.as_deref(), comment: comment.as_deref() };
            format_output(&mut output_buffer, opt, &global_mi, &ri, &pq, out_cfg);
            write!(handle, "{}", output_buffer)?;
            total_stats = total_stats + pq.stats;
        }
    }

    // Clean up temp files
    split_rm_tmp(prefix, n_parts);

    Ok(total_stats)
}
