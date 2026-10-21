//! Annotated splice junctions (`--junc-bed`), used by the spliced-alignment DP.
//!
//! Loaded from a BED file (BED12 transcripts: introns between blocks; BED6: one
//! interval per junction) and queried per DP segment, marking annotated donor and
//! acceptor sites so the splice kernel gives them `--junc-bonus` (and unannotated
//! ones `--junc-pen`). Not part of the index.

use std::collections::HashMap;
use std::io::{self, BufRead, BufReader};
use std::fs::File;

/// One annotated intron on a target sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Intron {
    start: i32,
    end: i32,
    strand: i32, // +1, -1, or 0
}

/// Annotated introns per target sequence (indexed by reference id), each list
/// sorted by (start, end) with duplicates removed.
#[derive(Debug, Default)]
pub struct JunctionDb {
    introns: Vec<Vec<Intron>>,
}

impl JunctionDb {
    /// Load junctions from a BED file. `names` maps reference names to ids; lines
    /// on other references are skipped.
    pub fn load_bed(path: &str, names: &HashMap<String, usize>, n_seq: usize) -> io::Result<JunctionDb> {
        let reader = BufReader::new(File::open(path)?);
        let mut introns: Vec<Vec<Intron>> = vec![Vec::new(); n_seq];

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 { continue; }
            let Some(&rid) = names.get(fields[0]) else { continue };
            let (Ok(st), Ok(en)) = (fields[1].parse::<i32>(), fields[2].parse::<i32>()) else { continue };
            if st < 0 || st >= en { continue; }
            let strand = match fields.get(5) {
                Some(&"+") => 1,
                Some(&"-") => -1,
                _ => 0,
            };

            if fields.len() >= 12 {
                // BED12: the introns are the gaps between consecutive blocks.
                let Ok(n_blk) = fields[9].parse::<usize>() else { continue };
                if n_blk < 2 { continue; }
                let parse = |s: &str| -> Vec<i32> { s.split(',').filter(|x| !x.is_empty()).filter_map(|x| x.parse().ok()).collect() };
                let (sizes, starts) = (parse(fields[10]), parse(fields[11]));
                if sizes.len() < n_blk || starts.len() < n_blk { continue; }
                let mut prev_en = st + starts[0] + sizes[0];
                for i in 1..n_blk {
                    let intron_en = st + starts[i];
                    if intron_en > prev_en {
                        introns[rid].push(Intron { start: prev_en, end: intron_en, strand });
                    }
                    prev_en = st + starts[i] + sizes[i];
                }
            } else {
                introns[rid].push(Intron { start: st, end: en, strand });
            }
        }

        let mut total = 0usize;
        for v in &mut introns {
            v.sort_by_key(|x| (x.start, x.end));
            v.dedup_by(|b, a| a.start == b.start && a.end == b.end);
            total += v.len();
        }
        eprintln!("[M::load_bed_junctions] loaded {} non-redundant junctions", total);
        Ok(JunctionDb { introns })
    }

    /// Load junctions from a BED file for the reference sequences of `mi`.
    pub fn load_bed_for_index(path: &str, mi: &crate::align::index::Index) -> io::Result<JunctionDb> {
        let names: HashMap<String, usize> = mi.seqs.iter().enumerate().map(|(i, s)| (s.name.clone(), i)).collect();
        Self::load_bed(path, &names, mi.seqs.len())
    }

    /// Fill `s[..en - st]` with the annotated splice sites of reference `rid` in
    /// `[st, en)` (introns entirely inside the window): bit 0 forward donor, bit 1
    /// forward acceptor, bit 2 reverse acceptor, bit 3 reverse donor.
    pub fn fill(&self, rid: usize, st: i32, en: i32, s: &mut [u8]) {
        let len = (en - st) as usize;
        s[..len].fill(0);
        let Some(v) = self.introns.get(rid) else { return };
        let first = v.partition_point(|x| x.start < st);
        for iv in &v[first..] {
            if iv.start >= en { break; }
            if iv.end > en || iv.strand == 0 { continue; }
            let (d, a) = ((iv.start - st) as usize, (iv.end - 1 - st) as usize);
            if iv.strand > 0 {
                s[d] |= 1;
                s[a] |= 2;
            } else {
                s[d] |= 8;
                s[a] |= 4;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed12_introns_and_fill() {
        let dir = std::env::temp_dir().join(format!("bf_junc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("j.bed");
        // Two blocks [100,110) and [200,220): intron [110,200) on +; a BED6 intron on -.
        std::fs::write(&p, "chr1\t100\t220\ttx\t0\t+\t100\t220\t0\t2\t10,20,\t0,100,\nchr1\t300\t350\tj\t0\t-\nchrX\t1\t5\tx\t0\t+\n").unwrap();
        let names: HashMap<String, usize> = [("chr1".to_string(), 0)].into_iter().collect();
        let db = JunctionDb::load_bed(p.to_str().unwrap(), &names, 1).unwrap();
        let mut s = vec![0u8; 400];
        db.fill(0, 0, 400, &mut s);
        assert_eq!((s[110], s[199], s[300], s[349]), (1, 2, 8, 4));
        assert_eq!(s.iter().filter(|&&b| b != 0).count(), 4);
        // A window that cuts an intron does not mark it.
        let mut s = vec![0u8; 100];
        db.fill(0, 150, 250, &mut s);
        assert!(s.iter().all(|&b| b == 0));
        std::fs::remove_dir_all(&dir).ok();
    }
}
