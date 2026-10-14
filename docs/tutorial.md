# Tutorial: from a dorado BAM to a sorted BAM and methylation calls in a few minutes

This walks through bonefish on ONT's public HG002 demo data (the small dataset ONT ships for its
`wf-human-variation` pipeline): 3,781 R10.4.1 reads over regions of chr6, chr20 and chrM, basecalled
with modified-base models, so every read carries `MM`/`ML` methylation tags and dorado's other tags.
Everything runs on a laptop.

## 1. Install

Download the binary for your platform from the
[latest release](https://github.com/DylanCosto/bonefish/releases/latest) and unpack it:

```bash
tar -xzf bonefish-v*-linux-x86_64.tar.gz    # or -linux-arm64, -macos-arm64
export PATH="$PWD/$(ls -d bonefish-v*/ | head -1):$PATH"
bonefish --version
```

(On macOS, a downloaded binary may need `xattr -d com.apple.quarantine bonefish` once.) Or build
from source: `cargo build --release --features jemalloc`.

You will also want [samtools](https://www.htslib.org/) and, for step 4,
[modkit](https://github.com/nanoporetech/modkit) (`conda install -c bioconda samtools ont-modkit`).

## 2. Get the demo data (142 MB)

```bash
curl -LO https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/hg002%2bmods.v1/wf-human-variation-demo.tar.gz
tar -xzf wf-human-variation-demo.tar.gz
cd wf-human-variation-demo
md5sum -c md5sums.txt
```

`demo.bam` holds the reads (aligned by ONT, which does not matter: bonefish restores each read's
original orientation and skips secondary and supplementary records), and `demo.fasta` the matching
reference (chr6, chr20, chrM of GRCh38).

## 3. Align to a sorted, indexed BAM

```bash
bonefish -t 8 -y -o demo.bonefish.bam demo.fasta demo.bam
```

```
[*] Mapped 3781 reads in 4.233s (1 index part)
[*] Sorted BAM: 5205 records, merged and written in 0.287s after mapping
```

That is the whole alignment step: the reference is indexed on the fly, the output is a
coordinate-sorted BAM with its `.bai` index (no `samtools sort` needed), and `-y` keeps each read's
tags. A dorado unaligned BAM (`dorado basecaller ... > calls.bam`) or a FASTQ works the same way.

```bash
samtools flagstat demo.bonefish.bam | head -4
samtools view -H demo.bonefish.bam | grep '^@RG'                 # dorado's read group, with the model
samtools view demo.bonefish.bam | head -1 | cut -f12- | tr '\t' '\n' | grep -E '^(MM|ML|RG|qs)'
```

The `@RG` line (which records the basecalling model) and the per-read `MM`/`ML`, `RG`, `qs`, `ch`,
... tags are carried over; old alignment tags in the input are replaced by the new alignment's.

## 4. Methylation calls

```bash
modkit pileup demo.bonefish.bam demo.cpg.bed --modified-bases 5mC 5hmC --cpg --combine-strands \
    --ref demo.fasta -t 8
```

`demo.cpg.bed` lists 5mC and 5hmC per CpG. On a full chromosome, CpG methylation from bonefish's
BAM agrees with that from minimap2's at r = 0.999
([details](validation/2026-10-14-modkit-methylation/README.md)).

## 5. Larger runs

- **Whole human genome.** Build the index once (about 1 minute), in the raw format that loads in a
  fraction of a second:

  ```bash
  bonefish -t 16 -d GRCh38.raw --idx-raw GRCh38.fa
  bonefish -t 16 -y -o sample.bam GRCh38.raw sample.bam     # dorado BAM(s) or FASTQ(s)
  ```

  A 30x genome takes about 35-40 minutes on a 16-thread desktop and about 20 GB of memory
  (about 23 GB for a sorted BAM with methylation tags; `--sort-mem 1G` lowers it).
- **Several input files** are mapped in turn: `bonefish -o sample.bam ref.raw run1.bam run2.bam`.
- **Streaming**: `dorado basecaller ... | bonefish -y -o sample.bam ref.raw -`.
- **SAM or PAF**: `bonefish -a ref reads > out.sam`, `bonefish -c --cs ref reads > out.paf`.
- **Variant calling**: the BAM is a drop-in replacement for minimap2's in Clair3, Sniffles2, modkit
  and similar tools ([accuracy](benchmarks.md)).
