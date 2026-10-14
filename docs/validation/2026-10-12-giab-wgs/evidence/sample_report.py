#!/usr/bin/env python3
"""Markdown report for one bench/sample_eval.sh run: speed and resources, BAM statistics, Clair3
small-variant accuracy (overall and stratified), SVs, and the primary-alignment comparison,
minimap2 vs bonefish side by side.

Usage: sample_report.py OUT_DIR > report.md
"""
import csv
import os
import sys
from collections import defaultdict

O = sys.argv[1]
AL = ["minimap2", "bonefish"]
GENOME = 3.1e9


def kv(path):
    d = {}
    if os.path.exists(path):
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                d[p[0]] = p[1]
    return d


def timefile(path):
    """GNU time "%e %U %S %M" -> (wall s, cpu s, peak RSS GB)."""
    if not os.path.exists(path):
        return None
    w, u, s, m = open(path).read().split()[-4:]
    return float(w), float(u) + float(s), int(m) / 1048576


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" if i == 0 else "---:" for i in range(len(header))) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def fmt_t(t):
    return f"{t / 60:.1f} min" if t is not None else "-"


def rows_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


sample = os.path.basename(O.rstrip("/"))
print(f"# {sample}: bonefish vs minimap2, whole genome\n")

# Input
s_mm = kv(f"{O}/minimap2/bam_summary.tsv")
bases = int(s_mm.get("primary_query_bases", 0))
dl = timefile(f"{O}/download.time")
print("## Input\n")
print(f"- Reads: {int(s_mm.get('reads', 0)):,}; primary query bases {bases / 1e9:.1f} Gbp (~{bases / GENOME:.1f}x, "
      "unmapped reads' bases not counted)")
cv = timefile(f"{O}/convert.time")
if dl:
    print(f"- Download of ONT's BAM from S3: {fmt_t(dl[0])}" + (f"; BAM -> FASTQ: {fmt_t(cv[0])}" if cv else ""))
print()

# Stage resources from the monitor
stages = defaultdict(lambda: dict(n=0, cpu=0.0, mem=0.0, swap=0.0, iow=0.0, rss=0.0))
swap0 = None
if os.path.exists(f"{O}/resources.csv"):
    for r in csv.DictReader(open(f"{O}/resources.csv")):
        if swap0 is None:
            swap0 = float(r["swap_used_gb"])
        st = stages[r.get("stage", "-")]
        st["n"] += 1
        st["cpu"] += float(r["cpu_pct"])
        st["iow"] += float(r.get("iowait_pct") or 0)
        st["mem"] = max(st["mem"], float(r["mem_used_gb"]))
        st["swap"] = max(st["swap"], float(r["swap_used_gb"]) - swap0)
        st["rss"] = max(st["rss"], float(r["proc_rss_gb"]))

print("## Speed and resources\n")
rows = []
for a in AL:
    al = timefile(f"{O}/{a}/align.time")
    srt = timefile(f"{O}/{a}/sort.time")
    e2e = float(open(f"{O}/{a}/to_sorted_bam.seconds").read()) if os.path.exists(f"{O}/{a}/to_sorted_bam.seconds") else None
    bam_gb = int(open(f"{O}/{a}/bam_bytes").read()) / 1e9 if os.path.exists(f"{O}/{a}/bam_bytes") else 0
    st = stages.get(f"align_{a}", {})
    cpu = al[1] + (srt[1] if srt else 0) if al else None
    rows.append([a, fmt_t(e2e), fmt_t(al[0]) if al else "-", f"{cpu / 3600:.2f} h" if cpu else "-",
                 f"{al[2]:.1f} GB" if al else "-", f"{srt[2]:.1f} GB" if srt else "(built in)",
                 f"{st.get('mem', 0):.1f} GB" if st else "-", f"{st['cpu'] / st['n']:.0f}%" if st and st.get("n") else "-",
                 f"{st['iow'] / st['n']:.1f}%" if st and st.get("n") else "-", f"{bam_gb:.1f} GB"])
print(table(["aligner", "reads -> sorted BAM", "aligner wall", "CPU (aligner + sort)", "aligner peak RSS",
             "sort peak RSS", "system mem peak", "avg CPU", "avg iowait", "BAM size"], rows))
t = [float(open(f"{O}/{a}/to_sorted_bam.seconds").read()) for a in AL if os.path.exists(f"{O}/{a}/to_sorted_bam.seconds")]
if len(t) == 2:
    print(f"\nbonefish is **{t[0] / t[1]:.2f}x** faster than minimap2 + samtools sort from reads to sorted, indexed BAM.\n")

print("### Every stage (from the 10 s resource samples)\n")
print("System memory includes other programs running on the machine; swap growth is relative to the "
      f"swap in use when the run started ({swap0 or 0:.1f} GB).\n")
rows = []
for name, st in stages.items():
    if st["n"] == 0 or name in ("-",):
        continue
    rows.append([name, fmt_t(st["n"] * 10), f"{st['cpu'] / st['n']:.0f}%", f"{st['iow'] / st['n']:.1f}%",
                 f"{st['mem']:.1f} GB", f"{st['rss']:.1f} GB", f"{st['swap']:.1f} GB"])
print(table(["stage", "duration", "avg CPU", "avg iowait", "system mem peak", "largest process RSS", "swap growth"], rows))
print()

# BAM statistics
print("## BAM statistics\n")
keys = ["reads", "primary_mapped", "unmapped", "pct_reads_mapped", "supplementary", "secondary",
        "pct_primary_bases_aligned", "primary_identity", "pct_primary_mapq_ge_20",
        "mapq_0-0", "mapq_1-9", "mapq_10-19", "mapq_20-29", "mapq_30-59", "mapq_60-60"]
summ = {a: kv(f"{O}/{a}/bam_summary.tsv") for a in AL}
print(table(["statistic"] + AL, [[k] + [summ[a].get(k, "-") for a in AL] for k in keys]))
print()
cmp_ = kv(f"{O}/compare_primary.tsv")
if cmp_:
    print("### Primary alignments, bonefish vs minimap2\n")
    print(table(["category", "reads"], [[k, v] for k, v in cmp_.items()]))
    print()

# Clair3
print("## Clair3 (sup v5.0.0 model, longphase intermediate phasing)\n")
rows = []
for a in AL:
    c3 = timefile(f"{O}/{a}/clair3.time")
    cnt = kv(f"{O}/{a}/clair3_counts.tsv")
    st = stages.get(f"clair3_{a}", {})
    rows.append([a, fmt_t(c3[0]) if c3 else "-", f"{c3[1] / 3600:.2f} h" if c3 else "-",
                 f"{st.get('mem', 0):.1f} GB" if st else "-", cnt.get("calls_total", "-"),
                 cnt.get("calls_pass_snps", "-"), cnt.get("calls_pass_indels", "-")])
print(table(["BAM from", "wall", "CPU", "system mem peak", "calls", "PASS SNPs", "PASS indels"], rows))
print()
print("### Small-variant accuracy vs GIAB v4.2.1 (rtg vcfeval, PASS calls)\n")
acc = {a: {(r["stratum"], r["type"]): r for r in rows_tsv(f"{O}/{a}/accuracy.tsv") + rows_tsv(f"{O}/{a}/strat.tsv")} for a in AL}
rows = []
for key in sorted(set(acc["minimap2"]) | set(acc["bonefish"]), key=lambda k: (k[0] != "all", k[0], k[1])):
    r = [f"{key[0]} {key[1]}"]
    for a in AL:
        x = acc[a].get(key)
        r += [x["f1"], f"{x['tp']}/{x['fp']}/{x['fn']}"] if x else ["-", "-"]
    m, b = acc["minimap2"].get(key), acc["bonefish"].get(key)
    r.append(f"{float(b['f1']) - float(m['f1']):+.4f}" if m and b else "-")
    rows.append(r)
print(table(["region / type", "minimap2 F1", "TP/FP/FN", "bonefish F1", "TP/FP/FN", "ΔF1"], rows))
print()

# SVs
print("## Structural variants (Sniffles2)\n")
rows = []
for a in AL:
    sn = timefile(f"{O}/{a}/sniffles.time")
    sv = {r["matching"]: r for r in rows_tsv(f"{O}/{a}/sv_accuracy.tsv")}
    r = [a, fmt_t(sn[0]) if sn else "-"]
    for mode in ("strict", "multi"):
        x = sv.get(mode)
        r += [f"{x['f1']} ({x['tp']}/{x['fp']}/{x['fn']})" if x else "-"]
    rows.append(r)
print(table(["BAM from", "Sniffles wall", "F1 vs GIAB v5.0q, strict (TP/FP/FN)", "F1, one call may match several truth SVs"], rows))
ag = kv(f"{O}/sv_agreement.tsv")
if ag:
    print("\nSV call agreement (Truvari, minimap2 calls as base): " + "; ".join(f"{k} {v}" for k, v in ag.items()))
print()
