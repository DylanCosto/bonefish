#!/usr/bin/env python3
"""Stream a public S3 object to stdout with parallel ranged GETs (a single stream tops out at
~10-15 MB/s; ranged parallel reads reach the multipart-download rate). Chunks are written in order;
at most WINDOW chunks are in memory. Retries each range on errors.

Usage: s3_parallel_cat.py s3://bucket/key [workers=16] [chunk_mb=32] [max_bytes] > out
(max_bytes: stream only the first max_bytes of the object, e.g. a prefix of a large unaligned BAM)
"""
import sys, time, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

url = sys.argv[1]
workers = int(sys.argv[2]) if len(sys.argv) > 2 else 16
chunk = (int(sys.argv[3]) if len(sys.argv) > 3 else 32) << 20
bucket, key = url[5:].split("/", 1)
http = f"https://{bucket}.s3.amazonaws.com/{key}"

def size():
    req = urllib.request.Request(http, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])

def fetch(start, end):
    for attempt in range(8):
        try:
            req = urllib.request.Request(http, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) == end - start + 1:
                return data
        except (urllib.error.URLError, OSError, TimeoutError):
            pass
        time.sleep(2 ** attempt)
    raise RuntimeError(f"range {start}-{end} failed")

n = size()
if len(sys.argv) > 4:
    n = min(n, int(float(sys.argv[4])))
ranges = [(s, min(s + chunk, n) - 1) for s in range(0, n, chunk)]
out = sys.stdout.buffer
window = 2 * workers
with ThreadPoolExecutor(workers) as ex:
    pending = {}
    nxt = 0
    for i in range(min(window, len(ranges))):
        pending[i] = ex.submit(fetch, *ranges[i])
    while nxt < len(ranges):
        out.write(pending.pop(nxt).result())
        j = nxt + window
        if j < len(ranges):
            pending[j] = ex.submit(fetch, *ranges[j])
        nxt += 1
out.flush()
