#!/usr/bin/env python3
"""Split one 10x 3' library (one lane) into N pseudo-samples by cell barcode.

Every read pair goes to the pseudo-sample chosen by its cell barcode (first 16 bases of R1),
so each cell ends up in exactly one pseudo-sample, as if they were separate libraries.

Usage:
  split_by_barcode.py <R1.fastq.gz> <R2.fastq.gz> <outdir> <lane_number> [n_samples] [keep_fraction]

keep_fraction (0-1, default 1) downsamples: a read pair is kept when the hash of its read name
falls below the fraction, so R1 and R2 are always kept or dropped together.

Output (facility-style names, so the pipeline's run/lane parsing works):
  <outdir>/SLX-00000.TEST<k>.NEURON1K.s_<lane>.r_{1,2}.fq.gz
"""
import subprocess
import sys
import zlib

r1_path, r2_path, outdir, lane = sys.argv[1:5]
n_samples = int(sys.argv[5]) if len(sys.argv) > 5 else 4
keep_fraction = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
keep_below = int(keep_fraction * 2**32)


def reader(path):
    return subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 20).stdout


def writer(path):
    out = open(path, "wb")
    proc = subprocess.Popen(["pigz", "-p", "2", "-c"], stdin=subprocess.PIPE, stdout=out, bufsize=1 << 20)
    return proc, out


r1_in, r2_in = reader(r1_path), reader(r2_path)
writers = {}
for k in range(1, n_samples + 1):
    base = f"{outdir}/SLX-00000.TEST{k}.NEURON1K.s_{lane}.r_"
    writers[k] = (writer(base + "1.fq.gz.part"), writer(base + "2.fq.gz.part"))

counts = {k: 0 for k in writers}
n_in = 0
while True:
    rec1 = [r1_in.readline() for _ in range(4)]
    rec2 = [r2_in.readline() for _ in range(4)]
    if not rec1[0]:
        break
    n_in += 1
    name = rec1[0].split()[0]
    if not rec2[0] or name != rec2[0].split()[0]:
        sys.exit(f"R1/R2 out of step at read {n_in}")
    # Downsampling: decided per read pair (adler32 of the name, independent of the barcode hash)
    if keep_fraction < 1 and zlib.adler32(name) * 2654435761 % 2**32 >= keep_below:
        continue
    k = zlib.crc32(rec1[1][:16]) % n_samples + 1
    (p1, _), (p2, _) = writers[k]
    p1.stdin.write(b"".join(rec1))
    p2.stdin.write(b"".join(rec2))
    counts[k] += 1

for k, ((p1, o1), (p2, o2)) in writers.items():
    for proc, out in ((p1, o1), (p2, o2)):
        proc.stdin.close()
        proc.wait()
        out.close()
        if proc.returncode != 0:
            sys.exit(f"pigz failed for TEST{k}")

import os
for k in writers:
    base = f"{outdir}/SLX-00000.TEST{k}.NEURON1K.s_{lane}.r_"
    for r in (1, 2):
        os.replace(base + f"{r}.fq.gz.part", base + f"{r}.fq.gz")

print(f"lane {lane}: {n_in:,} read pairs in, kept {sum(counts.values()):,} (fraction {keep_fraction}): "
      + ", ".join(f"TEST{k} {n:,}" for k, n in counts.items()))
