#!/usr/bin/env python3
"""
make_exons.py - turn a UCSC RefSeq BED12 (one record per transcript) into a
BED4 list of unique protein-coding exons (CDS part only).

Usage: python make_exons.py exon.bed exons.bed

Filters (state these in the paper):
  - curated protein-coding transcripts only (accession starts with NM_)
  - standard chromosomes only (chr1-22, X, Y); no alt/random/fix contigs
  - only the coding part of each exon (clipped to thickStart..thickEnd)
  - identical exon coordinates from different transcripts counted once
"""
import sys, re

STD = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

def main(inp, out):
    seen, n_tx, n_nm = set(), 0, 0
    rows = []
    for line in open(inp):
        if line.startswith(("track", "browser", "#")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            sys.exit("This is not a BED12 file (need 12 columns).")
        n_tx += 1
        chrom, name = f[0], f[3]
        if chrom not in STD or not name.startswith("NM_"):
            continue
        n_nm += 1
        start = int(f[1]); tS = int(f[6]); tE = int(f[7])
        sizes = [int(x) for x in f[10].strip(",").split(",")]
        offs = [int(x) for x in f[11].strip(",").split(",")]
        for s, o in zip(sizes, offs):
            a, b = start + o, start + o + s
            a, b = max(a, tS), min(b, tE)      # coding part only
            if b <= a:
                continue
            key = (chrom, a, b)
            if key in seen:
                continue
            seen.add(key)
            rows.append((chrom, a, b, name))
    with open(out, "w") as w:
        for c, a, b, n in rows:
            w.write(f"{c}\t{a}\t{b}\t{n}\n")
    ge = sum(1 for r in rows if r[2] - r[1] >= 237)
    print(f"transcripts read: {n_tx}; NM_ on standard chromosomes: {n_nm}")
    print(f"unique coding exons written: {len(rows)}; of these >=237 bp: {ge}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
