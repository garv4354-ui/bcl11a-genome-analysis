#!/usr/bin/env python3
"""
fetch_windows.py - download DNA windows from the UCSC sequence service and write a FASTA file.

Input: a text file with one region per line (tab or space separated):
        chrom  start  end  name        e.g.  chr2  60462160  60462397  BCL11A_amp
       start/end are 0-based half-open, exactly as in a UCSC BED file.
Each region longer than WINDOW bp is cut to its central WINDOW bp; shorter regions are skipped,
so every sequence has the same length (default 237 bp, the length of the BCL11A amplicon).

Usage:  python fetch_windows.py regions.txt exons.fasta [window_bp]
Needs an internet connection. Not run in the sandbox where this was written: test it on one line first.
"""
import sys, json, time, urllib.request

API = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={c};start={s};end={e}"

def fetch(chrom, start, end):
    with urllib.request.urlopen(API.format(c=chrom, s=start, e=end), timeout=30) as r:
        return json.load(r)["dna"].upper()

def main():
    regions, out = sys.argv[1], sys.argv[2]
    window = int(sys.argv[3]) if len(sys.argv) > 3 else 237
    kept = skipped = 0
    with open(regions) as fh, open(out, "w") as fo:
        for line in fh:
            parts = line.split()
            if len(parts) < 3 or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) > 3 else f"{chrom}_{start}_{end}"
            if end - start < window:
                skipped += 1; continue
            mid = (start + end) // 2
            s = mid - window // 2
            e = s + window
            seq = fetch(chrom, s, e)
            fo.write(f">{name}_{chrom}:{s}-{e}\n{seq}\n")
            kept += 1
            time.sleep(0.3)          # be polite to the server
    print(f"wrote {kept} windows to {out}; skipped {skipped} regions shorter than {window} bp")

if __name__ == "__main__":
    main()
