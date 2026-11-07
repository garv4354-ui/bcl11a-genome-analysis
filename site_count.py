#!/usr/bin/env python3
"""
site_count.py  -  Supplementary file S1
Rule-based counts of candidate nuclease sites in one or more DNA sequences.
These are design-rule counts, NOT tool outputs and NOT efficiency/specificity predictions.

Usage
  python site_count.py                  # runs on the BCL11A 237 bp amplicon (chr2:60462161-60462397, hg38)
  python site_count.py exons.fasta      # runs on every sequence in a (multi-)FASTA file
  python site_count.py exons.fasta -o results.csv

Rules (same as Section 2.6 of the paper)
  SpCas9 : 20-nt protospacer + NGG, either strand, whole 23 nt inside the window.
  TALEN  : top-strand string = T0 + left arm + spacer + right arm + terminal A.
           tool-like : arms 17+17, spacer 14-16 bp, each (T0, spacer) pair counted.
           broad     : arms 15-20 each, spacer 12-30 bp, counted as distinct T0 positions
                       that have at least one valid partner.
  ZFN    : left half-site on the top strand, right half-site on the bottom strand,
           spacer 5-7 bp. A half-site is a series of 3-bp fingers; a finger is "GNN" if its
           first base is G.
           strict  : 3 fingers, all GNN        (GNNx3)
           relaxed : 3 fingers, at least 2 of 3 GNN
           four    : 4 fingers, all GNN        (GNNx4)
"""
import sys, csv, statistics

# BCL11A amplicon, chr2:60462161-60462397 (hg38), same orientation as the CHOPCHOP output
BCL11A = ("CCTCTCTCGATACTGATCCTGGTATTCTTAGCAGGTTAAAGGGGTTATTGTCTGCAATATGAATCCCATGGAGAGGTGGCTGGGAAGGACATTCTG"
          "CACCTAGTCCTGAAGGGATACCAACCCGCGGGGTCAGGGGACTTCCGTGTTCGCTTTCTAAGTAGATTCTTAATCCATGAGTGTTCTGTGCGTGTTG"
          "CAAGAGAAACCATGCACTGGTGAATGGCTGTTTGCAAGTTGTAC")

def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

def count_cas9(seq):
    """Protospacer (20 nt) + NGG on either strand; the whole 23 nt must lie inside the window."""
    n, total = len(seq), 0
    for i in range(n - 22):
        w = seq[i:i + 23]
        if w[21:23] == "GG":            # top strand: protospacer = w[0:20], PAM = w[20:23]
            total += 1
        if revcomp(w)[21:23] == "GG":   # bottom strand
            total += 1
    return total

def count_talen_tool(seq, arm=17, spacers=(14, 15, 16)):
    """T0 + arm + spacer + arm + terminal A. Spacer variants are counted separately."""
    n, total = len(seq), 0
    for i in range(n):
        if seq[i] != "T":
            continue
        for s in spacers:
            j = i + 1 + arm + s + arm          # index of the terminal A
            if j < n and seq[j] == "A":
                total += 1
    return total

def count_talen_broad(seq, arms=range(15, 21), spacers=range(12, 31)):
    """Distinct T0 positions that have at least one valid (arm, spacer, arm, terminal A) partner."""
    n, total = len(seq), 0
    for i in range(n):
        if seq[i] != "T":
            continue
        ok = False
        for left in arms:
            for right in arms:
                for s in spacers:
                    j = i + 1 + left + s + right
                    if j < n and seq[j] == "A":
                        ok = True
                        break
                if ok: break
            if ok: break
        total += ok
    return total

def count_zfn(seq, fingers=3, min_gnn=3, spacers=(5, 6, 7)):
    """Pairs of half-sites (top strand + bottom strand) with a 5-7 bp spacer.
    A half-site passes if at least `min_gnn` of its fingers start with G."""
    n, L, total = len(seq), 3 * fingers, 0
    def passes(site):
        return sum(site[k] == "G" for k in range(0, L, 3)) >= min_gnn
    for i in range(n - L + 1):
        left = seq[i:i + L]
        if not passes(left):
            continue
        for sp in spacers:
            j = i + L + sp
            if j + L > n:
                continue
            right = revcomp(seq[j:j + L])       # read 5'->3' on the bottom strand
            if passes(right):
                total += 1
    return total

def count_all(seq):
    return {
        "length_bp": len(seq),
        "cas9_NGG": count_cas9(seq),
        "talen_tool_like": count_talen_tool(seq),
        "talen_broad": count_talen_broad(seq),
        "zfn_strict_GNNx3": count_zfn(seq, 3, 3),
        "zfn_relaxed_2of3": count_zfn(seq, 3, 2),
        "zfn_4_fingers_GNNx4": count_zfn(seq, 4, 4),
    }

def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()

def main(argv):
    out_csv = None
    if "-o" in argv:
        k = argv.index("-o"); out_csv = argv[k + 1]; del argv[k:k + 2]
    records = list(read_fasta(argv[1])) if len(argv) > 1 else [("BCL11A_amplicon", BCL11A)]
    rows = []
    for name, seq in records:
        if set(seq) - set("ACGT"):
            print(f"note: {name} contains characters other than A/C/G/T (e.g. N); windows containing them may be miscounted")
        r = count_all(seq); r["name"] = name
        rows.append(r)
    keys = [k for k in rows[0] if k not in ("name", "length_bp")]
    print(f"{'sequence':28s} {'bp':>5s} " + " ".join(f"{k:>20s}" for k in keys))
    for r in rows:
        print(f"{r['name'][:28]:28s} {r['length_bp']:5d} " + " ".join(f"{r[k]:20d}" for k in keys))
    if len(rows) > 1:
        print("\nMedian candidates per kb across", len(rows), "sequences")
        for k in keys:
            vals = [r[k] / r["length_bp"] * 1000 for r in rows]
            print(f"  {k:22s} {statistics.median(vals):8.1f}")
    else:
        r = rows[0]
        print("\nPer kb:", ", ".join(f"{k}={r[k] / r['length_bp'] * 1000:.1f}" for k in keys))
    if out_csv:
        with open(out_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["name", "length_bp"] + keys)
            w.writeheader(); w.writerows(rows)
        print("wrote", out_csv)

if __name__ == "__main__":
    main(sys.argv)
