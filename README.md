# ZFN vs TALEN vs CRISPR-Cas9: in-silico comparison at BCL11A

Consolidated comparison, shared code and the manuscript: *Design Architecture of ZFNs, TALENs and CRISPR-Cas9: An In-Silico Case Study at One BCL11A Coding Exon, with a 45-Exon Site-Count Extension.*

Related repositories: `crispr-cas9-bcl11a-design`, `talen-bcl11a-design`, `zfn-bcl11a-design`.

## Contents
| Path | What it is |
|---|---|
| `docs/` | Final manuscript PDF |
| `site_count.py` | Rule-based candidate-site counts for SpCas9, TALEN and ZFN (Supplement S1) |
| `fetch_windows.py` | Downloads the central 237 bp of listed regions from the UCSC sequence service (S2) |
| `make_exons.py` | Converts a UCSC RefSeq BED12 into unique coding exons, NM_ transcripts only (S3) |
| `data/` | Amplicon FASTA, test region, density table (Table 12) |
| `figures/` | Schematics and the script screenshots (Figures 1, 3, 20-22) |

## Reproduce the single-locus counts
```
python site_count.py
```
Expected: SpCas9 32, TALEN 47 (tool-like) and 58 (broad), ZFN 1 (strict), 19 (relaxed), 0 (4 fingers).

## Reproduce the 45-window run
1. UCSC Table Browser: hg38, NCBI RefSeq, BED output, save as `exon.bed` (one record per transcript).
2. `python make_exons.py exon.bed exons.bed`  (expect 216,361 unique coding exons; 23,754 of at least 237 bp)
3. `python -c "import random;L=[l for l in open('exons.bed')];random.seed(1);open('regions.txt','w').writelines(random.sample(L,500))"`
4. `python fetch_windows.py regions.txt exons.fasta`  (45 windows were kept in the reported run)
5. `python site_count.py exons.fasta -o results.csv`

Requires Python 3. If the download fails with a certificate error, give Python a trusted certificate store (for example `pip install certifi` and set `SSL_CERT_FILE`).

## Median candidate density (per kb, rule-based)
| Rule | BCL11A window | Median of 45 windows |
|---|---|---|
| SpCas9 NGG | 135.0 | 143.5 |
| TALEN tool-like | 198.3 | 126.6 |
| TALEN broad | 244.7 | 189.9 |
| ZFN strict GNNx3 | 4.2 | 0.0 |
| ZFN relaxed | 80.2 | 97.0 |

## Important limits
- Everything here is a computational prediction or rule-based count. Nothing was built or tested in the laboratory.
- Tools were accessed on 24 October 2026; the genome build is GRCh38/hg38.
- Tool screenshots show a faint "Activate Windows" watermark; recapture them before journal submission.

## Cite
See `CITATION.cff` (Zenodo DOI: add after the first release).
