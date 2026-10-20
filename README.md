# old_snps: van Heerwaarden 2011 markers on B73 v5

Positions on B73 v5 (Zm-B73-REFERENCE-NAM-5.0) for the 983 legacy Panzea SNP markers.

## Layout
- `extracted_vanH/vanHeerwaarden2011_markers_v5/`: unpacked `vanHeerwaarden2011_markers_v5_final.zip`
  - `vanHeerwaarden2011_markers_v5.tsv`: original table (MaizeGDB positions), unmodified
  - `flanks_to_align.fa`, `align_flanks.py`: published flanks for 49 markers and the aligner
  - `flanks_v5.tsv`: alignment output (this work)
  - `vanHeerwaarden2011_markers_v5.merged.tsv`: updated table (this work)
  - `README.md`: the original instructions
- `extracted_markers/maize_markers_v5/`: unpacked `maize_markers_v5.zip`, an earlier handoff for the same task. Not run.
- `environment.yml`: conda env `old-snps` (python 3.11, bioconda `mappy`)

## What was done
1. Created the `old-snps` conda env. `pip install mappy` fails here (no `zlib.h`), so the env uses bioconda.
2. Aligned the 49 flanks with `align_flanks.py` (minimap2 `sr` preset via mappy) against
   `/quobyte/jrigrp/DATA/RESOURCES/NAM/fasta/Zm-B73-REFERENCE-NAM-5.0.fa.gz`, on a compute node (4 CPUs, 24 GB):

       python -I align_flanks.py Zm-B73-REFERENCE-NAM-5.0.fa.gz flanks_to_align.fa > flanks_v5.tsv

3. Merged rows with status `ok` into the table as `precision=exact` (new chr/pos/base, plus-strand
   alleles, old precision and position recorded in `note`). Merge script: not kept in the repo.

## Results
| status | markers |
|---|---|
| ok, merged as exact | 41 (7 were amplicon, 33 unplaced, 1 was PZA00355.2) |
| low_mapq, not merged | 6: PZA00060.2, PZA00088.3, PZA00213.19, PZA00395.2, PZA03305.7, PZB00221.3 |
| unmapped, not merged | 2: PZA00225.8, PZA03001.15 |

Merged table: 548 exact, 312 amplicon, 17 gene, 106 unplaced (was 508 / 319 / 17 / 139).

Checks: all 7 amplicon-level markers landed inside their MaizeGDB amplicon. PZA00355.2 moved
112 bp (chr6:89415515 to chr6:89415403), as the original README suspected.

Not checked: the `alleles_v5_plus` strand flip for minus-strand hits against the other 508 exact rows,
and the low_mapq alignments.
