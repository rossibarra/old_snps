# van Heerwaarden 2011 markers on B73 v5: finishing step

`vanHeerwaarden2011_markers_v5.tsv` gives a v5 position for each of the 983 markers, taken from
MaizeGDB's "Panzea legacy SNPs on v1/v4/v5" table (built 2026-10-19).

| precision | markers | what it means |
|---|---|---|
| exact | 508 | exact SNP base on v5, with RefGen_v1 and RefGen_v4 positions |
| amplicon | 319 | SNP lies somewhere in v5_start–v5_end (median about 500 bp) |
| gene | 17 | only MaizeGDB's v5 gene model for the marker |
| unplaced | 139 | no position; 85 of these are not in MaizeGDB at all |

How well MaizeGDB's exact positions agree with independent sources:

| check | result |
|---|---|
| Published flanks from Yan 2010 Mol Breed S3 | 38 of 39 match MaizeGDB's v5 flank |
| Chromosome from the MR_0676 table (RefGen_v2) | 26 of 26 match |
| Chromosome from Mol Breed 2015 S6 (AGPv3) | 7 of 7 match |
| MaizeGDB marker-to-gene links | 44 of 44 fall within 20 kb |
| Yan 2010 S5 physical positions | 12 of 12 match RefGen_v1, within 0–2 bp |

The one flank mismatch is PZA00355.2. MaizeGDB also notes that its B73 base is not one of the
alleles, so its v1-based position is probably wrong.

## Remaining step
`flanks_to_align.fa` holds the published flanks for 49 markers. 48 of them are non-exact in the
table, and the 49th is PZA00355.2. Align them to v5 to get exact positions:

    pip install mappy
    python3 align_flanks.py Zm-B73-REFERENCE-NAM-5.0.fa.gz flanks_to_align.fa > flanks_v5.tsv

Then, for rows with status `ok`, replace those markers' rows in the main table with
precision=exact and the new chr/pos. That brings the table to about 556 exact positions.

`align_flanks.py` was tested on a synthetic genome with SNPs on both strands. It has not been run
on real data yet.
