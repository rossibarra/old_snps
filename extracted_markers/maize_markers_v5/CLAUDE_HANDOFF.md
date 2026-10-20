# Handoff: place maize SNP markers on B73 v5

## Goal
Give a bp position on the B73 v5 genome (Zm-B73-REFERENCE-NAM-5.0) for each of the
983 maize SNP markers in `names.txt`. The user is Jeff Ross-Ibarra (UC Davis, maize genetics),
so no background explanations are needed.

Marker types:
- Panzea SNPs: `PZA#####.#`, `PZB#####.#`, `PZD#####.#`
- Gene-named Panzea SNPs: `ae1.5`, `tb1.17`, `Ra2_promoter.1`, etc.

These are the legacy GoldenGate / NAM-era markers.

Deliverable: `markers_v5.tsv`, with one row per marker giving chr, 1-based bp, strand,
mapping QC, and the v4 position as a cross-check. Include a short summary of how many
markers placed cleanly and which ones need attention.

## Approach (agreed with the user)
1. Get each marker's context sequence (the SNP plus its flanks) from MaizeGDB.
2. Align the sequences to B73 v5 with minimap2, using the `mappy` Python package.
3. Convert each SNP's offset in the query into a reference coordinate by walking the CIGAR string.
4. Cross-check the result against the v4 positions in MaizeGDB's xref file.

Fallbacks for markers with no sequence:
- If a v4 position exists, lift it over to v5 with CrossMap and a v4→v5 chain file
  (from MaizeGDB or Ensembl Plants).
- Otherwise, report the v5 coordinates of the named gene (e.g. ae1 → its Zm00001eb gene
  model), flagged as gene-level only.

## Data to download (all from download.maizegdb.org)
- `SNPs_and_Maps/SNPs/MaizeGDB_genetic_SNPs.fa`: SNP context sequences
- `SNPs_and_Maps/SNPs/MaizeGDB_genetic_SNPs.xref`: SNP IDs, v3/v4 positions, gene models,
  and genetic map positions (including NAM)
- `SNPs_and_Maps/SNPs/MaizeGDB_genetic_SNPs.README`: read this first
- The B73 v5 genome FASTA. Look under `Zm-B73-REFERENCE-NAM-5.0/` or
  `Genomes/NAM_Founders/Zm-B73-REFERENCE-NAM-5.0/`; the file is
  `Zm-B73-REFERENCE-NAM-5.0.fa.gz`.

The previous session was blocked from reaching MaizeGDB, Ensembl and NCBI by network
policy, so these files were never inspected. Check the following before trusting the output:

- **FASTA headers:** confirm the marker names appear in them. The script matches any
  whitespace-, `|`-, `;`- or `,`-separated token in the header, case-insensitively.
  If the names live only in the xref file, join through the xref ID column.
- **SNP notation in the sequence:** the script handles `[A/G]` brackets and a single
  IUPAC ambiguity code. If the file uses something else, such as a separate offset
  column, adapt `parse_context()`.
- **xref column names:** the script stores the whole xref row as `v4_info`. You'll
  probably want to split it into chr and pos columns and compare v4 against v5:
  v4 and v5 positions are close on the same chromosome for most markers, but not identical.
- **Coverage:** count how many of the 983 names are in the files. Some gene-named or
  PZD markers may be missing.

## Code
- `place_markers.py` does the alignment. Usage: `pip install mappy`, then
  `python3 place_markers.py Zm-B73-REFERENCE-NAM-5.0.fa.gz`, run from the folder that
  holds `names.txt` and the two MaizeGDB files.
  - It builds and caches a minimap2 index, `B73v5.sr.mmi`. Building it takes several
    minutes and roughly 8–10 GB of RAM (an estimate).
  - Preset is `sr`. If flanks turn out to be long (>500 bp), consider switching to the
    default or `asm5`.
  - Output columns: marker, chr_v5, pos_v5 (1-based), strand, mapq, aln_identity, aln_len,
    query_len, n_hits, snp_in_seq, alleles, v5_ref_base, status, v4_info, fasta_header.
  - Values in `status`: `ok`, `no_sequence_found`, `unmapped`,
    `snp_offset_unknown(aln_start_given)`, `snp_outside_alignment`, `snp_in_insertion`.
    `;multi_or_lowmapq` and `;low_identity` can be appended to any of them.
- `test_synthetic.sh` runs a self-test on a random two-chromosome genome. It passed in the
  previous session: both a plus-strand and a minus-strand SNP were placed at the exact base.

## Sanity checks
- `v5_ref_base` should match one of the two SNP alleles; if it doesn't, check the strand.
  For minus-strand hits the alleles are given in the marker's orientation, so compare
  against the reverse complement.
- Gene-named markers should fall inside, or very near, the v5 gene model of the named gene.
- Markers with several SNPs in one amplicon (e.g. PZA03775.1–.11) should cluster within
  about 1 kb of each other.
- Report multi-mapping markers rather than silently picking one hit.

## Where to run
Run it on any Linux or Mac machine that can reach download.maizegdb.org, such as UC Davis
Farm or the user's Mac.
