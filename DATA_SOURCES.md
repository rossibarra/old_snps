# Data sources

Sources for the B73 v5 marker positions in `vanHeerwaarden2011_markers_v5.final.tsv`
(983 markers). The `method` column of that table says which source each position came from.

## References and data used

**van Heerwaarden J, Doebley J, Briggs WH, Glaubitz JC, Goodman MM, Sánchez González JJ,
Ross-Ibarra J (2011).** Genetic signals of origin, spread, and introgression in a large sample of
maize landraces. *PNAS* 108:1088–1092. doi:10.1073/pnas.1013011108
- The list of 983 SNP markers.

**MaizeGDB (2026).** *Panzea legacy SNPs on B73 RefGen_v1, RefGen_v4 and
Zm-B73-REFERENCE-NAM-5.0 (v5)*, Excel workbook (`panzea_legacy_snps_v1_v4_v5.xlsx`), produced by
MaizeGDB, 2026-10-19. https://maizegdb.org
- Exact v5 positions for 507 markers, with RefGen_v1/v4 positions, B73 base, alleles and v5 gene model.
- v5 amplicon intervals for 312 markers.

**MaizeGDB (2026).** Marker records exported from MaizeGDB, 2026-10-19
(`maizegdb_marker_list_20261019_110740.tsv`).
- v5 gene-model coordinates for 17 markers with no SNP- or amplicon-level position.

**Yan J, Yang X, Shah T, Sánchez-Villeda H, Li J, Warburton M, Zhou Y, Crouch JH, Xu Y (2010).**
High-throughput SNP genotyping with the GoldenGate assay in maize. *Molecular Breeding*
25:441–451. doi:10.1007/s11032-009-9343-2
- SNP flanking sequences (Supplementary Table 2, file `11032_2009_9343_MOESM3_ESM.xlsx`).
- Aligned to v5 to give exact positions for 39 markers.

**Tian HL, Wang FG, Zhao JR, et al. (2015).** Development of maizeSNP3072, a high-throughput compatible SNP array, for
DNA fingerprinting identification of Chinese maize varieties. *Molecular Breeding* 35:136.
doi:10.1007/s11032-015-0335-0
- SNP flanking sequences (file `11032_2015_335_MOESM6_ESM.xlsx`).
- Aligned to v5 to give exact positions for 4 markers: PZA00005.8, PZA00300.16, PZA03081.10, PZB00067.5.

**Hufford MB, et al. (2021).** De novo assembly, annotation, and comparative analysis of 26 diverse
maize genomes. *Science* 373:655–662. doi:10.1126/science.abg5289
- The B73 v5 reference genome (Zm-B73-REFERENCE-NAM-5.0) used for the flank alignments.

**Panzea germplasm-survey SNP list (J. Glaubitz, 2007, unpublished;
`Germplasm_survey-all_useful_genotyped_SNPs.xls`).**
- Assigns each SNP to its Panzea assay amplicon.
- Used to infer intervals for 17 markers (`amplicon_inferred`) from other SNPs in the same amplicon.
- Used to accept 2 low-confidence flank alignments (PZA03305.7, PZB00221.3) that fall inside a
  same-amplicon SNP's interval. Their flanks are among the 39 from Yan et al. (2010) above.

## Methods

The full table is rebuilt by three scripts, run in order:

1. **`align_flanks.py`**: aligns the published flanks (`flanks_to_align.fa`, with the source paper in
   each header) to Zm-B73-REFERENCE-NAM-5.0. It uses minimap2's `sr` preset (via mappy) and finds the
   SNP base by walking the alignment CIGAR. Output: `flanks_v5.tsv`.
2. **`merge_flanks.py`**: writes the hits with status `ok` (mapq ≥ 20, unique) into the MaizeGDB-based
   table as `exact`. Output: `vanHeerwaarden2011_markers_v5.merged.tsv`, 41 markers replaced.
3. **`finalize.py`**: uses `amplicons.tsv` (marker → Panzea amplicon, from the germplasm-survey list)
   to do two things. Output: `vanHeerwaarden2011_markers_v5.final.tsv`.
   - It accepts a low-mapq hit only when the hit lies inside the interval of a same-amplicon SNP.
     That applies to 2 markers.
   - For unplaced PZ markers, it infers an interval: the span of their placed same-amplicon SNPs,
     ±750 bp. That applies to 17 markers. These are not measured positions for the SNP itself.

Rerunning steps 2 and 3 on the committed inputs reproduces the committed tables byte for byte.

## Used only to check positions

These sources didn't contribute any positions:

- Yan et al. (2010), Supplementary file `11032_2009_9343_MOESM5_ESM.xlsx`: physical positions.
- Tian et al. (2015), AGPv3 coordinates.
- `MR_0676.xlsx`: RefGen_v2 positions for 44 markers, from a teosinte GBS genotype table.
