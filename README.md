# old_snps: van Heerwaarden 2011 markers on B73 v5

B73 v5 (Zm-B73-REFERENCE-NAM-5.0) positions for the 983 legacy Panzea SNP markers of
van Heerwaarden et al. (2011). Sources and references: [DATA_SOURCES.md](DATA_SOURCES.md).

**Result:** `vanHeerwaarden2011_markers_v5.final.tsv`. The `precision` column says how each
position was obtained:

| precision | markers | meaning |
|---|---|---|
| exact | 550 | exact SNP base on v5 |
| amplicon | 312 | SNP lies somewhere in v5_start–v5_end (MaizeGDB amplicon) |
| amplicon_inferred | 17 | interval inferred from placed SNPs in the same Panzea amplicon (±750 bp); not measured for the SNP itself |
| gene | 17 | only MaizeGDB's v5 gene model |
| unplaced | 87 | no position |

## Files
- `vanHeerwaarden2011_markers_v5.tsv`: input table, positions from MaizeGDB (unmodified)
- `flanks_to_align.fa`: published flanks for 49 markers
- `amplicons.tsv`: marker to Panzea assay amplicon
- `align_flanks.py`: aligns flanks to v5 with minimap2 `sr` (via mappy) and finds the SNP base from the CIGAR
- `merge_flanks.py`: writes flank hits with status `ok` into the table as `exact` (41 markers; an intermediate `*.merged.tsv`, gitignored)
- `finalize.py`: accepts 2 low-mapq hits that fall inside a same-amplicon SNP's interval, and infers intervals for 17 unplaced markers
- `flanks_v5.tsv`: alignment output
- `run_pipeline.sh`, `environment.yml`: rebuild the final table; conda env `old-snps` (python 3.11, bioconda `mappy`; `pip install mappy` fails here for lack of `zlib.h`)

## Rebuild
`./run_pipeline.sh` runs merge and finalize, which are light. Add `--align` to redo the alignment, which
indexes the genome and needs a compute node (it ran with 4 CPUs, 24 GB):

    HPC_CPUS=4 HPC_MEM=24G HPC_TIME=01:00:00 ~/.claude/bin/hpc_run './run_pipeline.sh --align'

Merge and finalize were rerun from the committed inputs and reproduce the final table byte for byte.
The alignment step was not rerun after the files were moved.

## Flank alignment results
41 of 49 flanks aligned cleanly (`ok`) and 2 low-mapq hits (PZA03305.7, PZB00221.3) were accepted by the
amplicon rule. The other 4 low-mapq hits (PZA00060.2, PZA00088.3, PZA00213.19, PZA00395.2) and the 2
unmapped (PZA00225.8, PZA03001.15) get no exact position; their precision in the final table is
unplaced (PZA00060.2, PZA00088.3, PZA00213.19, PZA00225.8), amplicon (PZA00395.2) or amplicon_inferred (PZA03001.15). PZA00355.2, which the MaizeGDB notes flagged as probably wrong, moved 112 bp
(chr6:89415515 to chr6:89415403).

Not checked: the `alleles_v5_plus` strand flip for minus-strand hits against the other exact rows, and the
low-mapq alignments beyond the amplicon rule.
