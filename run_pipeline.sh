#!/bin/bash
# Rebuild the final table. Step 1 needs a compute node (see README); steps 2-3 are light.
# Usage: ./run_pipeline.sh [--align]      (--align also redoes the flank alignment)
set -euo pipefail
GENOME=/quobyte/jrigrp/DATA/RESOURCES/NAM/fasta/Zm-B73-REFERENCE-NAM-5.0.fa.gz
if [[ "${1:-}" == "--align" ]]; then
  python -I align_flanks.py "$GENOME" flanks_to_align.fa > flanks_v5.tsv
fi
python -I merge_flanks.py .
python -I finalize.py . amplicons.tsv
