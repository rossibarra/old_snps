#!/usr/bin/env bash
# Self-test for place_markers.py on a random 2-chromosome genome (+ and - strand SNPs).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
T="$(mktemp -d)"; cd "$T"
python3 - << 'EOF'
import random
random.seed(1)
g = {f"chr{i}": "".join(random.choice("ACGT") for _ in range(200000)) for i in (1, 2)}
with open("Zm-B73-REFERENCE-NAM-5.0.fa", "w") as f:
    for k, v in g.items(): f.write(f">{k}\n{v}\n")
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
s, p = g["chr1"], 50000
alt = "A" if s[p] != "A" else "C"
a = f"{s[p-100:p]}[{s[p]}/{alt}]{s[p+1:p+101]}"
t, p2 = g["chr2"], 120000
r = rc(t[p2-100:p2+101]); b = r[:100] + "[" + r[100] + "/N]" + r[101:]
open("MaizeGDB_genetic_SNPs.fa", "w").write(f">PZA00003.11 x\n{a}\n>ae1.5\n{b}\n")
open("names.txt", "w").write("PZA00003.11\nae1.5\nmissing.1\n")
EOF
python3 "$HERE/place_markers.py" Zm-B73-REFERENCE-NAM-5.0.fa
cut -f1-4,12,13 markers_v5.tsv
echo "EXPECT: PZA00003.11 chr1 50001 + ; ae1.5 chr2 120001 - ; missing.1 no_sequence_found"
