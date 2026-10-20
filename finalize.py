"""Final step: rescue low-mapq flank hits and infer intervals from same-amplicon SNPs.

Usage: python -I finalize.py . amplicons.tsv
Reads  <dir>/vanHeerwaarden2011_markers_v5.merged.tsv and <dir>/flanks_v5.tsv
Writes <dir>/vanHeerwaarden2011_markers_v5.final.tsv

amplicons.tsv maps each marker to its Panzea assay amplicon (from the Panzea germplasm-survey
SNP list); markers with no entry use the name minus its ".N" suffix.
1. A flank hit with status low_mapq is accepted as exact only if it lies inside the v5 interval
   of another SNP from the same amplicon.
2. An unplaced PZ* marker gets precision=amplicon_inferred: the span of its same-amplicon SNPs
   that have a v5 position, +/-750 bp (skipped if those SNPs disagree on chromosome).
"""
import csv, sys, collections
d, ampfile = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(f"{d}/vanHeerwaarden2011_markers_v5.merged.tsv"), delimiter="\t"))
cols = list(rows[0].keys())
M = {r["marker"]: r for r in rows}
F = {r["marker"]: r for r in csv.DictReader(open(f"{d}/flanks_v5.tsv"), delimiter="\t")}
A = {r["marker"]: r["amplicon"] for r in csv.DictReader(open(ampfile), delimiter="\t")}
amp = lambda m: A.get(m) or m.rsplit(".", 1)[0]
G = collections.defaultdict(list)
for m in M: G[amp(m)].append(m)
comp = str.maketrans("ACGT", "TGCA")
snap = {m: dict(r) for m, r in M.items()}
placed = lambda x: snap[x]["precision"] in ("exact", "amplicon")
n = collections.Counter()

for m, f in F.items():
    if f["status"] != "low_mapq": continue
    p = int(f["v5_pos"])
    sib = [x for x in G[amp(m)] if x != m and placed(x) and snap[x]["v5_chr"] == f["v5_chr"]
           and int(snap[x]["v5_start"]) <= p <= int(snap[x]["v5_end"])]
    if not sib: continue
    r = M[m]
    al = f["alleles_marker_strand"]
    if f["strand"] == "-": al = al.translate(comp)
    old = f"{r['precision']} {r['v5_chr']}:{r['v5_start']}-{r['v5_end']}".strip()
    r.update(precision="exact", v5_chr=f["v5_chr"], v5_pos=f["v5_pos"], v5_start=f["v5_pos"], v5_end=f["v5_pos"],
             v5_B73_base=f["v5_B73_base"], alleles_v5_plus=al,
             method=f"published flank aligned to v5 (mapq {f['mapq']}; accepted because hit lies inside amplicon of sibling {sib[0]})")
    r["note"] = "; ".join(x for x in [r["note"], f"was {old}; LOW MAPQ rescued by sibling amplicon"] if x)
    n["rescued"] += 1

for m, r in M.items():
    if r["precision"] != "unplaced" or not m.startswith("PZ"): continue
    sib = [x for x in G[amp(m)] if x != m and placed(x)]
    if not sib or len({snap[x]["v5_chr"] for x in sib}) != 1: continue
    lo = min(int(snap[x]["v5_start"]) for x in sib) - 750
    hi = max(int(snap[x]["v5_end"]) for x in sib) + 750
    r.update(precision="amplicon_inferred", v5_chr=snap[sib[0]]["v5_chr"], v5_pos=str((lo + hi) // 2),
             v5_start=str(lo), v5_end=str(hi), method=f"inferred from same-amplicon SNPs {','.join(sib)} (±750 bp)")
    r["note"] = "; ".join(x for x in [r["note"], "position NOT measured for this SNP; same Panzea amplicon as siblings"] if x)
    n["inferred"] += 1

w = csv.DictWriter(open(f"{d}/vanHeerwaarden2011_markers_v5.final.tsv", "w", newline=""), cols, delimiter="\t")
w.writeheader(); w.writerows(rows)
print(dict(n), dict(collections.Counter(r["precision"] for r in rows)), file=sys.stderr)
