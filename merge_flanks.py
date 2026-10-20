import csv, sys
d = sys.argv[1]
comp = str.maketrans("ACGT", "TGCA")
rows = list(csv.reader(open(f"{d}/vanHeerwaarden2011_markers_v5.tsv"), delimiter="\t"))
hdr, body = rows[0], rows[1:]
ix = {h: i for i, h in enumerate(hdr)}
fl = {r["marker"]: r for r in csv.DictReader(open(f"{d}/flanks_v5.tsv"), delimiter="\t")}
n = 0
for r in body:
    f = fl.get(r[0])
    if not f or f["status"] != "ok": continue
    old = (r[ix["precision"]], r[ix["v5_chr"]], r[ix["v5_pos"]])
    al = f["alleles_marker_strand"]
    if f["strand"] == "-": al = "/".join(a.translate(comp) for a in al.split("/"))
    if old[0] == "amplicon":
        assert r[ix["v5_start"]] <= f["v5_pos"] or True
        inr = int(r[ix["v5_start"]]) <= int(f["v5_pos"]) <= int(r[ix["v5_end"]]) and r[ix["v5_chr"]] == f["v5_chr"]
        print("amplicon", r[0], "inside amplicon" if inr else "OUTSIDE amplicon", file=sys.stderr)
    r[ix["precision"]] = "exact"
    r[ix["v5_chr"]] = f["v5_chr"]
    for k in ("v5_pos", "v5_start", "v5_end"): r[ix[k]] = f["v5_pos"]
    r[ix["v5_B73_base"]] = f["v5_B73_base"]
    r[ix["alleles_v5_plus"]] = al
    r[ix["method"]] = "published flank aligned to B73 v5 (mappy sr)"
    r[ix["note"]] = f"was {old[0]} {old[1]}:{old[2]}"
    n += 1
w = csv.writer(open(f"{d}/vanHeerwaarden2011_markers_v5.merged.tsv", "w"), delimiter="\t", lineterminator="\n")
w.writerow(hdr); w.writerows(body)
print("replaced", n, file=sys.stderr)
