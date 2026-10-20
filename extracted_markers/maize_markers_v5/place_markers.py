#!/usr/bin/env python3
"""
Place maize SNP markers (Panzea PZA/PZB/PZD and gene-named, e.g. ae1.5) on B73 v5
by aligning their MaizeGDB context sequences with minimap2 (via mappy).

Inputs (same folder):
  names.txt                         one marker name per line
  MaizeGDB_genetic_SNPs.fa          from download.maizegdb.org/SNPs_and_Maps/SNPs/
  MaizeGDB_genetic_SNPs.xref        (optional) same place; v4 positions used as a cross-check
  Zm-B73-REFERENCE-NAM-5.0.fa(.gz)  B73 v5 genome from download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/

Output: markers_v5.tsv

  pip install mappy
  python3 place_markers.py [genome.fa.gz]
"""
import gzip, re, sys, os, csv
import mappy as mp

IUPAC = set("RYKMSWBDHVN")

def opn(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)

def read_fasta(fn):
    name, buf = None, []
    with opn(fn) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)

def parse_context(seq):
    """Return (clean_seq, snp_offset or None, alleles string)."""
    m = re.search(r"\[([^\]]+)\]", seq)
    if m:
        alleles = m.group(1)
        left, right = seq[:m.start()], seq[m.end():]
        a1 = alleles.split("/")[0].replace("-", "")
        a1 = a1 if a1 else ""
        clean = left + (a1 if a1 else "") + right
        return re.sub(r"[^ACGTNacgtn]", "N", clean).upper(), len(left), alleles
    up = seq.upper()
    amb = [i for i, c in enumerate(up) if c in IUPAC and c != "N"]
    if len(amb) == 1:
        return re.sub(r"[^ACGTN]", "N", up), amb[0], up[amb[0]]
    return re.sub(r"[^ACGTN]", "N", up), None, ""

def norm(s):
    return s.strip().lower()

def main():
    genome = sys.argv[1] if len(sys.argv) > 1 else next(
        (f for f in ("Zm-B73-REFERENCE-NAM-5.0.fa.gz", "Zm-B73-REFERENCE-NAM-5.0.fa") if os.path.exists(f)), None)
    if not genome:
        sys.exit("Genome FASTA not found; pass its path as the first argument.")
    names = [l.strip() for l in open("names.txt") if l.strip()]
    want = {norm(n): n for n in names}

    # marker sequences: match on any whitespace/|-separated token of the header
    seqs = {}
    for hdr, s in read_fasta("MaizeGDB_genetic_SNPs.fa"):
        for tok in re.split(r"[\s|;,]+", hdr):
            k = norm(tok)
            if k in want and want[k] not in seqs:
                seqs[want[k]] = (hdr, s)

    v4 = {}
    if os.path.exists("MaizeGDB_genetic_SNPs.xref"):
        with open("MaizeGDB_genetic_SNPs.xref") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for row in csv.reader(f, delimiter="\t"):
                for tok in row[:2]:
                    if norm(tok) in want:
                        v4[want[norm(tok)]] = dict(zip(hdr, row))

    idx_fn = "B73v5.sr.mmi"
    if os.path.exists(idx_fn):
        aligner = mp.Aligner(idx_fn, preset="sr")
    else:
        print("Building minimap2 index (one-time, several minutes)...", file=sys.stderr)
        aligner = mp.Aligner(genome, preset="sr", fn_idx_out=idx_fn)
    if not aligner:
        sys.exit("Failed to load/build index")

    out = open("markers_v5.tsv", "w")
    cols = ["marker", "chr_v5", "pos_v5", "strand", "mapq", "aln_identity", "aln_len", "query_len",
            "n_hits", "snp_in_seq", "alleles", "v5_ref_base", "status", "v4_info", "fasta_header"]
    out.write("\t".join(cols) + "\n")
    n_ok = 0
    for n in names:
        if n not in seqs:
            out.write("\t".join([n] + [""] * 11 + ["no_sequence_found", str(v4.get(n, "")), ""]) + "\n")
            continue
        hdr, raw = seqs[n]
        q, off, alleles = parse_context(raw)
        hits = [h for h in aligner.map(q)]
        prim = [h for h in hits if h.is_primary]
        if not prim:
            out.write("\t".join([n] + [""] * 8 + [str(off is not None), alleles, "", "unmapped",
                                   str(v4.get(n, "")), hdr]) + "\n")
            continue
        h = max(prim, key=lambda x: x.mlen)
        ident = h.mlen / h.blen if h.blen else 0
        pos, refb = "", ""
        status = "ok"
        if off is not None and h.q_st <= off < h.q_en:
            # walk CIGAR to translate query offset -> reference coordinate
            if h.strand == 1:
                qpos, rpos = h.q_st, h.r_st
                target = off
            else:
                # for reverse strand, mappy CIGAR is in reference order; query coords run backwards
                qpos, rpos = len(q) - h.q_en, h.r_st
                target = len(q) - 1 - off
            found = None
            for length, op in h.cigar:
                if op in (0, 7, 8):  # M,=,X
                    if qpos <= target < qpos + length:
                        found = rpos + (target - qpos); break
                    qpos += length; rpos += length
                elif op == 1:  # I
                    if qpos <= target < qpos + length:
                        found = rpos; status = "snp_in_insertion"; break
                    qpos += length
                elif op in (2, 3):  # D,N
                    rpos += length
            if found is not None:
                pos = found + 1  # 1-based
                refb = aligner.seq(h.ctg, found, found + 1) or ""
        elif off is None:
            status = "snp_offset_unknown(aln_start_given)"
            pos = h.r_st + 1
        else:
            status = "snp_outside_alignment"
        if len(prim) > 1 or (h.mapq < 20):
            status += ";multi_or_lowmapq"
        if ident < 0.9:
            status += ";low_identity"
        n_ok += status.startswith("ok")
        out.write("\t".join(map(str, [n, h.ctg, pos, "+" if h.strand == 1 else "-", h.mapq, f"{ident:.3f}",
                                      h.blen, len(q), len(hits), off is not None, alleles, refb, status,
                                      v4.get(n, ""), hdr])) + "\n")
    out.close()
    print(f"{len(names)} markers; {len(seqs)} with sequence; {n_ok} placed cleanly -> markers_v5.tsv",
          file=sys.stderr)

if __name__ == "__main__":
    main()
