#!/usr/bin/env python3
"""Align published SNP flanks ([A/B] notation) to B73 v5 and report the exact SNP base.
Usage: pip install mappy; python3 align_flanks.py Zm-B73-REFERENCE-NAM-5.0.fa.gz flanks_to_align.fa > flanks_v5.tsv"""
import re, sys
import mappy as mp

def snp_ref_pos(h, qlen, off):
    if not (h.q_st <= off < h.q_en): return None
    qpos, rpos = (h.q_st, h.r_st) if h.strand == 1 else (qlen - h.q_en, h.r_st)
    t = off if h.strand == 1 else qlen - 1 - off
    for L, op in h.cigar:
        if op in (0, 7, 8):
            if qpos <= t < qpos + L: return rpos + t - qpos
            qpos += L; rpos += L
        elif op == 1:
            if qpos <= t < qpos + L: return None
            qpos += L
        elif op in (2, 3): rpos += L
    return None

genome, fa = sys.argv[1], sys.argv[2]
al = mp.Aligner(genome, preset="sr")
print("marker\tv5_chr\tv5_pos\tstrand\tmapq\tidentity\tn_primary\tv5_B73_base\talleles_marker_strand\tstatus")
for name, seq, _ in mp.fastx_read(fa):
    m = re.search(r"\[([^\]]+)\]", seq)
    left, right, alleles = seq[:m.start()], seq[m.end():], m.group(1)
    q = re.sub(r"[^ACGT]", "N", (left + alleles.split("/")[0] + right).upper())
    hits = [h for h in al.map(q) if h.is_primary]
    if not hits:
        print(f"{name}\t\t\t\t\t\t0\t\t{alleles}\tunmapped"); continue
    h = max(hits, key=lambda x: x.mlen)
    p = snp_ref_pos(h, len(q), len(left))
    ident = h.mlen / max(h.blen, 1)
    status = []
    if p is None: status.append("snp_not_in_alignment")
    if h.mapq < 20: status.append("low_mapq")
    if len(hits) > 1: status.append("multiple_hits")
    if ident < 0.95: status.append("identity<0.95")
    base = al.seq(h.ctg, p, p + 1) if p is not None else ""
    bases = set(re.sub(r"[^ACGT]", "", alleles.upper()))
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    if base and base not in (bases if h.strand == 1 else {comp[b] for b in bases}):
        status.append("B73_base_not_an_allele")
    print(f"{name}\t{h.ctg}\t{'' if p is None else p + 1}\t{'+' if h.strand == 1 else '-'}\t{h.mapq}\t{ident:.3f}\t{len(hits)}\t{base}\t{alleles}\t{';'.join(status) or 'ok'}")
