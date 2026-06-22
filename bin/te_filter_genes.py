#!/usr/bin/env python3
"""te_filter_genes.py — drop transposable-element-derived genes from a gene set,
using RepeatMasker interspersed-repeat coordinates (matches what reference
annotation pipelines do, e.g. excluding TE-overlapping genes). A gene is removed
if >= THRESH of its CDS bases lie inside interspersed repeats (simple/low-
complexity/satellite repeats are NOT used — real genes carry those).

Usage: te_filter_genes.py rm.out in.gff3 in.pep.fa out.gff3 out.pep.fa [thresh=0.5]
"""
import sys, re
from bisect import bisect_right

rm_out, in_gff, in_pep, out_gff, out_pep = sys.argv[1:6]
THRESH = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
GID = re.compile(r"[A-Za-z][A-Za-z0-9]*_[AER]\d+")   # any <prefix>_[AR]<digits> gene ID
SKIP = ("Simple_repeat", "Low_complexity", "Satellite", "Simple", "Low_comp",
        "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA")

# 1. interspersed-repeat intervals per scaffold
ivs = {}
with open(rm_out) as fh:
    for line in fh:
        c = line.split()
        if len(c) < 11 or not c[0].replace(".", "").isdigit():
            continue                                  # header / blank
        cls = c[10]
        if cls.startswith(SKIP):
            continue
        scaf, b, e = c[4], int(c[5]), int(c[6])
        ivs.setdefault(scaf, []).append((b, e))

# merge per scaffold -> starts[], ends[] for bisect
merged = {}
for scaf, lst in ivs.items():
    lst.sort()
    ms = []
    for b, e in lst:
        if ms and b <= ms[-1][1]:
            ms[-1][1] = max(ms[-1][1], e)
        else:
            ms.append([b, e])
    merged[scaf] = ([m[0] for m in ms], [m[1] for m in ms])


def te_overlap(scaf, b, e):
    if scaf not in merged:
        return 0
    starts, ends = merged[scaf]
    i = bisect_right(starts, e)                       # intervals starting <= e
    ov = 0
    j = i - 1
    while j >= 0 and ends[j] >= b:
        ov += max(0, min(e, ends[j]) - max(b, starts[j]) + 1)
        j -= 1
    return ov

# 2. CDS per gene
cds = {}        # gene -> list of (scaf,b,e)
for line in open(in_gff):
    if line.startswith("#"):
        continue
    c = line.split("\t")
    if len(c) < 9 or c[2] != "CDS":
        continue
    m = GID.search(c[8])
    if not m:
        continue
    cds.setdefault(m.group(), []).append((c[0], int(c[3]), int(c[4])))

te_genes = set()
for g, segs in cds.items():
    tot = sum(e - b + 1 for _, b, e in segs)
    ov = sum(te_overlap(s, b, e) for s, b, e in segs)
    if tot and ov / tot >= THRESH:
        te_genes.add(g)

# 3. write filtered gff + pep
g_in = g_out = 0
with open(out_gff, "w") as o:
    for line in open(in_gff):
        if line.startswith("#"):
            o.write(line); continue
        c = line.split("\t")
        if len(c) > 2 and c[2] == "gene":
            g_in += 1
        m = GID.search(line)
        if m and m.group() in te_genes:
            continue
        if len(c) > 2 and c[2] == "gene":
            g_out += 1
        o.write(line)

keep = True
p_in = p_out = 0
with open(out_pep, "w") as o:
    for line in open(in_pep):
        if line.startswith(">"):
            p_in += 1
            m = GID.search(line)
            keep = not (m and m.group() in te_genes)
            if keep:
                p_out += 1
        if keep:
            o.write(line)

print(f"TE threshold (CDS fraction in interspersed repeats): {THRESH}")
print(f"TE-derived genes dropped: {len(te_genes)}")
print(f"genes:    {g_in} -> {g_out}")
print(f"proteins: {p_in} -> {p_out}")
