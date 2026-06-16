#!/usr/bin/env python3
"""filter_taxonomy.py — drop residual non-eukaryotic (bacterial) genes from a
gene set, using eggNOG-mapper taxonomy. Whole-body insect assemblies carry gut
microbiome; contig-level decontamination removes bacterial *contigs* but smaller
/low-GC ones leak through and produce bacterial gene models. A gene is kept only
if it has at least one eukaryotic-rooted eggNOG orthogroup.

Usage: filter_taxonomy.py emapper.annotations in.gff3 in.pep.fa out.gff3 out.pep.fa
"""
import sys, re

ann, in_gff, in_pep, out_gff, out_pep = sys.argv[1:6]
GID = re.compile(r"ZscG_[AR]\d+")                      # adapt to your gene-ID scheme
EUK = {"2759", "33208", "6656", "50557", "7147", "33392", "7214", "7227", "33340"}
BAC = {"2", "1224", "1236", "91347", "543", "1239", "201174", "976", "1760"}


def klass(ogs):
    if not ogs or ogs == "-":
        return "noOG"
    tids = {t.split("@", 1)[1].split("|", 1)[0] for t in ogs.split(",") if "@" in t}
    if tids & EUK:
        return "euk"
    if tids & BAC:
        return "bac"
    return "other"


# bacterial gene IDs from eggNOG
bac = set()
H = []
for line in open(ann):
    if line.startswith("#"):
        if line.startswith("#query"):
            H = line.lstrip("#").rstrip("\n").split("\t")
        continue
    r = dict(zip(H, line.rstrip("\n").split("\t")))
    m = GID.search(r.get("query", ""))
    if m and klass(r.get("eggNOG_OGs", "-")) == "bac":
        bac.add(m.group())

# filter GFF (drop every feature line of a bacterial gene)
g_in = g_out = 0
with open(out_gff, "w") as o:
    for line in open(in_gff):
        if line.startswith("#"):
            o.write(line); continue
        c = line.split("\t")
        if len(c) > 2 and c[2] == "gene":
            g_in += 1
        m = GID.search(line)
        if m and m.group() in bac:
            continue
        if len(c) > 2 and c[2] == "gene":
            g_out += 1
        o.write(line)

# filter proteins
keep = True
p_in = p_out = 0
with open(out_pep, "w") as o:
    for line in open(in_pep):
        if line.startswith(">"):
            p_in += 1
            m = GID.search(line)
            keep = not (m and m.group() in bac)
            if keep:
                p_out += 1
        if keep:
            o.write(line)

print(f"bacterial genes dropped: {len(bac)}")
print(f"genes:    {g_in} -> {g_out}")
print(f"proteins: {p_in} -> {p_out}")
