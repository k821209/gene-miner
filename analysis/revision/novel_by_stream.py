#!/usr/bin/env python3
"""novel_by_stream.py — novel / recovered loci split by the stream that founded
them (gene IDs are <prefix>_A / _R / _E). Usage: novel_by_stream.py gm.gff3 ref.gff3"""
import os, re, sys
from collections import Counter
sys.path.insert(0, os.path.expanduser("~/gene-miner-runs/gene-miner/bin"))
from gm_compare import parse_gff, index, overlaps
gm, _ = parse_gff(sys.argv[1]); ref, _ = parse_gff(sys.argv[2], want_protein_coding=True)
ridx = index(ref)
tag = {"A": "AUGUSTUS", "R": "RNA-seq", "E": "GeneMark-ETP"}
novel, rec = Counter(), Counter()
for gid, gv in gm.items():
    best = 0.0
    for rs, re_, rid in ridx.get((gv["chrom"], gv["strand"]), []):
        if re_ < gv["start"] or rs > gv["end"]:
            continue
        o = overlaps(gv, ref[rid])
        if o > best:
            best = o
    m = re.search(r"_([ARE])\d+$", gid)
    k = tag.get(m.group(1) if m else "?", "?")
    (rec if best >= 0.5 else novel)[k] += 1
print("novel   ", dict(novel), sum(novel.values()))
print("matching", dict(rec), sum(rec.values()))
