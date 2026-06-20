#!/usr/bin/env python3
"""dedup_gff_genes.py — collapse duplicate/overlapping gene models that arise
when AUGUSTUS is run over OVERLAPPING genome windows (a gene fully inside an
overlap is predicted in both adjacent windows). Genes on the same (chrom,strand)
whose CDS footprints overlap reciprocally >= THRESH are collapsed to the single
model with the longest total CDS. Preserves all child feature lines of kept genes.

Usage: dedup_gff_genes.py <in.gff3> <out.gff3> [min_reciprocal_overlap=0.8]
"""
import sys, re
from collections import defaultdict

inp, outp = sys.argv[1], sys.argv[2]
THRESH = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8

def attr(a, k):
    m = re.search(k + r'=([^;]+)', a)
    return m.group(1) if m else None

lines = open(inp).read().splitlines()
recs = []
for i, l in enumerate(lines):
    if not l or l.startswith('#'):
        continue
    f = l.split('\t')
    if len(f) >= 9:
        recs.append((i, f))

# map each line to its gene id
tx_gene = {}
for i, f in recs:
    if f[2] in ('mRNA', 'transcript'):
        tx_gene[attr(f[8], 'ID')] = attr(f[8], 'Parent')

line_gene = {}
cds = defaultdict(list)
gchrom, gstrand = {}, {}
for i, f in recs:
    t, a = f[2], f[8]
    if t == 'gene':
        g = attr(a, 'ID')
    elif t in ('mRNA', 'transcript'):
        g = tx_gene.get(attr(a, 'ID'))
    else:
        g = tx_gene.get(attr(a, 'Parent'))
    line_gene[i] = g
    if t == 'CDS' and g:
        cds[g].append((int(f[3]), int(f[4])))
        gchrom[g] = f[0]; gstrand[g] = f[6]

def clen(iv):
    return sum(e - s + 1 for s, e in iv)

genes = {}
for g, iv in cds.items():
    iv = sorted(set(iv))
    genes[g] = dict(chrom=gchrom[g], strand=gstrand[g], cds=iv,
                    start=min(s for s, e in iv), end=max(e for s, e in iv), L=clen(iv))

# group by (chrom,strand), sort by start, cluster overlapping, keep longest CDS
def recip(a, b):
    ov = 0
    for s, e in a['cds']:
        for bs, be in b['cds']:
            lo, hi = max(s, bs), min(e, be)
            if hi >= lo:
                ov += hi - lo + 1
    return ov / (min(a['L'], b['L']) or 1)

drop = set()
by_cs = defaultdict(list)
for g, v in genes.items():
    by_cs[(v['chrom'], v['strand'])].append(g)
for key, gl in by_cs.items():
    gl.sort(key=lambda g: genes[g]['start'])
    for i in range(len(gl)):
        if gl[i] in drop:
            continue
        a = genes[gl[i]]
        for j in range(i + 1, len(gl)):
            b = genes[gl[j]]
            if b['start'] > a['end']:
                break
            if gl[j] in drop:
                continue
            if recip(a, b) >= THRESH:
                # keep the longer; drop the shorter
                if b['L'] > a['L']:
                    drop.add(gl[i]); a = b
                else:
                    drop.add(gl[j])

kept = 0
with open(outp, 'w') as o:
    for i, l in enumerate(lines):
        if not l or l.startswith('#'):
            o.write(l + '\n'); continue
        g = line_gene.get(i)
        if g in drop:
            continue
        o.write(l + '\n')
n_total = len(genes); n_drop = len(drop)
sys.stderr.write(f"dedup: {n_total} genes, dropped {n_drop} window-overlap duplicates, kept {n_total - n_drop}\n")
