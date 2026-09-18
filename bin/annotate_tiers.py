#!/usr/bin/env python3
"""annotate_tiers.py — write the independent evidence behind every locus into the
catalogue itself, so downstream users can filter it (Reviewer 2, major point 1).

Three independent lines of evidence per locus:
  rna     CDS overlaps an RNA-seq-assembled (TransDecoder) transcript
  sprot   DIAMOND blastp hit to UniProtKB/Swiss-Prot at E < 1e-5
  eggnog  an eggNOG ortholog group was assigned
Tier: high (>= 2), medium (1), low (none, i.e. ab initio structure only).

Each gene and mRNA line gains  evidence=rna,sprot;tier=high  (evidence=none for
the low tier). With a reference annotation, loci are also tagged
status=novel|matched so a user can select "novel and high-confidence" directly.

Usage: annotate_tiers.py <in.gff3> <out.gff3> --td transdecoder.gff3
                         --eggnog emapper.annotations --sprot blastp.outfmt6
                         [--ref reference.gff3] [--summary out.tsv]
"""
import argparse
import os
import re
import sys
from collections import defaultdict

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), 'analysis'))
from gm_compare import parse_gff, index  # noqa: E402
from gm_overlap import ofrac, merge, shared  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('gff'); ap.add_argument('out')
ap.add_argument('--td', required=True)
ap.add_argument('--eggnog', required=True)
ap.add_argument('--sprot', required=True)
ap.add_argument('--ref')
ap.add_argument('--summary')
a = ap.parse_args()

GENE_ID = re.compile(r'ID=([^;]+)')
PARENT = re.compile(r'Parent=([^;]+)')


def gene_of_transcript(path):
    t2g = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] in ('mRNA', 'transcript'):
            t2g[GENE_ID.search(f[8]).group(1)] = PARENT.search(f[8]).group(1)
    return t2g


gm, _ = parse_gff(a.gff)
t2g = gene_of_transcript(a.gff)

# --- RNA-seq support: merged TransDecoder CDS per (chrom, strand) --------------
td = defaultdict(list)
for line in open(a.td):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) > 8 and f[2] == 'CDS':
        td[(f[0], f[6])].append((int(f[3]), int(f[4])))
td = {k: merge(v) for k, v in td.items()}

rna = set()
for gid, gv in gm.items():
    m = td.get((gv['chrom'], gv['strand']))
    if m and shared(merge(gv['cds']), m) > 0:
        rna.add(gid)

# --- homology and orthology ---------------------------------------------------
sprot, eggnog = set(), set()
for line in open(a.sprot):
    q = line.split('\t', 1)[0]
    sprot.add(t2g.get(q, q))
H = []
for line in open(a.eggnog):
    if line.startswith('#'):
        if line.startswith('#query'):
            H = line.lstrip('#').rstrip('\n').split('\t')
        continue
    r = dict(zip(H, line.rstrip('\n').split('\t')))
    ogs = r.get('eggNOG_OGs', '-')
    if ogs and ogs != '-':
        q = r.get('query', '')
        eggnog.add(t2g.get(q, q))

# --- novel vs matched ---------------------------------------------------------
status = {}
if a.ref:
    ref, _ = parse_gff(a.ref, want_protein_coding=True)
    ridx = index(ref)
    for gid, gv in gm.items():
        best = 0.0
        for rs, re_, rid in ridx.get((gv['chrom'], gv['strand']), []):
            if re_ < gv['start'] or rs > gv['end']:
                continue
            best = max(best, ofrac(gv, ref[rid]))
        status[gid] = 'matched' if best >= 0.5 else 'novel'

tiers, ev_of = {}, {}
for gid in gm:
    ev = [n for n, s in (('rna', rna), ('sprot', sprot), ('eggnog', eggnog)) if gid in s]
    ev_of[gid] = ev
    tiers[gid] = 'high' if len(ev) >= 2 else 'medium' if ev else 'low'

with open(a.out, 'w') as o:
    for line in open(a.gff):
        if line.startswith('#') or not line.strip():
            o.write(line); continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            o.write(line); continue
        if f[2] == 'gene':
            g = GENE_ID.search(f[8]).group(1)
        elif f[2] in ('mRNA', 'transcript'):
            g = PARENT.search(f[8]).group(1)
        else:
            o.write(line); continue
        if g in tiers:
            extra = f";evidence={','.join(ev_of[g]) or 'none'};tier={tiers[g]}"
            if g in status:
                extra += f";status={status[g]}"
            f[8] = f[8].rstrip(';') + extra
        o.write('\t'.join(f) + '\n')

# --- summary ------------------------------------------------------------------
rows = defaultdict(lambda: defaultdict(int))
for gid in gm:
    rows[status.get(gid, 'all')][tiers[gid]] += 1
    rows[status.get(gid, 'all')]['total'] += 1
    for e in ev_of[gid]:
        rows[status.get(gid, 'all')][e] += 1
out = open(a.summary, 'w') if a.summary else sys.stdout
out.write('set\tloci\thigh\tmedium\tlow\trna\tsprot\teggnog\n')
for k, c in sorted(rows.items()):
    out.write(f"{k}\t{c['total']}\t{c['high']}\t{c['medium']}\t{c['low']}\t"
              f"{c['rna']}\t{c['sprot']}\t{c['eggnog']}\n")
if a.summary:
    out.close()
    print(open(a.summary).read())
