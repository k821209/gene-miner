#!/usr/bin/env python3
"""split_fusion.py — gene-boundary errors of a query annotation against a reference.

Each gene is reduced to its merged CDS footprint (the union of coding bases over
all its isoforms, same strand). For a query gene q and reference gene r,
ov(q,r) = number of coding bases shared.

  split  : reference gene r with >= 2 query genes q such that ov(q,r) >= F * |q|
           (two or more query genes each lie mostly inside r -> r was broken up)
  fusion : query gene q with >= 2 reference genes r such that ov(q,r) >= F * |r|
           (two or more reference genes each lie mostly inside q -> q merges them)

F defaults to 0.5. Split rate is over reference genes, fusion rate over query genes.

Usage: split_fusion.py ref.gff3 query.gff3 [F]   -> one TSV line on stdout
"""
import re, sys
from collections import defaultdict


def load(path):
    t2g, cds = {}, defaultdict(list)
    info = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9:
            continue
        if c[2] in ('mRNA', 'transcript'):
            tid = re.search(r'ID=([^;]+)', c[8]); par = re.search(r'Parent=([^;]+)', c[8])
            if tid:
                t2g[tid.group(1)] = par.group(1) if par else tid.group(1)
        elif c[2] == 'CDS':
            par = re.search(r'Parent=([^;]+)', c[8])
            if not par:
                continue
            tid = par.group(1).split(',')[0]
            cds[tid].append((int(c[3]), int(c[4])))
            info[tid] = (c[0], c[6])
    genes = defaultdict(list)
    ginfo = {}
    for tid, iv in cds.items():
        g = t2g.get(tid, tid)
        genes[g].extend(iv)
        ginfo[g] = info[tid]
    out = {}
    for g, iv in genes.items():
        iv.sort()
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[g] = dict(chrom=ginfo[g][0], strand=ginfo[g][1], iv=merged,
                      start=merged[0][0], end=merged[-1][1],
                      len=sum(e - s + 1 for s, e in merged))
    return out


def shared(a, b):
    i = j = ov = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi >= lo:
            ov += hi - lo + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return ov


def pairs(ref, qry):
    idx = defaultdict(list)
    for g, v in ref.items():
        idx[(v['chrom'], v['strand'])].append((v['start'], v['end'], g))
    for k in idx:
        idx[k].sort()
    for q, qv in qry.items():
        for rs, re_, r in idx.get((qv['chrom'], qv['strand']), []):
            if rs > qv['end']:
                break
            if re_ < qv['start']:
                continue
            ov = shared(qv['iv'], ref[r]['iv'])
            if ov:
                yield q, r, ov


def main():
    ref, qry = load(sys.argv[1]), load(sys.argv[2])
    F = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    q_in_r = defaultdict(set)   # r -> query genes mostly inside r
    r_in_q = defaultdict(set)   # q -> ref genes mostly inside q
    for q, r, ov in pairs(ref, qry):
        if ov >= F * qry[q]['len']:
            q_in_r[r].add(q)
        if ov >= F * ref[r]['len']:
            r_in_q[q].add(r)
    split_ref = [r for r, s in q_in_r.items() if len(s) >= 2]
    fusion_q = [q for q, s in r_in_q.items() if len(s) >= 2]
    fused_ref = set(r for q in fusion_q for r in r_in_q[q])
    print('\t'.join(map(str, [
        len(ref), len(qry), F,
        len(split_ref), f"{100*len(split_ref)/len(ref):.2f}",
        sum(len(q_in_r[r]) for r in split_ref),
        len(fusion_q), f"{100*len(fusion_q)/len(qry):.2f}",
        len(fused_ref), f"{100*len(fused_ref)/len(ref):.2f}"])))


if __name__ == '__main__':
    main()
