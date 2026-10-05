#!/usr/bin/env python3
"""pick_boundary_examples.py — list compact gene-boundary errors worth drawing.

Same footprint rule as split_fusion.py. For each candidate prints coordinates,
the genes involved, whether one transcript alone spans the fused pair and which
input stream that transcript matches, and how many spliced-read junctions fall in
the window, so an example can be chosen for the per-locus figure.

Usage: pick_boundary_examples.py ref.gff3 union.gff3 junctions.tsv
                                 [max_span=30000] [augustus.gff3 transdecoder.gff3 genemark.gtf]
"""
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/analysis/revision'))
from split_fusion import load, shared, pairs  # noqa: E402

F = 0.5


def tx_cds(path, gtf=False):
    t2g, cds, loc = {}, defaultdict(list), {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9:
            continue
        if gtf:
            if c[2] != 'CDS':
                continue
            m = re.search(r'transcript_id "([^"]+)"', c[8])
            if not m:
                continue
            t = m.group(1)
        else:
            if c[2] in ('mRNA', 'transcript'):
                a = re.search(r'ID=([^;]+)', c[8])
                p = re.search(r'Parent=([^;]+)', c[8])
                if a:
                    t2g[a.group(1)] = p.group(1) if p else a.group(1)
                continue
            if c[2] != 'CDS':
                continue
            p = re.search(r'Parent=([^;]+)', c[8])
            if not p:
                continue
            t = p.group(1).split(',')[0]
        cds[t].append((int(c[3]), int(c[4])))
        loc[t] = (c[0], c[6])
    return {t: sorted(v) for t, v in cds.items()}, loc, t2g


def main():
    ref_gff, qry_gff, junc = sys.argv[1:4]
    max_span = int(sys.argv[4]) if len(sys.argv) > 4 else 30000
    streams = sys.argv[5:8]

    ref, qry = load(ref_gff), load(qry_gff)
    q_in_r, r_in_q = defaultdict(set), defaultdict(set)
    for q, r, ov in pairs(ref, qry):
        if ov >= F * qry[q]['len']:
            q_in_r[r].add(q)
        if ov >= F * ref[r]['len']:
            r_in_q[q].add(r)

    jun = defaultdict(list)
    for line in open(junc):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 4:
            jun[f[0]].append((int(f[1]), int(f[2]), int(f[3])))

    qtx, _, qt2g = tx_cds(qry_gff)
    stream_tx = {}
    for name, path in zip(('AUGUSTUS', 'RNA-seq', 'GeneMark-ETP'), streams):
        if path and os.path.exists(path):
            c, _, _ = tx_cds(path, gtf=path.endswith('.gtf'))
            stream_tx[name] = {tuple(v) for v in c.values()}

    def njunc(chrom, s, e, minreads=3):
        return sum(1 for a, b, n in jun.get(chrom, []) if n >= minreads and a >= s and b <= e)

    print('kind\tchrom\tstrand\tstart\tend\tspan\tquery\treference_genes\tchimeric_tx\tstream\tjunctions')
    rows = []
    for q, rs in r_in_q.items():
        if len(rs) != 2:
            continue
        v = qry[q]
        span = v['end'] - v['start']
        if span > max_span:
            continue
        chim, stream = 'no', '.'
        for t, segs in qtx.items():
            if qt2g.get(t, t) != q:
                continue
            if sum(1 for r in rs if shared(segs, ref[r]['iv']) >= F * ref[r]['len']) >= 2:
                chim = t
                stream = '+'.join(n for n, s in stream_tx.items() if tuple(segs) in s) or 'unmatched'
                break
        rows.append(('fusion', v, span, q, sorted(rs), chim, stream))
    for r, qs in q_in_r.items():
        if len(qs) != 2:
            continue
        v = ref[r]
        span = v['end'] - v['start']
        if span > max_span:
            continue
        rows.append(('split', v, span, r, sorted(qs), '.', '.'))

    for kind, v, span, g, others, chim, stream in sorted(
            rows, key=lambda x: (x[0], x[2])):
        print(f'{kind}\t{v["chrom"]}\t{v["strand"]}\t{v["start"]}\t{v["end"]}\t{span}\t{g}\t'
              f'{",".join(others)}\t{chim}\t{stream}\t{njunc(v["chrom"], v["start"], v["end"])}')


if __name__ == '__main__':
    main()
