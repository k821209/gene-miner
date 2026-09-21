#!/usr/bin/env python3
"""fusion_origin.py — decompose Gene-Miner fusion loci (see split_fusion.py) by cause.

For each Gene-Miner gene whose merged CDS footprint holds >= 2 reference genes
(each >= F of its CDS bases inside), ask whether any SINGLE transcript does so:
  chimeric_tx   : at least one transcript alone spans >= 2 reference genes
  isoform_merge : no single transcript does; different isoforms of the locus
                  cover different reference genes (fusion made by the union's
                  isoform folding)
Each chimeric transcript is traced to its stream by exact CDS-interval match to
the AUGUSTUS, TransDecoder and GeneMark-ETP inputs of build_union.py.

Usage: fusion_origin.py ref.gff3 union.final.gff3 augustus.gff3 transdecoder.gff3 genemark.gtf [F]
"""
import re, sys
from collections import defaultdict, Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
from split_fusion import load, shared, pairs


def tx_cds_gff(path):
    t2g, cds, loc = {}, defaultdict(list), {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9:
            continue
        if c[2] in ('mRNA', 'transcript'):
            m = re.search(r'ID=([^;]+)', c[8]); p = re.search(r'Parent=([^;]+)', c[8])
            if m:
                t2g[m.group(1)] = p.group(1) if p else m.group(1)
        elif c[2] == 'CDS':
            p = re.search(r'Parent=([^;]+)', c[8])
            if p:
                t = p.group(1).split(',')[0]
                cds[t].append((int(c[3]), int(c[4])))
                loc[t] = (c[0], c[6])
    return t2g, cds, loc


def keyset_gff(path):
    _, cds, loc = tx_cds_gff(path)
    return {(loc[t][0], loc[t][1], tuple(sorted(set(v)))) for t, v in cds.items()}


def keyset_gtf(path):
    cds, loc = defaultdict(list), {}
    for line in open(path):
        c = line.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'CDS':
            continue
        t = re.search(r'transcript_id "([^"]+)"', c[8]).group(1)
        cds[t].append((int(c[3]), int(c[4]))); loc[t] = (c[0], c[6])
    return {(loc[t][0], loc[t][1], tuple(sorted(set(v)))) for t, v in cds.items()}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main():
    refp, gmp, augp, tdp, gmtf = sys.argv[1:6]
    F = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
    ref, gm = load(refp), load(gmp)
    r_in_q = defaultdict(set)
    for q, r, ov in pairs(ref, gm):
        if ov >= F * ref[r]['len']:
            r_in_q[q].add(r)
    fused = {q: rs for q, rs in r_in_q.items() if len(rs) >= 2}

    t2g, tcds, tloc = tx_cds_gff(gmp)
    g2t = defaultdict(list)
    for t in tcds:
        g2t[t2g.get(t, t)].append(t)
    src = {'AUGUSTUS': keyset_gff(augp), 'RNA-seq': keyset_gff(tdp), 'GeneMark-ETP': keyset_gtf(gmtf)}

    kind, origin = Counter(), Counter()
    for q, rs in fused.items():
        chim = []
        for t in g2t[q]:
            iv = merge(tcds[t]); ln = sum(e - s + 1 for s, e in iv)
            hit = [r for r in rs if shared(iv, ref[r]['iv']) >= F * ref[r]['len']]
            if len(hit) >= 2:
                chim.append(t)
        if chim:
            kind['chimeric_tx'] += 1
            for t in chim:
                k = (tloc[t][0], tloc[t][1], tuple(sorted(set(tcds[t]))))
                names = [n for n, ks in src.items() if k in ks]
                origin['+'.join(names) if names else 'unmatched'] += 1
        else:
            kind['isoform_merge'] += 1
    print('fused_loci', len(fused), dict(kind))
    print('chimeric_transcript_origin', dict(origin))


if __name__ == '__main__':
    main()
