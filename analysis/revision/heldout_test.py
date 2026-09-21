#!/usr/bin/env python3
"""heldout_test.py — held-out-library check of a catalogue built from a subset of
the RNA-seq libraries (Reviewer 1, point 3).

The catalogue was built from the training libraries only, so its RNA-seq flags and
tiers know nothing about the held-out libraries. Every coding intron (the gap
between consecutive CDS segments of a transcript) is looked up in the splice
junctions of the training and of the held-out libraries; an intron is supported by
a set when a junction with the same two splice sites carries >= MINREADS reads.

Reported:
  1. distinct coding introns of the catalogue and of the reference, split by
     training support, with the fraction the held-out libraries confirm
  2. multi-exon loci by reference status (matched / novel) x tier and x RNA flag:
     fraction with at least one held-out-confirmed intron, and with all introns
     confirmed

Usage: heldout_test.py <tiers.gff3> <reference.gff3> <heldout.tsv> <minreads>
                       <train.tsv> [<train.tsv> ...]
"""
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/analysis'))
sys.path.insert(0, '/data2/k821209/gm/gene-miner/bin')
sys.path.insert(0, '/data2/k821209/gm/gene-miner/analysis')
from gm_compare import parse_gff, index  # noqa: E402
from gm_overlap import ofrac  # noqa: E402


def junctions(paths, minreads):
    c = defaultdict(int)
    for p in paths:
        for line in open(p):
            ch, s, e, n = line.split('\t')
            c[(ch, int(s), int(e))] += int(n)
    return {k for k, n in c.items() if n >= minreads}


def coding_introns(gff, want_protein_coding=False):
    """gene -> set of (chrom, first intron base, last intron base), from CDS gaps."""
    t2g, cds = {}, defaultdict(list)
    pc_genes = None
    if want_protein_coding:
        pc_genes = set(parse_gff(gff, want_protein_coding=True)[0])
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or line.startswith('#'):
            continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        if f[2] in ('mRNA', 'transcript') and 'Parent' in a:
            t2g[a['ID']] = a['Parent']
        elif f[2] == 'CDS' and 'Parent' in a:
            for p in a['Parent'].split(','):
                cds[p].append((f[0], int(f[3]), int(f[4])))
    out = defaultdict(set)
    for t, segs in cds.items():
        g = t2g.get(t)
        if g is None or (pc_genes is not None and g not in pc_genes):
            continue
        segs.sort(key=lambda x: x[1])
        for (c, _, e1), (_, s2, _) in zip(segs, segs[1:]):
            if s2 - e1 > 1:
                out[g].add((c, e1 + 1, s2 - 1))
    return out


def status(gff, ref):
    gm, _ = parse_gff(gff)
    rf, _ = parse_gff(ref, want_protein_coding=True)
    ridx = index(rf)
    st = {}
    for gid, gv in gm.items():
        hit = False
        for rs, re_, rid in ridx.get((gv['chrom'], gv['strand']), []):
            if re_ < gv['start'] or rs > gv['end']:
                continue
            if ofrac(gv, rf[rid]) >= 0.5:
                hit = True
                break
        st[gid] = 'matched' if hit else 'novel'
    return st


def tiers(gff):
    tier, ev = {}, {}
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'gene':
            g = re.search(r'ID=([^;]+)', f[8]).group(1)
            t = re.search(r'tier=([^;]+)', f[8])
            e = re.search(r'evidence=([^;]+)', f[8])
            tier[g] = t.group(1) if t else 'NA'
            ev[g] = set() if not e or e.group(1) == 'none' else set(e.group(1).split(','))
    return tier, ev


def pct(a, b):
    return f'{a:,} of {b:,} ({100 * a / b:.1f}%)' if b else '0 of 0'


def main():
    gff, ref, held, minreads = sys.argv[1:5]
    minreads = int(minreads)
    train = junctions(sys.argv[5:], minreads)
    heldj = junctions([held], minreads)
    print(f'junctions >= {minreads} reads: training {len(train):,}, held-out {len(heldj):,}, '
          f'held-out not in training {len(heldj - train):,}')

    cat = coding_introns(gff)
    refi = coding_introns(ref, want_protein_coding=True)
    print('\n== distinct coding introns')
    print('set\tintrons\ttraining-supported\theld-out confirmed (all)\t'
          'held-out confirmed | training-supported\theld-out confirmed | no training support')
    for name, d in (('catalogue', cat), ('reference', refi)):
        allI = set().union(*d.values()) if d else set()
        tr = allI & train
        notr = allI - train
        print(f'{name}\t{len(allI):,}\t{pct(len(tr), len(allI))}\t{pct(len(allI & heldj), len(allI))}\t'
              f'{pct(len(tr & heldj), len(tr))}\t{pct(len(notr & heldj), len(notr))}')

    st = status(gff, ref)
    tier, ev = tiers(gff)
    print('\n== multi-exon loci: held-out confirmation per locus')
    print('status\tgroup\tmulti-exon loci\t>=1 intron confirmed\tall introns confirmed')
    groups = defaultdict(list)
    for g, ints in cat.items():
        if not ints:
            continue
        s = st.get(g, 'NA')
        for key in (('all', 'all'), ('tier', tier.get(g, 'NA')),
                    ('rna', 'rna' if 'rna' in ev.get(g, set()) else 'no-rna')):
            groups[(s, key[0], key[1])].append(ints)
    order = {'all': 0, 'tier': 1, 'rna': 2}
    torder = {'all': 0, 'high': 1, 'medium': 2, 'low': 3, 'rna': 4, 'no-rna': 5}
    for (s, kind, lab) in sorted(groups, key=lambda k: (k[0], order[k[1]], torder.get(k[2], 9))):
        L = groups[(s, kind, lab)]
        one = sum(1 for ints in L if ints & heldj)
        full = sum(1 for ints in L if ints <= heldj)
        print(f'{s}\t{lab}\t{len(L):,}\t{pct(one, len(L))}\t{pct(full, len(L))}')


if __name__ == '__main__':
    main()
