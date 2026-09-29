#!/usr/bin/env python3
"""paired_codon.py — is the 'revised' class less codon-biased because the genes it
revises are themselves less biased, or because revising their boundaries diluted
the coding signal?

Pairs each catalogue locus with the reference gene it matches (gm_compare's rule:
same strand, shared coding bases >= 0.5 of the smaller merged CDS footprint), then
reports GC3, ENC and CAI on both sides of the same pairs, and the per-pair
difference. A class whose paired reference genes carry the same weak bias is a
property of those genes; a drop from reference to catalogue within the pair is
something the revision did.

Usage: paired_codon.py <genome.fa> <catalogue.gff3> <reference.gff3> <classes.tsv>
"""
import os
import statistics as stat
import sys
from collections import defaultdict

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/revision'))
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
from codon_usage import fasta, cds_by_gene, codons, gc3, enc, cai, cai_weights  # noqa: E402
from gm_compare import parse_gff, index, overlaps  # noqa: E402


def main():
    genome_fa, cat_gff, ref_gff, classes_tsv = sys.argv[1:5]
    genome = dict(fasta(genome_fa))
    cat_seq = cds_by_gene(cat_gff, genome)
    ref_seq = cds_by_gene(ref_gff, genome)
    w = cai_weights(list(ref_seq.values()))

    cls = {}
    for line in open(classes_tsv):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2:
            cls[f[0]] = f[1]

    cat_g, _ = parse_gff(cat_gff)
    ref_g, _ = parse_gff(ref_gff, want_protein_coding=True)
    ridx = index(ref_g)

    pairs = defaultdict(list)          # class -> [(cat gene, ref gene)]
    for g, v in cat_g.items():
        c = cls.get(g)
        if c not in ('identical', 'revised'):
            continue
        best, bov = None, 0.0
        for s, e, rid in ridx.get((v['chrom'], v['strand']), []):
            if e < v['start']:
                continue
            if s > v['end']:
                break
            ov = overlaps(v, ref_g[rid])
            if ov > bov:
                best, bov = rid, ov
        if best and bov >= 0.5:
            pairs[c].append((g, best))

    def metrics(seq):
        cs = codons(seq)
        if len(cs) < 20:
            return None
        v = (gc3(cs), enc(cs), cai(cs, w))
        return None if any(x != x for x in v) else v

    print('class\tpairs\tcatalogue GC3/ENC/CAI\tpaired reference GC3/ENC/CAI\t'
          'median per-pair difference (catalogue - reference)')
    for c in ('identical', 'revised'):
        rows = []
        for g, rid in pairs[c]:
            a = metrics(cat_seq.get(g, ''))
            b = metrics(ref_seq.get(rid, ''))
            if a and b:
                rows.append((a, b))
        if not rows:
            continue
        med = lambda xs: stat.median(xs)
        A = [med([r[0][i] for r in rows]) for i in range(3)]
        B = [med([r[1][i] for r in rows]) for i in range(3)]
        D = [med([r[0][i] - r[1][i] for r in rows]) for i in range(3)]
        print(f'{c}\t{len(rows):,}\t{A[0]:.3f} / {A[1]:.1f} / {A[2]:.3f}\t'
              f'{B[0]:.3f} / {B[1]:.1f} / {B[2]:.3f}\t'
              f'{D[0]:+.3f} / {D[1]:+.1f} / {D[2]:+.3f}')
        same = sum(1 for r in rows if abs(r[0][0] - r[1][0]) < 1e-9)
        print(f'  pairs with an identical CDS sequence on both sides: {same:,}')


if __name__ == '__main__':
    main()
