#!/usr/bin/env python3
"""wild_gene_overlap.py — do the rice loci land on genes that the wild relative's
own annotation already records?

The dN/dS test (Supplementary Table 5) aligns each protein to the wild genome and
never consults that genome's annotation, so it cannot say whether the aligned
region is an annotated gene there. This script answers that: it takes the
miniprot alignment of each sampled protein and asks whether the aligned interval
falls inside a protein-coding gene of the wild relative's annotation.

Usage: wild_gene_overlap.py <miniprot.gff> <wild.gff3> <out.tsv> [min_frac=0.5]
  miniprot.gff : miniprot --gff output for one query class
  wild.gff3    : the wild relative's Ensembl annotation
Writes one row per query: query, target seq, aligned span, overlapped gene, the
fraction of the aligned span inside it.
"""
import gzip
import re
import sys
from bisect import bisect_right
from collections import defaultdict


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def wild_genes(path):
    """chrom -> sorted [(start, end, gene_id)] for protein-coding genes."""
    genes = defaultdict(list)
    for line in opener(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene':
            continue
        if 'biotype=protein_coding' not in f[8] and 'gene_biotype=protein_coding' not in f[8]:
            continue
        m = re.search(r'ID=(?:gene:)?([^;]+)', f[8])
        genes[f[0]].append((int(f[3]), int(f[4]), m.group(1) if m else '.'))
    for c in genes:
        genes[c].sort()
    return genes


def best_alignments(path):
    """query -> (chrom, start, end) of its highest-scoring miniprot mRNA."""
    best = {}
    for line in opener(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'mRNA':
            continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        q = a.get('Target', '').split()[0]
        score = float(f[5]) if f[5] not in ('.', '') else 0.0
        if q and (q not in best or score > best[q][3]):
            best[q] = (f[0], int(f[3]), int(f[4]), score)
    return {q: v[:3] for q, v in best.items()}


def overlap_gene(genes, chrom, s, e):
    """The annotated gene sharing most bases with [s, e], and that share."""
    lst = genes.get(chrom) or []
    if not lst:
        return None, 0.0
    i = bisect_right([g[0] for g in lst], e)
    best, frac = None, 0.0
    for gs, ge, gid in reversed(lst[:i]):
        if ge < s:
            if gs < s - 2_000_000:
                break
            continue
        ov = min(e, ge) - max(s, gs) + 1
        f = ov / (e - s + 1)
        if f > frac:
            best, frac = gid, f
    return best, frac


def main():
    mp_gff, wild_gff, out = sys.argv[1:4]
    minf = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    genes = wild_genes(wild_gff)
    aln = best_alignments(mp_gff)
    hit = 0
    with open(out, 'w') as o:
        o.write('query\ttarget\tstart\tend\twild_gene\toverlap_frac\n')
        for q, (c, s, e) in sorted(aln.items()):
            g, f = overlap_gene(genes, c, s, e)
            hit += f >= minf
            o.write(f'{q}\t{c}\t{s}\t{e}\t{g or "."}\t{f:.3f}\n')
    n = len(aln)
    print(f'{out}: aligned {n}, on an annotated gene (>= {minf:g} of the aligned span) '
          f'{hit} ({100 * hit / n:.1f}%)' if n else f'{out}: no alignments')


if __name__ == '__main__':
    main()
