#!/usr/bin/env python3
"""Per-pair spread of the revised class' codon metrics against the reference gene
each locus revises: quartiles of the difference, and how often the CDS is the same
sequence on both sides."""
import os, statistics as stat, sys
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/revision'))
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
from codon_usage import fasta, cds_by_gene, codons, gc3, enc, cai, cai_weights
from gm_compare import parse_gff, index, overlaps

genome = dict(fasta(os.path.expanduser('~/gene-miner-runs/rice/genome/IRGSP.fa')))
R = os.path.expanduser('~/gene-miner-runs/revision/')
cat_seq = cds_by_gene(R + 'regen/rice_s/union.final.tiers.gff3', genome)
ref_seq = cds_by_gene(os.path.expanduser('~/gene-miner-runs/rice/genome/IRGSP.ref.gff3'), genome)
w = cai_weights(list(ref_seq.values()))
cls = dict(l.rstrip('\n').split('\t')[:2] for l in open(R + 'contamination/seqs_rice/classes.tsv') if '\t' in l)
cat_g, _ = parse_gff(R + 'regen/rice_s/union.final.tiers.gff3')
ref_g, _ = parse_gff(os.path.expanduser('~/gene-miner-runs/rice/genome/IRGSP.ref.gff3'), want_protein_coding=True)
ridx = index(ref_g)

def met(seq):
    cs = codons(seq)
    if len(cs) < 20:
        return None
    v = (gc3(cs), enc(cs), cai(cs, w))
    return None if any(x != x for x in v) else v

d, same_seq, longer, shorter = [], 0, 0, 0
for g, v in cat_g.items():
    if cls.get(g) != 'revised':
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
    if not best or bov < 0.5:
        continue
    a, b = met(cat_seq.get(g, '')), met(ref_seq.get(best, ''))
    if not (a and b):
        continue
    if cat_seq[g] == ref_seq[best]:
        same_seq += 1
    elif len(cat_seq[g]) > len(ref_seq[best]):
        longer += 1
    else:
        shorter += 1
    d.append(tuple(x - y for x, y in zip(a, b)))

q = lambda i, p: stat.quantiles([x[i] for x in d], n=100)[p - 1]
print(f'pairs {len(d):,}; representative CDS the same sequence on both sides {same_seq:,} '
      f'({100*same_seq/len(d):.1f}%), catalogue longer {longer:,}, catalogue shorter {shorter:,}')
for i, name in ((0, 'GC3'), (1, 'ENC'), (2, 'CAI')):
    vals = [x[i] for x in d]
    print(f'{name}: median {stat.median(vals):+.3f}  mean {stat.mean(vals):+.3f}  '
          f'IQR {q(i,25):+.3f} to {q(i,75):+.3f}  '
          f'fraction above the reference {100*sum(v>0 for v in vals)/len(vals):.1f}%')
