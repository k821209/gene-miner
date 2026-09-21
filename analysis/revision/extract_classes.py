#!/usr/bin/env python3
"""extract_classes.py — label each Gene-Miner locus against the reference with the
same rule as gm_compare.py (same-strand CDS-footprint overlap >= 0.5 of the
shorter locus) and write the primary transcript's CDS and protein per class.

Usage: extract_classes.py <genome.fa> <union.final.gff3> <reference.gff3> <outdir>
                          [sample_per_class] [extra.gff3=TE-removed models]
Writes <outdir>/<class>.cds.fa, <class>.pep.fa and classes.tsv
"""
import os
import random
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
from gm_compare import parse_gff, index, overlaps  # noqa: E402

CODON = {}
for i, b1 in enumerate('TCAG'):
    for j, b2 in enumerate('TCAG'):
        for k, b3 in enumerate('TCAG'):
            CODON[b1 + b2 + b3] = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[i * 16 + j * 4 + k]


def fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, ''.join(seq)


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


def translate(nt):
    return ''.join(CODON.get(nt[i:i + 3].upper(), 'X') for i in range(0, len(nt) - 2, 3))


def primary_transcripts(gff):
    """gene -> (transcript id, sorted CDS blocks, chrom, strand) for the longest CDS"""
    t2g, cds, loc = {}, defaultdict(list), {}
    for line in open(gff):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] in ('mRNA', 'transcript'):
            t2g[re.search(r'ID=([^;]+)', f[8]).group(1)] = re.search(r'Parent=([^;]+)', f[8]).group(1)
        elif f[2] == 'CDS':
            t = re.search(r'Parent=([^;]+)', f[8]).group(1).split(',')[0]
            cds[t].append((int(f[3]), int(f[4])))
            loc[t] = (f[0], f[6])
    best = {}
    for t, iv in cds.items():
        g = t2g.get(t, t)
        ln = sum(e - s + 1 for s, e in iv)
        if g not in best or ln > best[g][0]:
            best[g] = (ln, t, sorted(iv), loc[t][0], loc[t][1])
    return {g: v[1:] for g, v in best.items()}


def main():
    genome_fa, gm_gff, ref_gff, outdir = sys.argv[1:5]
    n_sample = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    extra_gff = sys.argv[6] if len(sys.argv) > 6 else None
    os.makedirs(outdir, exist_ok=True)

    gm, _ = parse_gff(gm_gff)
    ref, _ = parse_gff(ref_gff, want_protein_coding=True)
    ridx = index(ref)
    label = {}
    for gid, gv in gm.items():
        best, best_ro = None, 0.0
        for rs, re_, rid in ridx.get((gv['chrom'], gv['strand']), []):
            if re_ < gv['start'] or rs > gv['end']:
                continue
            ro = overlaps(gv, ref[rid])
            if ro > best_ro:
                best_ro, best = ro, rid
        if best and best_ro >= 0.5:
            label[gid] = 'identical' if gv['cds'] == ref[best]['cds'] else 'revised'
        else:
            label[gid] = 'novel'

    sets = {'gm': (gm_gff, label)}
    if extra_gff:                     # e.g. models dropped by the TE filter
        sets['te'] = (extra_gff, None)

    genome = dict(fasta(genome_fa))
    with open(os.path.join(outdir, 'classes.tsv'), 'w') as ct:
        for gid, lab in label.items():
            ct.write(f'{gid}\t{lab}\n')

    for key, (gff, lab) in sets.items():
        prim = primary_transcripts(gff)
        groups = defaultdict(list)
        for g, (t, iv, chrom, strand) in prim.items():
            cls = 'te_removed' if lab is None else lab.get(g)
            if cls:
                groups[cls].append((g, iv, chrom, strand))
        for cls, items in groups.items():
            if n_sample and len(items) > n_sample:
                random.seed(0)
                items = random.sample(items, n_sample)
            with open(f'{outdir}/{cls}.cds.fa', 'w') as fc, open(f'{outdir}/{cls}.pep.fa', 'w') as fp:
                n = 0
                for g, iv, chrom, strand in items:
                    seq = ''.join(genome[chrom][s - 1:e] for s, e in iv)
                    if strand == '-':
                        seq = revcomp(seq)
                    p = translate(seq).rstrip('*')
                    if len(p) < 30 or '*' in p:
                        continue
                    fc.write(f'>{g}\n{seq}\n')
                    fp.write(f'>{g}\n{p}\n')
                    n += 1
            print(cls, n)


if __name__ == '__main__':
    main()
