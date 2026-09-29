#!/usr/bin/env python3
"""codon_usage.py — codon usage of the rice catalogue by comparison class and
confidence tier, with an intergenic-ORF null.

Reviewer 2 offered codon usage bias as an alternative to Ka/Ks. It needs no
alignment to a relative, so unlike the dN/dS test it also covers the loci that do
not align, and it asks a different question: does the sequence read like coding
DNA of this species at all? Three measures per representative CDS:

  GC3   GC at third codon positions (rice coding sequence is GC3-rich; random
        ORFs sit near the genomic background)
  ENC   Wright's effective number of codons, 20 (extreme bias) to 61 (none)
  CAI   codon adaptation index against the codon usage of the reference
        annotation's own CDS, so 1 is "uses the codons rice genes use"

The null is ORFs of the same length distribution drawn from intergenic sequence:
whatever these three measures give there is what "no coding constraint" looks
like on this genome.

Usage: codon_usage.py <genome.fa> <catalogue.tiers.gff3> <classes.tsv>
                      <reference.gff3> <out.tsv> [n_null=3000] [n_replicates=10]
"""
import random
import re
import sys
from collections import Counter, defaultdict
from math import isnan

STOPS = {'TAA', 'TAG', 'TGA'}
BASES = 'TCAG'
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AA = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
CODON2AA = dict(zip(CODONS, AA))
SYN = defaultdict(list)
for c, a in CODON2AA.items():
    if a != '*':
        SYN[a].append(c)


def fasta(path):
    name, seq = None, []
    for line in open(path):
        if line[0] == '>':
            if name:
                yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, ''.join(seq)


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


def cds_by_gene(gff, genome):
    """representative (longest) CDS per gene, spliced and oriented."""
    t2g, parts = {}, defaultdict(list)
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or line[0] == '#':
            continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        if f[2] in ('mRNA', 'transcript'):
            t2g[a.get('ID', '')] = a.get('Parent', '')
        elif f[2] == 'CDS':
            for p in a.get('Parent', '').split(','):
                parts[p].append((f[0], int(f[3]), int(f[4]), f[6]))
    best = {}
    for t, segs in parts.items():
        g = t2g.get(t)
        if not g:
            continue
        segs.sort(key=lambda x: x[1])
        seq = ''.join(genome.get(c, '')[s - 1:e] for c, s, e, _ in segs)
        if segs[0][3] == '-':
            seq = revcomp(seq)
        if len(seq) > len(best.get(g, ('', ''))[1]):
            best[g] = (t, seq.upper())
    return {g: s for g, (_, s) in best.items()}


def codons(seq):
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]


def gc3(cs):
    third = [c[2] for c in cs if c in CODON2AA and CODON2AA[c] != '*']
    return sum(b in 'GC' for b in third) / len(third) if third else float('nan')


def enc(cs):
    """Wright's Nc, computed over the amino acids present."""
    counts = defaultdict(Counter)
    for c in cs:
        a = CODON2AA.get(c)
        if a and a != '*':
            counts[a][c] += 1
    by_family = defaultdict(list)
    for a, cnt in counts.items():
        n = sum(cnt.values())
        k = len(SYN[a])
        if n <= 1 or k == 1:
            continue
        f = (n * sum((v / n) ** 2 for v in cnt.values()) - 1) / (n - 1)
        if f > 0:
            by_family[k].append(f)
    if not by_family:
        return float('nan')
    fbar = {k: sum(v) / len(v) for k, v in by_family.items()}
    total = 2.0                                   # the two single-codon families
    for k, n_fam in ((2, 9), (3, 1), (4, 5), (6, 3)):
        if k in fbar and fbar[k] > 0:
            total += n_fam / fbar[k]
        else:
            total += n_fam                        # no data: assume no bias
    return min(total, 61.0)


def cai_weights(seqs):
    counts = Counter()
    for s in seqs:
        counts.update(c for c in codons(s) if c in CODON2AA and CODON2AA[c] != '*')
    w = {}
    for a, cods in SYN.items():
        mx = max(counts[c] for c in cods) or 1
        for c in cods:
            w[c] = max(counts[c] / mx, 0.005)
    return w


def cai(cs, w):
    vals = [w[c] for c in cs if c in w and len(SYN[CODON2AA[c]]) > 1]
    if not vals:
        return float('nan')
    from math import exp, log
    return exp(sum(log(v) for v in vals) / len(vals))


def intergenic_orfs(genome, ref_gff, lengths, n, seed=0):
    """ORFs of the given lengths drawn from sequence no reference gene covers."""
    covered = defaultdict(list)
    for line in open(ref_gff):
        f = line.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'gene':
            covered[f[0]].append((int(f[3]), int(f[4])))
    rng = random.Random(seed)
    free = {}
    for c, seq in genome.items():
        iv = sorted(covered.get(c, []))
        pos, gaps = 1, []
        for s, e in iv:
            if s - pos > 5000:
                gaps.append((pos, s - 1))
            pos = max(pos, e + 1)
        if len(seq) - pos > 5000:
            gaps.append((pos, len(seq)))
        if gaps:
            free[c] = gaps
    out, chroms = [], list(free)
    while len(out) < n and chroms:
        c = rng.choice(chroms)
        s, e = rng.choice(free[c])
        L = rng.choice(lengths)
        if e - s <= L:
            continue
        p = rng.randrange(s, e - L)
        seq = genome[c][p:p + L].upper()
        if 'N' in seq or len(seq) < 60:
            continue
        if rng.random() < 0.5:
            seq = revcomp(seq)
        out.append(seq)
    return out


def medians(seqs, w):
    rows = []
    for s in seqs:
        cs = codons(s)
        if len(cs) < 20:
            continue
        r = (gc3(cs), enc(cs), cai(cs, w))
        if not any(isnan(x) for x in r):
            rows.append(r)
    if not rows:
        return None
    med = lambda i: sorted(r[i] for r in rows)[len(rows) // 2]
    return len(rows), med(0), med(1), med(2)


def summarize(name, seqs, w, o):
    m = medians(seqs, w)
    if not m:
        return
    n, g, e, c = m
    o.write(f'{name}\t{n}\t{g:.3f}\t{e:.1f}\t{c:.3f}\n')
    print(f'{name:<28} n={n:>6}  GC3={g:.3f}  ENC={e:.1f}  CAI={c:.3f}')


def summarize_replicates(name, draws, w, o):
    """One row per measure: median over replicate draws, and their spread."""
    per = [medians(d, w) for d in draws]
    per = [p for p in per if p]
    if not per:
        return
    import statistics as stat
    n = int(stat.mean(p[0] for p in per))
    cells = []
    for i, fmt in ((1, '.3f'), (2, '.1f'), (3, '.3f')):
        vals = [p[i] for p in per]
        mean = stat.mean(vals)
        sd = stat.stdev(vals) if len(vals) > 1 else 0.0
        cells.append(f'{mean:{fmt}} ± {sd:{fmt}}')
    o.write(f'{name} ({len(per)} draws)\t{n}\t' + '\t'.join(cells) + '\n')
    print(f'{name:<28} n={n:>6}  GC3={cells[0]}  ENC={cells[1]}  CAI={cells[2]}  '
          f'({len(per)} draws)')


def main():
    genome_fa, cat_gff, classes_tsv, ref_gff, out = sys.argv[1:6]
    n_null = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
    n_rep = int(sys.argv[7]) if len(sys.argv) > 7 else 10
    genome = dict(fasta(genome_fa))
    cat = cds_by_gene(cat_gff, genome)
    ref = cds_by_gene(ref_gff, genome)
    w = cai_weights(list(ref.values()))

    cls = {}
    for line in open(classes_tsv):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2:
            cls[f[0]] = f[1]
    tier = {}
    for line in open(cat_gff):
        f = line.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'gene':
            g = re.search(r'ID=([^;]+)', f[8]).group(1)
            t = re.search(r'tier=([^;]+)', f[8])
            tier[g] = t.group(1) if t else 'NA'

    groups = defaultdict(list)
    for g, s in cat.items():
        c = cls.get(g, 'NA')
        groups[c].append(s)
        if c == 'novel':
            groups[f'novel {tier.get(g, "NA")} tier'].append(s)

    with open(out, 'w') as o:
        o.write('set\tn\tGC3_median\tENC_median\tCAI_median\n')
        summarize('reference annotation', list(ref.values()), w, o)
        for k in ('identical', 'revised', 'novel', 'novel high tier',
                  'novel medium tier', 'novel low tier'):
            if groups.get(k):
                summarize(k, groups[k], w, o)
        lengths = [len(s) for s in groups.get('novel', [])] or [300]
        draws = [intergenic_orfs(genome, ref_gff, lengths, n_null, seed=i) for i in range(n_rep)]
        summarize_replicates('intergenic ORF null', draws, w, o)


if __name__ == '__main__':
    main()
