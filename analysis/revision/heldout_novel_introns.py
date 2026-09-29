#!/usr/bin/env python3
"""heldout_novel_introns.py — of the coding introns a catalogue proposes that the
reference annotation does not contain, how many does an unseen tissue confirm?

The held-out test (Supplementary Table 6) splits introns by whether the libraries
the catalogue was built from support them. This splits them the other way: by
whether the reference annotation already has the intron. An intron the reference
lacks, confirmed by a library the catalogue never saw, is structure the reference
is missing and that independent data backs.

Usage: heldout_novel_introns.py <catalogue.gff3> <reference.gff3>
                                <heldout.junctions.tsv> <minreads>
                                <input.junctions.tsv> [<input.junctions.tsv> ...]
"""
import sys
from collections import defaultdict


def junctions(paths, minreads):
    c = defaultdict(int)
    for p in paths:
        for line in open(p):
            ch, s, e, n = line.split('\t')
            c[(ch, int(s), int(e))] += int(n)
    return {k for k, n in c.items() if n >= minreads}


def coding_introns(gff, want_protein_coding=False):
    """set of (chrom, first intron base, last intron base) over all transcripts."""
    import re
    keep = None
    if want_protein_coding:
        keep = set()
        for line in open(gff):
            f = line.rstrip('\n').split('\t')
            if len(f) > 8 and f[2] == 'mRNA' and (
                    'protein_coding' in f[8] or 'biotype=mRNA' in f[8] or True):
                m = re.search(r'ID=([^;]+)', f[8])
                if m:
                    keep.add(m.group(1))
    cds = defaultdict(list)
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        for p in a.get('Parent', '').split(','):
            if p and (keep is None or p in keep):
                cds[p].append((f[0], int(f[3]), int(f[4])))
    out = set()
    for segs in cds.values():
        segs.sort(key=lambda x: x[1])
        for (c, _, e1), (_, s2, _) in zip(segs, segs[1:]):
            if s2 - e1 > 1:
                out.add((c, e1 + 1, s2 - 1))
    return out


def pct(a, b):
    return f'{a:,} of {b:,} ({100 * a / b:.1f}%)' if b else '0'


def main():
    cat_gff, ref_gff, held, minreads = sys.argv[1:5]
    minreads = int(minreads)
    inp = junctions(sys.argv[5:], minreads)
    hout = junctions([held], minreads)
    cat = coding_introns(cat_gff)
    ref = coding_introns(ref_gff, want_protein_coding=True)

    added = cat - ref          # introns the reference does not contain
    shared = cat & ref
    print('set\tintrons\tsupported by the input libraries\tconfirmed by the held-out libraries\t'
          'confirmed among the input-supported\tconfirmed among the rest')
    for name, s in (('catalogue introns absent from the reference', added),
                    ('catalogue introns also in the reference', shared),
                    ('reference introns absent from the catalogue', ref - cat)):
        sup = s & inp
        rest = s - inp
        print(f'{name}\t{len(s):,}\t{pct(len(sup), len(s))}\t{pct(len(s & hout), len(s))}\t'
              f'{pct(len(sup & hout), len(sup))}\t{pct(len(rest & hout), len(rest))}')


if __name__ == '__main__':
    main()
