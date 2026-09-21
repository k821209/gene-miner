#!/usr/bin/env python3
"""revised_features.py — how each revised locus differs from the reference gene it
matches, counted two ways so the unit is explicit (Reviewer 1, point 9):

  representative : coding exons and CDS length are taken from the representative
                   (longest-CDS) transcript of each locus; isoform count is the
                   number of transcripts at the locus
  aggregate      : coding exons are the distinct CDS intervals over all isoforms
                   and CDS length is the merged coding footprint

Prints the counts of revised loci that add isoforms, add coding exons and extend
the CDS, plus the seven-way overlap used by the Figure 2 Venn diagram.

Usage: revised_features.py union.final.gff3 reference.gff3 [out.tsv]
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/analysis'))
from gm_compare import parse_gff, index  # noqa: E402
from gm_overlap import ofrac, merge, flen  # noqa: E402
import re  # noqa: E402


def loci_with_transcripts(path, protein_coding=False):
    """gene -> {tx: [(s,e)...]}, plus chrom/strand and reference biotype filter"""
    t2g, cds, loc, biotype = {}, defaultdict(list), {}, {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            gid = re.search(r'ID=([^;]+)', f[8])
            bt = re.search(r'biotype=([^;]+)', f[8])
            if gid:
                biotype[gid.group(1)] = bt.group(1) if bt else ''
        elif f[2] in ('mRNA', 'transcript'):
            tid = re.search(r'ID=([^;]+)', f[8]); par = re.search(r'Parent=([^;]+)', f[8])
            if tid:
                t2g[tid.group(1)] = par.group(1) if par else tid.group(1)
        elif f[2] == 'CDS':
            par = re.search(r'Parent=([^;]+)', f[8])
            if par:
                t = par.group(1).split(',')[0]
                cds[t].append((int(f[3]), int(f[4])))
                loc[t] = (f[0], f[6])
    genes = defaultdict(dict)
    info = {}
    for t, iv in cds.items():
        g = t2g.get(t, t)
        if protein_coding and biotype.get(g, '') not in ('', 'protein_coding'):
            continue
        genes[g][t] = sorted(iv)
        info[g] = loc[t]
    return genes, info


def features(txd):
    """(n_isoforms, rep_exons, rep_cds_len, agg_exons, agg_cds_len)"""
    rep = max(txd.values(), key=lambda iv: sum(e - s + 1 for s, e in iv))
    allc = sorted({iv for t in txd.values() for iv in t})
    m = merge([iv for t in txd.values() for iv in t])
    return (len(txd), len(rep), sum(e - s + 1 for s, e in rep), len(allc), flen(m))


def main():
    gmp, refp = sys.argv[1:3]
    out = sys.argv[3] if len(sys.argv) > 3 else None
    gm_tx, _ = loci_with_transcripts(gmp)
    ref_tx, _ = loci_with_transcripts(refp, protein_coding=True)
    gm, _ = parse_gff(gmp)
    ref, _ = parse_gff(refp, want_protein_coding=True)
    ridx = index(ref)

    counts = {k: defaultdict(int) for k in ('representative', 'aggregate')}
    venn = {k: defaultdict(int) for k in ('representative', 'aggregate')}
    n_revised = 0
    rows = []
    for gid, gv in gm.items():
        best, bro = None, 0.0
        for rs, re_, rid in ridx.get((gv['chrom'], gv['strand']), []):
            if re_ < gv['start'] or rs > gv['end']:
                continue
            r = ofrac(gv, ref[rid])
            if r > bro:
                bro, best = r, rid
        if not best or bro < 0.5 or gv['cds'] == ref[best]['cds']:
            continue                      # novel, unmatched, or identical
        n_revised += 1
        g_iso, g_rex, g_rlen, g_aex, g_alen = features(gm_tx[gid])
        r_iso, r_rex, r_rlen, r_aex, r_alen = features(ref_tx[best])
        for mode, (ge, gl, re_, rl) in (('representative', (g_rex, g_rlen, r_rex, r_rlen)),
                                        ('aggregate', (g_aex, g_alen, r_aex, r_alen))):
            iso, ex, cl = g_iso > r_iso, ge > re_, gl > rl
            c = counts[mode]
            c['adds_isoforms'] += iso; c['adds_exons'] += ex; c['longer_cds'] += cl
            c['all_three'] += iso and ex and cl
            c['none_of_three'] += not (iso or ex or cl)
            venn[mode][(iso, ex, cl)] += 1
        rows.append((gid, best, g_iso, r_iso, g_rex, r_rex, g_rlen, r_rlen))

    print(f'revised_loci\t{n_revised}')
    for mode in ('representative', 'aggregate'):
        c = counts[mode]
        print(f'-- {mode}')
        for k in ('adds_isoforms', 'adds_exons', 'longer_cds', 'all_three', 'none_of_three'):
            print(f'   {k}\t{c[k]}\t{100*c[k]/max(n_revised,1):.1f}%')
        print('   venn(iso,exon,cds)=' + ', '.join(
            f'{k}:{v}' for k, v in sorted(venn[mode].items(), reverse=True)))
    if out:
        with open(out, 'w') as o:
            o.write('locus\treference\tgm_isoforms\tref_isoforms\tgm_rep_exons\tref_rep_exons\t'
                    'gm_rep_cds\tref_rep_cds\n')
            for r in rows:
                o.write('\t'.join(map(str, r)) + '\n')


if __name__ == '__main__':
    main()
