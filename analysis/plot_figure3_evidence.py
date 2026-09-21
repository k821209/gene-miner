#!/usr/bin/env python3
"""plot_figure3_evidence.py — Figure 3: independent evidence behind each catalogue,
split by comparison with the reference (identical / revised / novel), for the two
crops. Everything is computed here from the released files:

  union.final.tiers.gff3   evidence flags written by annotate_tiers.py
  new.annotations          eggNOG-mapper output (functional annotation bar)
  reference GFF3           to label each locus identical / revised / novel

Usage: plot_figure3_evidence.py <out prefix> <rice dir> <rice ref> <soy dir> <soy ref>
"""
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/analysis'))
from gm_compare import parse_gff, index  # noqa: E402
from gm_overlap import ofrac  # noqa: E402

CATS = ('identical', 'revised', 'novel')
BARS = (('rna', 'RNA-seq transcript'), ('sprot', 'Swiss-Prot hit'),
        ('eggnog', 'eggNOG ortholog'), ('func', 'functional annotation'))
COLS = ('#4e9a5f', '#5b9bd5', '#e8983a', '#b07cc6')


def label(gff, ref):
    gm, _ = parse_gff(gff)
    rf, _ = parse_gff(ref, want_protein_coding=True)
    ridx = index(rf)
    out = {}
    for gid, gv in gm.items():
        best, bro = None, 0.0
        for rs, re_, rid in ridx.get((gv['chrom'], gv['strand']), []):
            if re_ < gv['start'] or rs > gv['end']:
                continue
            r = ofrac(gv, rf[rid])
            if r > bro:
                bro, best = r, rid
        if best and bro >= 0.5:
            out[gid] = 'identical' if gv['cds'] == rf[best]['cds'] else 'revised'
        else:
            out[gid] = 'novel'
    return out


def evidence(tiers_gff, ann):
    ev, t2g = {}, {}
    for line in open(tiers_gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            g = re.search(r'ID=([^;]+)', f[8]).group(1)
            e = re.search(r'evidence=([^;]+)', f[8])
            ev[g] = set() if not e or e.group(1) == 'none' else set(e.group(1).split(','))
        elif f[2] == 'mRNA':
            t2g[re.search(r'ID=([^;]+)', f[8]).group(1)] = re.search(r'Parent=([^;]+)', f[8]).group(1)
    H = []
    for line in open(ann):
        if line.startswith('#'):
            if line.startswith('#query'):
                H = line.lstrip('#').rstrip('\n').split('\t')
            continue
        r = dict(zip(H, line.rstrip('\n').split('\t')))
        g = t2g.get(r.get('query', ''))
        if g in ev and any((r.get(k) or '-') not in ('-', '') for k in ('COG_category', 'Description', 'GOs')):
            ev[g].add('func')
    return ev


def fractions(d, ref):
    lab = label(os.path.join(d, 'union.final.gff3'), ref)
    ev = evidence(os.path.join(d, 'union.final.tiers.gff3'), os.path.join(d, 'new.annotations'))
    n = defaultdict(int)
    hit = defaultdict(lambda: defaultdict(int))
    for g, c in lab.items():
        n[c] += 1
        for k, _ in BARS:
            hit[c][k] += k in ev.get(g, set())
    return n, {c: {k: 100 * hit[c][k] / max(n[c], 1) for k, _ in BARS} for c in CATS}


def main():
    out, rdir, rref, sdir, sref = sys.argv[1:6]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, (tag, d, ref, crop) in zip(axes, (('a', rdir, rref, 'rice'), ('b', sdir, sref, 'soybean'))):
        n, fr = fractions(d, ref)
        x = np.arange(len(CATS))
        w = 0.2
        for i, ((k, name), col) in enumerate(zip(BARS, COLS)):
            vals = [fr[c][k] for c in CATS]
            bars = ax.bar(x + (i - 1.5) * w, vals, w, color=col, label=name)
            for b, v in zip(bars, vals):
                ax.text(b.get_x() + b.get_width() / 2, v + 1, f'{v:.0f}', ha='center', fontsize=7.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f'{c}\n(n = {n[c]:,})' for c in CATS], fontsize=9.5)
        ax.set_ylim(0, 108)
        ax.set_title(tag, loc='left', fontweight='bold', fontsize=13)
        ax.text(0.99, 0.98, crop, transform=ax.transAxes, ha='right', va='top', fontsize=11)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        print(crop, dict(n), {c: {k: round(v, 1) for k, v in fr[c].items()} for c in CATS})
    axes[0].set_ylabel('Loci with the evidence (%)', fontsize=10)
    axes[0].legend(fontsize=8.5, loc='upper right', bbox_to_anchor=(1.0, 0.9), frameon=True)
    plt.tight_layout()
    plt.savefig(out + '.png', dpi=160)
    plt.savefig(out + '.svg')


if __name__ == '__main__':
    main()
