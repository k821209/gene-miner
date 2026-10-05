#!/usr/bin/env python3
"""plot_boundary_error.py — genome-browser view of one gene-boundary error.

Draws over one window: a subsample of RNA-seq read alignments (spliced reads
green, unspliced blue), the spliced-read junctions as arcs labelled with their
counts, then the Gene-Miner locus or loci and the reference gene or genes. Both
tracks take several genes, so a fusion (two reference genes inside one
Gene-Miner locus) and a split (one reference gene across two Gene-Miner loci)
are drawn by the same code. `--highlight-tx` singles out the one transcript
that spans both reference genes, the cause the Results section names.

Adapted from plot_splice_support.py (Supplementary Figures 1 and 2).

Usage:
  plot_boundary_error.py --chrom 3 --start 639689 --end 645897 --strand - \
      --gm-genes GM_R002604 --ref-genes gene:Os03g0111000,gene:Os03g0111100 \
      --highlight-tx GM_R002604.t4 --gm-gff union.final.gff3 \
      --ref-gff reference.gff3 --junctions rice.junctions.tsv \
      --bam s1.bam --bam s2.bam --out fusion.png
"""
import argparse
import re
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.patches import Rectangle             # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402
import pysam                                         # noqa: E402

GM_COL, REF_COL, UNSPL_COL, OTHER_COL, HL_COL = '#2e7d32', '#c0392b', '#4a90d9', '#b0b0b0', '#6a1b9a'


def attr(a, k):
    m = re.search(k + r'=([^;]+)', a)
    return m.group(1) if m else None


def gene_tx(path):
    """gene -> {transcript: sorted CDS blocks}"""
    txg, tc = {}, defaultdict(list)
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip().split('\t')
        if len(f) < 9:
            continue
        if f[2] in ('mRNA', 'transcript'):
            txg[attr(f[8], 'ID')] = attr(f[8], 'Parent')
        elif f[2] == 'CDS':
            tc[attr(f[8], 'Parent').split(',')[0]].append((int(f[3]), int(f[4])))
    out = defaultdict(dict)
    for tx, segs in tc.items():
        out[txg.get(tx, tx)][tx] = sorted(segs)
    return out


def introns_of(tx_map):
    s = set()
    for segs in tx_map.values():
        for i in range(len(segs) - 1):
            s.add((segs[i][1] + 1, segs[i + 1][0] - 1))
    return s


def subsample(lst, n):
    if len(lst) <= n:
        return lst
    return [lst[i] for i in np.linspace(0, len(lst) - 1, n).astype(int)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--end', type=int, required=True)
    ap.add_argument('--strand', default='+')
    ap.add_argument('--gm-genes', required=True)
    ap.add_argument('--ref-genes', required=True)
    ap.add_argument('--highlight-tx', default=None)
    ap.add_argument('--gm-gff', required=True)
    ap.add_argument('--ref-gff', required=True)
    ap.add_argument('--junctions', required=True)
    ap.add_argument('--bam', action='append', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-reads', type=int, default=3)
    ap.add_argument('--title', default='')
    ap.add_argument('--max-spliced', type=int, default=30)
    ap.add_argument('--max-unspliced', type=int, default=18)
    a = ap.parse_args()

    pad = int((a.end - a.start) * 0.05)
    L, R = a.start - pad, a.end + pad
    gm_all, ref_all = gene_tx(a.gm_gff), gene_tx(a.ref_gff)
    gm = {g: gm_all[g] for g in a.gm_genes.split(',')}
    ref = {g: ref_all[g] for g in a.ref_genes.split(',')}
    gm_introns = set()
    for m in gm.values():
        gm_introns |= introns_of(m)

    J = []
    for line in open(a.junctions):
        f = line.rstrip().split('\t')
        if len(f) >= 4 and f[0] == a.chrom and int(f[3]) >= a.min_reads \
                and int(f[1]) >= L and int(f[2]) <= R:
            J.append((int(f[1]), int(f[2]), int(f[3])))

    spliced, unspliced = [], []
    for b in a.bam:
        bam = pysam.AlignmentFile(b, 'rb')
        for r in bam.fetch(a.chrom, L, R):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            blk = r.get_blocks()
            if blk:
                (spliced if len(blk) > 1 else unspliced).append(blk)
        bam.close()
    reads = subsample(spliced, a.max_spliced) + subsample(unspliced, a.max_unspliced)

    rows, placed = [], []
    for blk in reads:
        s, e = blk[0][0], blk[-1][1]
        for ri, occ in enumerate(rows):
            if s > occ + 120:
                rows[ri] = e
                placed.append((ri, blk))
                break
        else:
            rows.append(e)
            placed.append((len(rows) - 1, blk))
    nrow = len(rows)

    n_models = sum(len(m) for m in gm.values()) + sum(len(m) for m in ref.values())
    n_labels = len(gm) + len(ref)
    fig, ax = plt.subplots(figsize=(11, min(max(4.2 + nrow * 0.08 + (n_models + n_labels) * 0.22, 5.5), 12)))
    base = 3.2
    for ri, blk in placed:
        y = base + ri * 0.11
        sp = len(blk) > 1
        ax.plot([blk[0][0], blk[-1][1]], [y, y],
                color=(GM_COL if sp else '#9bbcd8'), lw=0.4, zorder=2)
        for bs, be in blk:
            ax.add_patch(Rectangle((bs, y - 0.045), be - bs, 0.09,
                                   color=(GM_COL if sp else UNSPL_COL), lw=0, zorder=3))
    readtop = base + nrow * 0.11 + 0.15
    ax.text(L, readtop, 'RNA-seq read alignments (spliced reads green; thin line = intron gap)',
            fontsize=8, color='#2b6cb0')

    jmax = max((sc for _, _, sc in J), default=1)
    jy = base - 0.15
    for s, e, sc in J:
        col = '#1b5e20' if (s, e) in gm_introns else OTHER_COL
        mid, w = (s + e) / 2, e - s
        th = np.linspace(0, np.pi, 60)
        ax.plot(mid + (w / 2) * np.cos(th), jy - 0.9 * np.sin(th) * (0.4 + 0.6 * sc / jmax),
                color=col, lw=0.6 + 2.0 * sc / jmax, alpha=0.85, zorder=2)
        ax.text(mid, jy - 0.9 * (0.4 + 0.6 * sc / jmax) - 0.12, str(sc),
                ha='center', va='top', fontsize=6, color=col)

    y = jy - 1.3

    def draw(gene, tx_map, color, label):
        nonlocal y
        ax.text(L, y, label, fontsize=9, fontweight='bold', color=color, va='center')
        y -= 0.3
        for tx, segs in sorted(tx_map.items()):
            hl = (tx == a.highlight_tx)
            col = HL_COL if hl else color
            ax.plot([segs[0][0], segs[-1][1]], [y, y], color=col, lw=1.2 if hl else 1.0, zorder=3)
            for s, e in segs:
                ax.add_patch(Rectangle((s, y - (0.1 if hl else 0.075)), e - s,
                                       (0.2 if hl else 0.15), color=col, zorder=4))
            if hl:
                ax.text(segs[0][0] - (R - L) * 0.004, y,
                        'spans both reference genes  ', fontsize=7.5, color=HL_COL,
                        va='center', ha='right')
            y -= 0.3
        y -= 0.25

    for gene, tx_map in gm.items():
        draw(gene, tx_map, GM_COL, gene + ' (Gene-Miner)')
    for gene, tx_map in ref.items():
        draw(gene, tx_map, REF_COL, gene.replace('gene:', '') + ' (reference)')

    ax.set_xlim(L, R)
    ax.ticklabel_format(axis='x', style='plain', useOffset=False)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.set_ylim(y - 0.2, readtop + 0.3)
    ax.set_yticks([])
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel(f'{a.chrom}:{L:,}-{R:,} ({a.strand})', fontsize=9)
    if a.title:
        ax.set_title(a.title, fontsize=9)
    handles = [Line2D([0], [0], color=GM_COL, lw=3, label='spliced read / Gene-Miner junction'),
               Line2D([0], [0], color=UNSPL_COL, lw=3, label='unspliced read'),
               Line2D([0], [0], color=OTHER_COL, lw=2, label='other junction'),
               Line2D([0], [0], color=REF_COL, lw=3, label='reference gene')]
    if a.highlight_tx:
        handles.append(Line2D([0], [0], color=HL_COL, lw=3, label='spanning transcript'))
    ax.legend(handles=handles, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.11),
              frameon=False, ncol=len(handles))
    plt.subplots_adjust(bottom=0.16)
    plt.savefig(a.out, dpi=150, bbox_inches='tight')
    plt.savefig(a.out.rsplit('.', 1)[0] + '.svg', bbox_inches='tight')
    print(f'wrote {a.out} rows={nrow} spliced={len(spliced)} unspliced={len(unspliced)} junctions={len(J)}')


if __name__ == '__main__':
    main()
