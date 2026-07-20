#!/usr/bin/env python3
"""plot_splice_support.py — genome-browser-style view of RNA-seq splice-alignment
support for a gene model, comparing a Gene-Miner locus with a reference gene.

It draws, over one locus:
  * individual RNA-seq read alignments (a subsample) — reads that cross an intron
    are drawn as blocks joined by a thin line (a *spliced* alignment, green);
    unspliced reads are blue;
  * splice-junction arcs labelled with their spliced-read counts (green where the
    junction confirms a Gene-Miner intron);
  * the Gene-Miner and the reference gene models (CDS exons).

This makes visible whether a *revised* Gene-Miner model's extra exons are backed
by spliced reads (a real improvement) rather than merely differing from the
reference — the per-locus companion to the genome-wide junction analysis.

Dependencies: pysam, numpy, matplotlib (e.g. the `genemark`/any env, or
  conda create -n viz -c conda-forge -c bioconda pysam numpy matplotlib).

The junctions TSV is `chrom<TAB>intron_start<TAB>intron_end<TAB>read_count`
(1-based inclusive intron coordinates). Produce it from the pipeline's BAMs with
regtools, e.g.:
  regtools junctions extract -s XS -o lib.bed sample.bam
then convert each BED12 row to an intron: istart = chromStart+blockSizes[0]+1,
iend = chromEnd-blockSizes[1], summing read_count (BED col 5) across libraries.

Usage:
  plot_splice_support.py --chrom 10 --start 1093271 --end 1096798 --strand - \\
      --gm-gene Os_A001065 --ref-gene gene:Os10g0117000 \\
      --gm-gff union.final.gff3 --ref-gff reference.gff3 \\
      --junctions merged_junctions.tsv --out locus.png \\
      --bam s1.bam --bam s2.bam --bam s3.bam --bam s4.bam
Writes <out> (PNG) and the matching .svg.
"""
import argparse, re
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import pysam

GM_COL, REF_COL, UNSPL_COL, OTHER_COL = '#2e7d32', '#c0392b', '#4a90d9', '#b0b0b0'


def attr(a, k):
    m = re.search(k + r'=([^;]+)', a)
    return m.group(1) if m else None


def tx_cds(path, gene):
    """Return the CDS exon lists (one per transcript) of `gene` in a GFF3."""
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
            tc[attr(f[8], 'Parent')].append((int(f[3]), int(f[4])))
    return [sorted(segs) for tx, segs in tc.items() if txg.get(tx, tx) == gene]


def introns_of(tx_list):
    s = set()
    for segs in tx_list:
        for i in range(len(segs) - 1):
            s.add((segs[i][1] + 1, segs[i + 1][0] - 1))
    return s


def subsample(lst, n):
    if len(lst) <= n:
        return lst
    idx = np.linspace(0, len(lst) - 1, n).astype(int)
    return [lst[i] for i in idx]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--end', type=int, required=True)
    ap.add_argument('--strand', default='+')
    ap.add_argument('--gm-gene', required=True)
    ap.add_argument('--ref-gene', required=True)
    ap.add_argument('--gm-gff', required=True)
    ap.add_argument('--ref-gff', required=True)
    ap.add_argument('--junctions', required=True, help='chrom/istart/iend/count TSV')
    ap.add_argument('--bam', action='append', required=True, help='repeat for each library')
    ap.add_argument('--out', required=True, help='output PNG (a .svg is also written)')
    ap.add_argument('--min-reads', type=int, default=3, help='min spliced reads for a junction')
    ap.add_argument('--gm-label', default='Gene-Miner')
    ap.add_argument('--ref-label', default='Reference')
    ap.add_argument('--max-spliced', type=int, default=32)
    ap.add_argument('--max-unspliced', type=int, default=22)
    a = ap.parse_args()

    pad = int((a.end - a.start) * 0.06)
    L, R = a.start - pad, a.end + pad
    gm_tx = tx_cds(a.gm_gff, a.gm_gene)
    ref_tx = tx_cds(a.ref_gff, a.ref_gene)
    gm_introns = introns_of(gm_tx)

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
            blk = r.get_blocks()          # aligned segments, split at introns (N)
            if blk:
                (spliced if len(blk) > 1 else unspliced).append(blk)
        bam.close()
    reads = subsample(spliced, a.max_spliced) + subsample(unspliced, a.max_unspliced)

    # greedy interval packing into rows
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

    base = 3.2
    fig, ax = plt.subplots(figsize=(11, min(max(3.4 + nrow * 0.077, 4.5), 11)))
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

    def model(tx_list, y, color, label):
        ax.text(L, y + 0.16, label, fontsize=9, fontweight='bold', color=color)
        for k, segs in enumerate(tx_list):
            yy = y - 0.16 * k
            ax.plot([segs[0][0], segs[-1][1]], [yy, yy], color=color, lw=1.0, zorder=3)
            for s, e in segs:
                ax.add_patch(Rectangle((s, yy - 0.075), e - s, 0.15, color=color, zorder=4))

    mtop = jy - 1.15
    model(gm_tx, mtop, GM_COL, a.gm_label)
    model(ref_tx, mtop - 0.7, REF_COL, a.ref_label)
    ax.set_xlim(L, R)
    ax.set_ylim(mtop - 1.2, readtop + 0.3)
    ax.set_yticks([])
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel(f'{a.chrom}:{L:,}-{R:,} ({a.strand})', fontsize=9)
    ax.set_title(f'{a.gm_gene} vs {a.ref_gene}: spliced RNA-seq reads vs the revised gene model',
                 fontsize=9)
    ax.legend(handles=[
        Line2D([0], [0], color=GM_COL, lw=3, label=f'spliced read / {a.gm_label}-supported junction'),
        Line2D([0], [0], color=UNSPL_COL, lw=3, label='unspliced read'),
        Line2D([0], [0], color=OTHER_COL, lw=2, label='other junction')],
        fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.11), frameon=False, ncol=3)
    plt.subplots_adjust(bottom=0.16)
    plt.savefig(a.out, dpi=150)
    plt.savefig(a.out.rsplit('.', 1)[0] + '.svg')
    print(f"wrote {a.out}  rows={nrow} spliced={len(spliced)} unspliced={len(unspliced)} junctions={len(J)}")


if __name__ == '__main__':
    main()
