#!/usr/bin/env python3
"""collect_sub.py — one row per rice subsampling configuration (Reviewer 1, point 5):
input read pairs, catalogue size, BUSCO (representative and all isoforms), reference
recovery, boundary errors, gffcompare gene/transcript accuracy, and compute cost
from the Nextflow task records (trace.txt, or trace.recovered.txt where the run's
own trace could not be written).

Usage: collect_sub.py <sub dir> <config> [<config> ...]
"""
import glob
import os
import re
import sys

REF_GENES = 35804


def secs(s):
    s = s.strip()
    if s in ('', '-'):
        return 0.0
    t = 0.0
    for v, u in re.findall(r'([\d.]+)\s*(ms|h|m|s|d)', s):
        t += float(v) * {'ms': 1e-3, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}[u]
    return t


def gb(s):
    m = re.match(r'([\d.]+)\s*(KB|MB|GB|TB|B)', s.strip())
    if not m:
        return 0.0
    return float(m.group(1)) * {'B': 1e-9, 'KB': 1e-6, 'MB': 1e-3, 'GB': 1, 'TB': 1e3}[m.group(2)]


def kv(path):
    d = {}
    for line in open(path):
        f = line.strip().split('\t')
        if len(f) >= 2:
            d[f[0].strip()] = f[1]
    return d


def busco(path):
    for p in glob.glob(path):
        m = re.search(r'C:([\d.]+)%', open(p).read())
        if m:
            return float(m.group(1))
    return float('nan')


def gffcmp(path):
    out = {}
    for line in open(path):
        m = re.match(r'\s*(Locus|Transcript) level:\s*([\d.]+)\s*\|\s*([\d.]+)', line)
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return out


def cost(d):
    tr = os.path.join(d, 'trace.recovered.txt')
    if os.path.exists(tr):
        # task_id hash name status exit submit duration realtime pcpu peak_rss ...
        rows = [l.rstrip('\n').split('\t') for l in open(tr) if l.strip()]
        cols = dict(name=2, status=3, submit=5, duration=6, realtime=7, pcpu=8, rss=9)
    else:
        lines = [l.rstrip('\n').split('\t') for l in open(os.path.join(d, 'trace.txt'))]
        h = lines[0]
        rows = lines[1:]
        cols = dict(name=h.index('name'), status=h.index('status'), submit=h.index('submit'),
                    duration=h.index('duration'), realtime=h.index('realtime'),
                    pcpu=h.index('%cpu'), rss=h.index('peak_rss'))
    cpu_h, peak, stages = 0.0, 0.0, {}
    for r in rows:
        if r[cols['status']] != 'COMPLETED':
            continue
        rt = secs(r[cols['realtime']])
        pc = float(r[cols['pcpu']].rstrip('%') or 0) / 100
        cpu_h += rt * pc / 3600
        peak = max(peak, gb(r[cols['rss']]))
        st = r[cols['name']].split(' (')[0]
        stages[st] = stages.get(st, 0.0) + rt * pc / 3600
    return cpu_h, peak, stages


def main():
    sub = sys.argv[1]
    print('config\tread_pairs\tgenes\ttx_per_gene\tbusco_rep\tbusco_alliso\trecovered\trecovered_pct\t'
          'novel\tsplit_pct\tfusion_pct\tlocus_R\tlocus_P\tlocus_F1\ttx_R\ttx_P\ttx_F1\tcpu_h\tpeak_rss_gb\t'
          'cpu_h_genemark')
    for c in sys.argv[2:]:
        d = os.path.join(sub, c)
        reads = int(open(os.path.join(d, 'reads.txt')).read().split()[1])
        g = kv(os.path.join(d, 'gm_compare.tsv'))
        rec = REF_GENES - int(g['reference_loci_missed_by_gm'])
        sf = open(os.path.join(d, 'split_fusion.txt')).read().split()
        cm = gffcmp(os.path.join(d, 'cmp.stats'))
        f1 = lambda rp: 2 * rp[0] * rp[1] / (rp[0] + rp[1])
        cpu_h, peak, stages = cost(d)
        print('\t'.join(str(x) for x in (
            c, reads, g['gene_miner_genes'], g['gene_miner_isoforms_per_gene'],
            busco(os.path.join(d, 'busco_primary', 'short_summary.specific*.txt')),
            busco(os.path.join(d, 'gm_out', 'busco', 'busco_union', 'short_summary.specific*.txt')),
            rec, f'{100 * rec / REF_GENES:.1f}', g['gm_novel_loci_absent_from_reference'],
            sf[4], sf[7], cm['Locus'][0], cm['Locus'][1], f'{f1(cm["Locus"]):.1f}',
            cm['Transcript'][0], cm['Transcript'][1], f'{f1(cm["Transcript"]):.1f}',
            f'{cpu_h:.1f}', f'{peak:.1f}', f'{stages.get("GENEMARK_ETP", 0):.1f}')))


if __name__ == '__main__':
    main()
