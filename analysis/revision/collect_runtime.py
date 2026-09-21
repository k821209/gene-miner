#!/usr/bin/env python3
"""collect_runtime.py — per-stage cost of each full Gene-Miner run from its Nextflow
trace (Reviewer 1 point 5, Reviewer 2 minor 2): summed task real time, CPU-hours
(real time x CPU utilisation) and peak resident memory per stage.
Usage: collect_runtime.py <runtime dir> <genome> [<genome> ...]"""
import os, re, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_sub import secs, gb

STAGE = [('MASK_GENOME', 'Repeat masking (RepeatMasker, supplied library)'),
         ('HISAT2', 'Read alignment (HISAT2 index + align)'),
         ('STRINGTIE', 'Transcript assembly (StringTie + merge)'),
         ('TRANSDECODER', 'ORF prediction (TransDecoder + DIAMOND)'),
         ('BUSCO_TRAIN', 'AUGUSTUS training (BUSCO)'),
         ('AUGUSTUS', 'AUGUSTUS prediction'),
         ('GENEMARK', 'GeneMark-ETP'),
         ('EGGNOG', 'eggNOG-mapper'),
         ('SPROT', 'Swiss-Prot search (tiers)'),
         ('BUSCO', 'BUSCO (catalogue)'),
         ('', 'Union, filters and tiers')]

def stage(name):
    for k, lab in STAGE:
        if k and name.startswith(k):
            return lab
    return STAGE[-1][1]

print('genome\tstage\ttasks\treal_h\tcpu_h\tpeak_rss_gb')
for g in sys.argv[2:]:
    L = [l.rstrip('\n').split('\t') for l in open(os.path.join(sys.argv[1], g, 'trace.txt'))]
    h = L[0]; ix = {k: h.index(k) for k in ('name', 'status', 'realtime', '%cpu', 'peak_rss')}
    agg = defaultdict(lambda: [0, 0.0, 0.0, 0.0])
    for r in L[1:]:
        if r[ix['status']] != 'COMPLETED':
            continue
        rt = secs(r[ix['realtime']]); pc = float(r[ix['%cpu']].rstrip('%') or 0) / 100
        for s in (stage(r[ix['name']]), 'TOTAL'):
            a = agg[s]; a[0] += 1; a[1] += rt / 3600; a[2] += rt * pc / 3600; a[3] = max(a[3], gb(r[ix['peak_rss']]))
    for k, lab in STAGE + [('', 'TOTAL')]:
        if lab in agg:
            a = agg[lab]; print(f'{g}\t{lab}\t{a[0]}\t{a[1]:.2f}\t{a[2]:.1f}\t{a[3]:.1f}')
