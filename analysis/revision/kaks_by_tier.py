#!/usr/bin/env python3
"""kaks_by_tier.py — selection on the novel loci, split by the confidence tier the
catalogue assigns them. Joins the per-locus dN/dS table to the tier attribute in
the annotated GFF3.

Usage: kaks_by_tier.py union.final.tiers.gff3 <kaks dir>
"""
import os
import re
import statistics
import sys
from collections import defaultdict

gff, kdir = sys.argv[1], sys.argv[2]
tier = {}
for line in open(gff):
    f = line.split('\t')
    if len(f) > 8 and f[2] == 'gene':
        g = re.search(r'ID=([^;]+)', f[8]).group(1)
        t = re.search(r'tier=([^;\s]+)', f[8])
        if t:
            tier[g] = t.group(1)

hdr = ('set', 'n', 'median dN/dS', 'frac<1', 'frac<0.5', 'median dS')
print('%-34s %5s %12s %7s %9s %10s' % hdr)
for sp in ('glaberrima', 'punctata'):
    for cls in ('novel', 'identical', 'revised', 'te_removed'):
        p = os.path.join(kdir, '%s.%s' % (sp, cls), 'kaks.tsv')
        if not os.path.exists(p):
            continue
        rows = defaultdict(list)
        for line in open(p):
            if line.startswith('query'):
                continue
            f = line.rstrip('\n').split('\t')
            try:
                w, ds = float(f[5]), float(f[4])
            except (ValueError, IndexError):
                continue
            if ds <= 0 or ds > 2:        # zero-divergence and saturated pairs carry no signal
                continue
            key = tier.get(f[0], 'n/a') if cls == 'novel' else 'all'
            rows[key].append((w, ds))
        for key in ('high', 'medium', 'low', 'all', 'n/a'):
            v = rows.get(key)
            if not v:
                continue
            w = [x[0] for x in v]
            ds = [x[1] for x in v]
            name = '%s %s%s' % (sp, cls, '' if key in ('all', 'n/a') else ' ' + key)
            print('%-34s %5d %12.3f %7.2f %9.2f %10.3f' % (
                name, len(v), statistics.median(w),
                sum(x < 1 for x in w) / len(w), sum(x < 0.5 for x in w) / len(w),
                statistics.median(ds)))
