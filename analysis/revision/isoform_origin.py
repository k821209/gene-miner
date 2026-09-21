#!/usr/bin/env python3
"""isoform_origin.py — trace every transcript of a Gene-Miner catalogue to its
stream by exact CDS-interval match, and count the transcripts beyond one per gene
("added isoforms") by stream.

Usage: isoform_origin.py union.final.gff3 augustus.gff3 transdecoder.gff3 genemark.gtf
"""
import sys
from collections import defaultdict, Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
from fusion_origin import tx_cds_gff, keyset_gff, keyset_gtf

gm, aug, td, gmk = sys.argv[1:5]
t2g, cds, loc = tx_cds_gff(gm)
src = {'AUGUSTUS': keyset_gff(aug), 'RNA-seq': keyset_gff(td), 'GeneMark-ETP': keyset_gtf(gmk)}
g2t = defaultdict(list)
for t in cds:
    g2t[t2g.get(t, t)].append(t)
all_c, added = Counter(), Counter()
for g, ts in g2t.items():
    lab = []
    for t in ts:
        k = (loc[t][0], loc[t][1], tuple(sorted(set(cds[t]))))
        n = [s for s, ks in src.items() if k in ks]
        # precedence when a structure is emitted by several streams: AUGUSTUS > GeneMark > RNA
        lab.append("RNA-seq" if "RNA-seq" in n else "GeneMark-ETP" if "GeneMark-ETP" in n
                   else "AUGUSTUS" if n else "unmatched")
        all_c['+'.join(n) or 'unmatched'] += 1
    # the first transcript founds the locus; the rest are added isoforms
    for l in lab[1:]:
        added[l] += 1
print('transcripts_by_exact_source', dict(all_c))
print('added_isoforms_by_stream', dict(added), 'total', sum(added.values()))
