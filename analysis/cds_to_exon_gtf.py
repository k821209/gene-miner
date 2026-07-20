#!/usr/bin/env python3
# Convert a GFF3 to a CDS-only GTF where each CDS is relabeled "exon", so
# gffcompare compares protein-CODING structures (drops ncRNA, ignores UTRs).
import sys, re
inf, outf = sys.argv[1], sys.argv[2]
def attr(s, k):
    m = re.search(rf'(?:^|;)\s*{k}=([^;]+)', s)
    return m.group(1) if m else None
mrna2gene = {}
with open(inf) as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9: continue
        if c[2] in ('mRNA','transcript') or c[2].endswith('RNA'):
            tid = attr(c[8],'ID'); gid = attr(c[8],'Parent') or tid
            if tid: mrna2gene[tid] = gid
n=0
with open(inf) as f, open(outf,'w') as o:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'CDS': continue
        p = attr(c[8],'Parent')
        if not p: continue
        tid = p.split(',')[0]; gid = mrna2gene.get(tid, tid)
        c[2] = 'exon'
        c[8] = f'transcript_id "{tid}"; gene_id "{gid}";'
        o.write('\t'.join(c)+'\n'); n+=1
sys.stderr.write(f"{outf}: {n} CDS-exon lines, {len(mrna2gene)} transcripts mapped\n")
