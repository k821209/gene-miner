#!/usr/bin/env python3
# classify_missed.py raw_union.gff3 final_union.gff3 ref.gff3 [min_ro=0.5]
# For each reference gene NOT recovered by the final catalogue, classify why:
#   qc_removed / wrong_strand / structural_mismatch / prediction_failure
import sys, re
from collections import defaultdict
raw_f, fin_f, ref_f = sys.argv[1:4]
RO = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
def attr(a,k):
    m=re.search(k+r'=([^;]+)',a); return m.group(1) if m else None
def parse(path, pc=False):
    txg={}; bt={}; cds=defaultdict(list); ch={}; st={}
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9: continue
        t,a=f[2],f[8]
        if t=='gene': bt[attr(a,'ID')]=attr(a,'biotype') or ''
        elif t in ('mRNA','transcript'): txg[attr(a,'ID')]=attr(a,'Parent')
        elif t=='CDS':
            g=txg.get(attr(a,'Parent'),attr(a,'Parent'))
            cds[g].append((int(f[3]),int(f[4]))); ch[g]=f[0]; st[g]=f[6]
    G={}
    for g,iv in cds.items():
        if pc and bt.get(g,'') not in ('','protein_coding'): continue
        iv=sorted(set(iv)); G[g]=dict(chrom=ch[g],strand=st[g],cds=iv,start=min(s for s,e in iv),end=max(e for s,e in iv))
    return G
def clen(iv): return sum(e-s+1 for s,e in iv)
def ofrac(a,b):
    ov=0
    for s,e in a['cds']:
        for bs,be in b['cds']:
            lo,hi=max(s,bs),min(e,be)
            if hi>=lo: ov+=hi-lo+1
    return ov/(min(clen(a['cds']),clen(b['cds'])) or 1)
raw=parse(raw_f); fin=parse(fin_f); ref=parse(ref_f, pc=True)
# index final and raw by chrom
def idx_by_chrom(G):
    d=defaultdict(list)
    for gid,v in G.items(): d[v['chrom']].append(v)
    for k in d: d[k].sort(key=lambda x:x['start'])
    return d
fin_c=idx_by_chrom(fin); raw_c=idx_by_chrom(raw)
removed_ids=set(raw)-set(fin)
removed_c=idx_by_chrom({g:raw[g] for g in removed_ids})
# ref index by (chrom,strand) -> list of (gid, v)
ref_cs=defaultdict(list)
for gid,v in ref.items(): ref_cs[(v['chrom'],v['strand'])].append((gid,v))
# GM-centric matching (same as gm_compare): each final locus -> best same-strand ref >=RO
matched_ref=set()
for gid,gv in fin.items():
    best=None; bro=0.0
    for rgid,rv in ref_cs.get((gv['chrom'],gv['strand']),[]):
        if rv['end']<gv['start'] or rv['start']>gv['end']: continue
        r=ofrac(gv,rv)
        if r>bro: bro=r; best=rgid
    if best is not None and bro>=RO: matched_ref.add(best)
total_ref=len(ref); recovered=len(matched_ref); nr=total_ref-recovered
def best_opp_same(R):
    bs=bo=0.0
    for v in fin_c.get(R['chrom'],[]):
        if v['end']<R['start'] or v['start']>R['end']: continue
        r=ofrac(R,v)
        if v['strand']==R['strand']: bs=max(bs,r)
        else: bo=max(bo,r)
    return bs,bo
cnt=defaultdict(int)
for gid,R in ref.items():
    if gid in matched_ref: continue
    # QC-removed
    qc=False
    for v in removed_c.get(R['chrom'],[]):
        if v['strand']!=R['strand']: continue
        if v['end']<R['start'] or v['start']>R['end']: continue
        if ofrac(R,v)>=RO: qc=True; break
    if qc: cnt['qc_removed']+=1; continue
    bs,bo=best_opp_same(R)
    if bo>=RO: cnt['wrong_strand']+=1; continue
    if bs>0: cnt['structural_mismatch']+=1; continue
    cnt['prediction_failure']+=1
print(f"reference_genes\t{total_ref}")
print(f"recovered\t{recovered}\t{100*recovered/total_ref:.1f}%")
print(f"not_recovered\t{nr}\t{100*nr/total_ref:.1f}%")
for k in ('qc_removed','wrong_strand','structural_mismatch','prediction_failure'):
    print(f"  {k}\t{cnt[k]}\t{100*cnt[k]/max(nr,1):.1f}%")
