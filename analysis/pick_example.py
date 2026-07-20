#!/usr/bin/env python3
# pick_example.py gm.gff3 ref.gff3 junctions.tsv min_reads
# Find revised loci where GM adds RNA-confirmed introns the reference lacks (clear improvements).
import sys, re
from collections import defaultdict
gm_f, ref_f, jtsv, minr = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv)>4 else 3
def attr(a,k):
    m=re.search(k+r'=([^;]+)',a); return m.group(1) if m else None
def parse(path, pc=False):
    txg={}; bt={}; tx_cds=defaultdict(list)
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9: continue
        t,a=f[2],f[8]
        if t=='gene': bt[attr(a,'ID')]=attr(a,'biotype') or ''
        elif t in ('mRNA','transcript'): txg[attr(a,'ID')]=attr(a,'Parent')
        elif t=='CDS': tx_cds[attr(a,'Parent')].append((int(f[3]),int(f[4]),f[0],f[6]))
    genes={}
    for tx,segs in tx_cds.items():
        g=txg.get(tx,tx); chrom=segs[0][2]; strand=segs[0][3]
        iv=sorted((s,e) for s,e,_,_ in segs)
        d=genes.setdefault(g,dict(chrom=chrom,strand=strand,cds=set(),introns=set()))
        d['cds'].update(iv)
        for i in range(len(iv)-1): d['introns'].add((iv[i][1]+1,iv[i+1][0]-1))
    out={}
    for g,d in genes.items():
        if pc and bt.get(g,'') not in ('','protein_coding'): continue
        cds=sorted(d['cds']); out[g]=dict(chrom=d['chrom'],strand=d['strand'],cds=cds,introns=d['introns'],start=min(s for s,e in cds),end=max(e for s,e in cds))
    return out
def clen(iv): return sum(e-s+1 for s,e in iv)
def ofrac(a,b):
    ov=0
    for s,e in a:
        for bs,be in b:
            lo,hi=max(s,bs),min(e,be)
            if hi>=lo: ov+=hi-lo+1
    return ov/(min(clen(a),clen(b)) or 1)
gm=parse(gm_f); ref=parse(ref_f,pc=True)
J=set()
for line in open(jtsv):
    f=line.rstrip('\n').split('\t')
    if len(f)>=4 and int(f[3])>=minr: J.add((f[0],int(f[1]),int(f[2])))
ridx=defaultdict(list)
for gid,v in ref.items(): ridx[(v['chrom'],v['strand'])].append((gid,v))
cand=[]
for gid,gv in gm.items():
    best=None;bro=0.0
    for rgid,rv in ridx.get((gv['chrom'],gv['strand']),[]):
        if rv['end']<gv['start'] or rv['start']>gv['end']: continue
        r=ofrac(gv['cds'],rv['cds'])
        if r>bro: bro=r;best=(rgid,rv)
    if not best or bro<0.5 or gv['cds']==best[1]['cds']: continue
    rv=best[1]
    gmonly=[it for it in gv['introns']-rv['introns'] if (gv['chrom'],it[0],it[1]) in J]
    refonly_sup=[it for it in rv['introns']-gv['introns'] if (gv['chrom'],it[0],it[1]) in J]
    # good example: GM adds >=2 confirmed introns, ref-only has 0 confirmed, compact span
    span=gv['end']-gv['start']
    if len(gmonly)>=2 and len(refonly_sup)==0 and span<12000:
        cand.append((len(gmonly), gid, best[0], gv['chrom'], gv['start'], gv['end'], gv['strand'], span))
cand.sort(reverse=True)
print("n_gmonly_confirmed\tgm_locus\tref_gene\tchrom\tstart\tend\tstrand\tspan")
for c in cand[:15]: print("\t".join(map(str,c)))
