#!/usr/bin/env python3
# novel_tiers.py gm.gff3 ref.gff3 transdecoder.gff3 emapper.annotations blastp.outfmt6
# Classifies GM loci vs reference (same-strand CDS overlap >=0.5 -> novel if none),
# then for NOVEL loci computes 3 independent-evidence flags and confidence tiers.
import sys, re, bisect
from collections import defaultdict
gm_f, ref_f, td_f, emap_f, blast_f = sys.argv[1:6]
def attr(a,k):
    m=re.search(k+r'=([^;]+)',a); return m.group(1) if m else None
def parse_cds(path, pc=False):
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
gm=parse_cds(gm_f); ref=parse_cds(ref_f,pc=True)
ridx=defaultdict(list)
for v in ref.values(): ridx[(v['chrom'],v['strand'])].append(v)
novel=[]; recovered=set()
for gid,gv in gm.items():
    best=0.0
    for rv in ridx.get((gv['chrom'],gv['strand']),[]):
        if rv['end']<gv['start'] or rv['start']>gv['end']: continue
        r=ofrac(gv,rv)
        if r>best: best=r
    (recovered.add(gid) if best>=0.5 else novel.append(gid))
# RNA overlap via merged transdecoder exon intervals
tdi=defaultdict(list)
for line in open(td_f):
    if line.startswith('#') or not line.strip(): continue
    f=line.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='CDS': continue
    tdi[(f[0],f[6])].append((int(f[3]),int(f[4])))
merged={}
for k,ivs in tdi.items():
    ivs.sort(); m=[]
    for s,e in ivs:
        if m and s<=m[-1][1]+1: m[-1]=[m[-1][0],max(m[-1][1],e)]
        else: m.append([s,e])
    merged[k]=(([x[0] for x in m]),m)
def has_rna(gv):
    kv=merged.get((gv['chrom'],gv['strand']))
    if not kv: return False
    starts,m=kv
    for s,e in gv['cds']:
        i=bisect.bisect_right(starts,e)-1
        if i>=0 and m[i][1]>=s: return True
    return False
def gene_of(q): return re.sub(r'\.t\d+$','',q)
egg=set()
for line in open(emap_f):
    if line.startswith('#') or not line.strip(): continue
    egg.add(gene_of(line.split('\t',1)[0]))
sp=set()
for line in open(blast_f):
    if not line.strip(): continue
    sp.add(gene_of(line.split('\t',1)[0]))
def marg(ids):
    n=len(ids) or 1
    r=sum(has_rna(gm[g]) for g in ids); s=sum(g in sp for g in ids); e=sum(g in egg for g in ids)
    return len(ids), 100*r/n, 100*s/n, 100*e/n
# verify against published: identical/revised/novel marginals
allids=list(gm); rec=[g for g in gm if g in recovered]
print("category\tloci\tRNA%\tSwissP%\teggOrth%")
for name,ids in [("novel",novel),("recovered",rec),("all",allids)]:
    n,r,s,e=marg(ids); print(f"{name}\t{n}\t{r:.1f}\t{s:.1f}\t{e:.1f}")
# tiers on novel
cnt=defaultdict(int)
for g in novel:
    k=int(has_rna(gm[g]))+int(g in sp)+int(g in egg); cnt[k]+=1
N=len(novel) or 1
print("\n#evidence\tloci\tpct")
for k in (3,2,1,0): print(f"{k}\t{cnt[k]}\t{100*cnt[k]/N:.1f}")
hi=cnt[3]+cnt[2]; me=cnt[1]; lo=cnt[0]
print(f"\nHIGH(>=2)\t{hi}\t{100*hi/N:.1f}")
print(f"MED(1)\t{me}\t{100*me/N:.1f}")
print(f"LOW(0)\t{lo}\t{100*lo/N:.1f}")

# --- optional per-locus TSV for novel loci ---
if len(sys.argv) > 7:
    label, outtsv = sys.argv[6], sys.argv[7]
    with open(outtsv, 'w') as o:
        o.write("genome\tlocus_id\tstream\trna_seq\tswissprot\teggnog\tn_evidence\tconfidence\n")
        for g in novel:
            r=int(has_rna(gm[g])); s=int(g in sp); e=int(g in egg); n=r+s+e
            tier = "high" if n>=2 else ("medium" if n==1 else "low")
            stream = ("AUGUSTUS" if "_A" in g else "GeneMark-ETP" if "_E" in g else "RNA-seq" if "_R" in g else "?")
            o.write(f"{label}\t{g}\t{stream}\t{r}\t{s}\t{e}\t{n}\t{tier}\n")
    sys.stderr.write(f"wrote {outtsv}\n")
