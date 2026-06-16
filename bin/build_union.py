import re
from collections import defaultdict
def load(gff, kind):
    cds=defaultdict(list); m2g={}; gscore={}; gstr={}; gchr={}
    for l in open(gff):
        if l.startswith('#') or not l.strip(): continue
        f=l.rstrip('\n').split('\t')
        if len(f)<9: continue
        at=f[8]
        def A(k):
            m=re.search(k+r'=([^;]+)',at); return m.group(1) if m else None
        if f[2]=='gene':
            try: gscore[A('ID')]=float(f[5])
            except: gscore[A('ID')]=0.0
        elif f[2] in ('mRNA','transcript'):
            m2g[A('ID')]=A('Parent')
        elif f[2]=='CDS':
            g=m2g.get(A('Parent'),A('Parent')); cds[g].append((int(f[3]),int(f[4]),f[7]))
            gstr[g]=f[6]; gchr[g]=f[0]
    genes={}
    for g,iv in cds.items():
        iv.sort()
        aa=sum(b-a+1 for a,b,p in iv)//3
        genes[g]=dict(chrom=gchr[g], start=min(a for a,b,p in iv), end=max(b for a,b,p in iv),
                      strand=gstr[g], aa=aa, score=gscore.get(g,1.0), cds=iv)
    return genes
aug=load('augustus_scaffold.gff3','aug')
td =load('annot/genome.transdecoder.gff3','td')
SC=0.8; MINAA=100
passing={g:v for g,v in aug.items() if v['score']>=SC and v['aa']>=MINAA}
# overlap index of passing augustus
idx=defaultdict(list)
for g,v in passing.items(): idx[v['chrom']].append((v['start'],v['end'],v['strand']))
def ov(c,s,e,st):
    for as_,ae,ast in idx[c]:
        if ast==st and not(e<as_ or s>ae): return True
    return False
rna_only={g:v for g,v in td.items() if not ov(v['chrom'],v['start'],v['end'],v['strand'])}
# write union gff3 (gene/mRNA/CDS, clean IDs)
def write(o, src, genes, prefix):
    n=0
    for g,v in genes.items():
        n+=1; gid="%s%06d"%(prefix,n); mid=gid+".t1"
        c,s,e,st=v['chrom'],v['start'],v['end'],v['strand']
        o.write(f"{c}\t{src}\tgene\t{s}\t{e}\t.\t{st}\t.\tID={gid}\n")
        o.write(f"{c}\t{src}\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID={mid};Parent={gid}\n")
        for a,b,ph in v['cds']:
            o.write(f"{c}\t{src}\tCDS\t{a}\t{b}\t.\t{st}\t{ph}\tID={mid}.cds;Parent={mid}\n")
    return n
with open('union.gff3','w') as o:
    na=write(o,'AUGUSTUS',passing,'ZscG_A')
    nr=write(o,'RNAseq',rna_only,'ZscG_R')
print(f"usable AUGUSTUS (score>={SC},>={MINAA}aa): {na}")
print(f"RNA-only (AUGUSTUS-missed): {nr}")
print(f"=== UNION final genes: {na+nr} ===")
