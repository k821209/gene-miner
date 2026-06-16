import sys, re
from collections import defaultdict
def load_aug(gff):
    genes={}; cds=defaultdict(list); m2g={}; score={}
    for l in open(gff):
        if l.startswith('#') or not l.strip(): continue
        f=l.rstrip('\n').split('\t')
        if len(f)<9: continue
        at=f[8]; A=lambda k:(re.search(k+r'=([^;]+)',at) or [None,None])[1] if re.search(k+r'=([^;]+)',at) else None
        if f[2]=='gene':
            gid=A('ID'); 
            try: score[gid]=float(f[5])
            except: score[gid]=0.0
        elif f[2] in('mRNA','transcript'):
            m2g[A('ID')]=A('Parent')
        elif f[2]=='CDS':
            g=m2g.get(A('Parent'),A('Parent')); cds[g].append((int(f[3]),int(f[4]),f[0],f[6]))
    for g,iv in cds.items():
        genes[g]=(iv[0][2],min(a for a,b,c,d in iv),max(b for a,b,c,d in iv),iv[0][3],sum(b-a+1 for a,b,c,d in iv)//3,score.get(g,0.0))
    return genes
def load_td(gff):
    genes={}; cds=defaultdict(list); m2g={}
    for l in open(gff):
        if l.startswith('#') or not l.strip(): continue
        f=l.rstrip('\n').split('\t')
        if len(f)<9: continue
        at=f[8]; A=lambda k:(re.search(k+r'=([^;]+)',at).group(1) if re.search(k+r'=([^;]+)',at) else None)
        if f[2]=='mRNA': m2g[A('ID')]=A('Parent')
        elif f[2]=='CDS':
            g=m2g.get(A('Parent'),A('Parent')); cds[g].append((int(f[3]),int(f[4]),f[0],f[6]))
    for g,iv in cds.items():
        genes[g]=(iv[0][2],min(a for a,b,c,d in iv),max(b for a,b,c,d in iv),iv[0][3])
    return genes
aug=load_aug(sys.argv[1]); td=load_td(sys.argv[2])
tidx=defaultdict(list)
for g,(c,s,e,st) in td.items(): tidx[c].append((s,e,st))
def ovTD(c,s,e,st):
    for ts,te,tst in tidx[c]:
        if tst==st and not(e<ts or s>te): return True
    return False
for SC in [0.0,0.3,0.5,0.8]:
    usable={g:v for g,v in aug.items() if v[4]>=100 and v[5]>=SC}
    uidx=defaultdict(list)
    for g,(c,s,e,st,aa,sc) in usable.items(): uidx[c].append((s,e,st))
    def ovAUG(c,s,e,st):
        for as_,ae,ast in uidx[c]:
            if ast==st and not(e<as_ or s>ae): return True
        return False
    rna_only=[g for g,(c,s,e,st) in td.items() if not ovAUG(c,s,e,st)]
    union=len(usable)+len(rna_only)
    print("score>=%.1f & >=100aa: usable_AUG=%d  RNA-only=%d  UNION=%d"%(SC,len(usable),len(rna_only),union))
