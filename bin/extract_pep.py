import re,sys
from collections import defaultdict
gen=sys.argv[1]; gff=sys.argv[2]; out=sys.argv[3]
# load genome
seq={}; name=None; buf=[]
for l in open(gen):
    if l[0]=='>':
        if name: seq[name]="".join(buf)
        name=l[1:].split()[0]; buf=[]
    else: buf.append(l.strip())
if name: seq[name]="".join(buf)
comp=str.maketrans("ACGTNacgtn","TGCANtgcan")
def rc(s): return s.translate(comp)[::-1]
code={}
bases="TCAG"
aas="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i=0
for a in bases:
    for b in bases:
        for c in bases:
            code[a+b+c]=aas[i]; i+=1
def transl(nt):
    p=[]
    for i in range(0,len(nt)-2,3):
        p.append(code.get(nt[i:i+3].upper(),'X'))
    return "".join(p)
# parse gff: mRNA -> CDS list (start,end,strand,phase,chrom)
cds=defaultdict(list); chrom={}; strand={}
for l in open(gff):
    if l.startswith('#') or not l.strip(): continue
    f=l.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='CDS': continue
    par=re.search(r'Parent=([^;]+)',f[8]).group(1)
    cds[par].append((int(f[3]),int(f[4]),f[6],f[7])); chrom[par]=f[0]; strand[par]=f[6]
n=0; bad=0
with open(out,'w') as o:
    for mid,iv in cds.items():
        c=chrom[mid]; st=strand[mid]
        if c not in seq: bad+=1; continue
        iv.sort()
        try:
            if st=='+':
                nt="".join(seq[c][a-1:b] for a,b,s,p in iv)
                ph=int(iv[0][3]) if iv[0][3].isdigit() else 0
            else:
                nt="".join(rc(seq[c][a-1:b]) for a,b,s,p in reversed(iv))
                ph=int(iv[-1][3]) if iv[-1][3].isdigit() else 0
            nt=nt[ph:]
            pep=transl(nt).rstrip('*')
            if len(pep)>=1:
                o.write(">%s\n%s\n"%(mid,pep)); n+=1
        except Exception: bad+=1
print("proteins written:",n," skipped:",bad)
