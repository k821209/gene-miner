import sys
agp, gff, out = sys.argv[1], sys.argv[2], sys.argv[3]
# contig -> (scaffold, scaf_start, orientation, contig_len)
m={}
for l in open(agp):
    f=l.rstrip("\n").split("\t")
    if len(f)<9 or f[4]!="W": continue
    scaf, sstart, send, part, typ, cid, cstart, cend, orient = f[:9]
    m[cid]=(scaf, int(sstart), orient, int(cend))  # cend = contig length (cstart=1)
no=0; lifted=0
with open(out,"w") as o:
    for l in open(gff):
        if l.startswith("#") or not l.strip(): 
            o.write(l); continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9: continue
        cid=f[0]
        if cid not in m: no+=1; continue
        scaf, S, orient, L = m[cid]
        a,b=int(f[3]),int(f[4]); strand=f[6]
        if orient=="+":
            na, nb = S+(a-1), S+(b-1); ns=strand
        else:
            na, nb = S+(L-b), S+(L-a); ns={"+":"-","-":"+","." :"."}.get(strand,strand)
        f[0]=scaf; f[3]=str(na); f[4]=str(nb); f[6]=ns
        o.write("\t".join(f)+"\n"); lifted+=1
print("lifted features:",lifted,"; contigs-not-in-agp skipped:",no)
