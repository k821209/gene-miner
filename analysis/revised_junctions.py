#!/usr/bin/env python3
# revised_junctions.py gm_final.gff3 ref.gff3 junctions.tsv min_reads
# junctions.tsv = chrom<TAB>intron_start<TAB>intron_end<TAB>read_count (1-based, from regtools; see README)
# For REVISED loci (GM matched a ref gene same-strand >=0.5 but CDS structure differs),
# compare how well each side's CDS introns are supported by RNA-seq splice junctions.
import sys, re
from collections import defaultdict
gm_f, ref_f, jbed, minr = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv)>4 else 1
def attr(a,k):
    m=re.search(k+r'=([^;]+)',a); return m.group(1) if m else None
def parse(path, pc=False):
    # returns genes: gid -> dict(chrom,strand,cds=[(s,e)],introns=set((s,e)),start,end)
    txg={}; bt={}; tx_cds=defaultdict(list); tx_ord={}
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9: continue
        t,a=f[2],f[8]
        if t=='gene': bt[attr(a,'ID')]=attr(a,'biotype') or ''
        elif t in ('mRNA','transcript'): txg[attr(a,'ID')]=attr(a,'Parent')
        elif t=='CDS':
            par=attr(a,'Parent'); tx_cds[par].append((int(f[3]),int(f[4]),f[0],f[6]))
    genes={}
    for tx,segs in tx_cds.items():
        g=txg.get(tx,tx)
        chrom=segs[0][2]; strand=segs[0][3]
        iv=sorted((s,e) for s,e,_,_ in segs)
        if g not in genes:
            genes[g]=dict(chrom=chrom,strand=strand,cds=set(),introns=set())
        genes[g]['cds'].update(iv)
        for i in range(len(iv)-1):
            genes[g]['introns'].add((iv[i][1]+1, iv[i+1][0]-1))
    out={}
    for g,d in genes.items():
        if pc and bt.get(g,'') not in ('','protein_coding'): continue
        cds=sorted(d['cds'])
        out[g]=dict(chrom=d['chrom'],strand=d['strand'],cds=cds,introns=d['introns'],
                    start=min(s for s,e in cds),end=max(e for s,e in cds))
    return out
def clen(iv): return sum(e-s+1 for s,e in iv)
def ofrac(a,b):
    ov=0
    for s,e in a:
        for bs,be in b:
            lo,hi=max(s,bs),min(e,be)
            if hi>=lo: ov+=hi-lo+1
    return ov/(min(clen(a),clen(b)) or 1)
gm=parse(gm_f); ref=parse(ref_f, pc=True)
# junction set: regtools BED12 -> intron (chrom, s1based, e1based) with score>=minr
J=set()
for line in open(jbed):
    f=line.rstrip('\n').split('\t')
    if len(f)<4: continue
    if int(f[3])<minr: continue
    J.add((f[0],int(f[1]),int(f[2])))
def supp(chrom, introns):
    if not introns: return None
    n=sum(1 for (s,e) in introns if (chrom,s,e) in J)
    return n/len(introns)
# ref index by (chrom,strand)
ridx=defaultdict(list)
for gid,v in ref.items(): ridx[(v['chrom'],v['strand'])].append((gid,v))
# GM-centric matching, find revised pairs
gm_better=gm_eq=gm_worse=0; sum_gm=sum_ref=0.0; nrev=0; both_multi=0
gmonly=gmonly_sup=refonly=refonly_sup=shared=shared_sup=0
for gid,gv in gm.items():
    best=None; bro=0.0
    for rgid,rv in ridx.get((gv['chrom'],gv['strand']),[]):
        if rv['end']<gv['start'] or rv['start']>gv['end']: continue
        r=ofrac(gv['cds'],rv['cds'])
        if r>bro: bro=r; best=(rgid,rv)
    if not best or bro<0.5: continue
    rv=best[1]
    if gv['cds']==rv['cds']: continue   # identical, skip
    # revised
    nrev+=1
    gs=supp(gv['chrom'], gv['introns']); rs=supp(rv['chrom'], rv['introns'])
    if gs is None or rs is None: continue   # need both multi-exon to compare junctions
    both_multi+=1; sum_gm+=gs; sum_ref+=rs
    gi=gv['introns']; ri=rv['introns']; c=gv['chrom']
    for it in gi-ri:
        gmonly+=1; gmonly_sup+= (c,it[0],it[1]) in J
    for it in ri-gi:
        refonly+=1; refonly_sup+= (c,it[0],it[1]) in J
    for it in gi&ri:
        shared+=1; shared_sup+= (c,it[0],it[1]) in J
    if gs>rs+1e-9: gm_better+=1
    elif gs<rs-1e-9: gm_worse+=1
    else: gm_eq+=1
print(f"revised_loci\t{nrev}")
print(f"comparable_both_multiexon\t{both_multi}")
print(f"mean_GM_intron_support\t{sum_gm/max(both_multi,1):.3f}")
print(f"mean_REF_intron_support\t{sum_ref/max(both_multi,1):.3f}")
print(f"GM_better\t{gm_better}\t{100*gm_better/max(both_multi,1):.1f}%")
print(f"equal\t{gm_eq}\t{100*gm_eq/max(both_multi,1):.1f}%")
print(f"GM_worse\t{gm_worse}\t{100*gm_worse/max(both_multi,1):.1f}%")

print(f"--- intron-level differences (revised, both multi-exon) ---")
print(f"shared_introns\t{shared}\tsupported\t{100*shared_sup/max(shared,1):.1f}%")
print(f"GM_only_introns\t{gmonly}\tsupported\t{100*gmonly_sup/max(gmonly,1):.1f}%")
print(f"REF_only_introns\t{refonly}\tsupported\t{100*refonly_sup/max(refonly,1):.1f}%")
