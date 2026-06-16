#!/usr/bin/env python3
"""
filter_models.py — tiered evidence reporting for RNA-seq/TransDecoder gene models.

Philosophy: the upstream StringTie+TransDecoder pipeline is intentionally
HIGH-SENSITIVITY (capture as many genes/mRNAs/isoforms as possible). This
script does NOT throw that away — it keeps the comprehensive set and, in
addition, derives an evidence-supported HIGH-CONFIDENCE subset and reports
detailed tiers so the spurious "protein-shaped" ORFs can be quantified and
excluded when a clean count is needed.

Evidence = protein homology (DIAMOND/BLAST blastp) OR a Pfam domain (hmmscan).

Outputs:
  <out>                 high-confidence gff3 (genes with >=1 evidence-supported mRNA)
  <out>.comprehensive   passthrough of all models (= input; the max-capture set)
  <summary>             tiered statistics (genes, mRNAs, isoforms, evidence types)
"""
import argparse, re, sys
from collections import defaultdict

def ids_blastp(p):
    s=set()
    for l in open(p):
        if l.strip(): s.add(l.split('\t',1)[0])
    return s

def ids_pfam(p):
    s=set()
    for l in open(p):
        if l.startswith('#') or not l.strip(): continue
        f=l.split()
        if len(f)>3: s.add(f[3])   # domtblout query name = ORF/pep id
    return s

def attr(s,k):
    m=re.search(k+r'=([^;]+)',s); return m.group(1) if m else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--gff',required=True)
    ap.add_argument('--blastp',required=True)
    ap.add_argument('--pfam')
    ap.add_argument('--out',required=True)
    ap.add_argument('--summary',required=True)
    a=ap.parse_args()

    bp=ids_blastp(a.blastp)
    pf=ids_pfam(a.pfam) if a.pfam else set()
    ev=bp|pf

    lines=[l for l in open(a.gff) if l.strip() and not l.startswith('#')]
    genes={}; gene_lines=defaultdict(list); gene_order=[]
    mrna_of_gene=defaultdict(list); gene_of_mrna={}
    for l in lines:
        f=l.rstrip('\n').split('\t')
        if len(f)<9: continue
        t=f[2]; at=f[8]
        if t=='gene':
            gid=attr(at,'ID'); gene_order.append(gid); gene_lines[gid].append(l)
        elif t=='mRNA':
            mid=attr(at,'ID'); par=attr(at,'Parent')
            gene_of_mrna[mid]=par; mrna_of_gene[par].append(mid); gene_lines[par].append(l)
        else:
            par=attr(at,'Parent'); gid=gene_of_mrna.get(par)
            if gid: gene_lines[gid].append(l)

    def mrna_supported(mid):
        if mid in ev: return True
        base=mid.split('.mRNA')[0]
        return base in ev or any(mid in e for e in ()) # exact/base match only

    sup_b={m:(m in bp or m.split('.mRNA')[0] in bp) for m in gene_of_mrna}
    sup_p={m:(m in pf or m.split('.mRNA')[0] in pf) for m in gene_of_mrna}
    sup  ={m:(sup_b[m] or sup_p[m]) for m in gene_of_mrna}

    hc_genes=[g for g in gene_order if any(sup[m] for m in mrna_of_gene[g])]
    hc_set=set(hc_genes)

    # write high-confidence gff
    with open(a.out,'w') as o:
        for g in gene_order:
            if g in hc_set:
                for l in gene_lines[g]: o.write(l)
    # comprehensive passthrough
    with open(a.out+'.comprehensive','w') as o:
        for g in gene_order:
            for l in gene_lines[g]: o.write(l)

    n_gene=len(gene_order); n_mrna=len(gene_of_mrna)
    iso=[len(mrna_of_gene[g]) for g in gene_order]
    multi=sum(1 for x in iso if x>1)
    nb=sum(sup_b.values()); npf=sum(sup_p.values()); ne=sum(sup.values())
    with open(a.summary,'w') as s:
        s.write("=== RNA-seq annotation: tiered summary ===\n\n")
        s.write("[ COMPREHENSIVE  (max-capture set) ]\n")
        s.write(f"  genes                 : {n_gene}\n")
        s.write(f"  mRNAs (isoforms)      : {n_mrna}\n")
        s.write(f"  mean mRNAs/gene       : {n_mrna/max(n_gene,1):.2f}\n")
        s.write(f"  multi-isoform genes   : {multi}\n\n")
        s.write("[ EVIDENCE per mRNA ]\n")
        s.write(f"  blastp (UniProt) hit  : {nb} ({100*nb/max(n_mrna,1):.1f}%)\n")
        s.write(f"  Pfam domain           : {npf} ({100*npf/max(n_mrna,1):.1f}%)\n")
        s.write(f"  either (supported)    : {ne} ({100*ne/max(n_mrna,1):.1f}%)\n\n")
        s.write("[ HIGH-CONFIDENCE  (>=1 supported mRNA) ]\n")
        s.write(f"  genes                 : {len(hc_genes)}\n")
        s.write(f"  (spurious/unsupported genes removed: {n_gene-len(hc_genes)})\n")
    print(open(a.summary).read())

if __name__=='__main__':
    main()
