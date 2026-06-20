#!/usr/bin/env python3
"""gm_compare.py — compare a Gene-Miner catalogue against a published reference
annotation to quantify a gene-model UPDATE: how many reference genes are
recovered, how many Gene-Miner loci are novel (absent from the reference), how
many reference genes are missed, and (among recovered) how many have an
identical vs a revised CDS structure.

Loci are matched by same-strand CDS-footprint overlap (reciprocal >= --min_ro of
the shorter locus). Structure is "identical" when the sorted CDS interval set is
exactly equal.

Usage:
  gm_compare.py --gm union.final.gff3 --ref reference.gff3 [--min_ro 0.5] [--out report.tsv]
"""
import argparse, re, sys
from collections import defaultdict

def parse_gff(path, want_protein_coding=False):
    """Return {gene_id: {chrom,strand,cds:[(s,e)...]}} keyed by the mRNA's gene.
    Robust to Ensembl (ID=transcript:..,Parent=gene:..) and Gene-Miner GFF3."""
    tx_gene = {}            # mRNA id -> gene id
    gene_biotype = {}       # gene id -> biotype (Ensembl) ; '' if unknown
    cds = defaultdict(list) # gene id -> [(s,e)]
    g_chrom, g_strand = {}, {}
    def attr(a, k):
        m = re.search(k + r'=([^;]+)', a)
        return m.group(1) if m else None
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        t, a = f[2], f[8]
        if t == 'gene':
            gid = attr(a, 'ID')
            gene_biotype[gid] = attr(a, 'biotype') or ''
        elif t in ('mRNA', 'transcript'):
            tx_gene[attr(a, 'ID')] = attr(a, 'Parent')
        elif t == 'CDS':
            par = attr(a, 'Parent')
            g = tx_gene.get(par, par)
            cds[g].append((int(f[3]), int(f[4])))
            g_chrom[g] = f[0]; g_strand[g] = f[6]
    genes = {}
    for g, iv in cds.items():
        if want_protein_coding:
            bt = gene_biotype.get(g, '')
            if bt and bt != 'protein_coding':
                continue
        iv = sorted(set(iv))
        genes[g] = dict(chrom=g_chrom[g], strand=g_strand[g], cds=iv,
                        start=min(s for s, e in iv), end=max(e for s, e in iv))
    # transcripts (isoforms) belonging to kept genes
    n_tx = sum(1 for tid, g in tx_gene.items() if g in genes)
    return genes, n_tx

def index(genes):
    idx = defaultdict(list)
    for gid, v in genes.items():
        idx[(v['chrom'], v['strand'])].append((v['start'], v['end'], gid))
    for k in idx:
        idx[k].sort()
    return idx

def cds_len(iv):
    return sum(e - s + 1 for s, e in iv)

def overlaps(a, b):
    """reciprocal CDS-footprint overlap fraction of the shorter locus"""
    ov = 0
    j = 0
    bi = b['cds']
    for s, e in a['cds']:
        for bs, be in bi:
            lo, hi = max(s, bs), min(e, be)
            if hi >= lo:
                ov += hi - lo + 1
    shorter = min(cds_len(a['cds']), cds_len(b['cds'])) or 1
    return ov / shorter

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gm', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--min_ro', type=float, default=0.5)
    ap.add_argument('--out', default='-')
    a = ap.parse_args()

    gm, gm_tx = parse_gff(a.gm)
    ref, ref_tx = parse_gff(a.ref, want_protein_coding=True)
    ridx = index(ref)

    matched_ref = set()
    gm_recovered = gm_novel = identical = revised = 0
    for gid, gv in gm.items():
        cand = ridx.get((gv['chrom'], gv['strand']), [])
        best = None; best_ro = 0.0
        for rs, re_, rid in cand:
            if re_ < gv['start'] or rs > gv['end']:
                continue
            ro = overlaps(gv, ref[rid])
            if ro > best_ro:
                best_ro, best = ro, rid
        if best and best_ro >= a.min_ro:
            gm_recovered += 1
            matched_ref.add(best)
            if gv['cds'] == ref[best]['cds']:
                identical += 1
            else:
                revised += 1
        else:
            gm_novel += 1
    ref_missed = len(ref) - len(matched_ref)

    rows = [
        ('reference_protein_coding_genes', len(ref)),
        ('reference_transcripts', ref_tx),
        ('reference_isoforms_per_gene', round(ref_tx / max(len(ref), 1), 2)),
        ('gene_miner_genes', len(gm)),
        ('gene_miner_transcripts', gm_tx),
        ('gene_miner_isoforms_per_gene', round(gm_tx / max(len(gm), 1), 2)),
        ('gm_recovered_reference_loci', gm_recovered),
        ('  recovered_identical_structure', identical),
        ('  recovered_revised_structure', revised),
        ('gm_novel_loci_absent_from_reference', gm_novel),
        ('reference_loci_missed_by_gm', ref_missed),
        ('min_reciprocal_overlap', a.min_ro),
    ]
    out = sys.stdout if a.out == '-' else open(a.out, 'w')
    for k, v in rows:
        out.write(f"{k}\t{v}\n")
    if out is not sys.stdout:
        out.close()
    sys.stderr.write("done\n")

if __name__ == '__main__':
    main()
