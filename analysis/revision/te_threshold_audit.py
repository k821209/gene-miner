#!/usr/bin/env python3
"""te_threshold_audit.py — what the TE-overlap threshold keeps and what it throws
away (Reviewer 2, major point 5), for each threshold setting:

  kept_TE_annotated     loci kept whose eggNOG description or Pfam domains name a
                        transposable-element function (transposase, reverse
                        transcriptase, integrase, gag/pol, ...): TE models the
                        filter let through
  removed_ref_matching  loci the filter removed that match a reference gene: genes
                        the reference treats as host genes that the filter discarded
  removed_ref_nonTE     the subset of those whose own annotation is not a TE
                        function, i.e. the plausible domesticated/host genes lost

Usage: te_threshold_audit.py <union.gff3> <new.annotations> <reference.gff3>
                             <thr>=<union.final.gff3> [<thr>=<...> ...]
"""
import os
import re
import sys

sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/bin'))
sys.path.insert(0, os.path.expanduser('~/gene-miner-runs/gene-miner/analysis'))
from gm_compare import parse_gff, index  # noqa: E402
from gm_overlap import ofrac  # noqa: E402

TE_RE = re.compile(r'transpos|retrotrans|reverse transcriptase|integrase|\bgag\b|\bpol\b|'
                   r'ribonuclease h|rnase h|\bty[13]\b|copia|gypsy|helitron|mutator|\bmule\b|'
                   r'\bhat\b|cacta|en/spm|line-1|non-ltr|\bltr\b|rve|zf-cchc|rvt_|transposon',
                   re.I)


def gene_ids(gff):
    return {re.search(r'ID=([^;\s]+)', l.split('\t')[8]).group(1)
            for l in open(gff) if len(l.split('\t')) > 8 and l.split('\t')[2] == 'gene'}


def te_annotated(union_gff, ann):
    t2g = {}
    for l in open(union_gff):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'mRNA':
            t2g[re.search(r'ID=([^;\s]+)', f[8]).group(1)] = re.search(r'Parent=([^;\s]+)', f[8]).group(1)
    H, te = [], set()
    for line in open(ann):
        if line.startswith('#'):
            if line.startswith('#query'):
                H = line.lstrip('#').rstrip('\n').split('\t')
            continue
        r = dict(zip(H, line.rstrip('\n').split('\t')))
        text = ' '.join(r.get(k, '') for k in ('Description', 'PFAMs', 'Preferred_name'))
        if TE_RE.search(text):
            g = t2g.get(r.get('query', ''))
            if g:
                te.add(g)
    return te


def main():
    union_gff, ann, ref_gff = sys.argv[1:4]
    runs = [a.split('=', 1) for a in sys.argv[4:]]
    te = te_annotated(union_gff, ann)
    gm, _ = parse_gff(union_gff)
    ref, _ = parse_gff(ref_gff, want_protein_coding=True)
    ridx = index(ref)
    matched = set()
    for gid, gv in gm.items():
        for rs, re_, rid in ridx.get((gv['chrom'], gv['strand']), []):
            if re_ < gv['start'] or rs > gv['end']:
                continue
            if ofrac(gv, ref[rid]) >= 0.5:
                matched.add(gid)
                break
    all_ids = set(gm)
    print('threshold\tkept\tkept_TE_annotated\tremoved\tremoved_ref_matching\tremoved_ref_nonTE')
    for thr, final in runs:
        kept = gene_ids(final) & all_ids
        removed = all_ids - kept
        rm_ref = removed & matched
        print(f'{thr}\t{len(kept)}\t{len(kept & te)} ({100*len(kept & te)/max(len(kept),1):.1f}%)\t'
              f'{len(removed)}\t{len(rm_ref)}\t{len(rm_ref - te)}')


if __name__ == '__main__':
    main()
