#!/usr/bin/env python3
"""check_merge_placement.py -- audit a merged annotation from the FILES.

merge_into_reference.py reports its own placement counts, and those come from
its own bookkeeping: if its placement logic is wrong, so are its counts. This
reads the reference and the merged output and asks the same questions
independently, so a delivered file can be checked by someone who did not run
the merge.

Three questions, all of them about an ADDED transcript (column 2 == --source):

  1. is it under the reference gene whose CDS it overlaps most (same strand)?
  2. does it overlap the CDS of its host gene at all?
  3. does its CDS span two or more reference genes while sitting under one?

Any of these going wrong is silent in the file: the GFF3 stays valid, nothing is
missing, and only a downstream "which gene holds this feature" lookup is wrong.
That is why the audit exists separately from the merge.

Exit status: 0 clean, 1 mis-attachments found, 2 COULD NOT ASSESS -- an empty
answer from a file with no added transcripts, or no reference genes, is not a
pass and does not report as one.

Usage:
  check_merge_placement.py --ref reference.gff3 --merged reference_updated.gff3
      [--source gene-miner] [--min-cds-overlap 30] [--list 10]
"""
import argparse, collections, re, sys

def attrs(s):
    d = {}
    for kv in s.rstrip(';').split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k.strip()] = v.strip()
    return d

def strip_prefix(i):
    return re.sub(r'^(gene|transcript|mRNA|CDS|exon):', '', i)

def read(path):
    for line in open(path, errors='replace'):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9: yield f, attrs(f[8])

def merge_iv(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return [(s, e) for s, e in out]

def iv_overlap(a, b):
    n = i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if e >= s: n += e - s + 1
        if a[i][1] < b[j][1]: i += 1
        else: j += 1
    return n

def gene_types(recs):
    """Same structural rule merge_into_reference.py uses: an ID, no Parent, and
    the Parent of something. Enumerating type names would miss the SO subtypes
    references actually use (transposable_element_gene, ncRNA_gene)."""
    has_child, any_parent = set(), False
    for f, d in recs:
        if 'Parent' in d:
            any_parent = True
            for p in d['Parent'].split(','): has_child.add(strip_prefix(p.strip()))
    t = set()
    for f, d in recs:
        if 'Parent' in d or 'ID' not in d: continue
        if not any_parent or strip_prefix(d['ID']) in has_child: t.add(f[2])
    return t | {'gene'}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True)
    ap.add_argument('--merged', required=True)
    ap.add_argument('--source', default='gene-miner')
    ap.add_argument('--min-cds-overlap', type=int, default=30)
    ap.add_argument('--list', type=int, default=5, help='how many examples to print')
    a = ap.parse_args()

    ref = list(read(a.ref))
    types = gene_types(ref)
    span, tx_gene = {}, {}
    for f, d in ref:
        if f[2] in types and 'ID' in d:
            span[strip_prefix(d['ID'])] = (f[0], f[6], int(f[3]), int(f[4]))
    for f, d in ref:
        if f[2] in types or 'Parent' not in d or 'ID' not in d: continue
        p = strip_prefix(d['Parent'].split(',')[0])
        if p in span: tx_gene[strip_prefix(d['ID'])] = p
    cds = collections.defaultdict(list)
    for f, d in ref:
        if f[2] != 'CDS' or 'Parent' not in d: continue
        g = tx_gene.get(strip_prefix(d['Parent'].split(',')[0]))
        if g: cds[g].append((int(f[3]), int(f[4])))
    cds = {g: merge_iv(v) for g, v in cds.items()}
    idx = collections.defaultdict(list)
    for g, (c, st, s, e) in span.items(): idx[(c, st)].append((s, e, g))
    for k in idx: idx[k].sort()

    tx, kid = {}, collections.defaultdict(list)
    for f, d in read(a.merged):
        if f[1] != a.source: continue
        if f[2] in ('mRNA', 'transcript') and 'ID' in d: tx[strip_prefix(d['ID'])] = (f, d)
        elif f[2] == 'CDS' and 'Parent' in d:
            kid[strip_prefix(d['Parent'].split(',')[0])].append((int(f[3]), int(f[4])))

    # COULD NOT ASSESS comes before any verdict. "0 wrong" out of nothing is the
    # permissive reading of a check that never ran, and it is the one that gets
    # believed.
    if not span:
        print(f"COULD NOT ASSESS: no gene records found in {a.ref}")
        return 2
    if not tx:
        print(f"COULD NOT ASSESS: {a.merged} has no transcripts with "
              f"'{a.source}' in column 2 -- nothing was added, or --source is wrong")
        return 2

    under_ref, own, no_cds = 0, 0, 0
    wrong, zero, multi = [], [], []
    for tid, (f, d) in sorted(tx.items()):
        host = strip_prefix(d.get('Parent', ''))
        if host not in span:
            own += 1
            continue
        under_ref += 1
        c = merge_iv(kid.get(tid, []))
        if not c:
            no_cds += 1
            continue
        lo, hi = c[0][0], c[-1][1]
        cand = []
        for s, e, g in idx[(f[0], f[6])]:
            if s > hi: break
            if e >= lo:
                n = iv_overlap(c, cds.get(g, []))
                if n >= a.min_cds_overlap: cand.append((n, g))
        cand.sort(reverse=True)
        if len(cand) >= 2: multi.append((tid, host, [g for _, g in cand]))
        if cand and cand[0][1] != host: wrong.append((tid, host, cand[0][1], cand[0][0]))
        if iv_overlap(c, cds.get(host, [])) == 0: zero.append((tid, host))

    print(f"reference genes                                 {len(span):,}")
    print(f"added transcripts under a reference gene        {under_ref:,}")
    print(f"added transcripts under their own gene record   {own:,}")
    if no_cds:
        print(f"  (of those under a reference gene, {no_cds:,} carry no CDS row "
              f"and cannot be judged)")
    print(f"NOT under their largest-CDS-overlap gene        {len(wrong):,}")
    print(f"no CDS overlap with their host at all          {len(zero):,}")
    print(f"CDS spans >=2 reference genes, filed under one  {len(multi):,}")
    for label, rows in (('wrong host', wrong), ('zero overlap', zero),
                        ('read-through under one gene', multi)):
        for r in rows[:a.list]: print(f"  {label}: {r}")

    if wrong or multi:
        print("FAIL: an added transcript is under a gene its own CDS does not "
              "place it in; a 'which gene holds this' lookup will disagree with "
              "the annotation", file=sys.stderr)
        return 1
    print("OK: every added transcript with a CDS is under the reference gene it "
          "overlaps most, and no read-through model is filed under one gene.")
    return 0

if __name__ == '__main__':
    sys.exit(main())
