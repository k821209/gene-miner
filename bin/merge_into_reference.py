#!/usr/bin/env python3
"""merge_into_reference.py -- deliver an UPDATED reference annotation rather than
a parallel catalogue.

The reference GFF3 is the base and is preserved: every reference line is written
back unchanged, keeping its IDs, its source, its attributes and its curation.
Onto it go two things, both tagged so they can be filtered out in one pass:

  * for a reference locus where our model has a DIFFERENT CDS structure, our
    model is added as an ADDITIONAL TRANSCRIPT of that reference gene. The
    reference transcript stays. Nothing curated is overwritten; a reviewer sees
    both side by side and decides.
  * for a locus of ours with no reference counterpart, a NEW GENE is added.

Identical-structure matches add nothing -- that would duplicate the reference.
Reference loci we did not recover are left exactly as they are.

Everything added carries `gene-miner` in column 2 (the source field), so
`awk -F'\t' '$2=="gene-miner"'` separates the additions from the reference, and
`awk -F'\t' '$2!="gene-miner"'` reconstructs the input.

The one place this modifies a reference line: GFF3 requires a parent to span its
children, so attaching a transcript that extends past the reference gene's
boundary forces that gene's start/end to widen. Those are counted and reported
rather than left for a reader to discover, and the original span is kept in a
`gm_original_span` attribute so the change is reversible.

Usage:
  merge_into_reference.py --ref reference.gff3 --gm union.final.gff3 \
      --pairs pairs.tsv --out reference_updated.gff3 [--report merge_report.tsv]
      [--tiers union.final.tiers.gff3]
"""
import argparse, sys, re, collections

def attrs(s):
    d = {}
    for kv in s.rstrip(';').split(';'):
        if not kv or '=' not in kv: continue
        k, v = kv.split('=', 1)
        d[k.strip()] = v.strip()
    return d

def fmt_attrs(d):
    # ';' and '=' are GFF3's attribute separators. A value carrying either
    # silently becomes extra attributes: a note of ours did exactly that and
    # left a valueless key in all 53,432 records it touched. Percent-encode here
    # rather than trusting every call site to remember.
    def esc(v):
        return (str(v).replace('%', '%25').replace(';', '%3B')
                      .replace('=', '%3D').replace('\t', '%09'))
    return ';'.join(f"{k}={esc(v)}" for k, v in d.items() if v is not None and v != '')

def strip_prefix(i):
    # Reference files vary: ITAG4.0 writes ID=gene:Solyc00g500001.1 while the
    # others write the bare id. gm_compare strips these, so this must match or
    # no pair would resolve.
    return re.sub(r'^(gene|transcript|mRNA|CDS|exon):', '', i)

def read_gff(path):
    """-> (header_lines, records) where a record is (fields, attrs, lineno)."""
    header, recs = [], []
    with open(path, errors='replace') as fh:
        for n, line in enumerate(fh):
            if line.startswith('#'):
                if not recs: header.append(line.rstrip('\n'))
                continue
            if not line.strip(): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9: continue
            recs.append((f, attrs(f[8]), n))
    return header, recs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True)
    ap.add_argument('--gm', required=True)
    ap.add_argument('--pairs', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--report')
    ap.add_argument('--tiers', help='union.final.tiers.gff3, to carry the evidence tier onto each addition')
    ap.add_argument('--source', default='gene-miner')
    a = ap.parse_args()

    # --- pairing -------------------------------------------------------------
    cls, refof = {}, {}
    with open(a.pairs) as fh:
        head = fh.readline()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4: continue
            cls[p[0]] = p[3]
            if p[1]: refof[p[0]] = (p[1], p[2])
    if not cls:
        sys.exit("FATAL: --pairs has no rows; nothing could be merged")

    # --- evidence tier per gene (optional) -----------------------------------
    tier = {}
    if a.tiers:
        _, trecs = read_gff(a.tiers)
        for f, d, _ in trecs:
            if f[2] == 'gene':
                gid = strip_prefix(d.get('ID', ''))
                for k in ('tier', 'confidence', 'gm_tier'):
                    if k in d: tier[gid] = d[k]; break

    # --- our annotation, grouped by gene -------------------------------------
    _, grecs = read_gff(a.gm)
    gm_children = collections.defaultdict(list)   # gene id -> non-gene records
    gm_gene = {}
    tx_parent = {}                                 # transcript id -> gene id
    for f, d, _ in grecs:
        if f[2] == 'gene':
            gm_gene[strip_prefix(d['ID'])] = (f, d)
        elif 'Parent' in d:
            par = strip_prefix(d['Parent'].split(',')[0])
            if f[2] in ('mRNA', 'transcript'):
                tx_parent[strip_prefix(d.get('ID', ''))] = par
                gm_children[par].append((f, d))
            else:
                gm_children[tx_parent.get(par, par)].append((f, d))

    # --- reference -----------------------------------------------------------
    ref_header, ref_recs = read_gff(a.ref)
    ref_gene_line = {}        # ref gene id -> index into ref_recs
    for i, (f, d, _) in enumerate(ref_recs):
        if f[2] == 'gene':
            ref_gene_line[strip_prefix(d.get('ID', ''))] = i

    # our gene -> the reference gene it will be attached to
    attach = {}
    for gid, c in cls.items():
        if c == 'revised' and gid in refof:
            rid = refof[gid][0]
            if rid in ref_gene_line: attach[gid] = rid
    novel = [g for g, c in cls.items() if c == 'novel' and g in gm_gene]

    # --- widen reference genes that would not span an attached transcript ----
    widened = []
    for gid, rid in attach.items():
        kids = gm_children.get(gid, [])
        if not kids: continue
        lo = min(int(f[3]) for f, _ in kids); hi = max(int(f[4]) for f, _ in kids)
        i = ref_gene_line[rid]
        f, d, n = ref_recs[i]
        rs, re_ = int(f[3]), int(f[4])
        if lo < rs or hi > re_:
            nf = list(f); nd = dict(d)
            # A reference gene can receive transcripts from MORE THAN ONE of our
            # loci (measured: 842 of our loci shared a best-match reference
            # locus on one crop). On the second pass the record being read has
            # already been widened, so writing gm_original_span again would
            # record the first widening's result and lose the true original --
            # which is the whole point of keeping it. Set it once; accumulate
            # the attribution.
            if 'gm_original_span' not in nd:
                nd['gm_original_span'] = f"{rs}-{re_}"
            prev = nd.get('gm_widened_for')
            nd['gm_widened_for'] = f"{prev},{gid}" if prev else gid
            nf[3] = str(min(rs, lo)); nf[4] = str(max(re_, hi))
            nf[8] = fmt_attrs(nd)
            ref_recs[i] = (nf, nd, n)
            widened.append((rid, gid, nd['gm_original_span'], f"{nf[3]}-{nf[4]}"))

    # --- emit ----------------------------------------------------------------
    added_tx = added_genes = 0
    out_extra = collections.defaultdict(list)     # ref gene id -> lines to place after it
    for gid, rid in sorted(attach.items()):
        ov = refof[gid][1]
        for f, d in gm_children.get(gid, []):
            nf = list(f); nd = dict(d)
            nf[1] = a.source
            if f[2] in ('mRNA', 'transcript'):
                nd['Parent'] = rid                    # child of the REFERENCE gene
                nd['gm_class'] = 'revised_structure_candidate'
                nd['gm_overlap'] = ov
                nd['gm_source_locus'] = gid
                if gid in tier: nd['gm_tier'] = tier[gid]
                nd['gm_note'] = 'alternative_structure_proposed_by_gene-miner,reference_transcript_retained'
                added_tx += 1
            nf[8] = fmt_attrs(nd)
            out_extra[rid].append('\t'.join(nf))

    novel_lines = []
    for gid in sorted(novel):
        f, d = gm_gene[gid]
        nf = list(f); nd = dict(d)
        nf[1] = a.source
        nd['gm_class'] = 'novel_locus'
        if gid in tier: nd['gm_tier'] = tier[gid]
        nd['gm_note'] = 'absent_from_reference,candidate_for_expert_review'
        nf[8] = fmt_attrs(nd)
        novel_lines.append('\t'.join(nf))
        added_genes += 1
        for cf, cd in gm_children.get(gid, []):
            xf = list(cf); xd = dict(cd); xf[1] = a.source; xf[8] = fmt_attrs(xd)
            novel_lines.append('\t'.join(xf))

    with open(a.out, 'w') as out:
        out.write('\n'.join(ref_header) + '\n' if ref_header else '##gff-version 3\n')
        out.write(f"# merged by merge_into_reference.py: reference preserved; "
                  f"{added_tx} alternative transcript(s) added to {len(attach)} reference gene(s), "
                  f"{added_genes} novel gene(s) appended; "
                  f"{len({r[0] for r in widened})} reference gene line(s) widened\n")
        out.write(f"# additions carry '{a.source}' in column 2; "
                  f"awk -F'\\t' '$2!=\"{a.source}\"' reconstructs the input\n")
        cur_gene = None
        for f, d, _ in ref_recs:
            if f[2] == 'gene':
                if cur_gene is not None:
                    for l in out_extra.get(cur_gene, []): out.write(l + '\n')
                cur_gene = strip_prefix(d.get('ID', ''))
            out.write('\t'.join(f) + '\n')
        if cur_gene is not None:
            for l in out_extra.get(cur_gene, []): out.write(l + '\n')
        for l in novel_lines: out.write(l + '\n')

    # --- checks that would otherwise be left to a reader --------------------
    ref_genes_in = len(ref_gene_line)
    _, chk = read_gff(a.out)
    genes_out = sum(1 for f, _, _ in chk if f[2] == 'gene')
    ref_genes_out = sum(1 for f, _, _ in chk if f[2] == 'gene' and f[1] != a.source)
    # Only gene and mRNA IDs must be unique. A multi-segment CDS or exon shares
    # one ID across its parts -- that is GFF3's discontinuous-feature form, not
    # an error (measured on one input: 320,833 CDS records carrying 68,451
    # distinct IDs, i.e. one per transcript). Counting every ID reported 46,320
    # "duplicates" on a merge whose gene and mRNA IDs were in fact all unique.
    ids = [strip_prefix(d['ID']) for f, d, _ in chk
           if 'ID' in d and f[2] in ('gene', 'mRNA', 'transcript')]
    dup = [k for k, v in collections.Counter(ids).items() if v > 1]
    parents = {strip_prefix(d['ID']) for f, d, _ in chk if 'ID' in d}
    orphan = [strip_prefix(d['Parent'].split(',')[0]) for f, d, _ in chk
              if 'Parent' in d and strip_prefix(d['Parent'].split(',')[0]) not in parents]

    print(f"reference genes in            {ref_genes_in}")
    print(f"reference genes out           {ref_genes_out}   (must equal the line above)")
    print(f"alternative transcripts added {added_tx}  to {len(attach)} reference gene(s)")
    print(f"novel genes appended          {added_genes}")
    widened_genes = len({r[0] for r in widened})
    print(f"reference gene lines widened  {widened_genes}"
          + (f"   ({len(widened)} widening events; {len(widened)-widened_genes} gene(s) widened twice)"
             if len(widened) != widened_genes else ""))
    print(f"total genes out               {genes_out}")
    print(f"duplicate IDs                 {len(dup)}" + (f"  e.g. {dup[:3]}" if dup else ""))
    print(f"orphan Parent references      {len(orphan)}" + (f"  e.g. {orphan[:3]}" if orphan else ""))

    if a.report:
        with open(a.report, 'w') as fh:
            fh.write("metric\tvalue\n")
            for k, v in (('reference_genes_in', ref_genes_in),
                         ('reference_genes_out', ref_genes_out),
                         ('alternative_transcripts_added', added_tx),
                         ('reference_genes_given_an_alternative', len(attach)),
                         ('novel_genes_appended', added_genes),
                         ('reference_gene_lines_widened', len({r[0] for r in widened})),
                         ('widening_events', len(widened)),
                         ('total_genes_out', genes_out),
                         ('duplicate_ids', len(dup)),
                         ('orphan_parent_references', len(orphan))):
                fh.write(f"{k}\t{v}\n")
            fh.write("\n# reference gene lines widened to span an added transcript\n")
            fh.write("ref_gene\tgm_locus\toriginal_span\tnew_span\n")
            for r in widened: fh.write('\t'.join(r) + '\n')

    bad = 0
    if ref_genes_out != ref_genes_in:
        print(f"FAIL: {ref_genes_in} reference genes in, {ref_genes_out} out -- the reference was not preserved", file=sys.stderr); bad = 1
    if dup:
        print(f"FAIL: {len(dup)} duplicate ID(s) in the output", file=sys.stderr); bad = 1
    if orphan:
        print(f"FAIL: {len(orphan)} record(s) name a Parent that is not in the file", file=sys.stderr); bad = 1
    sys.exit(bad)

if __name__ == '__main__':
    main()
