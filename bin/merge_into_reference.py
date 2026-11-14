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
import argparse, sys, re, collections, subprocess, shutil

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
    ap.add_argument('--gff3validator',
                    help="path to genometools `gt`; found on PATH if omitted. "
                         "Its absence is reported, never treated as a pass.")
    a = ap.parse_args()

    # --- pairing -------------------------------------------------------------
    cls, refof = {}, {}
    with open(a.pairs) as fh:
        head = fh.readline()
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4: continue
            # Strip here too. gm_compare writes the ID verbatim, and some
            # references prefix it by feature type -- ITAG4.0 writes
            # ID=gene:Solyc00g500001.1. The reference index below is built from
            # stripped IDs, so an unstripped id from the pairs file matches
            # nothing: every attachment would be skipped, no alternative
            # transcript would be added, and NOTHING would report it. Lettuce
            # and watermelon have no prefix, so this stayed invisible on them.
            cls[strip_prefix(p[0])] = p[3]
            if p[1]: refof[strip_prefix(p[0])] = (strip_prefix(p[1]), p[2])
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
    ref_gene_line = {}        # stripped ref gene id -> index into ref_recs
    ref_gene_id_raw = {}      # stripped -> the ID EXACTLY as the file writes it
    for i, (f, d, _) in enumerate(ref_recs):
        if f[2] == 'gene':
            raw = d.get('ID', '')
            ref_gene_line[strip_prefix(raw)] = i
            ref_gene_id_raw[strip_prefix(raw)] = raw

    # our gene -> the reference gene it will be attached to
    attach = {}
    unresolved = []
    for gid, c in cls.items():
        if c == 'revised' and gid in refof:
            rid = refof[gid][0]
            if rid in ref_gene_line: attach[gid] = rid
            else: unresolved.append((gid, rid))
    if unresolved:
        # A pairing that names a reference gene absent from the reference file
        # means the two sides disagree about identifiers. Skipping quietly would
        # produce a merge that looks finished and attached nothing.
        print(f"FAIL: {len(unresolved)} pairing(s) name a reference gene not found in --ref,"
              f" e.g. {unresolved[:3]}", file=sys.stderr)
        print("      The pairs file and the reference disagree about gene IDs.", file=sys.stderr)
        sys.exit(1)
    novel = [g for g, c in cls.items() if c == 'novel' and g in gm_gene]

    # What does THIS reference mean by a length-like attribute? ITAG4.0 writes
    # length = end - start; others may write end - start + 1. Assuming either
    # one puts a wrong number on every record we touch, and no GFF3 validator
    # would say so because the key is not in the spec. So it is derived from the
    # records as they arrive, and if they do not agree the attribute is dropped
    # rather than guessed -- a missing attribute is honest, a wrong one is not.
    SPAN_KEYS = ('length', 'len', 'size')
    span_off = collections.Counter()
    for f, d, _ in ref_recs:
        if f[2] != 'gene': continue
        for k, v in d.items():
            if k.lower() in SPAN_KEYS and str(v).isdigit():
                span_off[(k.lower(), int(v) - (int(f[4]) - int(f[3])))] += 1
    span_rule = {}
    for (k, off), n in span_off.items():
        tot = sum(c for (kk, _), c in span_off.items() if kk == k)
        # one convention has to account for essentially all of them
        if n / tot >= 0.99:
            span_rule[k] = off
    for k in {kk for kk, _ in span_off}:
        if k in span_rule:
            print(f"reference span convention: {k} = end - start + {span_rule[k]} "
                  f"({span_off[(k, span_rule[k])]:,} record(s))")
        else:
            print(f"reference span convention: {k} is NOT consistent "
                  f"({len(set(o for kk, o in span_off if kk == k))} different offsets); "
                  f"it will be dropped where a gene is widened rather than guessed")

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
            # The reference may carry its own derived span attribute -- ITAG4.0
            # writes length=9193 on the gene. Widening the coordinates and
            # leaving that behind puts a contradiction inside one record, and no
            # GFF3 validator catches it: `length` is non-standard, so the format
            # is still valid while the file disagrees with itself. Measured:
            # 5,945 of 5,945 widened genes on one crop. Recompute it, and keep
            # the original beside gm_original_span so nothing is lost.
            base = int(nf[4]) - int(nf[3])
            for lk in [k for k in nd if k.lower() in SPAN_KEYS]:
                if not str(nd[lk]).isdigit(): continue
                if 'gm_original_' + lk not in nd:
                    nd['gm_original_' + lk] = nd[lk]
                if lk.lower() in span_rule:
                    nd[lk] = str(base + span_rule[lk.lower()])
                else:
                    del nd[lk]        # no single convention to honour
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
                # The ID exactly as the reference writes it, not the stripped
                # form. ITAG4.0 writes ID=gene:Solyc00g500001.1, so a stripped
                # Parent points at nothing and the added transcript is orphaned
                # by GFF3's own rule -- 32,013 of them in one delivered file,
                # which a strict parser drops.
                nd['Parent'] = ref_gene_id_raw.get(rid, rid)   # child of the REFERENCE gene
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
    # VERBATIM, deliberately. Stripping both sides made this check blind to the
    # error it exists to catch: a Parent written in the stripped form against an
    # ID written with a prefix compared equal, so 32,013 genuinely orphaned
    # records were reported as zero orphans. GFF3 resolves Parent against ID as
    # written, so the check must too.
    parents = {d['ID'] for f, d, _ in chk if 'ID' in d}
    orphan = [d['Parent'].split(',')[0] for f, d, _ in chk
              if 'Parent' in d and d['Parent'].split(',')[0] not in parents]

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

    # Attributes that restate a coordinate must agree with it. A validator
    # cannot do this -- these keys are not in the spec -- so it is done here.
    # Judged against the convention derived from the reference, not against an
    # assumed one. Checking end-start+1 on a file that writes end-start flagged
    # 28,130 untouched reference genes as defective and would have driven a
    # "fix" that put a wrong value on every record this tool writes.
    inconsistent = []
    for f, d, _ in chk:
        if f[2] != 'gene': continue
        base = int(f[4]) - int(f[3])
        for k, v in d.items():
            kl = k.lower()
            if kl in SPAN_KEYS and str(v).isdigit() and kl in span_rule:
                if int(v) != base + span_rule[kl]:
                    inconsistent.append((d.get('ID', '?'), k, v, base + span_rule[kl]))
    print(f"span attributes disagreeing with the reference's own convention  {len(inconsistent)}"
          + (f"  e.g. {inconsistent[:2]}" if inconsistent else ""))

    # The external validator, because the checks above are mine and this one is
    # not. It caught a defect mine missed: 32,013 records whose Parent was
    # written in a stripped form against a prefixed ID, which my orphan check
    # normalised away and reported as zero.
    gt_state, gt_msg = 'not run', ''
    gt = a.gff3validator or shutil.which('gt')
    if gt:
        try:
            r = subprocess.run([gt, 'gff3validator', a.out],
                               capture_output=True, text=True, timeout=1800)
            err = (r.stderr or '') + (r.stdout or '')
            gt_state = 'valid' if r.returncode == 0 else 'INVALID'
            gt_msg = next((l for l in err.splitlines() if 'error' in l.lower()), '')
        except Exception as e:
            gt_state, gt_msg = 'could not run', str(e)
    else:
        gt_msg = "no `gt` on PATH; pass --gff3validator /path/to/gt to enable"
    print(f"gt gff3validator: {gt_state}" + (f"  {gt_msg}" if gt_msg else ""))

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
                         ('orphan_parent_references', len(orphan)),
                         ('span_attributes_inconsistent', len(inconsistent)),
                         ('gt_gff3validator', gt_state)):
                fh.write(f"{k}\t{v}\n")
            fh.write("\n# reference gene lines widened to span an added transcript\n")
            fh.write("ref_gene\tgm_locus\toriginal_span\tnew_span\n")
            for r in widened: fh.write('\t'.join(r) + '\n')

    bad = 0
    if gt_state == 'INVALID':
        print(f"FAIL: gt gff3validator rejected {a.out}: {gt_msg}", file=sys.stderr); bad = 1
    elif gt_state in ('not run', 'could not run'):
        # Not a pass. An absent validator and a satisfied one must not read the
        # same, which is the failure mode this whole file keeps running into.
        print(f"NOTE: the GFF3 validator did not run ({gt_msg}). The checks above "
              f"are self-written and one of them has been wrong before.", file=sys.stderr)
    if inconsistent:
        print(f"FAIL: {len(inconsistent)} gene record(s) carry a span attribute that "
              f"contradicts their coordinates", file=sys.stderr); bad = 1
    if ref_genes_out != ref_genes_in:
        print(f"FAIL: {ref_genes_in} reference genes in, {ref_genes_out} out -- the reference was not preserved", file=sys.stderr); bad = 1
    if dup:
        print(f"FAIL: {len(dup)} duplicate ID(s) in the output", file=sys.stderr); bad = 1
    if orphan:
        print(f"FAIL: {len(orphan)} record(s) name a Parent that is not in the file", file=sys.stderr); bad = 1
    sys.exit(bad)

if __name__ == '__main__':
    main()
