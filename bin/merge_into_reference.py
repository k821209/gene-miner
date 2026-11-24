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

WHICH reference gene a transcript goes under is decided per TRANSCRIPT, on its
own CDS overlap -- not per locus. gm_compare.py pairs a whole locus with one
best reference gene, and a locus can hold a read-through model, or isoforms
lying in different neighbouring genes; attaching all of them to the locus's
pairing put transcripts under genes they do not overlap. A transcript that
belongs under no single reference gene gets a gene record of its own, labelled
with why:

    gm_class=readthrough      its CDS overlaps two or more reference genes
    gm_class=no_overlap       its CDS overlaps no reference gene
    gm_class=exceeds_span     it would have required widening its host gene
                              into a neighbouring reference gene

Identical-structure matches add nothing -- that would duplicate the reference.
Reference loci we did not recover are left exactly as they are.

Everything added carries `gene-miner` in column 2 (the source field), so
`awk -F'\t' '$2=="gene-miner"'` separates the additions from the reference, and
`awk -F'\t' '$2!="gene-miner"'` reconstructs the input.

The one place this modifies a reference line: GFF3 requires a parent to span its
children, so attaching a transcript that extends past the reference gene's
boundary forces that gene's start/end to widen. Those are counted and reported
rather than left for a reader to discover, and the original span is kept in a
`gm_original_span` attribute so the change is reversible. A widening that would
reach into another reference gene on the same strand is REFUSED -- a gene that
grows over its neighbour then contains the neighbour's features, and every
"which gene holds this?" lookup downstream answers wrongly with no error
anywhere.

Usage:
  merge_into_reference.py --ref reference.gff3 --gm union.final.gff3 \
      --pairs pairs.tsv --out reference_updated.gff3 [--report merge_report.tsv]
      [--tiers union.final.tiers.gff3]
"""
import argparse, sys, re, collections, subprocess, shutil, os

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

def merge_iv(ivs):
    """Sorted, non-overlapping [(start, end)] from possibly overlapping ones.
    CDS blocks of several isoforms are compared as one footprint; counting them
    raw would double-count a shared exon."""
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return [(s, e) for s, e in out]

def iv_overlap(a, b):
    """Bases shared by two merged interval lists."""
    n = i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0]); e = min(a[i][1], b[j][1])
        if e >= s: n += e - s + 1
        if a[i][1] < b[j][1]: i += 1
        else: j += 1
    return n

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

def gene_record_types(recs):
    """Which feature types open a gene record IN THIS FILE -- derived from the
    file, not enumerated.

    References do not agree on how to type a gene. Phytozome C. reinhardtii
    CC-4532 v6.1 writes 810 of its 17,693 genes as `transposable_element_gene`;
    Ensembl writes `ncRNA_gene` and `pseudogene`. gm_compare.py pairs those
    correctly -- it keys a gene through its mRNAs' Parent, not through the gene
    row's type -- so an index built from `gene` rows alone rejected their
    pairings with "pairing(s) name a reference gene not found in --ref" on a
    reference that was perfectly valid. The message named the pairs file for a
    defect that was in this index.

    A list of type names would have to be extended for the next reference that
    uses one we have not seen, and would report a clean run until someone
    noticed. So the rule is structural: a gene record has an ID, has no Parent,
    and is the Parent of something.

    Returns (types, top_level_counts, gene_record_counts). The counts are
    printed by the caller, so what this included and what it left out is in the
    log -- rather than surfacing later as a gene count that does not add up.
    """
    has_child, any_parent = set(), False
    for f, d, _ in recs:
        if 'Parent' in d:
            any_parent = True
            for p in d['Parent'].split(','):
                has_child.add(strip_prefix(p.strip()))
    top, genes = collections.Counter(), collections.Counter()
    for f, d, _ in recs:
        if 'Parent' in d or 'ID' not in d: continue
        top[f[2]] += 1
        # A file with no Parent anywhere is a flat list of genes -- which is
        # what annotate_tiers.py writes. There every top-level record is one,
        # and requiring a child would find nothing.
        if not any_parent or strip_prefix(d['ID']) in has_child:
            genes[f[2]] += 1
    return set(genes) | {'gene'}, top, genes

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
        tier_types, _, _ = gene_record_types(trecs)
        for f, d, _ in trecs:
            if f[2] in tier_types:
                gid = strip_prefix(d.get('ID', ''))
                for k in ('tier', 'confidence', 'gm_tier'):
                    if k in d: tier[gid] = d[k]; break

    # --- our annotation, grouped by gene -------------------------------------
    _, grecs = read_gff(a.gm)
    gm_types, _, _ = gene_record_types(grecs)
    gm_children = collections.defaultdict(list)   # gene id -> non-gene records
    gm_gene = {}
    tx_parent = {}                                 # transcript id -> gene id
    for f, d, _ in grecs:
        if f[2] in gm_types:
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
    ref_types, ref_top, ref_gene_counts = gene_record_types(ref_recs)
    for t, n in ref_top.most_common():
        print(f"reference top-level {t:<26} {n:>8,}   "
              + (f"({ref_gene_counts[t]:,} with children -> gene record)"
                 if t in ref_gene_counts else "(no children -> not a gene record)"))
    ref_gene_line = {}        # stripped ref gene id -> index into ref_recs
    ref_gene_id_raw = {}      # stripped -> the ID EXACTLY as the file writes it
    for i, (f, d, _) in enumerate(ref_recs):
        if f[2] in ref_types:
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
        if f[2] not in ref_types: continue
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

    # --- where does each ADDED TRANSCRIPT belong? ---------------------------
    #
    # gm_compare.py pairs at the LOCUS level: one best reference gene for the
    # whole locus. Attaching every transcript of the locus to that gene, and
    # widening the gene to span them, puts a transcript under a gene it may
    # barely touch. Measured on a C. reinhardtii CC-4532 v6.1 build: 595 of
    # 10,156 added transcripts sat under a gene other than the one their CDS
    # overlaps most and 320 of those overlapped the CDS of their host gene not
    # at all; 2,272 reference gene spans were widened, some of them over a
    # neighbouring gene, so that the widened gene then contained every intron of
    # its neighbour. Nothing failed -- a downstream tool asking "which gene
    # contains this feature" simply got the wrong answer, for 10% of the introns
    # in one junction analysis.
    #
    # So each transcript is placed on its OWN CDS overlap, and a transcript that
    # does not belong under one reference gene is not forced under one.
    MIN_CDS_OV = 30        # bp; below this an overlap is not evidence of identity

    ref_tx_gene = {}       # reference transcript id -> reference gene id
    for f, d, _ in ref_recs:
        if f[2] in ref_types or 'Parent' not in d or 'ID' not in d: continue
        par = strip_prefix(d['Parent'].split(',')[0])
        if par in ref_gene_line:
            ref_tx_gene[strip_prefix(d['ID'])] = par
    ref_cds = collections.defaultdict(list)        # ref gene id -> [(s,e)]
    for f, d, _ in ref_recs:
        if f[2] != 'CDS' or 'Parent' not in d: continue
        g = ref_tx_gene.get(strip_prefix(d['Parent'].split(',')[0]))
        if g: ref_cds[g].append((int(f[3]), int(f[4])))
    ref_cds = {g: merge_iv(v) for g, v in ref_cds.items()}

    # Gene spans as the reference SUPPLIED them, and an index per sequence and
    # strand. Opposite-strand overlap is ordinary genome structure and the
    # reference already contains it, so it is not treated as a collision here;
    # the failure measured above is between same-strand neighbours.
    ref_span, by_strand = {}, collections.defaultdict(list)
    for rid, i in ref_gene_line.items():
        f = ref_recs[i][0]
        ref_span[rid] = (f[0], f[6], int(f[3]), int(f[4]))
        by_strand[(f[0], f[6])].append((int(f[3]), int(f[4]), rid))
    for k in by_strand: by_strand[k].sort()
    cur_span = {rid: (s, e) for rid, (_, _, s, e) in ref_span.items()}

    def genes_touching(seqid, strand, lo, hi, exclude=None):
        """Reference genes on this sequence and strand whose SUPPLIED span
        overlaps [lo, hi]."""
        out = []
        for s, e, rid in by_strand.get((seqid, strand), ()):
            if s > hi: break
            if e >= lo and rid != exclude: out.append(rid)
        return out

    # Regroup our records per transcript: gm_children holds a locus's
    # transcripts and their CDS/exon rows in one list, and a transcript's CDS is
    # what places it.
    tx_of = collections.defaultdict(list)          # locus -> [(tx_f, tx_d)]
    kids_of = collections.defaultdict(list)        # tx id -> [(f, d)]
    for gid in list(gm_children):
        for f, d in gm_children[gid]:
            if f[2] in ('mRNA', 'transcript'):
                tx_of[gid].append((f, d))
            elif 'Parent' in d:
                kids_of[strip_prefix(d['Parent'].split(',')[0])].append((f, d))

    place = {}             # tx id -> (kind, [ref gene ids], cds_overlap_bp, by_span)
    stats = collections.Counter()
    for gid, rid in attach.items():
        for tf, td in tx_of.get(gid, []):
            tid = strip_prefix(td.get('ID', ''))
            cds = merge_iv([(int(f[3]), int(f[4])) for f, _ in kids_of.get(tid, [])
                            if f[2] == 'CDS'])
            lo = min((s for s, _ in cds), default=int(tf[3]))
            hi = max((e for _, e in cds), default=int(tf[4]))
            ov = []
            for cand in genes_touching(tf[0], tf[6], lo, hi):
                n = iv_overlap(cds, ref_cds.get(cand, []))
                if n >= MIN_CDS_OV: ov.append((n, cand))
            ov.sort(reverse=True)
            if len(ov) >= 2:
                # A fusion / read-through model. It is not a transcript OF
                # either gene, and putting it under one of them is the error
                # this placement exists to stop.
                place[tid] = ('readthrough', [c for _, c in ov], ov[0][0], False)
            elif len(ov) == 1:
                place[tid] = ('gene', [ov[0][1]], ov[0][0], False)
                if ov[0][1] != rid: stats['reassigned'] += 1
            else:
                # No CDS overlap with any reference gene. It can still sit
                # inside one -- an intron, a UTR -- and that is a placement; if
                # it sits inside none, the locus's pairing came from somewhere
                # its CDS does not reach, and the transcript gets its own
                # record rather than a host picked for it.
                inside = genes_touching(tf[0], tf[6], int(tf[3]), int(tf[4]))
                if len(inside) == 1:
                    # The real overlap, not 0. It is below MIN_CDS_OV -- which
                    # is why this branch runs -- but 17 transcripts here had
                    # between 1 and 29 bp, and writing gm_cds_overlap_bp=0 on
                    # those made the record disagree with the file it describes.
                    place[tid] = ('gene', inside,
                                  iv_overlap(cds, ref_cds.get(inside[0], [])), True)
                    if inside[0] != rid: stats['reassigned'] += 1
                elif len(inside) >= 2:
                    place[tid] = ('readthrough', inside, 0, True)
                else:
                    place[tid] = ('no_overlap', [rid], 0, True)
            if place[tid][0] == 'gene' and place[tid][1][0] != rid:
                if iv_overlap(cds, ref_cds.get(rid, [])) == 0:
                    stats['reassigned_from_zero_overlap'] += 1

    # --- widen a host gene, but never into a neighbour ----------------------
    # Widening is what made the mis-attachment destructive: the gene grew over
    # the gene next door and then contained its features. A widening that would
    # enter another reference gene's span is refused and reported, and the
    # transcript that needed it is written as its own record instead -- GFF3
    # requires a parent to span its children, so the alternative is a file that
    # lies about one or the other.
    want = collections.defaultdict(list)           # ref gene -> [(lo, hi, tx id, locus)]
    for gid, rid in attach.items():
        for tf, td in tx_of.get(gid, []):
            tid = strip_prefix(td.get('ID', ''))
            k, hosts = place.get(tid, (None, [], 0, False))[:2]
            if k != 'gene': continue
            lo, hi = int(tf[3]), int(tf[4])
            for f, _ in kids_of.get(tid, []):
                lo = min(lo, int(f[3])); hi = max(hi, int(f[4]))
            want[hosts[0]].append((lo, hi, tid, gid))

    widened, refused = [], []
    for rid in sorted(want, key=lambda r: (ref_span[r][0], ref_span[r][2])):
        i = ref_gene_line[rid]
        seqid, strand, _, _ = ref_span[rid]
        for lo, hi, tid, gid in sorted(want[rid]):
            rs, re_ = cur_span[rid]
            if lo >= rs and hi <= re_: continue
            nlo, nhi = min(rs, lo), max(re_, hi)
            clash = [c for c in genes_touching(seqid, strand, nlo, nhi, exclude=rid)
                     if cur_span[c][0] <= nhi and cur_span[c][1] >= nlo]
            if clash:
                refused.append((rid, tid, gid, f"{rs}-{re_}", f"{nlo}-{nhi}", ','.join(sorted(clash))))
                place[tid] = ('exceeds_span', [rid], place[tid][2], place[tid][3])
                continue
            f, d, n = ref_recs[i]
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
            nf[3], nf[4] = str(nlo), str(nhi)
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
            cur_span[rid] = (nlo, nhi)
            widened.append((rid, gid, nd['gm_original_span'], f"{nf[3]}-{nf[4]}"))

    # Counted off the final placements, not accumulated while they were being
    # decided: a transcript placed under a gene and then refused its widening
    # would otherwise be counted in both.
    for k, _, _, by_span in place.values():
        # A transcript placed inside a gene's span without reaching the CDS
        # threshold -- intronic, UTR-only, or a few bases of coincidence -- is a
        # different claim from one placed on shared coding sequence, and a
        # reader checking "is every addition under its largest-CDS-overlap gene"
        # will find these and need the answer.
        stats['placed_' + (k if not (k == 'gene' and by_span) else 'gene_by_span')] += 1

    # --- emit ----------------------------------------------------------------
    added_tx = added_genes = 0
    out_extra = collections.defaultdict(list)     # ref gene id -> lines to place after it
    split_lines = []                              # additions that get their own gene
    split_groups = collections.OrderedDict()       # (locus, kind, hosts) -> [tx ids]
    for gid, rid in sorted(attach.items()):
        ov = refof[gid][1]
        for tf, td in tx_of.get(gid, []):
            tid = strip_prefix(td.get('ID', ''))
            kind, hosts, cds_bp, _ = place.get(tid, ('no_overlap', [rid], 0, True))
            if kind != 'gene':
                split_groups.setdefault((gid, kind, tuple(hosts)), []).append(tid)
                continue
            host = hosts[0]
            for f, d in [(tf, td)] + kids_of.get(tid, []):
                nf = list(f); nd = dict(d)
                nf[1] = a.source
                if f[2] in ('mRNA', 'transcript'):
                    # The ID exactly as the reference writes it, not the stripped
                    # form. ITAG4.0 writes ID=gene:Solyc00g500001.1, so a stripped
                    # Parent points at nothing and the added transcript is orphaned
                    # by GFF3's own rule -- 32,013 of them in one delivered file,
                    # which a strict parser drops.
                    nd['Parent'] = ref_gene_id_raw.get(host, host)  # child of the REFERENCE gene
                    nd['gm_class'] = 'revised_structure_candidate'
                    nd['gm_overlap'] = ov
                    nd['gm_source_locus'] = gid
                    nd['gm_cds_overlap_bp'] = cds_bp
                    if host != rid:
                        # The locus paired with one gene and this transcript
                        # went to another. Say so on the record, or the merge
                        # report is the only place it exists.
                        nd['gm_locus_paired_with'] = ref_gene_id_raw.get(rid, rid)
                    if gid in tier: nd['gm_tier'] = tier[gid]
                    nd['gm_note'] = 'alternative_structure_proposed_by_gene-miner,reference_transcript_retained'
                    added_tx += 1
                nf[8] = fmt_attrs(nd)
                out_extra[host].append('\t'.join(nf))

    # Additions that belong under no single reference gene get a gene record of
    # their own, grouped by the locus and the reference genes involved -- so a
    # locus whose isoforms disagree is split rather than averaged.
    NOTE = {
        'readthrough': 'cds_overlaps_two_or_more_reference_genes,'
                       'not_attached_to_either,candidate_for_expert_review',
        'no_overlap': 'cds_overlaps_no_reference_gene,'
                      'locus_pairing_retained_in_gm_locus_paired_with',
        'exceeds_span': 'would_require_widening_the_reference_gene_over_a_neighbour,'
                        'reference_span_left_unchanged',
    }
    tag = {'readthrough': 'rt', 'no_overlap': 'nx', 'exceeds_span': 'xs'}
    seen = collections.Counter()
    for (gid, kind, hosts), tids in split_groups.items():
        seen[gid] += 1
        new_id = f"{gid}_{tag[kind]}{seen[gid]}"
        rows = [(f, d) for t in tids for f, d in [next(x for x in tx_of[gid]
                 if strip_prefix(x[1].get('ID', '')) == t)] + kids_of.get(t, [])]
        if gid in gm_gene:
            gf = list(gm_gene[gid][0]); gd = dict(gm_gene[gid][1])
        else:
            # A locus whose own gene row is absent from --gm. Build one from the
            # transcripts, rather than failing on a file we can still describe.
            gf = list(rows[0][0]); gd = {}
            gf[2] = 'gene'
        gf[1] = a.source
        gf[3] = str(min(int(f[3]) for f, _ in rows))
        gf[4] = str(max(int(f[4]) for f, _ in rows))
        gd['ID'] = new_id
        gd['gm_class'] = kind
        gd['gm_source_locus'] = gid
        gd[{'readthrough': 'gm_readthrough_of',
            'no_overlap': 'gm_locus_paired_with',
            'exceeds_span': 'gm_intended_host'}[kind]] = \
            ','.join(ref_gene_id_raw.get(h, h) for h in hosts)
        if gid in tier: gd['gm_tier'] = tier[gid]
        gd['gm_note'] = NOTE[kind]
        for k in ('gm_original_span', 'gm_widened_for'): gd.pop(k, None)
        gf[8] = fmt_attrs(gd)
        split_lines.append('\t'.join(gf))
        added_genes += 1
        for f, d in rows:
            nf = list(f); nd = dict(d)
            nf[1] = a.source
            if f[2] in ('mRNA', 'transcript'):
                nd['Parent'] = new_id
                nd['gm_class'] = kind
                nd['gm_source_locus'] = gid
                if gid in tier: nd.setdefault('gm_tier', tier[gid])
                added_tx += 1
            nf[8] = fmt_attrs(nd)
            split_lines.append('\t'.join(nf))

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
            xf = list(cf); xd = dict(cd); xf[1] = a.source
            # Repeat the gene's tier on its transcripts. An added transcript
            # (above) carries gm_tier on the mRNA row because the transcript is
            # what is being proposed; a novel locus carries it on the gene row
            # because the locus is. Without this line the two kinds of addition
            # answer `$3=="mRNA" && $9 !~ /gm_tier=low/` differently -- the
            # obvious one-pass filter, which then keeps every transcript of
            # every low-tier novel locus and reports no error. Measured on one
            # delivered annotation: 94,904 transcripts kept where 89,806 was
            # meant, 5,098 of them from loci with no evidence at all.
            if xf[2] in ('mRNA', 'transcript') and gid in tier:
                xd.setdefault('gm_tier', tier[gid])
            xf[8] = fmt_attrs(xd)
            novel_lines.append('\t'.join(xf))

    with open(a.out, 'w') as out:
        out.write('\n'.join(ref_header) + '\n' if ref_header else '##gff-version 3\n')
        out.write(f"# merged by merge_into_reference.py: reference preserved; "
                  f"{added_tx} transcript(s) added, {added_genes} gene record(s) appended; "
                  f"{len({r[0] for r in widened})} reference gene line(s) widened, "
                  f"{len(refused)} widening(s) refused to protect a neighbour\n")
        out.write(f"# additions carry '{a.source}' in column 2; "
                  f"awk -F'\\t' '$2!=\"{a.source}\"' reconstructs the input\n")
        cur_gene = None
        for f, d, _ in ref_recs:
            if f[2] in ref_types:
                if cur_gene is not None:
                    for l in out_extra.get(cur_gene, []): out.write(l + '\n')
                cur_gene = strip_prefix(d.get('ID', ''))
            out.write('\t'.join(f) + '\n')
        if cur_gene is not None:
            for l in out_extra.get(cur_gene, []): out.write(l + '\n')
        for l in novel_lines: out.write(l + '\n')
        for l in split_lines: out.write(l + '\n')

    # --- checks that would otherwise be left to a reader --------------------
    ref_genes_in = len(ref_gene_line)
    _, chk = read_gff(a.out)
    # The output holds both sides, so it is counted against both sides' types.
    # Counting it by a different rule than `ref_genes_in` was built with would
    # make the "must equal the line above" check compare two different things.
    out_types = ref_types | gm_types
    genes_out = sum(1 for f, _, _ in chk if f[2] in out_types)
    ref_genes_out = sum(1 for f, _, _ in chk if f[2] in ref_types and f[1] != a.source)
    # Only gene and mRNA IDs must be unique. A multi-segment CDS or exon shares
    # one ID across its parts -- that is GFF3's discontinuous-feature form, not
    # an error (measured on one input: 320,833 CDS records carrying 68,451
    # distinct IDs, i.e. one per transcript). Counting every ID reported 46,320
    # "duplicates" on a merge whose gene and mRNA IDs were in fact all unique.
    ids = [strip_prefix(d['ID']) for f, d, _ in chk
           if 'ID' in d and (f[2] in out_types or f[2] in ('mRNA', 'transcript'))]
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
    print(f"transcripts added             {added_tx}   "
          f"({stats['placed_gene'] + stats['placed_gene_by_span']} under "
          f"{len({h[0] for k, h, _, _ in place.values() if k == 'gene'})} reference gene(s), "
          f"the rest under gene records of their own)")
    print(f"gene records appended         {added_genes}"
          + (f"   ({len(split_groups)} of them for additions that belong under no single "
             f"reference gene)" if split_groups else ""))
    print(f"added transcripts placed      "
          + ', '.join(f"{k.split('placed_')[1]}={v:,}" for k, v in sorted(stats.items())
                      if k.startswith('placed_')))
    print(f"  moved off the locus's own pairing   {stats['reassigned']:,}"
          f"   (of these, no CDS overlap with the paired gene at all: "
          f"{stats['reassigned_from_zero_overlap']:,})")
    widened_genes = len({r[0] for r in widened})
    print(f"widenings refused (neighbour) {len(refused)}"
          + (f"  e.g. {refused[0][0]} would have reached {refused[0][4]} over "
             f"{refused[0][5]}" if refused else ""))
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
        if f[2] not in out_types: continue
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
    #
    # The REFERENCE is validated too, and that is the point. Three of four
    # reference annotations in hand are not valid GFF3 as supplied -- two lack
    # the ##gff-version line, one declares every sequence-region twice, one
    # writes `.` as the phase on 98,453 CDS records. A merged file that inherits
    # such a defect is not the same event as one this tool introduced, and
    # failing identically on both would mean no reference with a pre-existing
    # flaw could ever be supplemented.
    def validate(path):
        gt = a.gff3validator or shutil.which('gt')
        if not gt:
            return 'not run', "no `gt` on PATH; pass --gff3validator /path/to/gt to enable"
        try:
            r = subprocess.run([gt, 'gff3validator', path],
                               capture_output=True, text=True, timeout=1800)
            err = (r.stderr or '') + (r.stdout or '')
            msg = next((l for l in err.splitlines() if 'error' in l.lower()), '')
            return ('valid' if r.returncode == 0 else 'INVALID'), msg
        except Exception as e:
            return 'could not run', str(e)

    def err_kind(msg):
        # Compare the KIND of error, not the text: it carries the filename and
        # a line number, which differ between the two files by construction.
        m = re.sub(r'in file "[^"]*"', '', msg)
        m = re.sub(r'\bline \d+\b', 'line N', m)
        m = re.sub(r'"[^"]*"', '"X"', m)
        return re.sub(r'\s+', ' ', m).strip()

    # Compare like with like. The validator stops at the FIRST error, so a
    # reference missing its ##gff-version line reports only that, and every
    # deeper defect it has stays hidden -- while the output, to which the merge
    # adds that line, reports the next one down. Measured: pepper's reference
    # reported only the missing version line, and the same file with the line
    # prepended reported the CDS phase error the output was being failed for.
    # Comparing the raw reference would therefore call an inherited defect new.
    ref_for_cmp = a.ref
    tmp_ref = None
    try:
        with open(a.ref, errors='replace') as fh:
            first = fh.readline()
        if not first.startswith('##gff-version'):
            tmp_ref = a.out + '.refcmp.gff3'
            with open(tmp_ref, 'w') as out, open(a.ref, errors='replace') as fh:
                out.write('##gff-version 3\n')
                for line in fh: out.write(line)
            ref_for_cmp = tmp_ref
    except OSError:
        pass
    ref_state, ref_msg = validate(ref_for_cmp)
    gt_state, gt_msg = validate(a.out)
    if tmp_ref:
        try: os.remove(tmp_ref)
        except OSError: pass
        if ref_state != 'valid':
            print("  (the reference was validated with a ##gff-version line prepended, "
                  "which it lacks; without it the validator stops there and hides the rest)")
    inherited = (ref_state == 'INVALID' and gt_state == 'INVALID'
                 and err_kind(ref_msg) == err_kind(gt_msg))
    print(f"gt gff3validator: reference {ref_state}"
          + (f"  {ref_msg}" if ref_msg else ""))
    print(f"gt gff3validator: output    {gt_state}"
          + (f"  {gt_msg}" if gt_msg else ""))
    if inherited:
        print("  ^ the same error is present in the reference as supplied: INHERITED, not introduced")
    elif ref_state == 'INVALID' and gt_state == 'valid':
        print("  ^ the reference was invalid and the output is valid (the merge writes a ##gff-version line)")

    if a.report:
        with open(a.report, 'w') as fh:
            fh.write("metric\tvalue\n")
            for k, v in (('reference_genes_in', ref_genes_in),
                         ('reference_genes_out', ref_genes_out),
                         # Renamed from alternative_transcripts_added, which used
                         # to mean the same thing and no longer does: a transcript
                         # can now be added under a record of its own, so the total
                         # and the under-a-reference-gene count differ. Reusing the
                         # old key would have changed what a number means without
                         # changing its name.
                         ('transcripts_added_total', added_tx),
                         ('reference_genes_paired_by_gm_compare', len(attach)),
                         ('reference_genes_that_received_a_transcript',
                          len({h[0] for k, h, _, _ in place.values() if k == 'gene'})),
                         ('gene_records_appended', added_genes),
                         ('transcripts_placed_on_cds_overlap', stats['placed_gene']),
                         ('transcripts_placed_inside_a_span_below_the_cds_threshold',
                          stats['placed_gene_by_span']),
                         ('transcripts_moved_off_the_locus_pairing', stats['reassigned']),
                         ('transcripts_moved_off_a_zero_overlap_pairing',
                          stats['reassigned_from_zero_overlap']),
                         ('transcripts_readthrough_own_record', stats['placed_readthrough']),
                         ('transcripts_no_reference_cds_overlap', stats['placed_no_overlap']),
                         ('transcripts_exceeding_a_reference_gene_span',
                          stats['placed_exceeds_span']),
                         ('reference_gene_lines_widened', len({r[0] for r in widened})),
                         ('widening_events', len(widened)),
                         ('widenings_refused_to_protect_a_neighbour', len(refused)),
                         ('total_genes_out', genes_out),
                         ('duplicate_ids', len(dup)),
                         ('orphan_parent_references', len(orphan)),
                         ('span_attributes_inconsistent', len(inconsistent)),
                         ('gt_gff3validator_output', gt_state),
                         ('gt_gff3validator_reference', ref_state),
                         ('output_invalidity_inherited', 'yes' if inherited else 'no')):
                fh.write(f"{k}\t{v}\n")
            fh.write("\n# reference gene lines widened to span an added transcript\n")
            fh.write("ref_gene\tgm_locus\toriginal_span\tnew_span\n")
            for r in widened: fh.write('\t'.join(r) + '\n')
            fh.write("\n# widenings refused because they would have reached into another "
                     "reference gene; the transcript was written as its own record instead\n")
            fh.write("ref_gene\ttranscript\tgm_locus\tspan_kept\tspan_refused\tblocked_by\n")
            for r in refused: fh.write('\t'.join(r) + '\n')

    bad = 0
    if gt_state == 'INVALID' and not inherited:
        print(f"FAIL: gt gff3validator rejected {a.out} with an error the reference does "
              f"not have: {gt_msg}", file=sys.stderr)
        print(f"      reference was: {ref_state}" + (f" ({ref_msg})" if ref_msg else ""), file=sys.stderr)
        bad = 1
    elif inherited:
        print(f"NOTE: the output is not valid GFF3, and neither is the reference it was built "
              f"from -- the same error. Not failing on it, because refusing would mean no "
              f"reference with a pre-existing flaw could be supplemented. The delivered file "
              f"carries the reference's defect; say so to whoever receives it.", file=sys.stderr)
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
