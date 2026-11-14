#!/usr/bin/env python3
"""annotate_rna_support.py -- per-TRANSCRIPT RNA-seq structural support.

The tier annotation already on these models is locus-level and categorical: it
says which KINDS of evidence a gene has (RNA / SwissProt / eggNOG), copied onto
each transcript. It does not say whether THIS transcript's structure is the one
the reads support, which is the question a curator asks of a revised model.

So for each mRNA this writes, from the pipeline's own alignments:

  rna_introns             introns in the transcript
  rna_introns_supported   how many carry >= --min-reads spliced reads
  rna_intron_support      the fraction, 0.00-1.00
  rna_min_junction_reads  reads on its WEAKEST intron -- the limiting evidence,
                          0 if any intron has none
  rna_support             full | partial | none | single_exon

`single_exon` is its own value and not `none`. A transcript with no introns has
no intron support to measure; calling that unsupported would condemn every
single-exon model on a quantity that does not exist for it.

Both sides of a merged annotation are measured, from the same junctions and on
the same quantity. Annotating only the added transcripts would tilt every
comparison toward whichever side carries numbers: a curator deciding between a
reference model and a proposed revision needs the reference's support too.

Introns come from CDS blocks on BOTH sides. The reference files here also carry
exon and UTR features, so a reference measured on exons would include UTR
introns and would not be comparable with our models, which have no exon
features. UTR introns are therefore outside this measure for everyone, which is
the price of the two numbers meaning the same thing.

The gene-level roll-up takes the BEST of a gene's transcripts
(rna_support_best); one well-supported isoform does not vouch for its siblings.

Usage:
  annotate_rna_support.py --gff in.gff3 --junctions j.tsv --out out.gff3
      [--min-reads 3] [--report r.tsv] [--selftest-locus ID]
"""
import argparse, sys, collections

def attrs(s):
    d = {}
    for kv in s.rstrip(';').split(';'):
        if kv and '=' in kv:
            k, v = kv.split('=', 1); d[k.strip()] = v.strip()
    return d

def fmt(d):
    return ';'.join(f"{k}={v}" for k, v in d.items() if v is not None and v != '')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', required=True)
    ap.add_argument('--junctions', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-reads', type=int, default=3)
    ap.add_argument('--report')
    ap.add_argument('--compare', help='TSV comparing the reference transcript with ours, per gene')
    ap.add_argument('--our-source', default='gene-miner',
                    help="column-2 value marking the added records")
    a = ap.parse_args()

    # junctions: (chrom, start, end) -> reads
    J = {}
    with open(a.junctions) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4: continue
            try: J[(p[0], int(p[1]), int(p[2]))] = int(p[3])
            except ValueError: continue
    if not J:
        sys.exit(f"FATAL: no junctions read from {a.junctions}")
    print(f"junctions loaded: {len(J):,}")

    # CDS blocks per transcript
    rows = []
    cds = collections.defaultdict(list)
    chrom_of = {}
    side_of = {}
    with open(a.gff, errors='replace') as fh:
        for line in fh:
            if line.startswith('#'): rows.append(('#', line.rstrip('\n'))); continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9: rows.append(('#', line.rstrip('\n'))); continue
            d = attrs(f[8]); rows.append((f[2], f, d))
            if f[2] == 'CDS' and 'Parent' in d:
                tx = d['Parent'].split(',')[0]
                cds[tx].append((int(f[3]), int(f[4])))
                chrom_of[tx] = f[0]
                side_of[tx] = 'ours' if f[1] == a.our_source else 'ref'

    # per-transcript support
    sup = {}
    for tx, blocks in cds.items():
        blocks.sort()
        introns = [(blocks[i][1] + 1, blocks[i+1][0] - 1) for i in range(len(blocks) - 1)]
        introns = [(s, e) for s, e in introns if e >= s]
        c = chrom_of[tx]
        if not introns:
            sup[tx] = dict(n=0, k=0, frac=None, weakest=None, cls='single_exon')
            continue
        counts = [J.get((c, s, e), 0) for s, e in introns]
        k = sum(1 for x in counts if x >= a.min_reads)
        cls = 'full' if k == len(introns) else ('none' if k == 0 else 'partial')
        sup[tx] = dict(n=len(introns), k=k, frac=k / len(introns),
                       weakest=min(counts), cls=cls)

    # gene roll-up: the best of its transcripts, so a strong isoform cannot
    # vouch for a weak sibling
    rank = {'full': 3, 'partial': 2, 'none': 1, 'single_exon': 0}
    best = {}
    for ftype, *rest in rows:
        if ftype != 'mRNA': continue
        f, d = rest
        tx, par = d.get('ID'), d.get('Parent', '').split(',')[0]
        if tx in sup and par:
            s = sup[tx]
            cur = best.get(par)
            if cur is None or rank[s['cls']] > rank[cur]: best[par] = s['cls']

    written = collections.Counter()
    with open(a.out, 'w') as out:
        for item in rows:
            if item[0] == '#': out.write(item[1] + '\n'); continue
            ftype, f, d = item
            if ftype == 'mRNA' and d.get('ID') in sup:
                s = sup[d['ID']]
                d['rna_introns'] = s['n']
                d['rna_introns_supported'] = s['k']
                if s['frac'] is not None:
                    d['rna_intron_support'] = f"{s['frac']:.2f}"
                    d['rna_min_junction_reads'] = s['weakest']
                d['rna_support'] = s['cls']
                written[s['cls']] += 1
                f = list(f); f[8] = fmt(d)
            elif ftype == 'gene' and d.get('ID') in best:
                d['rna_support_best'] = best[d['ID']]
                f = list(f); f[8] = fmt(d)
            out.write('\t'.join(f) + '\n')

    tot = sum(written.values())
    print(f"transcripts annotated: {tot:,}  (min-reads={a.min_reads})")
    by_side = collections.defaultdict(collections.Counter)
    for tx, sdat in sup.items():
        by_side[side_of.get(tx, 'ref')][sdat['cls']] += 1
    for sd in ('ref', 'ours'):
        n = sum(by_side[sd].values())
        if not n: continue
        label = '기준 전사체' if sd == 'ref' else '우리 전사체'
        print(f"  {label} ({n:,})")
        for k in ('full', 'partial', 'none', 'single_exon'):
            if by_side[sd][k]:
                print(f"    {k:<12} {by_side[sd][k]:>7,}  {100*by_side[sd][k]/n:5.1f}%")
    if a.report:
        with open(a.report, 'w') as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"min_reads\t{a.min_reads}\n")
            fh.write(f"junctions_in_rnaseq\t{len(J)}\n")
            fh.write(f"transcripts\t{tot}\n")
            for k in ('full', 'partial', 'none', 'single_exon'):
                fh.write(f"transcripts_{k}\t{written[k]}\n")
            for sd in ('ref', 'ours'):
                for k in ('full', 'partial', 'none', 'single_exon'):
                    fh.write(f"{sd}_transcripts_{k}\t{by_side[sd][k]}\n")
            mi = [s for s in sup.values() if s['n'] > 0]
            fh.write(f"multi_exon_transcripts\t{len(mi)}\n")
            if mi:
                fh.write(f"introns_total\t{sum(s['n'] for s in mi)}\n")
                fh.write(f"introns_supported\t{sum(s['k'] for s in mi)}\n")
    if a.compare:
        # Per gene, the best-supported reference transcript beside the
        # best-supported transcript of ours. This is the table a curator works
        # from: it says which proposed revisions the reads back BETTER than the
        # model they would revise.
        gene_tx = collections.defaultdict(lambda: {'ref': [], 'ours': []})
        for item in rows:
            if item[0] != 'mRNA': continue
            f, d = item[1], item[2]
            tx, par = d.get('ID'), d.get('Parent', '').split(',')[0]
            if tx in sup and par:
                gene_tx[par][side_of.get(tx, 'ref')].append(tx)
        def pick(txs):
            best = None
            for t in txs:
                s2 = sup[t]
                key = (s2['frac'] if s2['frac'] is not None else -1, s2['n'])
                if best is None or key > best[0]: best = (key, t)
            return best[1] if best else None
        n_rows = n_better = n_same = n_worse = 0
        with open(a.compare, 'w') as fh:
            fh.write("gene\tref_tx\tref_support\tref_introns\tref_supported\t"
                     "our_tx\tour_support\tour_introns\tour_supported\tcomparison\n")
            for g, d2 in sorted(gene_tx.items()):
                if not d2['ref'] or not d2['ours']: continue
                rt, ot = pick(d2['ref']), pick(d2['ours'])
                rs, os_ = sup[rt], sup[ot]
                # Compared on the supported fraction. A fraction over different
                # intron counts is not a like-for-like ratio, so both counts are
                # printed and this column is a ranking aid, not a verdict.
                rf = rs['frac'] if rs['frac'] is not None else -1
                of = os_['frac'] if os_['frac'] is not None else -1
                cmp_ = 'ours_better' if of > rf else ('reference_better' if rf > of else 'equal')
                n_rows += 1
                n_better += cmp_ == 'ours_better'; n_same += cmp_ == 'equal'; n_worse += cmp_ == 'reference_better'
                fh.write(f"{g}\t{rt}\t{rs['cls']}\t{rs['n']}\t{rs['k']}\t"
                         f"{ot}\t{os_['cls']}\t{os_['n']}\t{os_['k']}\t{cmp_}\n")
        print(f"comparison table: {n_rows:,} gene(s) carrying both a reference and an added transcript")
        if n_rows:
            print(f"  ours better        {n_better:>7,}  {100*n_better/n_rows:5.1f}%")
            print(f"  equal              {n_same:>7,}  {100*n_same/n_rows:5.1f}%")
            print(f"  reference better   {n_worse:>7,}  {100*n_worse/n_rows:5.1f}%")

    # A run where nothing is supported almost always means the coordinate
    # conventions disagree, not that the models are wrong. Say so loudly.
    if tot and written['none'] + written['single_exon'] == tot:
        print("WARNING: not one multi-exon transcript had a supported intron.", file=sys.stderr)
        print("         Suspect the intron coordinate convention before the models.", file=sys.stderr)
        sys.exit(2)

if __name__ == '__main__':
    main()
