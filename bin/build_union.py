import argparse
import re
from collections import defaultdict, OrderedDict
# build_union.py — gene set = usable AUGUSTUS loci (score>=SC, >=MINAA aa) UNION
#   RNA-only TransDecoder loci  UNION  GeneMark-ETP loci (optional, --genemark).
#   ISOFORM-AWARE: each usable AUGUSTUS locus is AUGMENTED with the RNA-seq /
#   GeneMark transcripts that overlap it on the same strand — AUGUSTUS provides
#   the clean primary model, the others add alternative splicing. Transcripts
#   with an identical CDS structure are de-duplicated. Evidence streams that
#   overlap no existing locus become their own loci:
#     <prefix>_A ab-initio AUGUSTUS, <prefix>_R RNA-only, <prefix>_E GeneMark-only.
#
#   Folding is per TRANSCRIPT, and overlap is measured on coding BASES shared
#   with the locus' merged CDS footprint (the largest-overlap locus wins).
#   Folding a whole StringTie locus on any span overlap, as earlier versions did,
#   moved every transcript of a read-through assembly onto one gene and so merged
#   neighbouring genes; per-transcript folding keeps each transcript with the
#   locus it actually shares coding sequence with.
#
# Inputs (fixed names, CWD): augustus_scaffold.gff3 , annot/genome.transdecoder.gff3
#   optional: GeneMark-ETP genemark.gtf via --genemark
# Output: union.gff3

SC, MINAA = 0.8, 100


def merge(iv):
    """merge [(start, end, ...)] into disjoint [start, end] blocks"""
    out = []
    for s, e, *_ in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def blen(m):
    return sum(e - s + 1 for s, e in m)


def shared(a, b):
    """coding bases shared by two merged-block lists"""
    i = j = ov = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi >= lo:
            ov += hi - lo + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return ov


def load(gff, stream, scored):
    """GFF3 -> one record per transcript; `scored` applies the AUGUSTUS gene score"""
    m2g, gscore, tx, loc = {}, {}, defaultdict(list), {}
    for l in open(gff):
        if l.startswith('#') or not l.strip():
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = f[8]

        def A(k):
            m = re.search(k + r'=([^;]+)', a)
            return m.group(1) if m else None
        if f[2] == 'gene':
            try:
                gscore[A('ID')] = float(f[5])
            except (TypeError, ValueError):
                gscore[A('ID')] = 0.0
        elif f[2] in ('mRNA', 'transcript'):
            m2g[A('ID')] = A('Parent')
        elif f[2] == 'CDS':
            t = A('Parent').split(',')[0]
            tx[t].append((int(f[3]), int(f[4]), f[7]))
            loc[t] = (f[0], f[6])
    out = []
    for t, cds in tx.items():
        g = m2g.get(t, t)
        if scored and gscore.get(g, 1.0) < SC:
            continue
        m = merge(cds)
        out.append(dict(id=f'{stream}:{t}', gene=f'{stream}:{g}', chrom=loc[t][0],
                        strand=loc[t][1], cds=cds, m=m, start=m[0][0], end=m[-1][1],
                        aa=blen(m) // 3))
    return out


def load_gtf(gtf):
    """GeneMark-ETP genemark.gtf: CDS rows carry gene_id/transcript_id/status.
    Keep only status 'complete' transcripts; GTF frame carried in col 8."""
    tx, loc, gene = defaultdict(list), {}, {}
    for l in open(gtf):
        if l.startswith('#') or not l.strip():
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = f[8]

        def A(k):
            m = re.search(k + r'\s+"([^"]+)"', a)
            return m.group(1) if m else None
        if (A('status') or 'complete') != 'complete':
            continue
        t, g = A('transcript_id'), A('gene_id')
        if not t or not g:
            continue
        tx[t].append((int(f[3]), int(f[4]), f[7] if f[7] in ('0', '1', '2') else '0'))
        loc[t] = (f[0], f[6])
        gene[t] = g
    out = []
    for t, cds in tx.items():
        m = merge(cds)
        out.append(dict(id=f'E:{t}', gene=f'E:{gene[t]}', chrom=loc[t][0], strand=loc[t][1],
                        cds=cds, m=m, start=m[0][0], end=m[-1][1], aa=blen(m) // 3))
    return out


ap = argparse.ArgumentParser()
ap.add_argument('--prefix', default='GMG',
                help='gene-ID prefix (<prefix>_A ab-initio, <prefix>_R RNA-only, <prefix>_E GeneMark-only)')
ap.add_argument('--genemark', default=None,
                help='optional GeneMark-ETP genemark.gtf to add as a third evidence stream')
ap.add_argument('--out', default='union.gff3')
args = ap.parse_args()

aug = [t for t in load('augustus_scaffold.gff3', 'A', True) if t['aa'] >= MINAA]
rna = load('annot/genome.transdecoder.gff3', 'R', False)
gmk = [t for t in load_gtf(args.genemark) if t['aa'] >= MINAA] if args.genemark else []

loci = []                        # dict(chrom, strand, m, tx, src)
index = defaultdict(list)        # (chrom, strand) -> locus indices


def new_locus(t, src):
    loci.append(dict(chrom=t['chrom'], strand=t['strand'], m=[list(b) for b in t['m']],
                     tx=OrderedDict([(t['id'], t['cds'])]), src=src))
    index[(t['chrom'], t['strand'])].append(len(loci) - 1)
    return loci[-1]


def add_tx(L, t):
    L['tx'][t['id']] = t['cds']
    L['m'] = merge([tuple(b) for b in L['m']] + [tuple(b) for b in t['m']])


def hosts(t, src=None):
    """loci on the same strand sharing coding bases with t, best overlap first"""
    res = []
    for i in index.get((t['chrom'], t['strand']), []):
        L = loci[i]
        if L['m'][-1][1] < t['start'] or L['m'][0][0] > t['end']:
            continue
        if src and L['src'] != src:
            continue
        ov = shared(t['m'], L['m'])
        if ov:
            res.append((ov, i))
    return sorted(res, reverse=True)


# the usable AUGUSTUS models seed the loci (isoforms of one AUGUSTUS gene together)
agene = OrderedDict()
for t in aug:
    agene.setdefault(t['gene'], []).append(t)
for ts in agene.values():
    L = new_locus(ts[0], 'A')
    for t in ts[1:]:
        add_tx(L, t)


def fold(stream, src):
    """attach each transcript to the locus it shares most coding sequence with;
    transcripts that overlap none found their own loci of this stream"""
    pending = []
    for t in stream:
        h = hosts(t)
        if h:
            add_tx(loci[h[0][1]], t)
        else:
            pending.append(t)
    for t in pending:
        h = hosts(t, src)
        if h and h[0][0] >= 0.5 * blen(t['m']):
            add_tx(loci[h[0][1]], t)
        else:
            new_locus(t, src)


fold(rna, 'R')
fold(gmk, 'E')


def dedup(txd):
    seen, out = set(), OrderedDict()
    for tid, cds in txd.items():
        key = tuple(sorted((s, e) for s, e, p in cds))
        if key and key not in seen:
            seen.add(key)
            out[tid] = cds
    return out


counts = {}
with open(args.out, 'w') as o:
    for src, source in (('A', 'GeneMiner'), ('R', 'RNAseq'), ('E', 'GeneMarkETP')):
        ng = nt = 0
        for L in loci:
            if L['src'] != src:
                continue
            txd = dedup(L['tx'])
            if not txd:
                continue
            ng += 1
            gid = "%s_%s%06d" % (args.prefix, src, ng)
            allc = [iv for c in txd.values() for iv in c]
            o.write(f"{L['chrom']}\t{source}\tgene\t{min(s for s, e, p in allc)}\t"
                    f"{max(e for s, e, p in allc)}\t.\t{L['strand']}\t.\tID={gid}\n")
            for i, (tid, cds) in enumerate(txd.items(), 1):
                nt += 1
                mid, cds = f"{gid}.t{i}", sorted(cds)
                o.write(f"{L['chrom']}\t{source}\tmRNA\t{cds[0][0]}\t{max(e for s, e, p in cds)}"
                        f"\t.\t{L['strand']}\t.\tID={mid};Parent={gid};source_tx={tid}\n")
                for s, e, ph in cds:
                    o.write(f"{L['chrom']}\t{source}\tCDS\t{s}\t{e}\t.\t{L['strand']}\t{ph}"
                            f"\tID={mid}.cds;Parent={mid}\n")
        counts[src] = (ng, nt)

print(f"AUGUSTUS loci (+RNA/GeneMark isoforms): {counts['A'][0]} genes / {counts['A'][1]} transcripts")
print(f"RNA-only loci: {counts['R'][0]} genes / {counts['R'][1]} transcripts")
print(f"GeneMark-only loci: {counts['E'][0]} genes / {counts['E'][1]} transcripts")
tot_g = sum(c[0] for c in counts.values())
tot_t = sum(c[1] for c in counts.values())
print(f"=== UNION: {tot_g} genes / {tot_t} transcripts ({tot_t / tot_g:.2f} iso/gene) ===")
