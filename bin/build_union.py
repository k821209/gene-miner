import re
from collections import defaultdict, OrderedDict
# build_union.py — gene set = usable AUGUSTUS loci (score>=SC, >=MINAA aa) UNION
#   RNA-only TransDecoder loci  UNION  GeneMark-ETP loci (optional, --genemark).
#   ISOFORM-AWARE: each usable AUGUSTUS locus is AUGMENTED with the isoforms of
#   any RNA-seq / GeneMark transcript that overlaps it (same strand) — AUGUSTUS
#   provides the clean primary model, the others add alternative splicing.
#   Transcripts with an identical CDS structure are de-duplicated. Evidence
#   streams that overlap no existing locus become their own loci:
#     <prefix>_A ab-initio AUGUSTUS, <prefix>_R RNA-only, <prefix>_E GeneMark-only.
# Inputs (fixed names, CWD): augustus_scaffold.gff3 , annot/genome.transdecoder.gff3
#   optional: GeneMark-ETP genemark.gtf via --genemark
# Output: union.gff3

def load(gff):
    m2g = {}
    gscore = {}
    tx = defaultdict(OrderedDict)
    gchr, gstr = {}, {}
    for l in open(gff):
        if l.startswith('#') or not l.strip():
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        at = f[8]
        def A(k):
            m = re.search(k + r'=([^;]+)', at)
            return m.group(1) if m else None
        if f[2] == 'gene':
            try:
                gscore[A('ID')] = float(f[5])
            except (TypeError, ValueError):
                gscore[A('ID')] = 0.0
        elif f[2] in ('mRNA', 'transcript'):
            m2g[A('ID')] = A('Parent')
        elif f[2] == 'CDS':
            tid = A('Parent')
            g = m2g.get(tid, tid)
            tx[g].setdefault(tid, []).append((int(f[3]), int(f[4]), f[7]))
            gchr[g] = f[0]
            gstr[g] = f[6]
    genes = {}
    for g, txd in tx.items():
        allcds = [iv for t in txd.values() for iv in t]
        maxaa = max(sum(e - s + 1 for s, e, p in t) // 3 for t in txd.values())
        genes[g] = dict(chrom=gchr[g], strand=gstr[g],
                        start=min(s for s, e, p in allcds),
                        end=max(e for s, e, p in allcds),
                        score=gscore.get(g, 1.0), aa=maxaa, tx=txd)
    return genes


def load_gtf(gtf):
    """GeneMark-ETP genemark.gtf: CDS rows carry gene_id/transcript_id/status.
    Keep only status 'complete' transcripts; GTF frame carried in col 8."""
    tx = defaultdict(OrderedDict)
    gchr, gstr = {}, {}
    for l in open(gtf):
        if l.startswith('#') or not l.strip():
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        at = f[8]
        def A(k):
            m = re.search(k + r'\s+"([^"]+)"', at)
            return m.group(1) if m else None
        st = A('status')
        if st is not None and st != 'complete':
            continue
        tid = A('transcript_id')
        g = A('gene_id')
        if not tid or not g:
            continue
        ph = f[7] if f[7] in ('0', '1', '2') else '0'
        tx[g].setdefault(tid, []).append((int(f[3]), int(f[4]), ph))
        gchr[g] = f[0]
        gstr[g] = f[6]
    genes = {}
    for g, txd in tx.items():
        allcds = [iv for t in txd.values() for iv in t]
        if not allcds:
            continue
        maxaa = max(sum(e - s + 1 for s, e, p in t) // 3 for t in txd.values())
        genes[g] = dict(chrom=gchr[g], strand=gstr[g],
                        start=min(s for s, e, p in allcds),
                        end=max(e for s, e, p in allcds),
                        score=1.0, aa=maxaa, tx=txd)
    return genes


import argparse as _argparse
_ap = _argparse.ArgumentParser()
_ap.add_argument('--prefix', default='GMG',
                 help='gene-ID prefix (<prefix>_A ab-initio, <prefix>_R RNA-only, <prefix>_E GeneMark-only)')
_ap.add_argument('--genemark', default=None,
                 help='optional GeneMark-ETP genemark.gtf to add as a third evidence stream')
_args = _ap.parse_known_args()[0]
PREFIX = _args.prefix

aug = load('augustus_scaffold.gff3')
td = load('annot/genome.transdecoder.gff3')
SC, MINAA = 0.8, 100

passing = {g: v for g, v in aug.items() if v['score'] >= SC and v['aa'] >= MINAA}


# index usable loci for overlap lookup; tag = which dict the host lives in
def build_index(*dicts):
    idx = defaultdict(list)
    for tag, d in dicts:
        for g, v in d.items():
            idx[(v['chrom'], v['strand'])].append((v['start'], v['end'], g, tag))
    for k in idx:
        idx[k].sort()
    return idx


def best_host(idx, c, st, s, e):
    best, bov, btag = None, 0, None
    for as_, ae, g, tag in idx.get((c, st), []):
        if as_ > e:
            break
        lo, hi = max(s, as_), min(e, ae)
        if hi >= lo and (hi - lo) > bov:
            bov, best, btag = hi - lo, g, tag
    return best, btag

# assign each RNA-seq locus to an overlapping AUGUSTUS locus (-> extra isoforms)
# or keep it as an RNA-only locus
idx = build_index(('A', passing))
rna_only = {}
for g, v in td.items():
    host, _ = best_host(idx, v['chrom'], v['strand'], v['start'], v['end'])
    if host:
        for tid, cds in v['tx'].items():
            passing[host]['tx']['rna_' + tid] = cds
    else:
        rna_only[g] = v

# assign each GeneMark-ETP locus to an overlapping AUGUSTUS *or* RNA-only locus
# (-> extra isoforms) or keep it as a GeneMark-only locus
gm_only = {}
if _args.genemark:
    gm = load_gtf(_args.genemark)
    gm = {g: v for g, v in gm.items() if v['aa'] >= MINAA}
    idx2 = build_index(('A', passing), ('R', rna_only))
    for g, v in gm.items():
        host, tag = best_host(idx2, v['chrom'], v['strand'], v['start'], v['end'])
        if host:
            tgt = passing if tag == 'A' else rna_only
            for tid, cds in v['tx'].items():
                tgt[host]['tx']['gm_' + tid] = cds
        else:
            gm_only[g] = v


def dedup_tx(txd):
    seen, out = set(), OrderedDict()
    for tid, cds in txd.items():
        key = tuple(sorted((s, e) for s, e, p in cds))
        if not key or key in seen:
            continue
        seen.add(key)
        out[tid] = cds
    return out

for d in (passing, rna_only, gm_only):
    for g in d:
        d[g]['tx'] = dedup_tx(d[g]['tx'])
        allcds = [iv for t in d[g]['tx'].values() for iv in t]
        if allcds:
            d[g]['start'] = min(s for s, e, p in allcds)
            d[g]['end'] = max(e for s, e, p in allcds)


def write(o, src, genes, prefix):
    ng = nt = 0
    for g, v in genes.items():
        if not v['tx']:
            continue
        ng += 1
        gid = "%s%06d" % (prefix, ng)
        c, st = v['chrom'], v['strand']
        o.write(f"{c}\t{src}\tgene\t{v['start']}\t{v['end']}\t.\t{st}\t.\tID={gid}\n")
        for ti, (tid, cdss) in enumerate(v['tx'].items(), 1):
            nt += 1
            mid = f"{gid}.t{ti}"
            cdss = sorted(cdss)
            o.write(f"{c}\t{src}\tmRNA\t{min(a for a,b,p in cdss)}\t{max(b for a,b,p in cdss)}\t.\t{st}\t.\tID={mid};Parent={gid}\n")
            for a, b, ph in cdss:
                o.write(f"{c}\t{src}\tCDS\t{a}\t{b}\t.\t{st}\t{ph}\tID={mid}.cds;Parent={mid}\n")
    return ng, nt

with open('union.gff3', 'w') as o:
    nag, nat = write(o, 'GeneMiner', passing, PREFIX + '_A')
    nrg, nrt = write(o, 'RNAseq', rna_only, PREFIX + '_R')
    neg, net = write(o, 'GeneMarkETP', gm_only, PREFIX + '_E')

print(f"AUGUSTUS loci (+RNA/GeneMark isoforms): {nag} genes / {nat} transcripts")
print(f"RNA-only loci: {nrg} genes / {nrt} transcripts")
print(f"GeneMark-only loci: {neg} genes / {net} transcripts")
print(f"=== UNION: {nag + nrg + neg} genes / {nat + nrt + net} transcripts "
      f"({(nat+nrt+net)/(nag+nrg+neg):.2f} iso/gene) ===")
