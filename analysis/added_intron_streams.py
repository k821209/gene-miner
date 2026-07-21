#!/usr/bin/env python3
# added_intron_streams.py gm_final.gff3 ref.gff3 junctions.tsv transdecoder.gff3 genemark.gtf augustus.gff3 min_reads
# For REVISED loci (GM matched a ref gene same-strand CDS>=0.5 but structure differs),
# take the GM-ADDED CDS introns (in GM, not in the reference) and classify each by which
# evidence stream contributes it (RNA-seq TransDecoder / GeneMark-ETP / AUGUSTUS), then
# report the RNA-seq splice-junction confirmation rate of each class.
# Tests the hypothesis: the ADDED introns that lack RNA-seq support are the ab-initio ones,
# while RNA-seq-derived additions are self-confirming.
import sys, re
from collections import defaultdict

gm_f, ref_f, jbed, td_f, gm_gtf, aug_f = sys.argv[1:7]
minr = int(sys.argv[7]) if len(sys.argv) > 7 else 3

def attr_gff(a, k):
    m = re.search(k + r'=([^;]+)', a); return m.group(1) if m else None
def tx_key(a):
    # GFF3 Parent= or GTF transcript_id "X"
    p = attr_gff(a, 'Parent')
    if p: return p
    m = re.search(r'transcript_id "([^"]+)"', a)
    return m.group(1) if m else None

def cds_introns(path, pc=False):
    """Return (gene_introns dict gid->set, all_introns set) from CDS features. gid via mRNA Parent when present, else tx key."""
    txg = {}; bt = {}; tx_cds = defaultdict(list)
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        t, a = f[2], f[8]
        if t == 'gene': bt[attr_gff(a, 'ID')] = attr_gff(a, 'biotype') or ''
        elif t in ('mRNA', 'transcript'):
            i = attr_gff(a, 'ID')
            if i: txg[i] = attr_gff(a, 'Parent')
        elif t == 'CDS':
            k = tx_key(a)
            if k: tx_cds[k].append((int(f[3]), int(f[4]), f[0]))
    genes = {}; allint = set()
    for tx, segs in tx_cds.items():
        g = txg.get(tx, tx); chrom = segs[0][2]
        iv = sorted((s, e) for s, e, _ in segs)
        d = genes.setdefault(g, {'chrom': chrom, 'cds': set(), 'introns': set(), 'bt': bt.get(g, '')})
        d['cds'].update(iv)
        for i in range(len(iv) - 1):
            it = (chrom, iv[i][1] + 1, iv[i + 1][0] - 1)
            d['introns'].add(it); allint.add(it)
    if pc:
        genes = {g: d for g, d in genes.items() if d['bt'] in ('', 'protein_coding')}
    return genes, allint

def clen(iv): return sum(e - s + 1 for s, e in iv)
def ofrac(a, b):
    ov = 0
    for s, e in a:
        for bs, be in b:
            lo, hi = max(s, bs), min(e, be)
            if hi >= lo: ov += hi - lo + 1
    return ov / (min(clen(a), clen(b)) or 1)

gm, _ = cds_introns(gm_f)
ref, _ = cds_introns(ref_f, pc=True)
_, td_int = cds_introns(td_f)
_, gmk_int = cds_introns(gm_gtf)
_, aug_int = cds_introns(aug_f)

# strand/cds for matching
def gmeta(d):
    cds = sorted(d['cds']);
    return dict(cds=cds, introns=d['introns'], start=min(s for s, e in cds), end=max(e for s, e in cds))
# need strand: re-parse quickly for strand per gene (first CDS strand)
def strands(path, use_parent=True):
    txg = {}; st = {}
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        t, a = f[2], f[8]
        if t in ('mRNA', 'transcript'):
            i = attr_gff(a, 'ID')
            if i: txg[i] = attr_gff(a, 'Parent')
        elif t == 'CDS':
            k = tx_key(a); g = txg.get(k, k)
            st.setdefault(g, f[6])
    return st
gm_st = strands(gm_f); ref_st = strands(ref_f)

J = set()
for line in open(jbed):
    f = line.rstrip('\n').split('\t')
    if len(f) < 4: continue
    if int(f[3]) < minr: continue
    J.add((f[0], int(f[1]), int(f[2])))

# reference index by (chrom, strand)
ridx = defaultdict(list)
for gid, d in ref.items():
    cds = sorted(d['cds']); chrom = d['chrom']; strand = ref_st.get(gid, '+')
    ridx[(chrom, strand)].append((gid, cds, d['introns']))

# tallies for GM-added introns, by stream membership
classes = ['rna', 'genemark', 'augustus', 'other']
added_total = 0; added_conf = 0
cnt = {c: 0 for c in classes}; conf = {c: 0 for c in classes}
# also cross-tab: rna vs abinitio-only (not in TD)
abinit_total = 0; abinit_conf = 0; rna_total = 0; rna_conf = 0

for gid, d in gm.items():
    cds = sorted(d['cds']); chrom = d['chrom']; strand = gm_st.get(gid, '+')
    gstart = min(s for s, e in cds); gend = max(e for s, e in cds)
    best = None; bro = 0.0
    for rgid, rcds, rint in ridx.get((chrom, strand), []):
        re_ = max(e for s, e in rcds); rs_ = min(s for s, e in rcds)
        if re_ < gstart or rs_ > gend: continue
        r = ofrac(cds, rcds)
        if r > bro: bro = r; best = (rcds, rint)
    if not best or bro < 0.5: continue
    rcds, rint = best
    if set(cds) == set(rcds): continue          # identical
    gi = d['introns']; ri = rint
    if not gi or not ri: continue                # both multi-exon
    for it in gi - ri:                           # GM-added introns
        added_total += 1
        c = (it in J)
        added_conf += c
        in_rna = it in td_int
        if in_rna:
            rna_total += 1; rna_conf += c
            cnt['rna'] += 1; conf['rna'] += c
        else:
            abinit_total += 1; abinit_conf += c
            if it in gmk_int:
                cnt['genemark'] += 1; conf['genemark'] += c
            elif it in aug_int:
                cnt['augustus'] += 1; conf['augustus'] += c
            else:
                cnt['other'] += 1; conf['other'] += c

def pct(a, b): return f"{100*a/max(b,1):.1f}%"
print(f"# min_reads={minr}")
print(f"GM_added_introns_total\t{added_total}\tconfirmed\t{pct(added_conf, added_total)}")
print(f"--- split: RNA-derived (in TransDecoder) vs ab-initio-only ---")
print(f"RNA_derived\t{rna_total}\t({pct(rna_total, added_total)} of added)\tconfirmed\t{pct(rna_conf, rna_total)}")
print(f"ab_initio_only\t{abinit_total}\t({pct(abinit_total, added_total)} of added)\tconfirmed\t{pct(abinit_conf, abinit_total)}")
print(f"--- ab-initio-only, by finder (mutually exclusive: GeneMark checked first) ---")
print(f"  GeneMark\t{cnt['genemark']}\tconfirmed\t{pct(conf['genemark'], cnt['genemark'])}")
print(f"  AUGUSTUS\t{cnt['augustus']}\tconfirmed\t{pct(conf['augustus'], cnt['augustus'])}")
print(f"  neither(merge)\t{cnt['other']}\tconfirmed\t{pct(conf['other'], cnt['other'])}")
