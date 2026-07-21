#!/usr/bin/env python3
# pick_worse_example.py gm_final.gff3 ref.gff3 junctions.tsv transdecoder.gff3 genemark.gtf augustus.gff3 min_reads
# Find clean "GM worse" example loci: revised, shared introns all RNA-confirmed, but GM adds
# ab-initio-only (not in TransDecoder) introns that are NOT confirmed -> GM's confirmed FRACTION
# drops below the reference's. Prints candidates for a per-locus figure (panel b).
import sys, re
from collections import defaultdict
gm_f, ref_f, jbed, td_f, gm_gtf, aug_f = sys.argv[1:7]
minr = int(sys.argv[7]) if len(sys.argv) > 7 else 3

def attr_gff(a, k):
    m = re.search(k + r'=([^;]+)', a); return m.group(1) if m else None
def tx_key(a):
    p = attr_gff(a, 'Parent')
    if p: return p
    m = re.search(r'transcript_id "([^"]+)"', a); return m.group(1) if m else None
def cds_introns(path, pc=False):
    txg = {}; bt = {}; tx_cds = defaultdict(list); st = {}
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
            if k: tx_cds[k].append((int(f[3]), int(f[4]), f[0], f[6]))
    genes = {}; allint = set()
    for tx, segs in tx_cds.items():
        g = txg.get(tx, tx); chrom = segs[0][2]; strand = segs[0][3]
        iv = sorted((s, e) for s, e, _, _ in segs)
        d = genes.setdefault(g, {'chrom': chrom, 'strand': strand, 'cds': set(), 'introns': set(), 'bt': bt.get(g, '')})
        d['cds'].update(iv)
        for i in range(len(iv) - 1):
            it = (chrom, iv[i][1] + 1, iv[i + 1][0] - 1)
            d['introns'].add(it); allint.add(it)
    if pc: genes = {g: d for g, d in genes.items() if d['bt'] in ('', 'protein_coding')}
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
_, td = cds_introns(td_f)
_, gmk = cds_introns(gm_gtf)
_, aug = cds_introns(aug_f)
J = set()
for line in open(jbed):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 4 and int(f[3]) >= minr: J.add((f[0], int(f[1]), int(f[2])))
ridx = defaultdict(list)
for gid, d in ref.items():
    cds = sorted(d['cds'])
    ridx[(d['chrom'], d['strand'])].append((gid, cds, d['introns'], min(s for s,e in cds), max(e for s,e in cds)))

def frac(introns):
    return sum(1 for it in introns if it in J) / len(introns) if introns else None

cands = []
for gid, d in gm.items():
    cds = sorted(d['cds']); chrom = d['chrom']; strand = d['strand']
    gs = min(s for s, e in cds); ge = max(e for s, e in cds)
    best = None; bro = 0.0
    for rgid, rcds, rint, rs, re_ in ridx.get((chrom, strand), []):
        if re_ < gs or rs > ge: continue
        r = ofrac(cds, rcds)
        if r > bro: bro = r; best = (rgid, rcds, rint)
    if not best or bro < 0.5: continue
    rgid, rcds, rint = best
    if set(cds) == set(rcds): continue
    gi = d['introns']; ri = rint
    if not gi or not ri: continue
    gf = frac(gi); rf = frac(ri)
    if gf is None or rf is None: continue
    if not (gf < rf - 1e-9): continue                       # want GM worse
    if rf < 0.999: continue                                  # reference fully confirmed
    added = gi - ri; shared = gi & ri
    if not shared or not added: continue
    if any(it not in J for it in shared): continue           # shared all confirmed
    # added introns should be ab-initio-only + unconfirmed
    added_unconf = [it for it in added if it not in J]
    if not added_unconf: continue
    abinit = [it for it in added_unconf if it not in td and (it in gmk or it in aug)]
    if len(abinit) != len(added_unconf): continue            # all unconfirmed-added are ab initio
    if len(shared) != len(ri): continue                      # GM RETAINS every reference intron
    if not (2 <= len(added) <= 3): continue                  # adds only 2-3 unconfirmed introns
    ni = len(gi)
    width = ge - gs
    if width > 8000: continue
    src = {'gmk': sum(1 for it in abinit if it in gmk), 'aug': sum(1 for it in abinit if it in aug and it not in gmk)}
    cands.append((len(abinit), width, chrom, gs, ge, strand, gid, rgid, ni, len(ri), len(shared), len(added), gf, rf, src))

cands.sort(key=lambda x: (-x[0], x[1]))
print(f"# {len(cands)} clean GM-worse candidates (ref fully confirmed, shared confirmed, added-unconf all ab initio)")
print("addU\twidth\tchrom\tstart\tend\tstr\tGM_gene\tref_gene\tgmI\trefI\tshar\tadd\tGMfrac\treffrac\tgmk/aug")
for c in cands[:15]:
    print(f"{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t{c[4]}\t{c[5]}\t{c[6]}\t{c[7]}\t{c[8]}\t{c[9]}\t{c[10]}\t{c[11]}\t{c[12]:.2f}\t{c[13]:.2f}\t{c[14]['gmk']}/{c[14]['aug']}")
