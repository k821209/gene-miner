#!/usr/bin/env python3
"""kaks_pipeline.py — selection test for Gene-Miner loci against a wild relative.

For each query protein, the best miniprot spliced alignment in the wild genome
gives the homologous coding sequence; the two proteins are aligned with MAFFT,
the alignment is back-translated to codons, and yn00 (PAML) estimates dN, dS and
dN/dS. Run once per query class (novel / reference-matching / TE-removed) so the
distributions can be compared.

Usage: kaks_pipeline.py <query.cds.fa> <query.pep.fa> <wild.fa> <outdir> [threads]
  query.cds.fa / query.pep.fa : same ids, CDS nucleotides and their protein
Writes <outdir>/kaks.tsv : query  target_locus  identity  dN  dS  dNdS  len_codons
"""
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

CODON_STOPS = {'TAA', 'TAG', 'TGA'}


def fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, ''.join(seq)


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


CODON = {}
_bases = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, b1 in enumerate(_bases):
    for j, b2 in enumerate(_bases):
        for k, b3 in enumerate(_bases):
            CODON[b1 + b2 + b3] = _aa[i * 16 + j * 4 + k]


def translate(nt):
    return ''.join(CODON.get(nt[i:i + 3].upper(), 'X') for i in range(0, len(nt) - 2, 3))


def load_genome(path):
    g = {}
    for n, s in fasta(path):
        g[n] = s
    return g


def miniprot(pep, wild, threads, outdir):
    gff = os.path.join(outdir, 'miniprot.gff')
    if not os.path.exists(gff) or os.path.getsize(gff) == 0:
        with open(gff, 'w') as o:
            subprocess.run(['miniprot', '-t', str(threads), '--gff', '-u', '--outn', '1',
                            wild, pep], stdout=o, check=True)
    return gff


def parse_miniprot(gff):
    """best alignment per query -> (chrom, strand, [(start,end)], identity)"""
    best, cur = {}, None
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'mRNA':
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            q = a.get('Target', a.get('ID', '')).split()[0]
            cur = dict(chrom=f[0], strand=f[6], cds=[], ident=float(a.get('Identity', 0)),
                       score=float(f[5]) if f[5] not in ('.', '') else 0.0, q=q)
            if q not in best or cur['score'] > best[q]['score']:
                best[q] = cur
            else:
                cur = None
        elif f[2] == 'CDS' and cur is not None:
            cur['cds'].append((int(f[3]), int(f[4])))
    return best


def target_cds(rec, genome):
    seq = ''.join(genome[rec['chrom']][s - 1:e] for s, e in sorted(rec['cds']))
    if rec['strand'] == '-':
        seq = revcomp(seq)
    return seq


def mafft_pair(a, b, tmp):
    p = os.path.join(tmp, 'pair.fa')
    with open(p, 'w') as o:
        o.write(f'>q\n{a}\n>t\n{b}\n')
    out = subprocess.run(['mafft', '--quiet', '--auto', p], capture_output=True, text=True)
    aln = dict(fasta_str(out.stdout))
    return aln.get('q'), aln.get('t')


def fasta_str(text):
    name, seq = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, ''.join(seq)


def back_translate(paln, taln, qnt, tnt):
    """protein alignment + ungapped CDS -> codon alignment (gap columns dropped)"""
    qi = ti = 0
    qo, to = [], []
    for pa, pb in zip(paln, taln):
        if pa != '-' and pb != '-':
            qc, tc = qnt[qi * 3:qi * 3 + 3], tnt[ti * 3:ti * 3 + 3]
            if len(qc) == 3 and len(tc) == 3 and qc.upper() not in CODON_STOPS \
               and tc.upper() not in CODON_STOPS and 'N' not in (qc + tc).upper():
                qo.append(qc)
                to.append(tc)
        if pa != '-':
            qi += 1
        if pb != '-':
            ti += 1
    return ''.join(qo), ''.join(to)


YN_CTL = """      seqfile = pair.phy
      outfile = yn.out
      verbose = 0
        icode = 0
    weighting = 0
   commonf3x4 = 0
"""


def yn00(qc, tc, tmp):
    n = len(qc) // 3
    with open(os.path.join(tmp, 'pair.phy'), 'w') as o:
        o.write(f' 2 {len(qc)}\nq         {qc}\nt         {tc}\n')
    with open(os.path.join(tmp, 'yn00.ctl'), 'w') as o:
        o.write(YN_CTL)
    r = subprocess.run(['yn00', 'yn00.ctl'], cwd=tmp, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    out = os.path.join(tmp, 'yn.out')
    if not os.path.exists(out):
        return None
    lines = open(out).read().splitlines()
    for i, l in enumerate(lines):
        if l.strip().startswith('seq. seq.'):
            for row in lines[i + 1:]:
                p = row.split()
                if len(p) >= 11 and p[0] == '2' and p[1] == '1':
                    try:
                        dnds, dn, ds = float(p[6]), float(p[7]), float(p[10])
                    except ValueError:
                        return None
                    return dn, ds, dnds, n
    return None


def main():
    qcds, qpep, wild, outdir = sys.argv[1:5]
    threads = sys.argv[5] if len(sys.argv) > 5 else '16'
    os.makedirs(outdir, exist_ok=True)
    gff = miniprot(qpep, wild, threads, outdir)
    best = parse_miniprot(gff)
    genome = load_genome(wild)
    cds = dict(fasta(qcds))
    pep = dict(fasta(qpep))
    done = 0
    with open(os.path.join(outdir, 'kaks.tsv'), 'w') as o, tempfile.TemporaryDirectory() as tmp:
        o.write('query\tchrom\tidentity\tdN\tdS\tdNdS\tcodons\n')
        for q, rec in best.items():
            if q not in cds or q not in pep:
                continue
            tnt = target_cds(rec, genome)
            tp = translate(tnt).rstrip('*')
            qp = pep[q].rstrip('*')
            if len(tp) < 30 or len(qp) < 30:
                continue
            pa, pb = mafft_pair(qp, tp, tmp)
            if not pa:
                continue
            qc, tc = back_translate(pa, pb, cds[q], tnt)
            if len(qc) < 90:
                continue
            res = yn00(qc, tc, tmp)
            if not res:
                continue
            dn, ds, dnds, n = res
            o.write(f'{q}\t{rec["chrom"]}\t{rec["ident"]:.3f}\t{dn}\t{ds}\t{dnds}\t{n}\n')
            done += 1
            if done % 500 == 0:
                print(done, 'pairs', flush=True)
    print('wrote', done, 'pairs to', outdir)


if __name__ == '__main__':
    main()
