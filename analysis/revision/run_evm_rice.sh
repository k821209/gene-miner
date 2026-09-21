#!/usr/bin/env bash
# EVidenceModeler on rice, from the same evidence Gene-Miner unions (R2-M4):
# AUGUSTUS + GeneMark-ETP predictions, StringTie transcripts, Swiss-Prot miniprot
# alignments. Weighted consensus vs detection-first union on identical inputs.
set -uo pipefail
R=~/gene-miner-runs; REV=$R/revision; W=$REV/evm/rice
E=$HOME/miniconda3/envs/augustus/opt/evidencemodeler-2.1.0
GENOME=$R/rice/genome/IRGSP.fa
MERGED=$R/rice/work/9e/af1e7a1500d63d8212d2f06e697fe4/merged.gtf
AUG=$REV/rice_aug_sample/augustus_scaffold.gff3
GMK=$R/rice/etp_run/genemark.gtf
SPROT=$R/rice/db/uniprot_sprot.fasta
mkdir -p $W; cd $W
source ~/miniconda3/etc/profile.d/conda.sh; conda activate augustus
export PATH=$HOME/miniconda3/envs/augustus/bin:$PATH

# 1. gene predictions -> EVM gff3 (gene/mRNA/exon/CDS, one source per stream)
python3 - "$AUG" augustus.evm.gff3 AUGUSTUS <<'PY'
import re, sys
inp, out, src = sys.argv[1:4]
rows = []
for l in open(inp):
    if l.startswith('#') or not l.strip():
        continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'mRNA', 'CDS'):
        continue
    f[1] = src
    rows.append(f)
with open(out, 'w') as o:
    for f in rows:
        o.write('\t'.join(f) + '\n')
        if f[2] == 'CDS':
            e = list(f); e[2] = 'exon'; e[7] = '.'
            e[8] = re.sub(r'ID=([^;]+)', r'ID=\1.exon', e[8])
            o.write('\t'.join(e) + '\n')
PY
python3 - "$GMK" genemark.evm.gff3 GeneMarkETP <<'PY'
import re, sys
from collections import defaultdict
inp, out, src = sys.argv[1:4]
tx = defaultdict(list); info = {}; gene_of = {}
for l in open(inp):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'CDS':
        continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    st = re.search(r'status "([^"]+)"', f[8])
    if st and st.group(1) != 'complete':
        continue
    tx[t].append((int(f[3]), int(f[4]), f[7] if f[7] in '012' else '0'))
    info[t] = (f[0], f[6]); gene_of[t] = g
with open(out, 'w') as o:
    genes = defaultdict(list)
    for t in tx:
        genes[gene_of[t]].append(t)
    for g, ts in genes.items():
        blocks = [b for t in ts for b in tx[t]]
        c, strand = info[ts[0]]
        o.write(f"{c}\t{src}\tgene\t{min(b[0] for b in blocks)}\t{max(b[1] for b in blocks)}\t.\t{strand}\t.\tID={g}\n")
        for t in ts:
            b = sorted(tx[t])
            o.write(f"{c}\t{src}\tmRNA\t{b[0][0]}\t{b[-1][1]}\t.\t{strand}\t.\tID={t};Parent={g}\n")
            for s, e, ph in b:
                o.write(f"{c}\t{src}\texon\t{s}\t{e}\t.\t{strand}\t.\tID={t}.ex{s};Parent={t}\n")
                o.write(f"{c}\t{src}\tCDS\t{s}\t{e}\t.\t{strand}\t{ph}\tID={t}.cds{s};Parent={t}\n")
PY
cat augustus.evm.gff3 genemark.evm.gff3 > predictions.gff3

# 2. StringTie transcripts -> EVM transcript alignments (cDNA_match)
python3 - "$MERGED" transcripts.gff3 <<'PY'
import re, sys
from collections import defaultdict
inp, out = sys.argv[1:3]
ex = defaultdict(list); info = {}
for l in open(inp):
    if l.startswith('#'):
        continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    ex[t].append((int(f[3]), int(f[4]))); info[t] = (f[0], f[6])
with open(out, 'w') as o:
    for t, blocks in ex.items():
        c, strand = info[t]
        for s, e in sorted(blocks):
            o.write(f"{c}\tStringTie\tcDNA_match\t{s}\t{e}\t.\t{strand}\t.\tID={t};Target={t} 1 {e-s+1} +\n")
PY

# 3. Swiss-Prot protein alignments (miniprot) -> EVM protein alignments
if [ ! -s miniprot.evm.gff3 ]; then
  miniprot -t 20 --gff -u --outn 1 $GENOME $SPROT > miniprot.gff 2> miniprot.log
  python3 - miniprot.gff miniprot.evm.gff3 <<'PY'
import re, sys
inp, out = sys.argv[1:3]
cur = None
with open(out, 'w') as o:
    for l in open(inp):
        if l.startswith('#'):
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'mRNA':
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            cur = a.get('Target', a.get('ID', 'prot')).split()[0]
        elif f[2] == 'CDS' and cur:
            o.write(f"{f[0]}\tminiprot_protAln\tnucleotide_to_protein_match\t{f[3]}\t{f[4]}\t{f[5]}\t{f[6]}\t.\tID={cur};Target={cur} 1 {int(f[4])-int(f[3])+1} +\n")
PY
fi

# 4. EVM
printf "ABINITIO_PREDICTION\tAUGUSTUS\t1\nABINITIO_PREDICTION\tGeneMarkETP\t1\nTRANSCRIPT\tStringTie\t10\nPROTEIN\tminiprot_protAln\t5\n" > weights.txt
$E/EVidenceModeler --sample_id rice_evm \
  --genome $GENOME --weights $PWD/weights.txt \
  --gene_predictions predictions.gff3 \
  --transcript_alignments transcripts.gff3 \
  --protein_alignments miniprot.evm.gff3 \
  --segmentSize 1000000 --overlapSize 100000 --CPU 20 > evm.run.log 2>&1
echo "[$(date +%T)] EVM genes=$(awk -F'\t' '$3=="gene"' rice_evm.EVM.gff3 2>/dev/null | wc -l)"

# 5. same measurements as the Gene-Miner catalogue
python3 $R/gene-miner/bin/gm_compare.py --gm rice_evm.EVM.gff3 --ref $R/rice/genome/IRGSP.ref.gff3 --min_ro 0.5 --out gm_compare.tsv 2>/dev/null
python3 $REV/split_fusion.py $R/rice/genome/IRGSP.ref.gff3 rice_evm.EVM.gff3 0.5 > split_fusion.txt
conda activate gffcmp
python3 $R/gene-miner/analysis/cds_to_exon_gtf.py rice_evm.EVM.gff3 evm.cds.gtf 2>/dev/null
gffcompare -r $REV/structural/rice.ref.cds.gtf -o cmp evm.cds.gtf > /dev/null 2>&1
grep -E "level:" cmp.stats
echo EVM_RICE_DONE
