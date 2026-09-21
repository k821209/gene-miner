#!/usr/bin/env bash
# EVidenceModeler on rice (R2-M4), on h100 with 32 threads (shared host).
# Same evidence Gene-Miner unions: AUGUSTUS (sampled, native model) + GeneMark-ETP
# predictions, StringTie transcripts, Swiss-Prot miniprot alignments.
# Then the same measurements as the Gene-Miner and BRAKER3 catalogues.
set -uo pipefail
B=/data/k821209/gm_evm
export PATH=$B/env/bin:$PATH
E=$B/env/opt/evidencemodeler-2.1.0
W=$B/run; mkdir -p $W; cd $W
GENOME=$B/rice/genome/IRGSP.fa; REF=$B/rice/genome/IRGSP.ref.gff3
T=32
log() { echo "[$(date +%T)] $*"; }

log "miniprot: Swiss-Prot onto the rice genome"
if [ ! -s miniprot.evm.gff3 ]; then
  miniprot -t $T -d genome.mpi $GENOME > /dev/null 2> miniprot_index.log
  miniprot -t $T --gff -u --outn 1 genome.mpi $B/uniprot_sprot.fasta > miniprot.gff 2> miniprot.log
  python3 - miniprot.gff miniprot.evm.gff3 <<'PY'
import sys
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
            cur = a.get('Target', a.get('ID', 'prot')).split()[0] + '_' + a.get('ID', '')
        elif f[2] == 'CDS' and cur:
            o.write(f"{f[0]}\tminiprot_protAln\tnucleotide_to_protein_match\t{f[3]}\t{f[4]}\t{f[5]}\t{f[6]}\t.\tID={cur};Target={cur} 1 {int(f[4])-int(f[3])+1} +\n")
PY
fi
log "protein alignments: $(cut -f9 miniprot.evm.gff3 | cut -d';' -f1 | sort -u | wc -l)"

printf "ABINITIO_PREDICTION\tAUGUSTUS\t1\nABINITIO_PREDICTION\tGeneMarkETP\t1\nTRANSCRIPT\tStringTie\t10\nPROTEIN\tminiprot_protAln\t5\n" > weights.txt
log "EVidenceModeler"
$E/EVidenceModeler --sample_id rice_evm --genome $GENOME --weights $W/weights.txt \
  --gene_predictions $B/revision/evm/rice/predictions.gff3 \
  --transcript_alignments $B/revision/evm/rice/transcripts.gff3 \
  --protein_alignments $W/miniprot.evm.gff3 \
  --segmentSize 1000000 --overlapSize 100000 --CPU $T > evm.run.log 2>&1
log "EVM genes: $(awk -F'\t' '$3=="gene"' rice_evm.EVM.gff3 2>/dev/null | wc -l)"

log "measurements"
python3 $B/gene-miner/bin/extract_pep.py $GENOME rice_evm.EVM.gff3 evm.pep.fa > pep.log 2>&1
cp $B/gene-miner/analysis/gm_overlap.py $B/gene-miner/bin/ 2>/dev/null
PYTHONPATH=$B/gene-miner/bin python3 $B/gene-miner/bin/gm_compare.py --gm rice_evm.EVM.gff3 --ref $REF --min_ro 0.5 --out gm_compare.tsv 2>/dev/null
python3 $B/revision/split_fusion.py $REF rice_evm.EVM.gff3 0.5 > split_fusion.txt
python3 $B/gene-miner/analysis/cds_to_exon_gtf.py rice_evm.EVM.gff3 evm.cds.gtf 2>/dev/null
gffcompare -r $B/revision/structural/rice.ref.cds.gtf -o cmp evm.cds.gtf > /dev/null 2>&1
cat gm_compare.tsv split_fusion.txt; grep -E "level:" cmp.stats
echo EVM_H100_DONE
