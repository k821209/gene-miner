#!/bin/bash
cd ~/hobak_comparative
export PATH=/home/k821209/miniconda3/envs/augustus/bin:$PATH
E=/home/k821209/miniconda3/envs/augustus/opt/evidencemodeler-2.1.0
rm -f evm3.DONE evm3.FAIL
trap 'echo "[$(date +%T)] FAIL line $LINENO"; touch evm3.FAIL' ERR
set -e
printf "ABINITIO_PREDICTION\tAUGUSTUS\t1\nTRANSCRIPT\tCufflinks\t10\nPROTEIN\tminiprot_protAln\t5\n" > weights3.txt
rm -rf Zsc3.partitions Zsc3.EVM.gff3
echo "[$(date +%T)] EVM (3 evidence: augustus+transcript+protein)"
$E/EVidenceModeler --sample_id Zsc3 \
  --genome ragtag_out/ragtag.scaffold.fasta \
  --weights weights3.txt \
  --gene_predictions augustus_evm.gff3 \
  --transcript_alignments annot/merged.gff3 \
  --protein_alignments miniprot.evm.gff3 \
  --segmentSize 1000000 --overlapSize 100000 --CPU 36 > evm3.run.out 2>&1
OUT=$(ls Zsc3.EVM.gff3 2>/dev/null)
echo "[$(date +%T)] EVM3 genes=$(awk -F'\t' '$3=="gene"' $OUT | wc -l) mRNAs=$(awk -F'\t' '$3=="mRNA"' $OUT | wc -l)"
touch evm3.DONE
