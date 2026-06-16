#!/bin/bash
cd ~/hobak_comparative
export PATH=/home/k821209/miniconda3/envs/augustus/bin:$PATH
rm -f miniprot.DONE miniprot.FAIL
trap 'touch miniprot.FAIL' ERR
set -e
echo "[$(date +%T)] miniprot align (24215 proteins -> genome)"
miniprot -t 36 --gff ragtag_out/ragtag.scaffold.fasta db/zcuc_protein.faa > miniprot.gff 2> miniprot.err
echo "[$(date +%T)] miniprot done: mRNA=$(grep -cP '\tmRNA\t' miniprot.gff)"
touch miniprot.DONE
