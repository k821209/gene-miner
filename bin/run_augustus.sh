#!/bin/bash
cd ~/hobak_comparative
AUG=/home/k821209/miniconda3/envs/augustus
export AUGUSTUS_CONFIG_PATH=$AUG/config
export PATH=/home/k821209/miniconda3/envs/annot/bin:$PATH
SP=BUSCO_Bs.genome.polished.fa.busco
rm -f augustus.DONE augustus.FAIL
trap 'echo "[$(date +%T)] FAIL line $LINENO"; touch augustus.FAIL' ERR
set -e
rm -rf aug_split aug_out; mkdir -p aug_split aug_out
echo "[$(date +%T)] split by contig"
awk '/^>/{if(f)close(f); name=substr($1,2); sub(/ .*/,"",name); f="aug_split/"name".fa"} {print > f}' clean_masked_contigs.fa
echo "  split files: $(ls aug_split | wc -l)"
echo "[$(date +%T)] AUGUSTUS parallel (-P 38, BUSCO model)"
ls aug_split/*.fa | xargs -P 38 -I {} bash -c 'b=$(basename "$1" .fa); '"$AUG"'/bin/augustus --strand=both --genemodel=complete --gff3=on --UTR=off --species='"$SP"' --AUGUSTUS_CONFIG_PATH='"$AUG"'/config "$1" > aug_out/"$b".gff3 2>/dev/null' _ {}
echo "[$(date +%T)] merge (prefix IDs by contig)"
: > augustus_contigs.gff3
for g in aug_out/*.gff3; do
  b=$(basename "$g" .gff3)
  awk -v c="$b" '!/^#/ && NF>=8 {gsub(/ID=g/,"ID="c"_g"); gsub(/Parent=g/,"Parent="c"_g"); gsub(/ID=t/,"ID="c"_t"); gsub(/Parent=t/,"Parent="c"_t"); print}' "$g" >> augustus_contigs.gff3
done
NG=$(grep -cP '\tgene\t' augustus_contigs.gff3)
NT=$(grep -cP '\ttranscript\t' augustus_contigs.gff3)
echo "[$(date +%T)] AUGUSTUS ab initio: genes=$NG transcripts=$NT"
touch augustus.DONE
echo "[$(date +%T)] DONE"
