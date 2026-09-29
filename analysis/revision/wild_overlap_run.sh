#!/usr/bin/env bash
# Do the sampled rice loci land on genes the wild relative already annotates?
set -uo pipefail
R=~/gene-miner-runs; REV=$R/revision; W=$REV/wildoverlap; mkdir -p $W
export PATH=$HOME/miniconda3/envs/busco/bin:$PATH
declare -A GEN=( [glaberrima]=$REV/kaks/genomes/Oryza_glaberrima.Oryza_glaberrima_V1.dna.toplevel.fa
                 [punctata]=$REV/kaks/genomes/Oryza_punctata.Oryza_punctata_v1.2.dna.toplevel.fa )
declare -A ANN=( [glaberrima]=$REV/wildann/Oryza_glaberrima.Oryza_glaberrima_V1.63.chr.gff3.gz
                 [punctata]=$REV/wildann/Oryza_punctata.Oryza_punctata_v1.2.63.gff3.gz )
for sp in glaberrima punctata; do
  [ -s $W/$sp.mpi ] || miniprot -t 12 -d $W/$sp.mpi ${GEN[$sp]} 2> $W/$sp.index.log
  for cls in novel identical revised te_removed; do
    PEP=$REV/kaks_s/seqs/$cls.pep.fa
    [ -s $PEP ] || continue
    [ -s $W/$sp.$cls.gff ] || miniprot -t 12 --gff -u --outn 1 $W/$sp.mpi $PEP > $W/$sp.$cls.gff 2> $W/$sp.$cls.log
    python3 $REV/wild_gene_overlap.py $W/$sp.$cls.gff ${ANN[$sp]} $W/$sp.$cls.tsv 0.5
  done
done
echo WILD_OVERLAP_DONE
