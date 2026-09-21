#!/usr/bin/env bash
# Rice RNA-seq subsampling (R1.5): how completeness, structural accuracy and
# compute cost change with the amount of RNA-seq evidence. Two axes are varied so
# that read count and tissue count can be told apart:
#   L4-D100  4 libraries, all reads   (~118 M)   <- the published input, re-run here
#                                                so every cost figure comes from one machine
#   L2-D100  2 libraries, all reads   (~59 M)
#   L1-D100  1 library,   all reads   (~30 M)
#   L4-D50   4 libraries, half reads  (~59 M)    pairs with L2-D100
#   L4-D25   4 libraries, quarter     (~30 M)    pairs with L1-D100
# Masking is skipped with the repeat library from the full run; everything else
# (HISAT2, StringTie, TransDecoder, GeneMark-ETP, eggNOG, BUSCO) is re-run.
# Usage: subsample_rice.sh [config...]   (default: all but L4-D100)
set -uo pipefail
B=/data2/k821209/gm
export PATH=$HOME/miniconda3/envs/seqtk/bin:$HOME/miniconda3/envs/nf/bin:$HOME/miniconda3/bin:$B/bin:$PATH
source $HOME/miniconda3/etc/profile.d/conda.sh
export EGGNOG_DB=$HOME/miniconda3/envs/eggnog/lib/python3.11/site-packages/data
export GENEMARK_ETP_DIR=$HOME/miniconda3/opt/GeneMark-ETP
LIBS_ALL="leaf panicle seedling endosperm"
LIBS_2="leaf panicle"
LIBS_1="leaf"

prep_reads() {   # $1=config $2=libs $3=fraction
  local cfg=$1 libs=$2 frac=$3 d=$B/sub/$cfg/rnaseq
  mkdir -p $d
  for s in $libs; do
    for m in 1 2; do
      local src=$B/rice/rnaseq/${s}_${m}.fastq.gz dst=$d/${s}_${m}.fastq.gz
      [ -s "$dst" ] && continue
      if [ "$frac" = "1.0" ]; then
        ln -sf $src $dst
      else
        seqtk sample -s100 $src $frac | gzip -1 > $dst
      fi
    done
  done
}

run_cfg() {      # $1=config
  local cfg=$1 W=$B/sub/$cfg
  mkdir -p $W; cd $W
  nextflow run $B/gene-miner/main.nf -c $B/nextflow.226.config \
    --genome $B/rice/genome/IRGSP.fa \
    --reads "$W/rnaseq/*_{1,2}.fastq.gz" \
    --proteome $B/rice/db/uniprot_sprot.fasta \
    --repeat_lib $B/revision/repeatlibs/rice.families.fa \
    --augustus_species rice --busco_lineage poales_odb10 \
    --gene_prefix RICE --run_genemark true \
    --outdir $W/gm_out -work-dir $W/work \
    -with-trace $W/trace.txt -with-report $W/report.html -with-timeline $W/timeline.html \
    > $W/nextflow.log 2>&1
  local rc=$?
  echo "[$(date +%F' '%T)] $cfg nextflow rc=$rc"
  if [ -s $W/gm_out/union.final.gff3 ]; then
    python3 $B/gene-miner/bin/gm_compare.py --gm $W/gm_out/union.final.gff3 \
      --ref $B/rice/genome/IRGSP.ref.gff3 --min_ro 0.5 --out $W/gm_compare.tsv 2>/dev/null
    python3 $B/split_fusion.py $B/rice/genome/IRGSP.ref.gff3 $W/gm_out/union.final.gff3 0.5 > $W/split_fusion.txt
    python3 $B/cds_to_exon_gtf.py $W/gm_out/union.final.gff3 $W/gm.cds.gtf 2>/dev/null
    PATH=$HOME/miniconda3/envs/gffcmp/bin:$PATH gffcompare -r $B/rice.ref.cds.gtf -o $W/cmp $W/gm.cds.gtf > /dev/null 2>&1
    zcat $W/rnaseq/*_1.fastq.gz | awk 'END{print "read_pairs\t" NR/4}' > $W/reads.txt
  fi
  # drop intermediates only after a clean run, so a failed run can be resumed
  [ $rc -eq 0 ] && [ -s $W/gm_out/union.final.gff3 ] && rm -rf $W/work
}

CFGS=${*:-"L4-D100 L2-D100 L1-D100 L4-D50 L4-D25"}
for cfg in $CFGS; do
  case $cfg in
    L4-D100) prep_reads $cfg "$LIBS_ALL" 1.0 ;;
    L2-D100) prep_reads $cfg "$LIBS_2"   1.0 ;;
    L1-D100) prep_reads $cfg "$LIBS_1"   1.0 ;;
    L4-D50)  prep_reads $cfg "$LIBS_ALL" 0.5 ;;
    L4-D25)  prep_reads $cfg "$LIBS_ALL" 0.25 ;;
    *) echo "unknown config $cfg"; continue ;;
  esac
  echo "[$(date +%F' '%T)] $cfg start"
  run_cfg $cfg
done
echo SUBSAMPLE_DONE
