#!/usr/bin/env bash
# Rebuild a genome's catalogue with the new build_union.py and recompute every
# headline number: QC filters, BUSCO (primary + all isoforms), reference
# recovery, gffcompare, split/fusion, isoform provenance.
# Usage: regen.sh <genome...>
set -uo pipefail
R=~/gene-miner-runs; BIN=$R/gene-miner/bin; REV=$R/revision; AN=$R/gene-miner/analysis
source ~/miniconda3/etc/profile.d/conda.sh
declare -A GEN=( [rice]=$R/rice/genome/IRGSP.fa [soybean]=$R/soybean/genome/Wm82.fa [dmel]=$R/dmel/genome/dmel.fa [cele]=$R/cele/genome/cele.fa )
declare -A REF=( [rice]=$R/rice/genome/IRGSP.ref.gff3 [soybean]=$R/soybean/genome/Wm82.ref.gff3 [dmel]=$R/dmel/genome/dmel.ref.gff [cele]=$R/cele/genome/cele.ref.gff )
declare -A RUN=( [rice]=$R/rice/native_etp [soybean]=$R/soybean/gm_etp [dmel]=$R/dmel/gm_etp [cele]=$R/cele/gm_etp )
declare -A GMK=( [rice]=$R/rice/etp_run/genemark.gtf [soybean]=$R/soybean/braker_soy/GeneMark-ETP/genemark.gtf [dmel]=$R/dmel/braker_dmel/GeneMark-ETP/genemark.gtf [cele]=$R/cele/braker_cele/GeneMark-ETP/genemark.gtf )
declare -A LIN=( [rice]=poales_odb10 [soybean]=fabales_odb10 [dmel]=diptera_odb10 [cele]=nematoda_odb10 )
declare -A ANN=( [rice]=$R/rice/native_etp/emap/union.emapper.annotations [soybean]=$R/soybean/gm_etp/combined.annotations [dmel]=$R/dmel/gm_etp/combined.annotations [cele]=$R/cele/gm_etp/combined.annotations )

for g in "$@"; do
  W=$REV/regen/$g; mkdir -p $W/annot; cd $W
  ln -sf ${RUN[$g]}/augustus_scaffold.gff3 augustus_scaffold.gff3
  ln -sf ${RUN[$g]}/annot/genome.transdecoder.gff3 annot/genome.transdecoder.gff3
  echo "[$(date +%T)] $g union"
  python3 $BIN/build_union.py --prefix GM --genemark ${GMK[$g]} > union.log 2>&1 || { echo "$g UNION FAILED"; continue; }
  python3 $BIN/extract_pep.py ${GEN[$g]} union.gff3 union.pep.fa > pep.log 2>&1

  echo "[$(date +%T)] $g TE filter"
  python3 $BIN/te_filter_genes.py $R/$g/gm_out/mask/rm.out union.gff3 union.pep.fa union.te.gff3 union.te.pep.fa > te.log 2>&1

  echo "[$(date +%T)] $g eggNOG transfer"
  python3 $REV/transfer_eggnog.py ${RUN[$g]}/union.pep.fa ${ANN[$g]} union.pep.fa new.annotations todo.faa > transfer.log 2>&1
  if [ -s todo.faa ]; then
    bash $BIN/run_eggnog.sh todo.faa todo emap_new > eggnog.log 2>&1
    grep -v '^#' emap_new/todo.emapper.annotations >> new.annotations 2>/dev/null
  fi
  echo "[$(date +%T)] $g taxonomy filter"
  python3 $BIN/filter_taxonomy.py new.annotations union.te.gff3 union.te.pep.fa union.final.gff3 union.final.pep.fa > tax.log 2>&1

  python3 - <<'PY'
import re
from collections import defaultdict
g2t=defaultdict(list); tl=defaultdict(int)
for l in open("union.final.gff3"):
    f=l.rstrip().split("\t")
    if len(f)<9: continue
    if f[2]=="mRNA": g2t[re.search(r"Parent=([^;]+)",f[8]).group(1)].append(re.search(r"ID=([^;]+)",f[8]).group(1))
    elif f[2]=="CDS": tl[re.search(r"Parent=([^;]+)",f[8]).group(1)]+=int(f[4])-int(f[3])+1
prim={max(ts,key=lambda t:tl[t]) for ts in g2t.values()}
keep=False; o=open("union.final.primary.pep.fa","w")
for l in open("union.final.pep.fa"):
    if l[0]==">": keep=l[1:].split()[0] in prim
    if keep: o.write(l)
PY

  echo "[$(date +%T)] $g comparisons"
  python3 $BIN/gm_compare.py --gm union.final.gff3 --ref ${REF[$g]} --min_ro 0.5 --out gm_compare.tsv 2>/dev/null
  python3 $REV/split_fusion.py ${REF[$g]} union.final.gff3 0.5 > split_fusion.txt
  python3 $REV/isoform_origin.py union.final.gff3 augustus_scaffold.gff3 annot/genome.transdecoder.gff3 ${GMK[$g]} > isoform_origin.txt 2>/dev/null
  python3 $AN/classify_missed.py union.gff3 union.final.gff3 ${REF[$g]} > classify_missed.txt 2>&1
  conda activate gffcmp
  python3 $AN/cds_to_exon_gtf.py union.final.gff3 gm.cds.gtf 2>/dev/null
  gffcompare -r $REV/structural/$g.ref.cds.gtf -o cmp gm.cds.gtf > /dev/null 2>&1
  conda activate busco
  for set in primary alliso; do
    IN=union.final.primary.pep.fa; [ $set = alliso ] && IN=union.final.pep.fa
    busco -i $IN -l ${LIN[$g]} -m proteins -o busco_$set -f --offline \
      --download_path $R/$g/gm_out/busco_downloads -c 20 > busco_$set.log 2>&1
  done
  conda deactivate
  echo "[$(date +%T)] $g DONE  $(grep -h 'C:' busco_primary/short_summary*.txt | tr -d '\t ')"
done
echo REGEN_ALL_DONE
