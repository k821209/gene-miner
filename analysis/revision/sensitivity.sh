#!/usr/bin/env bash
# Threshold sensitivity for the two parameters reviewers asked about:
#   (a) the AUGUSTUS score cut-off for "usable" models  (R1.4)
#   (b) the TE-overlap fraction at which a model is dropped (R2-M5)
# Each setting is rebuilt from the same three streams and measured the same way.
set -uo pipefail
R=~/gene-miner-runs; BIN=$R/gene-miner/bin; REV=$R/revision; AN=$R/gene-miner/analysis
source ~/miniconda3/etc/profile.d/conda.sh
declare -A GEN=( [rice]=$R/rice/genome/IRGSP.fa [soybean]=$R/soybean/genome/Wm82.fa [dmel]=$R/dmel/genome/dmel.fa [cele]=$R/cele/genome/cele.fa )
declare -A REF=( [rice]=$R/rice/genome/IRGSP.ref.gff3 [soybean]=$R/soybean/genome/Wm82.ref.gff3 [dmel]=$R/dmel/genome/dmel.ref.gff [cele]=$R/cele/genome/cele.ref.gff )
declare -A RUN=( [rice]=$REV/regen/rice_s [soybean]=$R/soybean/gm_etp [dmel]=$R/dmel/gm_etp [cele]=$R/cele/gm_etp )
declare -A GMK=( [rice]=$R/rice/etp_run/genemark.gtf [soybean]=$R/soybean/braker_soy/GeneMark-ETP/genemark.gtf [dmel]=$R/dmel/braker_dmel/GeneMark-ETP/genemark.gtf [cele]=$R/cele/braker_cele/GeneMark-ETP/genemark.gtf )
declare -A LIN=( [rice]=poales_odb10 [soybean]=fabales_odb10 [dmel]=diptera_odb10 [cele]=nematoda_odb10 )

primary() {  # union.final.gff3 + union.final.pep.fa -> union.final.primary.pep.fa
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
}

measure() {  # genome tag -> one TSV line
  local g=$1 tag=$2
  primary
  conda activate busco
  busco -i union.final.primary.pep.fa -l ${LIN[$g]} -m proteins -o busco_primary -f --offline \
    --download_path $R/$g/gm_out/busco_downloads -c 20 > busco.log 2>&1
  local C=$(grep -ho "C:[0-9.]*%" busco_primary/short_summary*.txt | head -1)
  conda activate base
  python3 $BIN/gm_compare.py --gm union.final.gff3 --ref ${REF[$g]} --min_ro 0.5 --out gm_compare.tsv 2>/dev/null
  local rec=$(grep -P "^gm_recovered" gm_compare.tsv | cut -f2)
  local nov=$(grep -P "^gm_novel" gm_compare.tsv | cut -f2)
  local mis=$(grep -P "^reference_loci_missed" gm_compare.tsv | cut -f2)
  local sf=$(python3 $REV/split_fusion.py ${REF[$g]} union.final.gff3 0.5 | cut -f5,8)
  echo -e "$g\t$tag\t$(awk -F'\t' '$3=="gene"' union.final.gff3 | wc -l)\t$C\t$rec\t$nov\t$mis\t$sf"
}

# (a) AUGUSTUS score cut-off ---------------------------------------------------
OUT=$REV/sensitivity; mkdir -p $OUT
mv $OUT/aug_score.tsv $OUT/aug_score.unsampled_rice.tsv 2>/dev/null
echo -e "genome\tsetting\tgenes\tBUSCO\trecovered_loci\tnovel\tref_missed\tsplit%\tfusion%" > $OUT/aug_score.tsv
for g in cele soybean rice; do
  for sc in 0.0 0.5 0.8 0.95 1.0; do
    W=$OUT/aug_$sc/$g; mkdir -p $W/annot; cd $W
    ln -sf ${RUN[$g]}/augustus_scaffold.gff3 augustus_scaffold.gff3
    ln -sf ${RUN[$g]}/annot/genome.transdecoder.gff3 annot/genome.transdecoder.gff3
    conda activate base
    AUG_SCORE=$sc python3 - $sc <<'PY' > union.log 2>&1
import re, subprocess, sys, pathlib
src = pathlib.Path.home()/'gene-miner-runs/gene-miner/bin/build_union.py'
code = src.read_text().replace('SC, MINAA = 0.8, 100', f'SC, MINAA = {sys.argv[1]}, 100')
pathlib.Path('bu.py').write_text(code)
PY
    python3 bu.py --prefix GM --genemark ${GMK[$g]} >> union.log 2>&1
    python3 $BIN/extract_pep.py ${GEN[$g]} union.gff3 union.pep.fa > /dev/null 2>&1
    python3 $BIN/te_filter_genes.py $R/$g/gm_out/mask/rm.out union.gff3 union.pep.fa union.final.gff3 union.final.pep.fa > te.log 2>&1
    measure $g "aug_score=$sc" >> $OUT/aug_score.tsv
  done
done

# (b) TE-overlap threshold -----------------------------------------------------
# cele TE-threshold rows were computed before the AUGUSTUS change and are unaffected by it
for g in rice; do
  for te in 0.2 0.35 0.5 0.7 0.9; do
    W=$OUT/te_$te/$g; mkdir -p $W; cd $W
    cp $REV/regen/${g}_s/union.gff3 $REV/regen/${g}_s/union.pep.fa . 2>/dev/null || cp $REV/regen/$g/union.gff3 $REV/regen/$g/union.pep.fa . 2>/dev/null
    conda activate base
    python3 $BIN/te_filter_genes.py $R/$g/gm_out/mask/rm.out union.gff3 union.pep.fa union.final.gff3 union.final.pep.fa $te > te.log 2>&1
    measure $g "te_thresh=$te" >> $OUT/te_thresh.tsv
  done
done
column -t $OUT/aug_score.tsv $OUT/te_thresh.tsv
echo SENSITIVITY_DONE
