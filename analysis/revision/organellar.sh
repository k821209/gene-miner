#!/usr/bin/env bash
# NUMT/NUPT screen of the final catalogues against each species' own organelle genomes.
set -uo pipefail
R=~/gene-miner-runs; REV=$R/revision; W=$REV/contamination; cd $W
source ~/miniconda3/etc/profile.d/conda.sh; conda activate busco
declare -A GEN=( [rice]=$R/rice/genome/IRGSP.fa [soybean]=$R/soybean/genome/Wm82.fa [dmel]=$R/dmel/genome/dmel.fa [cele]=$R/cele/genome/cele.fa )
declare -A REF=( [rice]=$R/rice/genome/IRGSP.ref.gff3 [soybean]=$R/soybean/genome/Wm82.ref.gff3 [dmel]=$R/dmel/genome/dmel.ref.gff [cele]=$R/cele/genome/cele.ref.gff )
declare -A CAT=( [rice]=$REV/regen/rice_s [soybean]=$REV/regen/soybean [dmel]=$REV/regen/dmel [cele]=$REV/regen/cele )
echo -e "genome\tclass\tloci_tested\tloci_50pct_CDS_on_organelle\tpct" > organellar.tsv
for g in rice soybean dmel cele; do
  python3 $REV/extract_classes.py ${GEN[$g]} ${CAT[$g]}/union.final.gff3 ${REF[$g]} seqs_$g 0 > seqs_$g.log 2>&1
  for cls in novel identical revised; do
    q=seqs_$g/$cls.cds.fa; [ -s $q ] || continue
    blastn -query $q -db org/$g.db -outfmt "6 qseqid qlen length pident" -evalue 1e-10 -num_threads 1 -max_target_seqs 5 > $g.$cls.blast 2>/dev/null
    python3 - $q $g.$cls.blast $g $cls >> organellar.tsv <<"PY"
import sys
from collections import defaultdict
q, b, g, cls = sys.argv[1:5]
n = sum(1 for l in open(q) if l.startswith(">"))
cov = defaultdict(int); ql = {}
for line in open(b):
    f = line.split("\t"); cov[f[0]] += int(f[2]); ql[f[0]] = int(f[1])
hit = sum(1 for k, v in cov.items() if v >= 0.5 * ql[k])
print(f"{g}\t{cls}\t{n}\t{hit}\t{100*hit/max(n,1):.2f}")
PY
  done
done
column -t organellar.tsv
echo ORGANELLAR_DONE
