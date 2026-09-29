#!/usr/bin/env bash
# Coding-level gffcompare (R1.1) + gene split/fusion counts (R2-M3) for
# Gene-Miner and BRAKER3 on the four benchmark genomes.
set -euo pipefail
R=~/gene-miner-runs
OUT=$R/revision/structural; mkdir -p $OUT; cd $OUT
source ~/miniconda3/etc/profile.d/conda.sh; conda activate gffcmp
CONV=$R/gene-miner/analysis/cds_to_exon_gtf.py
SF=$R/revision/split_fusion.py

declare -A REF=( [rice]=$R/rice/genome/IRGSP.ref.gff3 [soybean]=$R/soybean/genome/Wm82.ref.gff3
                 [dmel]=$R/dmel/genome/dmel.ref.gff [cele]=$R/cele/genome/cele.ref.gff )
declare -A GM=( [rice]=$R/rice/native_etp/union.final.gff3 [soybean]=$R/soybean/gm_etp/union.final.gff3
                [dmel]=$R/dmel/gm_etp/union.final.gff3 [cele]=$R/cele/gm_etp/union.final.gff3 )
declare -A BR=( [rice]=$R/rice/etp_run/braker_full/braker.gff3 [soybean]=$R/soybean/braker_soy/braker.gff3
                [dmel]=$R/dmel/braker_dmel/braker.gff3 [cele]=$R/cele/braker_cele/braker.gff3 )

echo -e "genome\tpipeline\tref_genes\tqry_genes\tF\tsplit_ref\tsplit_ref_pct\tsplit_qry_genes\tfusion_qry\tfusion_qry_pct\tfused_ref\tfused_ref_pct" > split_fusion.tsv
for g in rice soybean dmel cele; do
  python3 $CONV ${REF[$g]} $g.ref.cds.gtf
  for p in GM BR; do
    eval "q=\${$p[$g]}"
    python3 $CONV $q $g.$p.cds.gtf
    gffcompare -r $g.ref.cds.gtf -o $g.$p $g.$p.cds.gtf
    echo -e "$g\t$p\t$(python3 $SF ${REF[$g]} $q 0.5)" >> split_fusion.tsv
  done
done
echo done
