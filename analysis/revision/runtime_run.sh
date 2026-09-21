#!/usr/bin/env bash
# Full end-to-end Nextflow runs with resource tracing, for the runtime/memory table
# (R1.5, R2-minor2). De novo repeat modelling is skipped by supplying the library
# built for the original run (--repeat_lib); its own cost is reported separately
# from the original runs' records.
# Usage: runtime_run.sh <genome...>
set -uo pipefail
R=~/gene-miner-runs; REV=$R/revision
source ~/miniconda3/etc/profile.d/conda.sh; conda activate nf
export NXF_ANSI_LOG=false
export GENEMARK_ETP_DIR=$HOME/gene-miner-runs/GeneMark-ETP
declare -A GEN=( [rice]=$R/rice/genome/IRGSP.fa [soybean]=$R/soybean/genome/Wm82.fa [dmel]=$R/dmel/genome/dmel.fa [cele]=$R/cele/genome/cele.fa )
declare -A SP=( [rice]=rice [soybean]=auto [dmel]=fly [cele]=caenorhabditis )
declare -A LIN=( [rice]=poales_odb10 [soybean]=fabales_odb10 [dmel]=diptera_odb10 [cele]=nematoda_odb10 )
declare -A PFX=( [rice]=RICE [soybean]=SOY [dmel]=DMEL [cele]=CELE )

for g in "$@"; do
  W=$REV/runtime/$g; mkdir -p $W; cd $W
  SPEC=""; [ "${SP[$g]}" != auto ] && SPEC="--augustus_species ${SP[$g]}"
  echo "[$(date +%F' '%T)] $g start"
  nextflow run $R/gene-miner/main.nf -c $R/gene-miner/nextflow.config \
    --genome ${GEN[$g]} \
    --reads "$R/$g/rnaseq/*_{1,2}.fastq.gz" \
    --proteome $R/$g/db/uniprot_sprot.fasta \
    --repeat_lib $REV/repeatlibs/$g.families.fa \
    --busco_lineage ${LIN[$g]} \
    --gene_prefix ${PFX[$g]} \
    --run_genemark true \
    $SPEC \
    --outdir $W/gm_out \
    -work-dir $W/work \
    -with-trace $W/trace.txt -with-report $W/report.html -with-timeline $W/timeline.html \
    > nextflow.log 2>&1
  echo "[$(date +%F' '%T)] $g exit=$? $(grep -c COMPLETED trace.txt 2>/dev/null) tasks"
  # keep the outputs and the trace; drop the bulky work directory
  du -sh work 2>/dev/null
  rm -rf work
done
echo RUNTIME_ALL_DONE
