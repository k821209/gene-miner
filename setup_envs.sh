#!/usr/bin/env bash
# setup_envs.sh — create the conda environments the Gene-Miner Nextflow pipeline
# (main.nf) uses, in one command:  bash setup_envs.sh
#
# Requires conda (miniconda / mambaforge) on PATH. Uses mamba if available
# (much faster to solve); otherwise falls back to conda. The env names below are
# the defaults nextflow.config expects; override the location with
#   export GM_CONDA_BASE=/path/to/miniconda   (default: $HOME/miniconda3)
#
# GeneMark-ETP is bundled with the bioconda `braker3` package. If your braker3
# build does not include it, install GeneMark-ETP separately and point
# run_genemark_etp.sh at it via  export GENEMARK_PATH=/path/to/GeneMark-ETP/bin
set -euo pipefail

CH="-c bioconda -c conda-forge"
SOLVER=$(command -v mamba || command -v conda)
[ -n "$SOLVER" ] || { echo "ERROR: conda/mamba not found on PATH"; exit 1; }
echo "[setup] using: $SOLVER"

create () { echo "[setup] creating env '$1'"; "$SOLVER" create -y -n "$1" $CH "${@:2}"; }

create annot     hisat2 stringtie transdecoder samtools gffread
create augustus  augustus diamond
create rmod      repeatmasker repeatmodeler
create eggnog    eggnog-mapper diamond
create busco     busco
create braker3   braker3          # bundles GeneMark-ETP
# optional: `hmmer` in the annot env enables the Pfam ORF-retention step (--pfam)

cat <<'NOTE'

[setup] conda environments ready.

Databases are fetched on first use, not by this script:
  - eggNOG DB    : run_eggnog.sh downloads it on first run (or set $EGGNOG_DB to an
                   existing copy; ~50 GB).
  - BUSCO lineage: BUSCO downloads the lineage on first run, or pre-fetch with
                   `busco --download <lineage_odb10>`.
  - Pfam (opt.)  : download Pfam-A.hmm and `hmmpress` it; pass --pfam to use it.

nextflow.config resolves the env prefixes from $HOME/miniconda3/envs by default
(override with GM_CONDA_BASE). You are ready to run, e.g.:

  nextflow run main.nf -c nextflow.config -profile rice \
    --genome genome.fa --reads 'rnaseq/*_{1,2}.fastq.gz' \
    --proteome db/uniprot_sprot.fasta --outdir gm_out
NOTE
