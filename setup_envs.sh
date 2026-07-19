#!/usr/bin/env bash
# setup_envs.sh — create the conda environments the Gene-Miner Nextflow pipeline
# (main.nf) uses, in one command:  bash setup_envs.sh
#
# Requires conda (miniconda / mambaforge) on PATH. Uses mamba if available
# (much faster to solve); otherwise falls back to conda. The env names below are
# the defaults nextflow.config expects; override the location with
#   export GM_CONDA_BASE=/path/to/miniconda   (default: $HOME/miniconda3)
#
# This installs the five conda envs the pipeline needs for its default two-stream
# run (AUGUSTUS + RNA-seq + QC). The optional 3rd stream (GeneMark-ETP, enabled
# with --run_genemark) is NOT conda-installable on a clean machine: the bioconda
# `braker3` recipe requires genomethreader, which bioconda no longer ships. Use
# the official BRAKER container instead — see the note printed at the end.
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
# optional: `hmmer` in the annot env enables the Pfam ORF-retention step (--pfam)

cat <<'NOTE'

[setup] five conda environments ready (annot, augustus, rmod, eggnog, busco).

Databases are fetched on first use, not by this script:
  - eggNOG DB    : run_eggnog.sh downloads it on first run (or set $EGGNOG_DB to an
                   existing copy; ~50 GB).
  - BUSCO lineage: BUSCO downloads the lineage on first run, or pre-fetch with
                   `busco --download <lineage_odb10>`.
  - Pfam (opt.)  : download Pfam-A.hmm and `hmmpress` it; pass --pfam to use it.

3rd stream — GeneMark-ETP (only if you pass --run_genemark true):
  run_genemark_etp.sh drives GeneMark-ETP through braker.pl, which also needs
  GenomeThreader + several Perl modules. This stack is NOT conda-installable on
  a clean machine (bioconda `braker3` requires genomethreader, no longer
  shipped). Use the official BRAKER container, which bundles braker.pl +
  GeneMark-ETP + GenomeThreader + the Perl deps:
        singularity build braker3.sif docker://teambraker/braker3:latest
  then point run_genemark_etp.sh at it (BRAKER_ENV / GENEMARK_PATH). Cloning
  GeneMark-ETP alone (github.com/gatech-genemark/GeneMark-ETP) does not suffice.
  The two-stream default below needs none of this.

To run the bin/ scripts by hand (main.nf does this for you), prepend the tool's
env to PATH, e.g.  export PATH=$GM_CONDA_BASE/envs/annot/bin:$PATH

nextflow.config resolves the env prefixes from $HOME/miniconda3/envs by default
(override with GM_CONDA_BASE). Two-stream run (works with the five envs above):

  nextflow run main.nf -c nextflow.config \
    --genome genome.fa --reads 'rnaseq/*_{1,2}.fastq.gz' \
    --proteome db/uniprot_sprot.fasta --augustus_species rice \
    --busco_lineage poales_odb10 --outdir gm_out
NOTE
