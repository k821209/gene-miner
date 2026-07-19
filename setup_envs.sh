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
# run (AUGUSTUS + RNA-seq + QC), plus a sixth `genemark` env and a GeneMark-ETP
# checkout for the optional 3rd stream (--run_genemark). GeneMark-ETP is driven
# directly through its own `gmetp.pl` (NOT through BRAKER): it bundles GeneMark +
# ProtHint and ships static binaries of every third-party tool it needs
# (bedtools, samtools, hisat2, diamond, stringtie, gffread), so it needs NO
# GenomeThreader and NO container — only Perl (a few CPAN modules) + python3.
# Skip the 3rd-stream step with:  GM_SKIP_GENEMARK=1 bash setup_envs.sh
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

# --- RepeatMasker Dfam FamDB ---
# The bioconda repeatmasker package ships NO Dfam library, and RepeatMasker >=4.1.5
# aborts at startup ("FamDB data directory not found") even with --repeat_lib.
# Fetch the small Dfam root partition (~60 MB) and point famdb.conf at it. The
# FamDB format must match the RepeatMasker version (Dfam 4.0 = famdb 3.x = RM 4.2.x).
echo "[setup] configuring RepeatMasker Dfam library (~60 MB download)"
# Derive the conda base from $SOLVER's own location (works for conda or mamba,
# and avoids `mamba info --base`, which prints a labelled line rather than a path).
CB="$(dirname "$(dirname "$SOLVER")")"
FAMDIR="$CB/envs/rmod/share/RepeatMasker/Libraries/famdb"
mkdir -p "$FAMDIR"
if ! ls "$FAMDIR"/*.h5 >/dev/null 2>&1; then
  ( cd "$FAMDIR" \
    && curl -fSL -o dfam40.0.h5.gz https://www.dfam.org/releases/current/families/FamDB/dfam40.0.h5.gz \
    && gunzip -f dfam40.0.h5.gz ) \
    || echo "[setup] WARN: Dfam download failed — fetch dfam40.0.h5 into $FAMDIR manually"
fi
CONF="$(ls "$CB"/envs/rmod/share/famdb-*/famdb.conf 2>/dev/null | head -1 || true)"
if [ -n "$CONF" ]; then
  if grep -q '^FAMDB_DATA_DIR' "$CONF"; then
    sed -i "s|^FAMDB_DATA_DIR.*|FAMDB_DATA_DIR = $FAMDIR|" "$CONF"
  else
    printf 'FAMDB_DATA_DIR = %s\n' "$FAMDIR" >> "$CONF"
  fi
fi

# --- 3rd stream: GeneMark-ETP (conda-only, no BRAKER, no GenomeThreader) ---
# GeneMark-ETP's gmetp.pl bundles GeneMark + ProtHint and its own static
# bedtools/samtools/hisat2/diamond/stringtie/gffread under tools/, so the only
# external needs are Perl (+ a few CPAN modules) and python3. We create a
# `genemark` env for those and clone the GeneMark-ETP repo next to the envs.
# GeneMark-ETP is CC BY-NC-SA (academic / non-commercial; no licence key).
if [ "${GM_SKIP_GENEMARK:-0}" = "1" ]; then
  echo "[setup] GM_SKIP_GENEMARK=1 -> skipping the GeneMark-ETP 3rd stream"
else
  echo "[setup] setting up the GeneMark-ETP 3rd stream (env 'genemark' + repo clone)"
  # perl-app-cpanminus + make let us add the two CPAN modules bioconda lacks
  # (Statistics::LineFit, Math::Utils); the rest come from conda.
  # perl-yaml (YAML.pm) AND perl-yaml-libyaml (YAML::XS.pm) are both needed —
  # gmetp.pl uses YAML::XS, ProtHint's proteins_from_gtf.pl uses plain YAML.
  create genemark perl perl-yaml perl-yaml-libyaml perl-parallel-forkmanager \
                  perl-hash-merge perl-mce perl-app-cpanminus perl-list-moreutils \
                  perl-scalar-list-utils python make
  "$CB/envs/genemark/bin/cpanm" --notest Math::Utils Statistics::LineFit \
    || echo "[setup] WARN: cpanm of Math::Utils/Statistics::LineFit failed — install them by hand into envs/genemark"
  GM_ETP_DIR="${GENEMARK_ETP_DIR:-$CB/opt/GeneMark-ETP}"
  if [ ! -x "$GM_ETP_DIR/bin/gmetp.pl" ]; then
    mkdir -p "$(dirname "$GM_ETP_DIR")"
    git clone --depth 1 https://github.com/gatech-genemark/GeneMark-ETP "$GM_ETP_DIR" \
      || echo "[setup] WARN: git clone of GeneMark-ETP failed — clone it into $GM_ETP_DIR by hand"
  fi
fi

cat <<'NOTE'

[setup] five conda environments ready (annot, augustus, rmod, eggnog, busco).

Databases are fetched on first use, not by this script:
  - eggNOG DB    : run_eggnog.sh downloads it on first run (or set $EGGNOG_DB to an
                   existing copy; ~50 GB).
  - BUSCO lineage: BUSCO downloads the lineage on first run, or pre-fetch with
                   `busco --download <lineage_odb10>`.
  - Pfam (opt.)  : download Pfam-A.hmm and `hmmpress` it; pass --pfam to use it.

3rd stream — GeneMark-ETP (only if you pass --run_genemark true):
  Installed above unless GM_SKIP_GENEMARK=1: the `genemark` conda env (Perl +
  the required CPAN modules + python3) and a GeneMark-ETP checkout in
  <conda_base>/opt/GeneMark-ETP. run_genemark_etp.sh calls GeneMark-ETP's own
  gmetp.pl directly — NO BRAKER, NO GenomeThreader, NO container (GeneMark-ETP
  bundles GeneMark + ProtHint + static bedtools/samtools/hisat2/diamond/
  stringtie). Override the checkout location with GENEMARK_ETP_DIR. GeneMark-ETP
  is academic / non-commercial (CC BY-NC-SA; no licence key needed).

To run the bin/ scripts by hand (main.nf does this for you), prepend the tool's
env to PATH, e.g.  export PATH=$GM_CONDA_BASE/envs/annot/bin:$PATH

nextflow.config resolves the env prefixes from $HOME/miniconda3/envs by default
(override with GM_CONDA_BASE). Two-stream run (works with the five envs above):

  nextflow run main.nf -c nextflow.config \
    --genome genome.fa --reads 'rnaseq/*_{1,2}.fastq.gz' \
    --proteome db/uniprot_sprot.fasta --augustus_species rice \
    --busco_lineage poales_odb10 --outdir gm_out
NOTE
