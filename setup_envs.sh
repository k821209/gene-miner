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
# checkout for the optional 3rd stream (--run_genemark). That stream runs through
# GeneMark-ETP's own `gmetp.pl`, which bundles GeneMark + ProtHint and static
# builds of its third-party tools, so it needs only Perl (a few CPAN modules) +
# python3 from conda. Skip it with:  GM_SKIP_GENEMARK=1 bash setup_envs.sh
set -euo pipefail

# conda-forge FIRST and strict priority: mixing in the `defaults` channel pulls an
# older boost, and the bioconda augustus build then cannot load libboost_iostreams.
CH="--strict-channel-priority --override-channels -c conda-forge -c bioconda"
SOLVER=$(command -v mamba || command -v conda)
[ -n "$SOLVER" ] || { echo "ERROR: conda/mamba not found on PATH"; exit 1; }
echo "[setup] using: $SOLVER"

GM_CONDA_BASE="${GM_CONDA_BASE:-$HOME/miniconda3}"
ENVS="$GM_CONDA_BASE/envs"
mkdir -p "$ENVS"
echo "[setup] envs go to: $ENVS"

create () { echo "[setup] creating env '$1' in $ENVS"; "$SOLVER" create -y -p "$ENVS/$1" $CH "${@:2}"; }

create annot     hisat2 stringtie transdecoder samtools gffread
create augustus  augustus diamond
# Pin RepeatModeler: the unconstrained spec resolves to 1.0.8 (a 2015 release
# with no LTR pipeline), because that build has far fewer dependencies and the
# solver prefers it. It is not a conflict -- `repeatmodeler>=2.0.5 repeatmasker`
# resolves cleanly to 4.2.4 + 2.0.9 -- so asking is all that is needed. Left
# unpinned, de-novo masking silently runs a decade-old tool: on a 20 Mb test it
# reported "Discovery complete: 0 families found" after 1h27m.
create rmod      "repeatmodeler>=2.0.5" repeatmasker
# Pin the interpreter: eggnog-mapper is a noarch package and constrains no
# python version, so the solver takes the newest one -- 3.14 at the time of
# writing. But eggnogmapper/common.py does `from distutils.spawn import
# find_executable`, and distutils was REMOVED in python 3.12, so emapper.py
# cannot even print its version:
#     ModuleNotFoundError: No module named 'distutils'
# The install itself succeeds, which is what makes this expensive: the EGGNOG
# process runs after BUILD_UNION, so the failure lands with masking, alignment,
# assembly, both ab-initio streams and the union already computed. Four hosts
# here were in exactly that state. Asking for <3.12 costs nothing and the
# solver takes it without complaint.
create eggnog    "python<3.12" eggnog-mapper diamond
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
FAMDIR="$ENVS/rmod/share/RepeatMasker/Libraries/famdb"
DFAM_BASE=https://www.dfam.org/releases/current/families/FamDB
mkdir -p "$FAMDIR"
# The root partition alone leaves `Families : 0`: FamDB opens every partition the
# root's file_map declares, so RepeatMasker dies at startup until they are there.
# rmblast searches CONSENSUS sequences, so the hmm partitions are not needed for
# masking and are skipped. `famdb.py check <clade>` prints exactly which files a
# given lineage wants, if this set ever needs revisiting.
for f in dfam40.0.h5 dfam40.curated.consensus.0.h5 \
         dfam40.uncurated.consensus.0.h5 dfam40.uncurated.consensus.1.h5; do
  [ -s "$FAMDIR/$f" ] && continue
  echo "[setup] fetching $f"
  # -f matters: without it a 404 writes the error page under the partition's
  # name, and famdb then fails with "file signature not found" on a file that
  # looks present. Note the .gz suffix -- the bare .h5 paths do not exist.
  if curl -fSL --retry 3 -o "$FAMDIR/$f.gz" "$DFAM_BASE/$f.gz" && gunzip -f "$FAMDIR/$f.gz"; then
    sig=$(head -c 8 "$FAMDIR/$f" | od -An -tx1 | tr -d ' ')
    case "$sig" in
      89484446*) : ;;
      *) echo "[setup] WARN: $f is not HDF5 (sig=$sig) — removing"; rm -f "$FAMDIR/$f" ;;
    esac
  else
    rm -f "$FAMDIR/$f.gz"
    echo "[setup] WARN: could not fetch $f — see $DFAM_BASE/"
  fi
done
CONF="$(ls "$ENVS"/rmod/share/famdb-*/famdb.conf 2>/dev/null | head -1 || true)"
if [ -n "$CONF" ]; then
  if grep -q '^FAMDB_DATA_DIR' "$CONF"; then
    sed -i "s|^FAMDB_DATA_DIR.*|FAMDB_DATA_DIR = $FAMDIR|" "$CONF"
  else
    printf 'FAMDB_DATA_DIR = %s\n' "$FAMDIR" >> "$CONF"
  fi
fi

# Verify, rather than assume. Downloading the root partition and writing
# famdb.conf is not the same as having a usable library: the root partition can
# land intact and still report 0 families, in which case famdb.py raises and
# RepeatMasker dies at startup -- but setup has already printed success, so the
# failure surfaces hours later inside MASK_GENOME. Check it here, where the
# remedy is still in front of the person running setup.
FAMDB_PY="$(ls "$ENVS"/rmod/share/famdb-*/famdb.py 2>/dev/null | head -1 || true)"
DFAM_OK=0
if [ -n "$FAMDB_PY" ] && [ -x "$ENVS/rmod/bin/python3" ]; then
  FAM_INFO="$("$ENVS/rmod/bin/python3" "$FAMDB_PY" -i "$FAMDIR" info 2>&1 || true)"
  FAM_N="$(printf '%s' "$FAM_INFO" | sed -n 's/^ *Families *: *//p' | tr -d ', ' | head -1)"
  case "$FAM_N" in
    ''|0) DFAM_OK=0 ;;
    *)    DFAM_OK=1 ;;
  esac
fi
if [ "$DFAM_OK" = "1" ]; then
  echo "[setup] Dfam OK — $FAM_N families visible to RepeatMasker"
else
  DFAM_WARN="Dfam library in $FAMDIR is NOT usable (famdb reports no families).
           RepeatMasker will abort at startup, so de-novo masking and any run
           WITHOUT --masked_genome/--masked_rmout will fail in MASK_GENOME.
           The root partition alone is not always sufficient: the partitions the
           root declares must be present too. Fetch a matching set into
           $FAMDIR (see https://www.dfam.org/releases/ and
           FamDB/utils/download_famdb.py), or pin a repeatmasker build whose
           bundled library generation you already have."
  echo "[setup] WARN: $DFAM_WARN"
fi

# --- eggnog env: can emapper.py actually start? ---
# The pin above keeps a FRESH env off python 3.12+, but an env created before
# that pin existed is already broken and re-running this script will not fix it
# (conda sees the packages as satisfied). So check, and say so loudly: an
# emapper.py that cannot import distutils installs cleanly, reports nothing
# wrong, and fails at the EGGNOG process -- after the whole expensive part of
# the pipeline has run.
EGG_WARN=""
if [ -x "$ENVS/eggnog/bin/emapper.py" ]; then
  if ! PATH="$ENVS/eggnog/bin:$PATH" "$ENVS/eggnog/bin/emapper.py" --version >/dev/null 2>&1; then
    EGG_PY=$("$ENVS/eggnog/bin/python" -c 'import sys;print("%d.%d"%sys.version_info[:2])' 2>/dev/null || echo "?")
    if "$ENVS/eggnog/bin/python" -c "import distutils" >/dev/null 2>&1; then
      EGG_WARN="emapper.py in envs/eggnog does not run (python $EGG_PY). Run it by hand to see why:
             PATH=$ENVS/eggnog/bin:\$PATH emapper.py --version"
    else
      EGG_WARN="envs/eggnog has python $EGG_PY, which has no distutils, and
             eggnog-mapper imports it -- emapper.py cannot start. This env predates
             the python<3.12 pin. Fix it with either:
               $ENVS/eggnog/bin/python -m pip install setuptools     (adds a distutils replacement)
               $SOLVER remove -y -p $ENVS/eggnog && re-run this script  (clean, slower)
             Until then the EGGNOG process will fail AFTER BUILD_UNION."
    fi
    echo "[setup] WARN: $EGG_WARN"
  else
    echo "[setup] eggnog env ok: $(PATH="$ENVS/eggnog/bin:$PATH" "$ENVS/eggnog/bin/emapper.py" --version 2>/dev/null | grep -oE 'emapper-[0-9.]+' | head -1)"
  fi
fi

# --- 3rd stream: GeneMark-ETP ---
# GeneMark-ETP's gmetp.pl bundles GeneMark + ProtHint and its own static
# bedtools/samtools/hisat2/diamond/stringtie/gffread under tools/, so the only
# external needs are Perl (+ a few CPAN modules) and python3. We create a
# `genemark` env for those and clone the GeneMark-ETP repo next to the envs.
# GeneMark-ETP is CC BY-NC-SA (academic / non-commercial).
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
  # Run cpanm with envs/genemark FIRST on PATH: it needs that env's own `make`
  # (and perl) to build the modules — calling it by absolute path alone leaves
  # `make` off PATH and the build fails.
  PATH="$ENVS/genemark/bin:$PATH" "$ENVS/genemark/bin/cpanm" \
    --notest Math::Utils Statistics::LineFit \
    || echo "[setup] WARN: cpanm of Math::Utils/Statistics::LineFit failed — install them by hand into envs/genemark"
  GM_ETP_DIR="${GENEMARK_ETP_DIR:-$GM_CONDA_BASE/opt/GeneMark-ETP}"
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
  <conda_base>/opt/GeneMark-ETP. run_genemark_etp.sh runs GeneMark-ETP's own
  gmetp.pl, which bundles GeneMark + ProtHint + static bedtools/samtools/hisat2/
  diamond/stringtie. Override the checkout location with GENEMARK_ETP_DIR.
  GeneMark-ETP is academic / non-commercial (CC BY-NC-SA).

To run the bin/ scripts by hand (main.nf does this for you), prepend the tool's
env to PATH, e.g.  export PATH=$GM_CONDA_BASE/envs/annot/bin:$PATH

nextflow.config resolves the env prefixes from $HOME/miniconda3/envs by default
(override with GM_CONDA_BASE). Two-stream run (works with the five envs above):

  nextflow run main.nf -c nextflow.config \
    --genome genome.fa --reads 'rnaseq/*_{1,2}.fastq.gz' \
    --proteome db/uniprot_sprot.fasta --augustus_species rice \
    --busco_lineage poales_odb10 --outdir gm_out
NOTE

# Warnings printed above scroll off the top behind the forty-line note, which
# makes "we warned about it" weaker than it reads: the last thing on screen is
# a success-shaped summary. Re-print them here so an unusable install is the
# last thing seen, and exit non-zero so a caller that chains on this script
# stops instead of proceeding into a pipeline that cannot finish.
if [ -n "${DFAM_WARN:-}" ] || [ -n "${EGG_WARN:-}" ]; then
  echo
  echo "================= setup finished WITH PROBLEMS ================="
  [ -n "${DFAM_WARN:-}" ] && { echo "RepeatMasker / Dfam:"; echo "  $DFAM_WARN"; echo; }
  [ -n "${EGG_WARN:-}" ]  && { echo "eggNOG env:";          echo "  $EGG_WARN";  echo; }
  echo "The envs above exist; the parts named here will fail when the pipeline"
  echo "reaches them, which for eggNOG is after BUILD_UNION."
  echo "==============================================================="
  exit 1
fi
