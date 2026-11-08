#!/usr/bin/env bash
# run_repeatmasker.sh — mask interspersed repeats on ONE genome and emit the
# RepeatMasker .out, used to drive the TE QC filter (te_filter_genes.py).
# Builds a de-novo RepeatModeler library if one isn't supplied.
#
# Usage: bash run_repeatmasker.sh <genome.fa> [repeat_lib.fa] [out_dir=rm_out]
# Output: <out_dir>/<genome>.masked (soft) + <out_dir>/<genome>.out
set -uo pipefail

GENOME=${1:?usage: run_repeatmasker.sh genome.fa [repeat_lib.fa] [out_dir]}
LIB=${2:-}
OUT=${3:-rm_out}
PA=${PA:-9}                          # x4 threads each = 36
ENV=${REPEAT_ENV:-aleseq}            # conda env owning RepeatMasker/RepeatModeler;
                                     # a NAME under the conda base, or a full prefix path
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)              # absolutise so the RepeatModeler subshell cd is safe
GENOME=$(readlink -f "$GENOME")

# RepeatMasker (conda) needs its env ACTIVATED to set REPEATMASKER_DIR etc.;
# calling the binary directly fails with "REPEATMASKER_DIR does not exist".
#
# Finding conda.sh: `conda info --base` only works when conda is already on
# PATH, which it is not inside a Nextflow process, and the old fallback was the
# literal guess $HOME/miniconda3. That guess is wrong for every install the
# README actively recommends -- anything placed with GM_CONDA_BASE -- and the
# failure arrives as `No such file or directory` for a conda.sh nobody asked
# for. Derive the base from the env prefix instead: $ENV is normally
# <base>/envs/<name>, so the base is two levels up. Fall back in order.
CONDA_BASE=""
case "$ENV" in
  */envs/*) cand="${ENV%/envs/*}"; [ -r "$cand/etc/profile.d/conda.sh" ] && CONDA_BASE="$cand" ;;
esac
if [ -z "$CONDA_BASE" ] && [ -n "${GM_CONDA_BASE:-}" ] \
   && [ -r "$GM_CONDA_BASE/etc/profile.d/conda.sh" ]; then CONDA_BASE="$GM_CONDA_BASE"; fi
if [ -z "$CONDA_BASE" ]; then CONDA_BASE="$(conda info --base 2>/dev/null || true)"; fi
if [ -z "$CONDA_BASE" ] || [ ! -r "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  [ -r "$HOME/miniconda3/etc/profile.d/conda.sh" ] && CONDA_BASE="$HOME/miniconda3"
fi
if [ -z "$CONDA_BASE" ] || [ ! -r "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  echo "ERROR: cannot find conda.sh. Tried: the base implied by REPEAT_ENV" >&2
  echo "       ('$ENV'), \$GM_CONDA_BASE ('${GM_CONDA_BASE:-unset}'), 'conda info --base'," >&2
  echo "       and \$HOME/miniconda3. Set GM_CONDA_BASE or pass REPEAT_ENV as a" >&2
  echo "       full prefix path under <base>/envs/." >&2
  exit 1
fi
echo "[$(date +%T)] conda base: $CONDA_BASE"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
# $ENV may be a bare env name or a full prefix path. `conda activate` accepts
# both, but a bare name is only resolvable under the DEFAULT envs dir -- so an
# env created elsewhere (GM_CONDA_BASE, or any --env_rmod outside the base)
# must be activated by prefix or this silently fails at the masking step.
if ! conda activate "$ENV" 2>/dev/null; then
  echo "ERROR: could not 'conda activate $ENV'" >&2
  echo "       REPEAT_ENV must be an env name under $CONDA_BASE/envs, or a full" >&2
  echo "       prefix path to the env directory." >&2
  exit 1
fi

if [ -z "$LIB" ]; then
  echo "[$(date +%T)] no library given — de-novo RepeatModeler (slow, hours)"
  LIB=$OUT/$(basename "$GENOME").repeatlib.fa
  if [ ! -s "$LIB" ]; then
    BuildDatabase -name "$OUT/rmdb" "$GENOME" > "$OUT/builddb.log" 2>&1
    # RepeatModeler >= 2.0.3 takes -threads; older builds take -pa (in units of
    # 4 cores). Capture the help text FIRST: `RepeatModeler -h` prints its usage
    # and exits 1, and under `set -o pipefail` (line 8) that non-zero status wins
    # over grep's success -- so testing it inside a pipeline always picked the
    # -pa branch, no matter what the build supports. RepeatModeler 2.0.8 rejects
    # -pa outright ("The -pa parameter has been deprecated"), which made de-novo
    # masking impossible; on an older build the wrong branch happened to work,
    # which is why the detection could be dead code for this long.
    RM_HELP=$(RepeatModeler -h 2>&1 || true)
    if printf '%s' "$RM_HELP" | grep -q -- "-threads"; then
      RM_CPU="-threads $((PA*4))"
    else
      RM_CPU="-pa $PA"
    fi
    echo "[$(date +%T)] RepeatModeler parallelism: $RM_CPU"
    # Run inside $OUT so both the RM_* work dir AND rmdb-families.fa land there.
    # Keep the exit status: it is what separates "this genome has no repeat
    # families" from "RepeatModeler broke", and those need opposite handling.
    ( cd "$OUT" && RepeatModeler -database rmdb $RM_CPU > repeatmodeler.log 2>&1 )
    RM_RC=$?
    # RepeatModeler 2.x writes <db>-families.fa; 1.x writes RM_*/consensi.fa.classified.
    SRC=$(ls -t "$OUT"/rmdb-families.fa "$OUT"/RM_*/consensi.fa.classified 2>/dev/null | head -1)
    if [ -z "$SRC" ] || [ ! -s "$SRC" ]; then
      # A repeat-poor genome (yeast, many microbial assemblies) legitimately yields
      # no families. That is a result, not a failure: emit an unmasked genome and an
      # empty RepeatMasker .out, so the TE filter simply removes nothing downstream.
      # Decide on the EXIT STATUS, not on a log phrase. "0 families found" is
      # one wording among several (it varies with RepeatModeler version and with
      # how early the search gives up), so matching it made a successful run on
      # a repeat-poor genome exit 1 -- the opposite of what the branch is for.
      # A non-zero status is a real failure and must stay fatal.
      if [ "$RM_RC" -eq 0 ]; then
        echo "[$(date +%T)] WARNING: RepeatModeler exited 0 with no library — treating this genome as repeat-free, continuing unmasked"
        grep -i "families found\|no repeat" "$OUT/repeatmodeler.log" 2>/dev/null | tail -2
        G=$(basename "$GENOME")
        cp "$GENOME" "$OUT/$G.masked"
        printf '   SW   perc perc perc  query      position in query           matching       repeat              position in  repeat\n score   div. del. ins.  sequence   begin  end    (left)   repeat         class/family    begin  end (left)   ID\n\n' > "$OUT/$G.out"
        exit 0
      fi
      echo "ERROR: RepeatModeler failed (exit $RM_RC) — see $OUT/repeatmodeler.log"
      tail -5 "$OUT/repeatmodeler.log" 2>/dev/null
      exit 1
    fi
    cp "$SRC" "$LIB"
  fi
fi
# NB: to compare TE content BETWEEN genomes, run RepeatMasker on both with the
# SAME -lib so the difference reflects the genomes, not the library.

echo "[$(date +%T)] RepeatMasker -lib $(basename "$LIB") on $(basename "$GENOME")"
[ -s "$OUT/$(basename "$GENOME").tbl" ] || \
  RepeatMasker -lib "$LIB" -pa "$PA" -xsmall -gff -dir "$OUT" "$GENOME" > "$OUT/rm.log" 2>&1 || exit 1

echo "[$(date +%T)] DONE -> $OUT/$(basename "$GENOME").out (+ .masked soft-masked, .tbl summary)"
sed -n '1,40p' "$OUT/$(basename "$GENOME").tbl"
