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
ENV=${REPEAT_ENV:-aleseq}            # conda env that owns RepeatMasker/RepeatModeler
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)              # absolutise so the RepeatModeler subshell cd is safe
GENOME=$(readlink -f "$GENOME")

# RepeatMasker (conda) needs its env ACTIVATED to set REPEATMASKER_DIR etc.;
# calling the binary directly fails with "REPEATMASKER_DIR does not exist".
# Locate conda.sh portably (do not assume $HOME/miniconda3).
CONDA_BASE="$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV"

if [ -z "$LIB" ]; then
  echo "[$(date +%T)] no library given — de-novo RepeatModeler (slow, hours)"
  LIB=$OUT/$(basename "$GENOME").repeatlib.fa
  if [ ! -s "$LIB" ]; then
    BuildDatabase -name "$OUT/rmdb" "$GENOME" > "$OUT/builddb.log" 2>&1
    # RepeatModeler >= 2.0.3 takes -threads; older builds take -pa (in units of 4 cores).
    if RepeatModeler -h 2>&1 | grep -q -- "-threads"; then
      RM_CPU="-threads $((PA*4))"
    else
      RM_CPU="-pa $PA"
    fi
    # Run inside $OUT so both the RM_* work dir AND rmdb-families.fa land there.
    ( cd "$OUT" && RepeatModeler -database rmdb $RM_CPU > repeatmodeler.log 2>&1 )
    # RepeatModeler 2.x writes <db>-families.fa; 1.x writes RM_*/consensi.fa.classified.
    SRC=$(ls -t "$OUT"/rmdb-families.fa "$OUT"/RM_*/consensi.fa.classified 2>/dev/null | head -1)
    if [ -z "$SRC" ] || [ ! -s "$SRC" ]; then
      # A repeat-poor genome (yeast, many microbial assemblies) legitimately yields
      # no families. That is a result, not a failure: emit an unmasked genome and an
      # empty RepeatMasker .out, so the TE filter simply removes nothing downstream.
      if grep -q "0 families found" "$OUT/repeatmodeler.log" 2>/dev/null; then
        echo "[$(date +%T)] WARNING: RepeatModeler found no repeat families — continuing unmasked"
        G=$(basename "$GENOME")
        cp "$GENOME" "$OUT/$G.masked"
        printf '   SW   perc perc perc  query      position in query           matching       repeat              position in  repeat\n score   div. del. ins.  sequence   begin  end    (left)   repeat         class/family    begin  end (left)   ID\n\n' > "$OUT/$G.out"
        exit 0
      fi
      echo "ERROR: RepeatModeler produced no library — see $OUT/repeatmodeler.log"
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
