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
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$ENV"

if [ -z "$LIB" ]; then
  echo "[$(date +%T)] no library given — de-novo RepeatModeler (slow, hours)"
  LIB=$OUT/$(basename "$GENOME").repeatlib.fa
  if [ ! -s "$LIB" ]; then
    BuildDatabase -name "$OUT/rmdb" "$GENOME" > "$OUT/builddb.log" 2>&1
    # Run inside $OUT so both the RM_* work dir AND rmdb-families.fa land there.
    ( cd "$OUT" && RepeatModeler -database rmdb -threads $((PA*4)) > repeatmodeler.log 2>&1 )
    # RepeatModeler 2.x writes <db>-families.fa; 1.x writes RM_*/consensi.fa.classified.
    SRC=$(ls -t "$OUT"/rmdb-families.fa "$OUT"/RM_*/consensi.fa.classified 2>/dev/null | head -1)
    { [ -n "$SRC" ] && [ -s "$SRC" ]; } || { echo "ERROR: RepeatModeler produced no library"; exit 1; }
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
