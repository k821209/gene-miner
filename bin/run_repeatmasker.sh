#!/usr/bin/env bash
# run_repeatmasker.sh — mask interspersed repeats on ONE genome and emit the
# RepeatMasker .out, used to drive the TE QC filter (te_filter_genes.py).
# Builds a de-novo RepeatModeler library if one isn't supplied.
#
# Usage: bash run_repeatmasker.sh <genome.fa> [repeat_lib.fa] [out_dir=rm_out]
# Output: <out_dir>/<genome>.out  (feed to te_filter_genes.py)
set -uo pipefail

GENOME=${1:?usage: run_repeatmasker.sh genome.fa [repeat_lib.fa] [out_dir]}
LIB=${2:-}
OUT=${3:-rm_out}
PA=${PA:-9}                          # x4 threads each = 36
ENV=${REPEAT_ENV:-aleseq}            # conda env that owns RepeatMasker/RepeatModeler
mkdir -p "$OUT"

# RepeatMasker (conda) needs its env ACTIVATED to set REPEATMASKER_DIR etc.;
# calling the binary directly fails with "REPEATMASKER_DIR does not exist".
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "$ENV"

if [ -z "$LIB" ]; then
  echo "[$(date +%T)] no library given — de-novo RepeatModeler (slow, hours)"
  LIB=$OUT/$(basename "$GENOME").repeatlib.fa
  if [ ! -s "$LIB" ]; then
    BuildDatabase -name "$OUT/rmdb" "$GENOME" > "$OUT/builddb.log" 2>&1
    RepeatModeler -database "$OUT/rmdb" -threads $((PA*4)) > "$OUT/repeatmodeler.log" 2>&1
    cp "$OUT"/RM_*/consensi.fa.classified "$LIB"
  fi
fi
# NB: to compare TE content BETWEEN genomes, run RepeatMasker on both with the
# SAME -lib so the difference reflects the genomes, not the library.

echo "[$(date +%T)] RepeatMasker -lib $(basename "$LIB") on $(basename "$GENOME")"
[ -s "$OUT/$(basename "$GENOME").tbl" ] || \
  RepeatMasker -lib "$LIB" -pa "$PA" -xsmall -gff -dir "$OUT" "$GENOME" > "$OUT/rm.log" 2>&1 || exit 1

echo "[$(date +%T)] DONE -> $OUT/$(basename "$GENOME").out (+ .tbl summary, .out.gff coords)"
sed -n '1,40p' "$OUT/$(basename "$GENOME").tbl"
