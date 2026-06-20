#!/usr/bin/env bash
# run_augustus_windowed.sh — AUGUSTUS ab initio on a FINISHED (chromosome-scale)
# genome. Each sequence is scanned in overlapping windows via AUGUSTUS
# --predictionStart/--predictionEnd, so coordinates stay global to the sequence
# (no scaffold lift needed). Output: augustus_scaffold.gff3 in the CWD.
#
# Usage: run_augustus_windowed.sh <masked.fa> <species> <config_path> <win_mb> <overlap_bp> <threads>
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
MASKED=$1; SP=$2; CFG=$3; WINMB=$4; OV=$5; THREADS=$6
export AUGUSTUS_CONFIG_PATH=$CFG
# SP="auto" -> use the BUSCO-trained model registered in this config
if [ "$SP" = "auto" ]; then
  SP=$(ls -d "$CFG"/species/BUSCO_* 2>/dev/null | head -1 | xargs -r basename)
  [ -z "$SP" ] && { echo "ERROR: SP=auto but no BUSCO_* species in $CFG/species" >&2; exit 1; }
  echo "using BUSCO-trained model: $SP"
fi
WIN=$(( WINMB * 1000000 ))

rm -rf seqs out; mkdir -p seqs out
awk '/^>/{n=substr($1,2); sub(/ .*/,"",n); f="seqs/"n".fa"} {print > f}' "$MASKED"

: > jobs.txt
for fa in seqs/*.fa; do
  L=$(grep -v "^>" "$fa" | tr -d "\n" | wc -c)
  s=1
  while [ "$s" -le "$L" ]; do
    e=$(( s + WIN - 1 )); [ "$e" -gt "$L" ] && e=$L
    echo "$fa $s $e" >> jobs.txt
    [ "$e" -ge "$L" ] && break
    s=$(( e - OV + 1 ))
  done
done
echo "windows: $(wc -l < jobs.txt)"

run_one(){
  fa=$1; s=$2; e=$3; b=$(basename "$fa" .fa)
  augustus --species="$SP" --strand=both --genemodel=complete --gff3=on --UTR=off \
    --predictionStart="$s" --predictionEnd="$e" "$fa" > "out/${b}_${s}.gff3" 2>/dev/null
}
export -f run_one; export SP
xargs -P "$THREADS" -L1 bash -c 'run_one "$@"' _ < jobs.txt

# merge; prefix gene/transcript IDs per (sequence,window) to keep them unique
: > augustus_raw.gff3
for g in out/*.gff3; do
  tag=$(basename "$g" .gff3)
  awk -v c="$tag" '!/^#/ && NF>=8 {gsub(/ID=g/,"ID="c"_g"); gsub(/Parent=g/,"Parent="c"_g"); gsub(/ID=t/,"ID="c"_t"); gsub(/Parent=t/,"Parent="c"_t"); print}' "$g" >> augustus_raw.gff3
done
echo "augustus genes (pre-dedup): $(awk -F'\t' '$3=="gene"' augustus_raw.gff3 | wc -l)"
# collapse genes duplicated across overlapping windows -> augustus_scaffold.gff3
python3 "$SCRIPT_DIR/dedup_gff_genes.py" augustus_raw.gff3 augustus_scaffold.gff3 0.8
echo "augustus genes (post-dedup): $(awk -F'\t' '$3=="gene"' augustus_scaffold.gff3 | wc -l)"
