#!/usr/bin/env bash
# run_eggnog.sh — eggNOG-mapper functional annotation of ONE proteome, used to
# drive the taxonomy QC filter (filter_taxonomy.py). Self-contained: it creates
# the conda env and fetches the eggNOG DB on first run, then annotates.
#
# Usage: bash run_eggnog.sh <proteins.faa> <out_prefix> [out_dir=emap_out]
# Output: <out_dir>/<out_prefix>.emapper.annotations  (feed to filter_taxonomy.py)
set -uo pipefail

FAA=${1:?usage: run_eggnog.sh proteins.faa out_prefix [out_dir]}
PREFIX=${2:?need out_prefix}
OUT=${3:-emap_out}
DATA_DIR=${EGGNOG_DB:-$HOME/eggnog_db}
ENV=${EGGNOG_ENV:-$HOME/miniconda3/envs/eggnog}
CONDA=${CONDA_BIN:-$HOME/miniconda3/bin/conda}
THREADS=${THREADS:-36}
mkdir -p "$OUT" "$DATA_DIR"

echo "[$(date +%T)] env — eggnog-mapper"
if [ ! -x "$ENV/bin/emapper.py" ] || [ ! -x "$ENV/bin/diamond" ]; then
  "$CONDA" create -y -n eggnog -c bioconda -c conda-forge "eggnog-mapper>=2.1.12" diamond || exit 1
fi
EM=$ENV/bin
# emapper bundles diamond into its own bin in some builds; symlink if missing
[ -x "$EM/diamond" ] || ln -sf "$(command -v diamond || echo "$ENV/bin/diamond")" "$EM/diamond" 2>/dev/null || true

echo "[$(date +%T)] eggNOG DB (download_eggnog_data.py fetch is unreliable; wget -c the archives directly)"
# eggnog.db.gz 6.8G, eggnog_proteins.dmnd.gz 5.2G, eggnog.taxa.tar.gz (taxonomy)
EGGURL=http://eggnog6.embl.de/download/emapperdb-5.0.2
if [ ! -s "$DATA_DIR/eggnog.db" ] || [ ! -s "$DATA_DIR/eggnog_proteins.dmnd" ]; then
  ( cd "$DATA_DIR"
    find . -maxdepth 1 -name '*.gz' -size 0 -delete    # clear 0-byte stubs
    for f in eggnog.db.gz eggnog_proteins.dmnd.gz eggnog.taxa.tar.gz; do
      [ -s "${f%.gz}" ] && continue
      n=0
      until gzip -t "$f" 2>/dev/null; do               # embl.de drops long connections; resume until intact
        n=$((n+1)); [ "$n" -gt 80 ] && { echo "gave up on $f after $n tries"; exit 1; }
        wget -c --timeout=30 --read-timeout=45 --tries=3 --progress=dot:giga "$EGGURL/$f" -o "wget.$f.log" || true
        sleep 3
      done
      echo "  $f complete ($n resume(s))"
    done
    [ -s eggnog.db ]            || gunzip -kf eggnog.db.gz
    [ -s eggnog_proteins.dmnd ] || gunzip -kf eggnog_proteins.dmnd.gz
    [ -s eggnog.taxa.db ]       || tar xzf eggnog.taxa.tar.gz
  ) || exit 1
fi

echo "[$(date +%T)] emapper (diamond, single --sensitive pass) on $FAA"
if [ ! -s "$OUT/$PREFIX.emapper.annotations" ]; then
  # --dmnd_iterate no: one sensitive pass (the default escalating multi-pass
  # took >1h/proteome here for negligible extra annotation).
  "$EM/emapper.py" -i "$FAA" --itype proteins -m diamond --dmnd_iterate no \
    --cpu "$THREADS" --data_dir "$DATA_DIR" --output_dir "$OUT" -o "$PREFIX" \
    --override --temp_dir "$OUT" > "$OUT/$PREFIX.emap.log" 2>&1 || { echo "emapper failed"; exit 1; }
fi
echo "[$(date +%T)] DONE -> $OUT/$PREFIX.emapper.annotations ($(grep -vc '^#' "$OUT/$PREFIX.emapper.annotations") annotated)"
