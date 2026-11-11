#!/usr/bin/env bash
# busco_train_augustus.sh — BUSCO-train an AUGUSTUS species model (genome mode,
# --long self-training), the way Gene-Miner does for organisms without a native
# AUGUSTUS model. Produces a WRITABLE augustus config dir whose species/ holds the
# retrained model, so downstream `augustus --species=<NAME> --AUGUSTUS_CONFIG_PATH=<cfg>`
# works. Prints the trained species NAME on the last stdout line.
#
# Usage: busco_train_augustus.sh <genome.fa> <lineage> <threads> <src_augustus_config> <out_cfg_dir>
set -euo pipefail
GENOME=$1; LIN=$2; THREADS=$3; SRCCFG=$4; OUTCFG=$5

# writable copy of the augustus config (BUSCO must be able to write the new species)
rm -rf "$OUTCFG"; cp -r "$SRCCFG" "$OUTCFG"
export AUGUSTUS_CONFIG_PATH="$(cd "$OUTCFG" && pwd)"

# BUSCO contacts busco-data.ezlab.org on every start, even when the lineage is
# already on disk, to fetch file_versions.tsv -- and it retries for ~17 minutes
# before giving up. That turns a transient (or host-specific) network problem
# into a failed pipeline hours into a run: seen here when one instance could not
# reach ezlab while its siblings could, after the alignment and assembly stages
# had already completed.
#
# So: if $LIN names a directory, treat it as a pre-staged dataset and run
# offline. BUSCO's --offline skips the version check entirely and takes the
# dataset from --download_path, which must be the parent of lineages/.
OFFLINE=()
if [ -d "$LIN" ]; then
  LIN_DIR=$(cd "$LIN" && pwd)
  LIN_NAME=$(basename "$LIN_DIR")
  DL_PATH=$(dirname "$(dirname "$LIN_DIR")")      # .../busco_downloads
  if [ "$(basename "$(dirname "$LIN_DIR")")" = "lineages" ]; then
    OFFLINE=(--offline --download_path "$DL_PATH")
    LIN="$LIN_NAME"
    echo "[busco] offline: lineage '$LIN_NAME' from $DL_PATH" >&2
  else
    echo "[busco] WARNING: '$LIN' is a directory but not inside a lineages/ dir;" >&2
    echo "[busco]          --offline needs <download_path>/lineages/<name>. Running online." >&2
  fi
fi

busco -i "$GENOME" -l "$LIN" -m genome --long --augustus -o busco_train -c "$THREADS" -f \
  "${OFFLINE[@]}" \
  > busco_train.log 2>&1 || { echo "BUSCO training failed; see busco_train.log" >&2; tail -20 busco_train.log >&2; exit 1; }

# BUSCO registers the retrained species directly in AUGUSTUS_CONFIG_PATH/species/BUSCO_*
SP=$(ls -d "$AUGUSTUS_CONFIG_PATH"/species/BUSCO_* 2>/dev/null | head -1 || true)
if [ -z "$SP" ]; then
  # fallback: copy retraining_parameters into the config
  RP=$(find busco_train -type d -name "retraining_parameters*" 2>/dev/null | head -1)
  TRN=$(find "$RP" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | head -1)
  [ -z "${TRN:-}" ] && { echo "ERROR: no retrained AUGUSTUS model found" >&2; exit 1; }
  NAME=$(basename "$TRN")
  mkdir -p "$AUGUSTUS_CONFIG_PATH/species/$NAME"
  cp "$TRN"/* "$AUGUSTUS_CONFIG_PATH/species/$NAME/"
  SP="$AUGUSTUS_CONFIG_PATH/species/$NAME"
fi
echo "trained model: $(basename "$SP")" >&2
basename "$SP"
