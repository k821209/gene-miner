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

busco -i "$GENOME" -l "$LIN" -m genome --long --augustus -o busco_train -c "$THREADS" -f \
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
