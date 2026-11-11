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

# Resolving the lineage (and whether BUSCO can run offline) is shared with the
# proteins-mode call in main.nf -- see bin/busco_lineage.sh for why.
HELPER=$(command -v busco_lineage.sh 2>/dev/null || echo "$(dirname "$0")/busco_lineage.sh")
[ -r "$HELPER" ] || { echo "ERROR: cannot find busco_lineage.sh (looked at $HELPER)" >&2; exit 1; }
# shellcheck source=busco_lineage.sh
. "$HELPER"
busco_resolve_lineage "$LIN"

busco -i "$GENOME" -l "$BUSCO_LIN" -m genome --long --augustus -o busco_train -c "$THREADS" -f \
  ${BUSCO_OFFLINE[@]+"${BUSCO_OFFLINE[@]}"} \
  > busco_train.log 2>&1 || { echo "BUSCO training failed; see busco_train.log" >&2; tail -20 busco_train.log >&2; exit 1; }

# BUSCO registers the retrained species directly in AUGUSTUS_CONFIG_PATH/species/BUSCO_*
# Glob, not `ls | head -1`: see the note in run_augustus_windowed.sh -- head
# closing the pipe makes ls die of SIGPIPE and pipefail promotes it.
SP=""
for d in "$AUGUSTUS_CONFIG_PATH"/species/BUSCO_*; do [ -d "$d" ] || continue; SP=$d; break; done
if [ -z "$SP" ]; then
  # fallback: copy retraining_parameters into the config
  # -print -quit stops find itself, so nothing writes to a closed pipe. With
  # `| head -1` under pipefail this assignment aborts the script (set -e) on any
  # tree where find has more to list -- in the one path that only runs when the
  # normal model registration has already failed, i.e. the worst place for it.
  RP=$(find busco_train -type d -name "retraining_parameters*" -print -quit 2>/dev/null || true)
  TRN=$(find "$RP" -maxdepth 1 -mindepth 1 -type d -print -quit 2>/dev/null || true)
  [ -z "${TRN:-}" ] && { echo "ERROR: no retrained AUGUSTUS model found" >&2; exit 1; }
  NAME=$(basename "$TRN")
  mkdir -p "$AUGUSTUS_CONFIG_PATH/species/$NAME"
  cp "$TRN"/* "$AUGUSTUS_CONFIG_PATH/species/$NAME/"
  SP="$AUGUSTUS_CONFIG_PATH/species/$NAME"
fi
echo "trained model: $(basename "$SP")" >&2
basename "$SP"
