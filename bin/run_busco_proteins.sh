#!/usr/bin/env bash
# run_busco_proteins.sh -- BUSCO completeness on a protein set (proteins mode).
#
# A wrapper rather than an inline command in main.nf so that the lineage
# resolution in bin/busco_lineage.sh applies here too. It did not, once: the
# training call learned to run offline and this one kept contacting ezlab, so a
# host that cannot reach ezlab would have failed the run again at the last
# process -- after everything else had been computed.
#
# Usage: run_busco_proteins.sh <proteins.fa> <lineage> <threads> <outdir_name>
set -euo pipefail
PEP=$1; LIN=$2; THREADS=$3; OUT=$4

HELPER=$(command -v busco_lineage.sh 2>/dev/null || echo "$(dirname "$0")/busco_lineage.sh")
[ -r "$HELPER" ] || { echo "ERROR: cannot find busco_lineage.sh (looked at $HELPER)" >&2; exit 1; }
# shellcheck source=busco_lineage.sh
. "$HELPER"
busco_resolve_lineage "$LIN"

busco -i "$PEP" -l "$BUSCO_LIN" -m proteins -c "$THREADS" -o "$OUT" -f \
  ${BUSCO_OFFLINE[@]+"${BUSCO_OFFLINE[@]}"} \
  > busco.log 2>&1 || { echo "BUSCO (proteins) failed; see busco.log" >&2; tail -20 busco.log >&2; exit 1; }
