#!/usr/bin/env bash
# busco_lineage.sh -- resolve a --busco_lineage value into the arguments BUSCO needs.
# SOURCE this file; it defines one function and runs nothing on its own.
#
# BUSCO contacts busco-data.ezlab.org on every start to fetch file_versions.tsv,
# even when the lineage it needs is already on disk, and retries for ~17 minutes
# before failing the run. One instance here could not reach ezlab while two
# siblings in the same zone could, and the failure landed after HISAT2_ALIGN,
# STRINGTIE and STRINGTIE_MERGE had completed -- ~17 h of work lost to a version
# file. --offline skips that check and reads the dataset from --download_path.
#
# This lives in its own file because BUSCO is invoked from more than one place
# (genome/--long training, and proteins mode on the union set). Resolving the
# lineage in each caller means the next caller added will miss it.
#
# Usage:
#   . "$(dirname "$0")/busco_lineage.sh"   # or: . "$(command -v busco_lineage.sh)"
#   busco_resolve_lineage "$LIN"
#   busco -l "$BUSCO_LIN" ... "${BUSCO_OFFLINE[@]}"
#
# Sets:
#   BUSCO_LIN       the value to pass to -l (a bare name when offline)
#   BUSCO_OFFLINE   array: (--offline --download_path <dir>), or empty

busco_resolve_lineage() {
  local lin=$1
  BUSCO_LIN=$lin
  BUSCO_OFFLINE=()

  # A plain string is a lineage name: leave it alone and let BUSCO fetch it.
  [ -d "$lin" ] || return 0

  local lin_dir lin_name parent
  lin_dir=$(cd "$lin" && pwd)
  lin_name=$(basename "$lin_dir")
  parent=$(dirname "$lin_dir")

  # --offline wants <download_path>/lineages/<name>. Check the shape instead of
  # assuming it, and STOP if it is wrong rather than falling back to the online
  # path: passing a directory is an unambiguous request to use a staged dataset,
  # and the fallback cannot honour it -- it would hand BUSCO a filesystem path as
  # a lineage name, so the run dies either way, just with an error that does not
  # name the real problem. Failing here costs seconds and says what to fix.
  if [ "$(basename "$parent")" != "lineages" ]; then
    echo "[busco] ERROR: '$lin' is a directory, so a pre-staged lineage was intended," >&2
    echo "[busco]        but --offline needs it at <download_path>/lineages/<name>." >&2
    echo "[busco]        Found parent directory: $parent (expected basename 'lineages')" >&2
    echo "[busco]        Either restage it under a lineages/ dir, or pass the bare" >&2
    echo "[busco]        lineage name to download it instead." >&2
    return 1
  fi

  BUSCO_LIN=$lin_name
  BUSCO_OFFLINE=(--offline --download_path "$(dirname "$parent")")
  echo "[busco] offline: lineage '$lin_name' from $(dirname "$parent")" >&2
}
