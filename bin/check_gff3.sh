#!/usr/bin/env bash
# check_gff3.sh <file.gff3> [reference.gff3]
#
# The sanity check this repo did not have, and which found four defects in one
# week that every self-written check had passed.
#
# Validate a GFF3 with genometools, and -- given the reference it was built from
# -- say whether an error is INHERITED from that reference or INTRODUCED here.
# That distinction is the whole point: three of four reference annotations in
# hand are not valid GFF3 as supplied, and failing identically on inherited and
# introduced defects would mean no such reference could ever be supplemented.
#
# Exit: 0 valid, or invalid with the same error as the reference (inherited)
#       1 invalid with an error the reference does not have (introduced)
#       2 the validator could not be run -- NOT a pass
set -uo pipefail
F=${1:?usage: check_gff3.sh <file.gff3> [reference.gff3]}
REF=${2:-}
GT=${GT:-$(command -v gt || true)}

[ -s "$F" ] || { echo "check_gff3: $F missing or empty" >&2; exit 1; }
if [ -z "$GT" ]; then
  # An absent validator must not read like a satisfied one.
  echo "check_gff3: UNKNOWN -- genometools \`gt\` not found. Install it (conda: genometools-genometools)" >&2
  echo "            or set GT=/path/to/gt. This is not a pass." >&2
  exit 2
fi

# Strip the filename and line number so two files can be compared on the KIND
# of error rather than on text that differs by construction.
kind(){ sed -e 's/in file "[^"]*"//' -e 's/\bline [0-9]\+\b/line N/' -e 's/"[^"]*"/"X"/g' -e 's/[[:space:]]\+/ /g'; }

OE=$("$GT" gff3validator "$F" 2>&1 | grep -i 'error' | head -1 || true)
if [ -z "$OE" ]; then echo "check_gff3: $F is valid GFF3"; exit 0; fi
echo "check_gff3: $F is NOT valid GFF3"
echo "  $OE"

if [ -z "$REF" ]; then
  echo "  (no reference given, so inherited vs introduced cannot be decided)" >&2
  exit 1
fi

# Compare like with like. The validator stops at the FIRST error, so a reference
# without a ##gff-version line reports only that and hides every deeper defect,
# while the output -- which usually gains that line -- reports the next one down.
# Without this, an inherited defect looks introduced. Measured on a real pair:
# the raw reference reported only the missing version line; the same file with
# the line prepended reported the identical CDS phase error the output had.
REFCMP=$REF; TMP=""
if ! head -1 "$REF" | grep -q '^##gff-version'; then
  TMP=$(mktemp); { echo '##gff-version 3'; cat "$REF"; } > "$TMP"; REFCMP=$TMP
  echo "  (the reference lacks a ##gff-version line; validated with one prepended)"
fi
RE=$("$GT" gff3validator "$REFCMP" 2>&1 | grep -i 'error' | head -1 || true)
[ -n "$TMP" ] && rm -f "$TMP"

if [ -z "$RE" ]; then
  echo "  the reference IS valid, so this error was introduced here" >&2
  exit 1
fi
echo "  reference error: $RE"
if [ "$(printf '%s' "$OE" | kind)" = "$(printf '%s' "$RE" | kind)" ]; then
  echo "  INHERITED: the reference as supplied has the same error. The delivered file"
  echo "  carries the reference's defect -- say so to whoever receives it."
  exit 0
fi
echo "  INTRODUCED: this error is not in the reference" >&2
exit 1
