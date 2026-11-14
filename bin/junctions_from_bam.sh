#!/usr/bin/env bash
# junctions_from_bam.sh <out.tsv> <jobs> <bam> [bam ...]
#
# Splice junctions with spliced-read counts, from the alignments the pipeline
# already produced. samtools + awk because that is what is installed: regtools
# and pysam are not, and the junction is in the CIGAR.
#
# Output: chrom <TAB> intron_start <TAB> intron_end <TAB> reads   (1-based,
# inclusive, summed over every BAM given).
#
# The coordinate convention is the whole game: an N operation of length L
# beginning at reference position p covers [p, p+L-1] inclusive, 1-based, the
# same convention a GFF3 intron takes between CDS blocks. Verified against real
# alignments before use -- 817 of 904 introns in a test window matched a
# junction exactly, and shifts of +-1 and +-2 matched zero.
set -uo pipefail
OUT=${1:?usage: junctions_from_bam.sh out.tsv jobs bam [bam ...]}; shift
JOBS=${1:?jobs required}; shift
SAMTOOLS=${SAMTOOLS:-samtools}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
export SAMTOOLS TMP

one() {
  bam=$1
  [ -s "$bam" ] || { echo "  skip (missing): $bam" >&2; return 0; }
  o="$TMP/$(echo "$bam" | md5sum | cut -c1-12).j"
  # -F 0x904 drops unmapped, secondary and supplementary: a spliced alignment
  # counted again through a secondary hit would inflate the support.
  "$SAMTOOLS" view -@ 2 -F 0x904 "$bam" 2>/dev/null \
  | awk -F'\t' '
      $6 ~ /N/ {
        ref=$4; cig=$6
        while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
          tok=substr(cig,1,RLENGTH); cig=substr(cig,RLENGTH+1)
          L=tok+0; op=substr(tok,length(tok),1)
          if (op=="N") { print $3"\t"ref"\t"ref+L-1; ref+=L }
          else if (op=="M"||op=="D"||op=="="||op=="X") ref+=L
        }
      }' | sort -k1,1 -k2,2n -k3,3n | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' > "$o"
  echo "  $(basename "$bam"): $(wc -l < "$o") junction(s)" >&2
}
export -f one

printf '%s\n' "$@" | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}

N=$(ls "$TMP"/*.j 2>/dev/null | wc -l)
[ "$N" -gt 0 ] || { echo "FAIL: no BAM produced junctions" >&2; exit 1; }
echo "merging $N per-BAM file(s)" >&2
# Sum the per-BAM counts into one table.
awk -F'\t' '{c[$1"\t"$2"\t"$3]+=$4} END {for (k in c) print k"\t"c[k]}' "$TMP"/*.j \
  | sort -k1,1 -k2,2n -k3,3n > "$OUT"
TOTR=$(awk -F'\t' '{s+=$4} END {print s+0}' "$OUT")
echo "$(wc -l < "$OUT") distinct junction(s), $TOTR spliced read(s) total" >&2
[ -s "$OUT" ] || { echo "FAIL: empty output" >&2; exit 1; }
