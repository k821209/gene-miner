#!/usr/bin/env bash
# run_genemark_etp.sh — produce GeneMark-ETP predictions (the 3rd ab-initio
# stream for build_union.py --genemark) from the pipeline's own inputs: the
# soft-masked genome, the HISAT2 BAMs, and the protein evidence.
#
# GeneMark-ETP is driven DIRECTLY through its native entry point `gmetp.pl`
# (github.com/gatech-genemark/GeneMark-ETP) — NOT through BRAKER. GeneMark-ETP
# bundles GeneMark-ES/ET/EP+ and ProtHint, and ships static binaries of every
# third-party tool it needs (bedtools, samtools, hisat2, diamond, stringtie,
# gffread) under its own tools/ folder. It does NOT use GenomeThreader. So the
# only external dependencies are Perl (with a handful of CPAN modules) and
# python3, both installed by setup_envs.sh into the `genemark` conda env. There
# is no BRAKER container and no genomethreader in this path.
#
# gmetp.pl consumes a small YAML config (this script writes it) plus a directory
# of coordinate-sorted BAMs named sorted_1.bam..sorted_N.bam (passed with --bam,
# so no re-alignment) and produces <workdir>/genemark.gtf.
#
# Usage: run_genemark_etp.sh <masked_genome.fa> <proteins.fa> <cpus> "<bam1 bam2 ...>" <out.gtf>
set -uo pipefail

MASKED=${1:?masked genome}
PROT=${2:?protein fasta}
CPUS=${3:?cpus}
BAMS_RAW=${4:?space-separated bams}
OUT=${5:-genemark.gtf}

# GeneMark-ETP git checkout (cloned by setup_envs.sh). tools/ holds the bundled
# static third-party binaries; bin/ holds gmetp.pl + GeneMark + ProtHint.
GM_ETP=${GENEMARK_ETP_DIR:-$HOME/GeneMark-ETP}
[ -x "$GM_ETP/bin/gmetp.pl" ] || {
  echo "ERROR: GeneMark-ETP not found at '$GM_ETP' (bin/gmetp.pl missing)." >&2
  echo "       Set GENEMARK_ETP_DIR, or run setup_envs.sh which clones it." >&2
  exit 1; }

# Perl (+CPAN modules) and python3 come from the 'genemark' conda env if the
# caller put it on PATH (main.nf does). The bundled tools/ must come FIRST so
# gmetp.pl uses GeneMark-ETP's own static bedtools/samtools/hisat2/diamond.
export PATH="$GM_ETP/tools:$GM_ETP/bin:$PATH"

command -v perl >/dev/null   || { echo "ERROR: perl not found (genemark conda env not on PATH)" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found (genemark conda env not on PATH)" >&2; exit 1; }
[ -s "$PROT" ] || { echo "ERROR: protein DB '$PROT' is missing or empty" >&2; exit 1; }
SAMTOOLS="$GM_ETP/tools/samtools"

# GeneMark requires single-token FASTA headers; keep soft-masking (lower-case)
# intact — do NOT upper-case, gmetp.pl is invoked with --softmask below.
echo "[$(date +%T)] single-token headers -> genome.clean.fa"
awk '/^>/{print $1; next} {print}' "$MASKED" > genome.clean.fa
GENOME=$(readlink -f genome.clean.fa)
PROTABS=$(readlink -f "$PROT")

# Coordinate-sort the BAMs into bams/sorted_1.bam..N (gmetp.pl --bam symlinks
# <dir>/<set>.bam and skips its own alignment step).
mkdir -p bams
i=0; SETS=""
for b in $BAMS_RAW; do
  i=$((i+1)); sb="bams/sorted_${i}.bam"
  if "$SAMTOOLS" view -H "$b" 2>/dev/null | grep -q 'SO:coordinate'; then
    cp "$b" "$sb"
  else
    "$SAMTOOLS" sort -@ "$CPUS" -o "$sb" "$b"
  fi
  SETS="${SETS:+$SETS,}sorted_${i}"
done
BAMDIR=$(readlink -f bams)
[ -n "$SETS" ] || { echo "ERROR: no BAM files supplied" >&2; exit 1; }

# unique species tag (avoids model-dir clashes); $$ = PID, no clock needed
SP="gmetp_$$"
WD=$(readlink -f "etp_$$"); mkdir -p "$WD"

# The native GeneMark-ETP YAML config (this is exactly what BRAKER used to
# generate and hand to gmetp.pl). RepeatMasker_path/annot_path stay empty:
# the genome is already masked (--softmask) and we have no reference annotation.
cat > "$WD/etp_config.yaml" <<EOF
---
RepeatMasker_path: ''
annot_path: ''
genome_path: $GENOME
protdb_path: $PROTABS
rnaseq_sets: [$SETS]
species: $SP
EOF

echo "[$(date +%T)] gmetp.pl on $i BAM(s), $CPUS cores (softmask, pre-aligned BAMs)"
perl "$GM_ETP/bin/gmetp.pl" \
  --cfg "$WD/etp_config.yaml" \
  --workdir "$WD" \
  --cores "$CPUS" \
  --softmask \
  --bam "$BAMDIR" \
  --verbose > gmetp.log 2>&1
GM_RC=$?

GMG="$WD/genemark.gtf"
if [ ! -s "$GMG" ]; then
  echo "ERROR: gmetp.pl (rc=$GM_RC) did not produce genemark.gtf; see gmetp.log" >&2
  tail -40 gmetp.log >&2
  exit 1
fi
cp "$GMG" "$OUT"
echo "[$(date +%T)] DONE -> $OUT ($(grep -c $'\tCDS\t' "$OUT") CDS rows)"
