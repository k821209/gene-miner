#!/usr/bin/env bash
# run_genemark_etp.sh — produce GeneMark-ETP predictions (the 3rd ab-initio
# stream for build_union.py --genemark) from the pipeline's own inputs: the
# soft-masked genome, the HISAT2 BAMs, and the protein evidence. GeneMark-ETP is
# driven through BRAKER3 (which wraps GeneMark-ETP + ProtHint); only its
# GeneMark-ETP output (genemark.gtf) is exported — AUGUSTUS/TSEBRA steps are a
# by-product. Speed is not the point here: maximal de-novo gene recall is.
#
# Usage: run_genemark_etp.sh <masked_genome.fa> <proteins.fa> <cpus> "<bam1 bam2 ...>" <out.gtf>
set -uo pipefail

MASKED=${1:?masked genome}
PROT=${2:?protein fasta}
CPUS=${3:?cpus}
BAMS_RAW=${4:?space-separated bams}
OUT=${5:-genemark.gtf}

BRAKER_ENV=${BRAKER_ENV:-$HOME/miniconda3/envs/braker3}
GM_PATH=${GENEMARK_PATH:-$HOME/gene-miner-runs/GeneMark-ETP/bin}
PH_PATH=${PROTHINT_PATH:-$HOME/gene-miner-runs/GeneMark-ETP/bin/gmes/ProtHint/bin}

# GeneMark/BRAKER require single-token FASTA headers; keep soft-masking intact.
echo "[$(date +%T)] single-token headers -> genome.clean.fa"
awk '/^>/{print $1; next} {print}' "$MASKED" > genome.clean.fa

# coordinate-sort BAMs (GeneMark-ETP re-sorts, but provide tidy inputs); comma-join
ST="$BRAKER_ENV/bin/samtools"
BAMLIST=""
i=0
for b in $BAMS_RAW; do
  i=$((i+1))
  sb="sorted_${i}.bam"
  if "$ST" view -H "$b" 2>/dev/null | grep -q 'SO:coordinate'; then
    cp "$b" "$sb"
  else
    "$ST" sort -@ "$CPUS" -o "$sb" "$b"
  fi
  BAMLIST="${BAMLIST:+$BAMLIST,}$sb"
done

# unique species tag (avoid AUGUSTUS species-dir clashes); $$ = PID, no clock needed
SP="gmetp_$$"
WD="braker_etp_$$"

echo "[$(date +%T)] BRAKER3 (GeneMark-ETP) on $(echo "$BAMLIST" | tr ',' '\n' | wc -l) BAM(s)"
"$BRAKER_ENV/bin/braker.pl" \
  --genome=genome.clean.fa \
  --bam="$BAMLIST" \
  --prot_seq="$PROT" \
  --GENEMARK_PATH="$GM_PATH" \
  --PROTHINT_PATH="$PH_PATH" \
  --threads="$CPUS" \
  --species="$SP" \
  --softmasking \
  --gff3 \
  --workingdir="$WD" > braker_etp.log 2>&1
RC=$?

GMG="$WD/GeneMark-ETP/genemark.gtf"
if [ ! -s "$GMG" ]; then
  echo "ERROR: GeneMark-ETP did not produce genemark.gtf (braker rc=$RC); see braker_etp.log" >&2
  tail -30 braker_etp.log >&2
  exit 1
fi
cp "$GMG" "$OUT"
echo "[$(date +%T)] DONE -> $OUT ($(grep -c $'\tCDS\t' "$OUT") CDS rows)"
