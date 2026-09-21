#!/usr/bin/env bash
# Per-library rice splice junctions, so a catalogue built from a subset of the
# libraries can be checked against the libraries it never saw (the held-out test
# for Reviewer 1, point 3).
set -uo pipefail
R=~/gene-miner-runs; W=$R/revision/junctions; mkdir -p $W
source ~/miniconda3/etc/profile.d/conda.sh; conda activate regt
for lib in endosperm panicle seedling leaf; do
  bam=$(ls $R/rice/work/*/*/$lib.bam 2>/dev/null | head -1)
  [ -s "$bam" ] || { echo "no bam for $lib"; continue; }
  [ -s $W/rice.$lib.junctions.tsv ] && { echo "$lib cached"; continue; }
  regtools junctions extract -s XS -o $W/rice.$lib.bed $bam 2>/dev/null
  python3 - $W/rice.$lib.bed $W/rice.$lib.junctions.tsv <<"PY"
import sys
from collections import defaultdict
inp, out = sys.argv[1:3]
c = defaultdict(int)
for line in open(inp):
    f = line.rstrip("\n").split("\t")
    if len(f) < 11:
        continue
    b = [int(x) for x in f[10].rstrip(",").split(",")]
    c[(f[0], int(f[1]) + b[0] + 1, int(f[2]) - b[-1])] += int(f[4])
with open(out, "w") as o:
    for (ch, s, e), n in sorted(c.items()):
        o.write(f"{ch}\t{s}\t{e}\t{n}\n")
print(out, len(c), "junctions")
PY
done
# held-out set for the two-library configuration (built from leaf + panicle)
python3 - $W <<"PY"
import sys, os
from collections import defaultdict
W = sys.argv[1]
c = defaultdict(int)
for lib in ("seedling", "endosperm"):
    p = os.path.join(W, f"rice.{lib}.junctions.tsv")
    if not os.path.exists(p):
        continue
    for line in open(p):
        ch, s, e, n = line.split("\t")
        c[(ch, int(s), int(e))] += int(n)
out = os.path.join(W, "rice.heldout_seedling_endosperm.junctions.tsv")
with open(out, "w") as o:
    for (ch, s, e), n in sorted(c.items()):
        o.write(f"{ch}\t{s}\t{e}\t{n}\n")
print(out, len(c), "junctions")
PY
echo HELDOUT_JUNCTIONS_DONE
