#!/usr/bin/env bash
# representative (longest-CDS) transcript BUSCO for each subsampling configuration
source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate busco6
cd /data2/k821209/gm/sub
for c in L4-D100 L2-D100 L1-D100 L4-D50 L4-D25; do (
  cd $c
  python3 - <<"PY"
import re
from collections import defaultdict
g2t=defaultdict(list); tl=defaultdict(int)
for l in open("gm_out/union.final.gff3"):
    f=l.rstrip().split("\t")
    if len(f)<9: continue
    if f[2]=="mRNA": g2t[re.search(r"Parent=([^;]+)",f[8]).group(1)].append(re.search(r"ID=([^;]+)",f[8]).group(1))
    elif f[2]=="CDS": tl[re.search(r"Parent=([^;]+)",f[8]).group(1)]+=int(f[4])-int(f[3])+1
prim={max(ts,key=lambda t:tl[t]) for ts in g2t.values()}
keep=False; o=open("primary.pep.fa","w")
for l in open("gm_out/union.final.pep.fa"):
    if l[0]==">": keep=l[1:].split()[0] in prim
    if keep: o.write(l)
PY
  busco -i primary.pep.fa -l poales_odb10 -m proteins -o busco_primary -f --offline --download_path /data2/k821209/gm/busco_downloads -c 10 > busco_primary.log 2>&1
  echo "$c $(grep -h C: busco_primary/short_summary*.txt)"
) & done; wait; echo BUSCO_PRIMARY_DONE
