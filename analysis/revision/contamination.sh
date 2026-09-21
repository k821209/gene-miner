#!/usr/bin/env bash
# What a broader non-target screen would catch (R2-M5):
#   (a) eggNOG assignments outside the target kingdom — bacteria (what the filter
#       removes today), archaea, fungi, other eukaryotes, viruses
#   (b) organellar origin: catalogue CDS aligned to the species' own mitochondrial
#       and plastid genomes (NUMT / NUPT)
set -uo pipefail
R=~/gene-miner-runs; REV=$R/revision; W=$REV/contamination
mkdir -p $W/org; cd $W
source ~/miniconda3/etc/profile.d/conda.sh

# --- (a) taxonomic classes present in the eggNOG assignments ------------------
python3 - <<'PY' > eggnog_classes.tsv
import os, re
from collections import Counter
HOME = os.path.expanduser('~')
GID = re.compile(r"[A-Za-z][A-Za-z0-9]*_[AER]\d+")
CLADE = {'2': 'Bacteria', '2157': 'Archaea', '4751': 'Fungi', '33090': 'Viridiplantae',
         '33208': 'Metazoa', '2759': 'Eukaryota', '10239': 'Viruses', '554915': 'Amoebozoa',
         '33630': 'Alveolata-Stramenopiles', '5794': 'Apicomplexa', '5653': 'Kinetoplastida'}
TARGET = {'rice': '33090', 'soybean': '33090', 'dmel': '33208', 'cele': '33208'}
print('genome\tclass\tloci')
for g in ('rice', 'soybean', 'dmel', 'cele'):
    d = f'{HOME}/gene-miner-runs/revision/regen/{g}_s'
    if not os.path.isdir(d): d = f'{HOME}/gene-miner-runs/revision/regen/{g}'
    ann = f'{d}/new.annotations'
    if not os.path.exists(ann):
        continue
    H, per_gene = [], {}
    for line in open(ann):
        if line.startswith('#'):
            if line.startswith('#query'):
                H = line.lstrip('#').rstrip('\n').split('\t')
            continue
        r = dict(zip(H, line.rstrip('\n').split('\t')))
        m = GID.search(r.get('query', ''))
        if not m:
            continue
        ogs = r.get('eggNOG_OGs', '-')
        tids = {t.split('@', 1)[1].split('|', 1)[0] for t in ogs.split(',') if '@' in t}
        per_gene.setdefault(m.group(), set()).update(tids)
    c = Counter()
    tgt = TARGET[g]
    for gene, tids in per_gene.items():
        if tgt in tids:
            c['target kingdom'] += 1
        elif '2' in tids:
            c['Bacteria (removed today)'] += 1
        elif '2157' in tids:
            c['Archaea'] += 1
        elif '4751' in tids:
            c['Fungi'] += 1
        elif '10239' in tids:
            c['Viruses'] += 1
        elif tids & {'554915', '33630', '5794', '5653'}:
            c['other eukaryote (protist)'] += 1
        elif '2759' in tids:
            c['Eukaryota, no kingdom-level OG'] += 1
        else:
            c['other / unassigned OG'] += 1
    for k, v in c.most_common():
        print(f'{g}\t{k}\t{v}')
PY

# --- (b) organellar genomes ---------------------------------------------------
declare -A ORG=( [rice]="NC_011033.1 NC_001320.1" [soybean]="NC_020455.1 NC_007942.1"
                 [dmel]="NC_024511.2" [cele]="NC_001328.1" )
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="
for g in rice soybean dmel cele; do
  out=org/$g.organelle.fa; : > $out
  for acc in ${ORG[$g]}; do
    curl -s --max-time 120 "$E$acc" >> $out
    sleep 1
  done
  echo "$g organelle bp: $(grep -v '>' $out | tr -d '\n' | wc -c)"
done

conda activate busco   # provides blastn / makeblastdb
declare -A GEN=( [rice]=IRGSP [soybean]=Wm82 [dmel]=dmel [cele]=cele )
echo -e "genome\tclass\tloci_tested\tloci_hit_50pct_CDS\tpct" > organellar.tsv
for g in rice soybean dmel cele; do
  D=$REV/regen/$g; [ -d $REV/regen/${g}_s ] && D=$REV/regen/${g}_s
  [ -s $D/union.final.gff3 ] || continue
  python3 $REV/extract_classes.py $R/$g/genome/${GEN[$g]}.fa $D/union.final.gff3 \
      $R/$g/genome/${GEN[$g]}.ref.gff* $W/seqs_$g 0 > /dev/null 2>&1 || true
  makeblastdb -in org/$g.organelle.fa -dbtype nucl -out org/$g.db > /dev/null 2>&1
  for cls in novel identical revised; do
    q=$W/seqs_$g/$cls.cds.fa
    [ -s "$q" ] || continue
    blastn -query $q -db org/$g.db -outfmt '6 qseqid qlen length pident' \
           -evalue 1e-10 -num_threads 16 -max_target_seqs 5 > $W/$g.$cls.blast 2>/dev/null
    python3 - $q $W/$g.$cls.blast $g $cls <<'PY' >> organellar.tsv
import sys
from collections import defaultdict
q, b, g, cls = sys.argv[1:5]
n = sum(1 for l in open(q) if l.startswith('>'))
cov = defaultdict(int); qlen = {}
for line in open(b):
    f = line.split('\t')
    cov[f[0]] += int(f[2]); qlen[f[0]] = int(f[1])
hit = sum(1 for k, v in cov.items() if qlen[k] and v >= 0.5 * qlen[k])
print(f"{g}\t{cls}\t{n}\t{hit}\t{100*hit/max(n,1):.2f}")
PY
  done
done
column -t eggnog_classes.tsv; column -t organellar.tsv
echo CONTAM_DONE
