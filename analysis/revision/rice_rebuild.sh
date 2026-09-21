#!/usr/bin/env bash
# Rebuild the native-model rice catalogue after AUGUSTUS was re-run with
# --sample=100 (so its gene scores are real posteriors and the 0.8 filter acts),
# then recompute every rice-dependent number in the paper.
# Waits for revision/rice_aug_sample/augustus_scaffold.gff3.
set -uo pipefail
R=~/gene-miner-runs; BIN=$R/gene-miner/bin; AN=$R/gene-miner/analysis; REV=$R/revision
REF=$R/rice/genome/IRGSP.ref.gff3; GEN=$R/rice/genome/IRGSP.fa
GMK=$R/rice/etp_run/genemark.gtf
OLD=$REV/regen/rice                       # rebuild with the unsampled AUGUSTUS, kept for comparison
W=$REV/regen/rice_s; mkdir -p $W/annot
source ~/miniconda3/etc/profile.d/conda.sh
log() { echo "[$(date +%T)] $*"; }

until [ -s $REV/rice_aug_sample/augustus_scaffold.gff3 ] && ! ps -eo args | grep -q "[r]un_augustus_windowed"; do sleep 60; done
log "AUGUSTUS done: $(awk -F'\t' '$3=="gene"' $REV/rice_aug_sample/augustus_scaffold.gff3 | wc -l) genes"

cd $W
ln -sf $REV/rice_aug_sample/augustus_scaffold.gff3 augustus_scaffold.gff3
ln -sf $R/rice/native_etp/annot/genome.transdecoder.gff3 annot/genome.transdecoder.gff3
python3 $BIN/build_union.py --prefix GM --genemark $GMK > union.log 2>&1
python3 $BIN/extract_pep.py $GEN union.gff3 union.pep.fa > pep.log 2>&1
python3 $BIN/te_filter_genes.py $R/rice/gm_out/mask/rm.out union.gff3 union.pep.fa union.te.gff3 union.te.pep.fa > te.log 2>&1

log "eggNOG (transfer by sequence from both earlier rice rebuilds, emapper on the rest)"
cat $R/rice/native_etp/union.pep.fa $OLD/union.pep.fa > prev.pep.fa
{ head -n 5 $OLD/new.annotations | grep '^#'; grep -hv '^#' $R/rice/native_etp/emap/union.emapper.annotations $OLD/new.annotations; } > prev.annotations
python3 $REV/transfer_eggnog.py prev.pep.fa prev.annotations union.pep.fa new.annotations todo.faa > transfer.log 2>&1
cat transfer.log
[ -s todo.faa ] && { bash $BIN/run_eggnog.sh todo.faa todo emap_new > eggnog.log 2>&1; grep -v '^#' emap_new/todo.emapper.annotations >> new.annotations 2>/dev/null; }
python3 $BIN/filter_taxonomy.py new.annotations union.te.gff3 union.te.pep.fa union.final.gff3 union.final.pep.fa > tax.log 2>&1

python3 - <<'PY'
import re
from collections import defaultdict
g2t=defaultdict(list); tl=defaultdict(int)
for l in open("union.final.gff3"):
    f=l.rstrip().split("\t")
    if len(f)<9: continue
    if f[2]=="mRNA": g2t[re.search(r"Parent=([^;]+)",f[8]).group(1)].append(re.search(r"ID=([^;]+)",f[8]).group(1))
    elif f[2]=="CDS": tl[re.search(r"Parent=([^;]+)",f[8]).group(1)]+=int(f[4])-int(f[3])+1
prim={max(ts,key=lambda t:tl[t]) for ts in g2t.values()}
keep=False; o=open("union.final.primary.pep.fa","w")
for l in open("union.final.pep.fa"):
    if l[0]==">": keep=l[1:].split()[0] in prim
    if keep: o.write(l)
PY

log "reference comparison, boundary errors, provenance"
python3 $BIN/gm_compare.py --gm union.final.gff3 --ref $REF --min_ro 0.5 --out gm_compare.tsv 2>/dev/null
python3 $REV/split_fusion.py $REF union.final.gff3 0.5 > split_fusion.txt
python3 $REV/fusion_origin.py $REF union.final.gff3 augustus_scaffold.gff3 annot/genome.transdecoder.gff3 $GMK > fusion_origin.txt 2>&1
python3 $REV/isoform_origin.py union.final.gff3 augustus_scaffold.gff3 annot/genome.transdecoder.gff3 $GMK > isoform_origin.txt 2>/dev/null
python3 $REV/novel_by_stream.py union.final.gff3 $REF > novel_by_stream.txt 2>&1
python3 $AN/classify_missed.py union.gff3 union.final.gff3 $REF > classify_missed.txt 2>&1
python3 $REV/revised_features.py union.final.gff3 $REF revised_features.tsv > revised_features.txt 2>&1
python3 /tmp/novel_extra.py union.final.gff3 $REF new.annotations > novel_extra.txt 2>&1
python3 $AN/revised_junctions.py union.final.gff3 $REF $REV/junctions/rice.junctions.tsv 3 > revised_junctions.txt 2>&1
python3 $AN/added_intron_streams.py union.final.gff3 $REF $REV/junctions/rice.junctions.tsv \
    annot/genome.transdecoder.gff3 $GMK augustus_scaffold.gff3 3 > added_intron_streams.txt 2>&1
conda activate gffcmp
python3 $AN/cds_to_exon_gtf.py union.final.gff3 gm.cds.gtf 2>/dev/null
gffcompare -r $REV/structural/rice.ref.cds.gtf -o cmp gm.cds.gtf > /dev/null 2>&1
conda deactivate

log "Swiss-Prot hits and evidence tiers"
export PATH=$HOME/miniconda3/envs/augustus/bin:$PATH
diamond blastp -q union.final.pep.fa -d $REV/sprot.dmnd -p 12 -e 1e-5 -k1 --outfmt 6 -o sprot_hits.outfmt6 > diamond.log 2>&1
python3 $REV/annotate_tiers.py union.final.gff3 union.final.tiers.gff3 --td annot/genome.transdecoder.gff3 \
    --eggnog new.annotations --sprot sprot_hits.outfmt6 --ref $REF --summary tiers.tsv > tiers.log 2>&1

log "BUSCO"
conda activate busco
for set in primary alliso; do
  IN=union.final.primary.pep.fa; [ $set = alliso ] && IN=union.final.pep.fa
  busco -i $IN -l poales_odb10 -m proteins -o busco_$set -f --offline \
    --download_path $R/rice/gm_out/busco_downloads -c 12 > busco_$set.log 2>&1
done
conda deactivate

log "Ka/Ks against the two wild Oryza genomes"
K=$REV/kaks_s; mkdir -p $K/seqs; cd $K
ln -sf $REV/kaks/genomes genomes
python3 - $W/union.gff3 $W/union.te.gff3 $K/seqs/te_removed.gff3 <<'PY'
import re, sys
raw, kept, out = sys.argv[1:4]
keep = {re.search(r'ID=([^;]+)', l.split('\t')[8]).group(1)
        for l in open(kept) if len(l.split('\t')) > 8 and l.split('\t')[2] == 'gene'}
tx, w = set(), open(out, 'w')
for l in open(raw):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    if f[2] == 'gene':
        if re.search(r'ID=([^;]+)', f[8]).group(1) not in keep: w.write(l)
    elif f[2] == 'mRNA':
        if re.search(r'Parent=([^;]+)', f[8]).group(1) not in keep:
            tx.add(re.search(r'ID=([^;]+)', f[8]).group(1)); w.write(l)
    elif f[2] == 'CDS' and re.search(r'Parent=([^;]+)', f[8]).group(1) in tx:
        w.write(l)
PY
python3 $REV/extract_classes.py $GEN $W/union.final.gff3 $REF $K/seqs 3000 $K/seqs/te_removed.gff3
conda activate kaks
export PATH=$HOME/miniconda3/envs/busco/bin:$PATH
for sp in glaberrima punctata; do
  case $sp in
    glaberrima) FA=$REV/kaks/genomes/Oryza_glaberrima.Oryza_glaberrima_V1.dna.toplevel.fa;;
    punctata)   FA=$REV/kaks/genomes/Oryza_punctata.Oryza_punctata_v1.2.dna.toplevel.fa;;
  esac
  for cls in novel identical revised te_removed; do
    [ -s $K/seqs/$cls.pep.fa ] || continue
    python3 $REV/kaks_pipeline.py $K/seqs/$cls.cds.fa $K/seqs/$cls.pep.fa $FA $K/$sp.$cls 20 > $K/$sp.$cls.log 2>&1
    log "  $sp $cls: $(($(wc -l < $K/$sp.$cls/kaks.tsv)-1)) pairs"
  done
done
conda deactivate
python3 $REV/kaks_by_tier.py $W/union.final.tiers.gff3 $K > $K/by_tier.txt
cat $K/by_tier.txt
echo RICE_REBUILD_DONE
