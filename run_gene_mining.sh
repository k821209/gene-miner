#!/usr/bin/env bash
# run_gene_mining.sh — shell driver for MESSY / contaminated, fragmented
# assemblies: mine a clean gene catalog from a (non-model, possibly contaminated)
# genome using RNA-seq + ab-initio + protein evidence, then QC-filter out the junk.
# It uses per-contig AUGUSTUS + a RagTag AGP coordinate lift (see README "Lessons").
# Developed on real whole-body insect data and provided as-is; for finished
# chromosome-scale genomes the validated, recommended entry point is the
# one-command Nextflow pipeline `main.nf`.
#
#   Stages 1-4 (annotation):  final = usable AUGUSTUS (ab initio) ∪ RNA-only TransDecoder loci
#   Stages 5-7 (QC filters):  drop residual bacterial genes (eggNOG taxonomy) and
#                             TE-derived genes (RepeatMasker overlap) -> CLEAN catalog
#
# The QC pass is the point of "mining" vs "annotation": a raw union over-counts
# (here 14,311 -> clean 11,463; see README "Lessons"). Edit CONFIG, then:
#   bash run_gene_mining.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); BIN=$HERE/bin

############################## CONFIG ##############################
GENOME=ragtag_out/ragtag.scaffold.fasta      # annotation coordinate system (scaffold)
MASKED_CONTIGS=clean_masked_contigs.fa        # masked CONTIGS for AUGUSTUS
AGP=ragtag_out/ragtag.scaffold.agp            # contig->scaffold map
READS_DIR=rnaseq                              # holds <S>_1.fastq.gz / <S>_2.fastq.gz
SAMPLES="BS-1 BS-2 BS"
UNIPROT=db/uniprot-all.fasta                  # DIAMOND blastp DB
PFAM=db/Pfam-A.hmm                            # hmmscan DB (pressed)
AUG_SPECIES=BUSCO_Bs.genome.polished.fa.busco # BUSCO-retrained AUGUSTUS model (in AUGUSTUS_CONFIG_PATH/species)
AUG_SCORE=0.8; AUG_MINAA=100
REPEAT_LIB=                                   # de-novo RepeatModeler library (.fa); empty => build one (slow)
TE_THRESH=0.5                                 # drop a gene if >= this fraction of CDS is in interspersed repeats
RUN_QC=1                                      # 1 = run stages 5-7 (eggNOG + TE filters); 0 = stop at raw union
THREADS=36
# conda envs (adjust):  annot=hisat2,stringtie,transdecoder,samtools,hmmer,gffread ; nf=diamond ; augustus=augustus,evidencemodeler,miniprot
export PATH=$HOME/miniconda3/envs/annot/bin:$HOME/miniconda3/envs/nf/bin:$HOME/miniconda3/envs/augustus/bin:$PATH
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/augustus/config
UTIL=$HOME/miniconda3/envs/annot/opt/transdecoder/util
###################################################################

mkdir -p annot
echo "[1/5] RNA-seq evidence (HISAT2 -> StringTie -> TransDecoder + DIAMOND/Pfam)"
hisat2-build -p $THREADS "$GENOME" annot/idx >/dev/null 2>&1
for s in $SAMPLES; do
  hisat2 -p $THREADS --dta --max-intronlen 30000 -x annot/idx \
    -1 $READS_DIR/${s}_1.fastq.gz -2 $READS_DIR/${s}_2.fastq.gz 2>annot/$s.log \
    | samtools sort -@6 -o annot/$s.bam -; stringtie annot/$s.bam -p $THREADS -o annot/$s.gtf -l $s
done
ls annot/*.gtf > annot/gtf_list.txt
stringtie --merge -p $THREADS -o annot/merged.gtf annot/gtf_list.txt
( cd annot
  $UTIL/gtf_genome_to_cdna_fasta.pl merged.gtf ../$GENOME > transcripts.fa
  $UTIL/gtf_to_alignment_gff3.pl merged.gtf > merged.gff3
  TransDecoder.LongOrfs -t transcripts.fa >/dev/null 2>&1
  diamond makedb --in ../$UNIPROT -d ../db/uniprot -p $THREADS 2>/dev/null
  diamond blastp -q transcripts.fa.transdecoder_dir/longest_orfs.pep -d ../db/uniprot -p $THREADS -e 1e-5 -k1 --outfmt 6 -o blastp.outfmt6 2>/dev/null
  hmmscan --cpu $THREADS -E 1e-5 --domtblout pfam.domtblout ../$PFAM transcripts.fa.transdecoder_dir/longest_orfs.pep >/dev/null 2>&1
  TransDecoder.Predict -t transcripts.fa --retain_blastp_hits blastp.outfmt6 --retain_pfam_hits pfam.domtblout --cpu 16 >/dev/null 2>&1
  $UTIL/cdna_alignment_orf_to_genome_orf.pl transcripts.fa.transdecoder.gff3 merged.gff3 transcripts.fa > genome.transdecoder.gff3 2>/dev/null || true )

echo "[2/5] AUGUSTUS ab initio (BUSCO-trained model, per-contig parallel) on masked contigs"
rm -rf aug_split aug_out; mkdir -p aug_split aug_out
awk '/^>/{if(f)close(f); n=substr($1,2); sub(/ .*/,"",n); f="aug_split/"n".fa"} {print > f}' "$MASKED_CONTIGS"
ls aug_split/*.fa | xargs -P $THREADS -I {} bash -c 'b=$(basename "$1" .fa); augustus --strand=both --genemodel=complete --gff3=on --UTR=off --species='"$AUG_SPECIES"' "$1" > aug_out/"$b".gff3 2>/dev/null' _ {}
: > augustus_contigs.gff3
for g in aug_out/*.gff3; do b=$(basename "$g" .gff3); awk -v c="$b" '!/^#/&&NF>=8{gsub(/ID=g/,"ID="c"_g");gsub(/Parent=g/,"Parent="c"_g");gsub(/ID=t/,"ID="c"_t");gsub(/Parent=t/,"Parent="c"_t");print}' "$g" >> augustus_contigs.gff3; done
python3 $BIN/lift_agp.py "$AGP" augustus_contigs.gff3 augustus_scaffold.gff3

echo "[3/7] UNION: usable AUGUSTUS (score>=$AUG_SCORE,>=${AUG_MINAA}aa)  ∪  RNA-only loci"
python3 $BIN/build_union.py augustus_scaffold.gff3 annot/genome.transdecoder.gff3 $AUG_SCORE $AUG_MINAA | tee union_summary.txt

echo "[4/7] proteins (raw union)"
python3 $BIN/extract_pep.py "$GENOME" union.gff3 union.pep.fa

if [ "${RUN_QC:-1}" != 1 ]; then
  echo "DONE (RUN_QC=0). Raw union: union.gff3 / union.pep.fa  (see union_summary.txt)"
  echo "NB: a raw union OVER-COUNTS (residual bacterial + TE genes); set RUN_QC=1 for the clean catalog."
  exit 0
fi

echo "[5/7] QC-a: eggNOG taxonomy -> drop residual bacterial genes (gut microbiome surviving contig decontam)"
THREADS=$THREADS bash $BIN/run_eggnog.sh union.pep.fa union emap_out
python3 $BIN/filter_taxonomy.py emap_out/union.emapper.annotations \
  union.gff3 union.pep.fa union.clean.gff3 union.clean.pep.fa

echo "[6/7] QC-b: RepeatMasker overlap -> drop TE-derived genes (>= $TE_THRESH of CDS in interspersed repeats)"
PA=$((THREADS/4)) bash $BIN/run_repeatmasker.sh "$GENOME" "$REPEAT_LIB" rm_out
RM_OUT=rm_out/$(basename "$GENOME").out
python3 $BIN/te_filter_genes.py "$RM_OUT" \
  union.clean.gff3 union.clean.pep.fa union.final.gff3 union.final.pep.fa $TE_THRESH

echo "[7/7] BUSCO (protein) validation — run separately where busco+lineage live:"
echo "  busco -i union.final.pep.fa -l diptera_odb10 -m proteins -c $THREADS --offline --download_path <busco_downloads>"
echo "DONE. CLEAN catalog: union.final.gff3 / union.final.pep.fa"
echo "  raw union     -> union.gff3        ($(grep -c $'\tgene\t' union.gff3 2>/dev/null || echo '?') genes)"
echo "  -bacterial    -> union.clean.gff3  ($(grep -c $'\tgene\t' union.clean.gff3 2>/dev/null || echo '?') genes)"
echo "  -TE (final)   -> union.final.gff3  ($(grep -c $'\tgene\t' union.final.gff3 2>/dev/null || echo '?') genes)"
