#!/bin/bash
cd ~/hobak_comparative
export PATH=/home/k821209/miniconda3/envs/annot/bin:/home/k821209/miniconda3/envs/nf/bin:$PATH
UTIL=/home/k821209/miniconda3/envs/annot/opt/transdecoder/util
set -e
REF=ragtag_out/ragtag.scaffold.fasta
mkdir -p annot annot/hisat_idx
echo "[$(date +%T)] HISAT2 build"; [ -f annot/hisat_idx/idx.1.ht2 ] || hisat2-build -p 24 "$REF" annot/hisat_idx/idx >/dev/null 2>&1
for s in BS-1 BS-2 BS; do
  echo "[$(date +%T)] align $s"
  hisat2 -p 36 --dta --max-intronlen 30000 -x annot/hisat_idx/idx -1 rnaseq/${s}_1.fastq.gz -2 rnaseq/${s}_2.fastq.gz 2> annot/${s}.hisat.log | samtools sort -@ 6 -m 2G -o annot/${s}.bam -
  samtools index annot/${s}.bam; tail -1 annot/${s}.hisat.log
done
echo "[$(date +%T)] StringTie per-lib + merge (max-capture, keep isoforms)"
for s in BS-1 BS-2 BS; do stringtie annot/${s}.bam -p 24 --dta -o annot/${s}.gtf -l $s; done
ls annot/BS-1.gtf annot/BS-2.gtf annot/BS.gtf > annot/gtf_list.txt
stringtie --merge -p 24 -o annot/merged.gtf annot/gtf_list.txt
cd annot
echo "[$(date +%T)] transcripts + LongOrfs"
$UTIL/gtf_genome_to_cdna_fasta.pl merged.gtf ../$REF > transcripts.fa
$UTIL/gtf_to_alignment_gff3.pl merged.gtf > merged.gff3
TransDecoder.LongOrfs -t transcripts.fa >/dev/null 2>&1
NORF=$(grep -c '^>' transcripts.fa.transdecoder_dir/longest_orfs.pep)
echo "  candidate ORFs: $NORF"
echo "[$(date +%T)] DIAMOND blastp vs UniProt"
diamond makedb --in ../db/uniprot-all.fasta -d ../db/uniprot -p 36 2>/dev/null
diamond blastp -q transcripts.fa.transdecoder_dir/longest_orfs.pep -d ../db/uniprot -p 36 -e 1e-5 -k 1 --outfmt 6 -o blastp.outfmt6 2>/dev/null
echo "  blastp hits: $(cut -f1 blastp.outfmt6 | sort -u | wc -l)"
echo "[$(date +%T)] hmmscan vs Pfam-A (this is the slow step)"
hmmscan --cpu 38 -E 1e-5 --domtblout pfam.domtblout ../db/Pfam-A.hmm transcripts.fa.transdecoder_dir/longest_orfs.pep > /dev/null 2>&1
echo "  pfam hits: $(grep -v '^#' pfam.domtblout | awk '{print $4}' | sort -u | wc -l)"
echo "[$(date +%T)] TransDecoder.Predict (retain blastp + pfam; keep isoforms)"
TransDecoder.Predict -t transcripts.fa --retain_blastp_hits blastp.outfmt6 --retain_pfam_hits pfam.domtblout --cpu 16 >/dev/null 2>&1
$UTIL/cdna_alignment_orf_to_genome_orf.pl transcripts.fa.transdecoder.gff3 merged.gff3 transcripts.fa > genome.transdecoder.gff3 2>/dev/null || true
echo "[$(date +%T)] tiered filter (filter_models.py, blastp+pfam)"
python3 ../filter_models.py --gff genome.transdecoder.gff3 --blastp blastp.outfmt6 --pfam pfam.domtblout --out genes.highconf.gff3 --summary annotation_summary.txt
echo "[$(date +%T)] DONE"
