#!/usr/bin/env nextflow
/*
 * rnaseq_annotation.nf — RNA-seq-guided gene annotation (Nextflow DSL2)
 *
 * Converted & modernized from the lab's Snakemake `stringtie2genemodel` pipeline
 * (HISAT2 -> sambamba sort/merge -> StringTie -> TransDecoder[+Pfam]).
 *
 * Additions over the original:
 *   - StringTie per-sample assembly + --merge (cleaner than merge-then-assemble)
 *   - FINAL EVIDENCE FILTER: keep only ORFs supported by protein homology
 *     (DIAMOND blastp) OR a Pfam domain (hmmscan) — removes spurious
 *     "protein-shaped" ORFs (TE-derived, random frames, fragments) that
 *     otherwise inflate the gene count.
 *
 * Usage:
 *   nextflow run rnaseq_annotation.nf \
 *     --genome clean_asm.fa --reads 'rnaseq/*_{1,2}.fastq.gz' \
 *     --proteome swissprot.fasta [--pfam Pfam-A.hmm] --outdir annot
 */
nextflow.enable.dsl=2

params.genome   = null
params.reads    = null            // glob: '..._{1,2}.fastq.gz'
params.proteome = null            // protein FASTA for DIAMOND blastp homology filter
params.pfam     = null            // optional Pfam-A.hmm; if set, hmmscan adds domain evidence
params.outdir   = 'annot_out'
params.max_intronlen = 30000
params.min_orf_aa    = 100        // min aa length to keep an evidence-less complete ORF (not used if strict)
params.strict        = true       // strict = keep ONLY homology/Pfam-supported genes

log.info """
  R N A - s e q   A N N O T A T I O N  (Nextflow)
  genome   : ${params.genome}
  reads    : ${params.reads}
  proteome : ${params.proteome}
  pfam     : ${params.pfam ?: '(skipped)'}
  outdir   : ${params.outdir}   strict=${params.strict}
"""

process HISAT2_BUILD {
  cpus 24
  input:  path genome
  output: tuple path(genome), path('idx*.ht2')
  script: "hisat2-build -p ${task.cpus} ${genome} idx"
}

process HISAT2_ALIGN {
  tag "$sid"
  cpus 36
  input:
    tuple val(sid), path(r1), path(r2)
    tuple path(genome), path(index)
  output: path "${sid}.bam"
  script: """
    hisat2 --dta --max-intronlen ${params.max_intronlen} -p ${task.cpus} -x idx \
      -1 ${r1} -2 ${r2} 2> ${sid}.hisat.log \
      | samtools sort -@ 6 -m 2G -o ${sid}.bam -
    samtools index ${sid}.bam
  """
}

process STRINGTIE_MERGE {
  cpus 24
  publishDir params.outdir, mode: 'copy'
  input:  path bams
  output: path 'merged.gtf'
  script: """
    for b in ${bams}; do stringtie \$b -p ${task.cpus} --dta -o \${b%.bam}.gtf -l \${b%.bam}; done
    ls *.gtf > gtf_list.txt
    stringtie --merge -p ${task.cpus} -o merged.gtf gtf_list.txt
  """
}

process TRANSCRIPTS {
  input:
    path merged_gtf
    path genome
  output:
    path 'transcripts.fa'
    path 'merged.gff3'
  script: """
    UTIL=\$(dirname \$(readlink -f \$(which TransDecoder.LongOrfs)))/../opt/transdecoder/util
    [ -d "\$UTIL" ] || UTIL=\$CONDA_PREFIX/opt/transdecoder/util
    \$UTIL/gtf_genome_to_cdna_fasta.pl ${merged_gtf} ${genome} > transcripts.fa
    \$UTIL/gtf_to_alignment_gff3.pl ${merged_gtf} > merged.gff3
  """
}

process LONGORFS {
  cpus 8
  input:  path transcripts
  output: tuple path(transcripts), path('transcripts.fa.transdecoder_dir/longest_orfs.pep')
  script: "TransDecoder.LongOrfs -t ${transcripts}"
}

process DIAMOND {
  cpus 36
  input:
    tuple path(transcripts), path(orfs)
    path proteome
  output: path 'blastp.outfmt6'
  script: """
    diamond makedb --in ${proteome} -d prot -p ${task.cpus} 2>/dev/null
    diamond blastp -q ${orfs} -d prot -p ${task.cpus} -e 1e-5 -k 1 \
      --outfmt 6 -o blastp.outfmt6 2>/dev/null
  """
}

process PFAM {
  cpus 24
  when: params.pfam
  input:
    tuple path(transcripts), path(orfs)
    path pfamdb
  output: path 'pfam.domtblout'
  script: "hmmscan --cpu ${task.cpus} --domtblout pfam.domtblout ${pfamdb} ${orfs} > /dev/null"
}

process PREDICT {
  cpus 16
  publishDir params.outdir, mode: 'copy'
  input:
    tuple path(transcripts), path(orfs)
    path blastp
    path pfam      // may be a stub if Pfam skipped
    path merged_gff3
  output:
    path 'transcripts.fa.transdecoder.{gff3,pep,cds,bed}'
    path 'genome.transdecoder.gff3', emit: genome_gff
    path 'blastp.outfmt6', emit: blastp
  script:
    def pfam_opt = (params.pfam ? "--retain_pfam_hits ${pfam}" : "")
    """
    UTIL=\$CONDA_PREFIX/opt/transdecoder/util
    TransDecoder.Predict -t ${transcripts} --single_best_only \
      --retain_blastp_hits ${blastp} ${pfam_opt} --cpu ${task.cpus}
    \$UTIL/cdna_alignment_orf_to_genome_orf.pl \
      transcripts.fa.transdecoder.gff3 ${merged_gff3} ${transcripts} > genome.transdecoder.gff3
    cp ${blastp} blastp.outfmt6
    """
}

process FILTER {
  publishDir params.outdir, mode: 'copy'
  input:
    path genome_gff
    path blastp
    path pfam
  output:
    path 'genes.filtered.gff3'
    path 'annotation_summary.txt'
  script: """
    python3 ${projectDir}/filter_models.py \
      --gff ${genome_gff} --blastp ${blastp} \
      ${ params.pfam ? "--pfam ${pfam}" : "" } \
      --min_aa ${params.min_orf_aa} ${ params.strict ? "--strict" : "" } \
      --out genes.filtered.gff3 --summary annotation_summary.txt
  """
}

workflow {
  genome = file(params.genome)
  read_pairs = Channel.fromFilePairs(params.reads, flat:true)
                      .map { id, r1, r2 -> tuple(id, r1, r2) }
  proteome = file(params.proteome)
  pfam_ch  = params.pfam ? file(params.pfam) : file("${projectDir}/NO_PFAM")

  idx   = HISAT2_BUILD(genome)
  bams  = HISAT2_ALIGN(read_pairs, idx).collect()
  gtf   = STRINGTIE_MERGE(bams)
  (tx, mg3) = TRANSCRIPTS(gtf, genome)
  orfs  = LONGORFS(tx)
  bp    = DIAMOND(orfs, proteome)
  pf    = params.pfam ? PFAM(orfs, pfam_ch) : Channel.value(file("${projectDir}/NO_PFAM"))
  pred  = PREDICT(orfs, bp, pf, mg3)
  FILTER(pred.genome_gff, pred.blastp, pf)
}
