#!/usr/bin/env nextflow
/*
 * main.nf — Zeugodacus scutellatus genome annotation (union-based) — NEXTFLOW TEMPLATE
 *
 * This is a DSL2 skeleton that documents the workflow as Nextflow processes.
 * The PRODUCTION, battle-tested implementation is the shell driver
 * `run_annotation.sh` calling the scripts in bin/ (all verified on real data).
 * Use this .nf as a starting point to port the driver to a Nextflow cluster run.
 *
 * Final gene set = UNION (not evidence-only consensus):
 *   final = usable_AUGUSTUS(score>=0.8, >=100aa)  ∪  RNA-only TransDecoder loci
 * See ../README.md "Lessons" for why (AUGUSTUS and RNA-seq each catch genes the
 * other misses; EVM consensus alone under-counts by dropping ab-initio-only genes).
 */
nextflow.enable.dsl=2

params.genome           = null   // RagTag scaffold (annotation coordinate system)
params.masked_contigs   = null   // masked CONTIGS for AUGUSTUS
params.agp              = null   // ragtag.scaffold.agp (AUGUSTUS contig->scaffold lift)
params.reads            = null   // RNA-seq '..._{1,2}.fastq.gz'
params.uniprot          = null
params.pfam             = null
params.augustus_species = null   // BUSCO-retrained AUGUSTUS species model name
params.aug_score        = 0.8
params.aug_min_aa       = 100
params.outdir           = 'annotation_out'

process RNASEQ_TRANSDECODER {
  publishDir params.outdir, mode:'copy'
  input: path genome; path reads; path uniprot; path pfam
  output: path 'genome.transdecoder.gff3', emit: rna_gff; path 'merged.gtf'
  script: "bash ${projectDir}/bin/run_rnaseq_transdecoder.sh ${genome} '${reads}' ${uniprot} ${pfam}"
}
process AUGUSTUS_ABINITIO {
  publishDir params.outdir, mode:'copy'
  input: path masked_contigs; path agp; val species
  output: path 'augustus_scaffold.gff3', emit: aug_gff
  script: """
    bash ${projectDir}/bin/run_augustus.sh ${masked_contigs} ${species}
    python3 ${projectDir}/bin/lift_agp.py ${agp} augustus_contigs.gff3 augustus_scaffold.gff3
  """
}
process UNION {
  publishDir params.outdir, mode:'copy'
  input: path genome; path aug_gff; path rna_gff
  output: path 'union.gff3'; path 'union.pep.fa'; path 'union_summary.txt'
  script: """
    python3 ${projectDir}/bin/build_union.py ${aug_gff} ${rna_gff} ${params.aug_score} ${params.aug_min_aa} > union_summary.txt
    python3 ${projectDir}/bin/extract_pep.py ${genome} union.gff3 union.pep.fa
  """
}

workflow {
  g = file(params.genome)
  rna = RNASEQ_TRANSDECODER(g, Channel.fromPath(params.reads).collect(),
                            file(params.uniprot), params.pfam ? file(params.pfam):file('NO_PFAM'))
  aug = AUGUSTUS_ABINITIO(file(params.masked_contigs), file(params.agp), params.augustus_species)
  UNION(g, aug.aug_gff, rna.rna_gff)
}
