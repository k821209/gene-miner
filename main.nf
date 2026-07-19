#!/usr/bin/env nextflow
// main.nf — Gene-Miner end-to-end pipeline (DSL2) for finished chromosome-scale genomes.
// final = usable AUGUSTUS (score>=aug_score, >=aug_min_aa) U RNA-only TransDecoder loci,
// then QC: eggNOG taxonomy (-non-target) + RepeatMasker overlap (-TE).
// Tools are supplied per process by prepending the relevant conda env bin to PATH
// (deterministic; does not rely on Nextflow's conda integration).
nextflow.enable.dsl=2

def reqd(v,n){ if(!v) exit 1, "Missing required param --${n}"; return v }
// Coerce a param to boolean. A CLI flag like `--run_genemark false` arrives as the
// STRING "false", which is truthy in Groovy — so `if(params.x)` would wrongly fire.
// Treat only true/1/yes/on (any case) as true; everything else (incl. "false") false.
def asBool(v){
  if( v instanceof Boolean ) return v
  return v?.toString()?.toLowerCase() in ['true','1','yes','on']
}

process MASK_GENOME {
  tag "RepeatModeler/RepeatMasker"
  publishDir "${params.outdir}/mask", mode: "copy"
  cpus 12
  input:
    path genome
  output:
    path "genome.masked.fa", emit: fasta
    path "rm.out",           emit: rmout
  script:
  """
  REPEAT_ENV=${params.env_rmod.split('/')[-1]} PA=\$(( ${task.cpus} / 4 )) \
    bash ${projectDir}/bin/run_repeatmasker.sh ${genome} '${params.repeat_lib ?: ''}' rm_out
  cp rm_out/${genome}.masked genome.masked.fa
  cp rm_out/${genome}.out rm.out
  """
}

// BUSCO-train an AUGUSTUS model (genome --long --augustus) for organisms without a
// native model. Emits a writable augustus config dir holding the trained species.
process BUSCO_TRAIN {
  publishDir "${params.outdir}/augustus_model", mode: "copy"
  cpus 10
  input:
    path genome
  output:
    path "aug_config", emit: config
  script:
  """
  export PATH=${params.env_busco}/bin:\$PATH
  bash ${projectDir}/bin/busco_train_augustus.sh ${genome} ${params.busco_lineage} ${task.cpus} \
    ${params.env_augustus}/config aug_config > trained_species.txt
  """
}

process HISAT2_BUILD {
  input:
    path genome
  output:
    path "idx.*.ht2", emit: idx
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  hisat2-build -p ${task.cpus} ${genome} idx > build.log 2>&1
  """
}

process HISAT2_ALIGN {
  tag "${sid}"
  input:
    tuple val(sid), path(reads)
    path idx
  output:
    path "${sid}.bam", emit: bam
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  hisat2 -p ${task.cpus} --dta --max-intronlen 30000 -x idx -1 ${reads[0]} -2 ${reads[1]} 2> ${sid}.hisat.log \
    | samtools sort -@4 -m 2G -o ${sid}.bam -
  """
}

process STRINGTIE {
  tag "${bam.simpleName}"
  input:
    path bam
  output:
    path "${bam.simpleName}.gtf", emit: gtf
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  stringtie ${bam} -p ${task.cpus} -o ${bam.simpleName}.gtf -l ${bam.simpleName}
  """
}

process STRINGTIE_MERGE {
  input:
    path gtfs
    path genome
  output:
    path "transcripts.fa", emit: transcripts
    path "merged.gff3",    emit: aligngff
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  UTIL=${params.env_annot}/opt/transdecoder/util
  ls ${gtfs} > gtf_list.txt
  stringtie --merge -p ${task.cpus} -o merged.gtf gtf_list.txt
  \$UTIL/gtf_genome_to_cdna_fasta.pl merged.gtf ${genome} > transcripts.fa
  \$UTIL/gtf_to_alignment_gff3.pl merged.gtf > merged.gff3
  """
}

process TRANSDECODER {
  publishDir "${params.outdir}/rnaseq", mode: "copy"
  input:
    path transcripts
    path aligngff
    path proteome
  output:
    path "genome.transdecoder.gff3", emit: gff
  script:
  """
  UTIL=${params.env_annot}/opt/transdecoder/util
  export PATH=\$UTIL:${params.env_annot}/bin:${params.env_augustus}/bin:\$PATH
  # Fail loudly on a missing/empty proteome instead of hiding it (a dangling
  # --proteome symlink used to die inside diamond with the stderr swallowed).
  [ -s "${proteome}" ] || { echo "ERROR: --proteome '${proteome}' is missing or empty" >&2; exit 1; }
  TransDecoder.LongOrfs -t ${transcripts} > longorfs.log 2>&1
  diamond makedb --in ${proteome} -d prot -p ${task.cpus} > diamond_makedb.log 2>&1
  diamond blastp -q ${transcripts}.transdecoder_dir/longest_orfs.pep -d prot -p ${task.cpus} -e 1e-5 -k1 --outfmt 6 -o blastp.outfmt6 > diamond_blastp.log 2>&1
  TransDecoder.Predict -t ${transcripts} --retain_blastp_hits blastp.outfmt6 --single_best_only --cpu ${task.cpus} > predict.log 2>&1
  \$UTIL/cdna_alignment_orf_to_genome_orf.pl ${transcripts}.transdecoder.gff3 ${aligngff} ${transcripts} > genome.transdecoder.gff3 2> cdna_to_genome.log || true
  """
}

process AUGUSTUS_PREDICT {
  publishDir "${params.outdir}/augustus", mode: "copy"
  input:
    path masked
    path augconfig
    val species
  output:
    path "augustus_scaffold.gff3", emit: gff
  script:
  """
  export PATH=${params.env_augustus}/bin:\$PATH
  bash ${projectDir}/bin/run_augustus_windowed.sh ${masked} ${species} ${augconfig} \
    ${params.window_mb} ${params.window_overlap} ${task.cpus}
  """
}

process GENEMARK_ETP {
  publishDir "${params.outdir}/genemark", mode: "copy"
  input:
    path masked
    path bams
    path proteome
  output:
    path "genemark.gtf", emit: gtf
  script:
  """
  bash ${projectDir}/bin/run_genemark_etp.sh ${masked} ${proteome} ${task.cpus} '${bams}' genemark.gtf
  """
}

process BUILD_UNION {
  publishDir "${params.outdir}/union", mode: "copy"
  input:
    path genome
    path aug_gff
    path rna_gff
    path genemark_gtf
  output:
    path "union.gff3",   emit: gff
    path "union.pep.fa", emit: pep
  script:
  def gm_arg = genemark_gtf.name == 'NO_GENEMARK' ? '' : "--genemark ${genemark_gtf}"
  """
  export PATH=${params.env_annot}/bin:\$PATH
  [ -e augustus_scaffold.gff3 ] || cp ${aug_gff} augustus_scaffold.gff3
  mkdir -p annot && cp ${rna_gff} annot/genome.transdecoder.gff3
  python3 ${projectDir}/bin/build_union.py --prefix ${params.gene_prefix} ${gm_arg} > union_summary.txt 2>&1
  python3 ${projectDir}/bin/extract_pep.py ${genome} union.gff3 union.pep.fa
  """
}

process EGGNOG {
  publishDir "${params.outdir}/qc", mode: "copy"
  input:
    path pep
  output:
    path "union.emapper.annotations", emit: ann
  script:
  """
  export PATH=${params.env_eggnog}/bin:\$PATH
  emapper.py -i ${pep} --itype proteins -m diamond --dmnd_iterate no \
    --cpu ${task.cpus} --data_dir ${params.eggnog_db} --output_dir . -o union \
    --override --temp_dir . > emap.log 2>&1
  """
}

process FILTER_TAXONOMY {
  publishDir "${params.outdir}/qc", mode: "copy"
  input:
    path ann
    path gff
    path pep
  output:
    path "union.clean.gff3",   emit: gff
    path "union.clean.pep.fa", emit: pep
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  python3 ${projectDir}/bin/filter_taxonomy.py ${ann} ${gff} ${pep} union.clean.gff3 union.clean.pep.fa
  """
}

process TE_FILTER {
  publishDir "${params.outdir}", mode: "copy"
  input:
    path rmout
    path gff
    path pep
  output:
    path "union.final.gff3",   emit: gff
    path "union.final.pep.fa", emit: pep
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  python3 ${projectDir}/bin/te_filter_genes.py ${rmout} ${gff} ${pep} union.final.gff3 union.final.pep.fa ${params.te_thresh}
  """
}

process BUSCO {
  publishDir "${params.outdir}/busco", mode: "copy"
  input:
    path pep
  output:
    path "busco_union", emit: out
  script:
  """
  export PATH=${params.env_busco}/bin:\$PATH
  busco -i ${pep} -l ${params.busco_lineage} -m proteins -c ${task.cpus} -o busco_union -f > busco.log 2>&1
  """
}

workflow {
  genome   = file(reqd(params.genome,"genome"))
  proteome = file(reqd(params.proteome,"proteome"))
  reads_ch = Channel.fromFilePairs(reqd(params.reads,"reads"))

  MASK_GENOME(genome)
  HISAT2_BUILD(genome)
  HISAT2_ALIGN(reads_ch, HISAT2_BUILD.out.idx.collect())
  STRINGTIE(HISAT2_ALIGN.out.bam)
  STRINGTIE_MERGE(STRINGTIE.out.gtf.collect(), genome)
  TRANSDECODER(STRINGTIE_MERGE.out.transcripts, STRINGTIE_MERGE.out.aligngff, proteome)

  // AUGUSTUS model: native species model if given (e.g. rice), else BUSCO-trained
  def use_busco = !params.augustus_species || params.augustus_species == 'auto'
  aug_cfg = use_busco ? BUSCO_TRAIN(genome).config : file("${params.env_augustus}/config")
  aug_sp  = use_busco ? 'auto' : params.augustus_species
  AUGUSTUS_PREDICT(MASK_GENOME.out.fasta, aug_cfg, aug_sp)

  // optional 3rd ab-initio stream: GeneMark-ETP (precomputed gtf, or run in-pipeline)
  no_gm = file("${projectDir}/assets/NO_GENEMARK")
  if( params.genemark_gtf ) {
    gm_ch = Channel.value(file(params.genemark_gtf))
  } else if( asBool(params.run_genemark) ) {
    GENEMARK_ETP(MASK_GENOME.out.fasta, HISAT2_ALIGN.out.bam.collect(), proteome)
    gm_ch = GENEMARK_ETP.out.gtf
  } else {
    gm_ch = Channel.value(no_gm)
  }
  BUILD_UNION(genome, AUGUSTUS_PREDICT.out.gff, TRANSDECODER.out.gff, gm_ch)

  if( asBool(params.run_qc) ) {
    EGGNOG(BUILD_UNION.out.pep)
    FILTER_TAXONOMY(EGGNOG.out.ann, BUILD_UNION.out.gff, BUILD_UNION.out.pep)
    TE_FILTER(MASK_GENOME.out.rmout, FILTER_TAXONOMY.out.gff, FILTER_TAXONOMY.out.pep)
    BUSCO(TE_FILTER.out.pep)
  }
}
