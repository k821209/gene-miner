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

// Derive the injection pair from a PUBLISHED repeat GFF (--repeat_gff), so a
// caller can hand over the file a genome project already distributes instead of
// converting it themselves. The conversion has three silent-wrong-answer traps
// (the (left) column, the class column, GFF2 vs GFF3 attributes) documented in
// bin/repeat_gff_to_rmout.py; doing it here means the pipeline's own
// VALIDATE_MASKED_INPUT gate then runs over the result, so the caller gets the
// same guarantees as a hand-supplied pair.
process REPEAT_GFF_TO_RMOUT {
  tag "convert published repeat GFF"
  publishDir "${params.outdir}/mask", mode: "copy", pattern: "repeat_gff_convert.log"
  input:
    path genome
    path repeat_gff
  output:
    path "from_gff.softmasked.fa", emit: fasta
    path "from_gff.out",           emit: rmout
    path "repeat_gff_convert.log"
  script:
  """
  python3 ${projectDir}/bin/repeat_gff_to_rmout.py ${genome} ${repeat_gff} from_gff \
    > repeat_gff_convert.log 2>&1 || { echo "ERROR: repeat GFF conversion failed" >&2; cat repeat_gff_convert.log >&2; exit 1; }
  cat repeat_gff_convert.log
  ROWS=\$(awk 'NF>=15 && \$1+0>0' from_gff.out | wc -l)
  echo "converted .out data rows: \$ROWS"
  # An empty .out would pass VALIDATE_MASKED_INPUT -- a repeat-poor genome is
  # legitimate -- but here it means the GFF's seqids did not match the genome's.
  # That is the one outcome that looks like success and is not, so it stops here.
  [ "\$ROWS" -gt 0 ] || { echo "ERROR: conversion produced no intervals -- do the GFF seqids match the genome? See the dropped-row warning above." >&2; exit 1; }
  """
}

// Record WHICH versions produced the output. The gene count and the TE filter
// move with the RepeatMasker build and the Dfam release, the AUGUSTUS model
// source and the eggNOG DB version, so "annotated with Gene-Miner" is not a
// reproducible statement and a reader asking "which Dfam?" later has no other
// way to find out. Runs last, taking the final GFF so it cannot run early.
// Deliver an UPDATED reference annotation, not only a parallel catalogue.
// --ref_gff names the annotation the caller already has; the reference is
// preserved and our results are merged onto it (alternative transcripts where
// the structure differs, new genes where there is no counterpart). The locus
// pairing comes from gm_compare's --pairs so the matching rule is implemented
// once. See bin/merge_into_reference.py for what it does to the reference and
// what it refuses to do.
process MERGE_REFERENCE {
  tag "supplement the reference annotation"
  publishDir "${params.outdir}", mode: "copy"
  input:
    path final_gff
    path tiers_gff
    path ref_gff
  output:
    path "reference_updated.gff3"
    path "merge_report.tsv"
    path "gm_compare_vs_reference.tsv"
    path "locus_pairs.tsv"
  script:
  """
  python3 ${projectDir}/bin/gm_compare.py --gm ${final_gff} --ref ${ref_gff} \
    --out gm_compare_vs_reference.tsv --pairs locus_pairs.tsv
  python3 ${projectDir}/bin/merge_into_reference.py --ref ${ref_gff} --gm ${final_gff} \
    --pairs locus_pairs.tsv --tiers ${tiers_gff} \
    --out reference_updated.gff3 --report merge_report.tsv
  # The claim the deliverable rests on: everything that is not ours reproduces
  # the input. Checked here, not left to the reader -- the only reference rows
  # that may differ are gene lines widened to span an added transcript.
  awk -F'\t' '\$2!="gene-miner"' reference_updated.gff3 | grep -v '^#' > .recon
  grep -v '^#' ${ref_gff} | awk -F'\t' 'NF>=9' > .orig
  CH=\$(diff .orig .recon | grep -c '^>' || true)
  WD=\$(diff .orig .recon | grep '^>' | grep -c gm_widened_for || true)
  echo "reference rows differing: \$CH, of which widened: \$WD"
  [ "\$CH" = "\$WD" ] || { echo "ERROR: \$((CH-WD)) reference row(s) changed for a reason other than widening" >&2; exit 1; }
  rm -f .recon .orig
  """
}

// Splice junctions with read counts, from the alignments this run already
// produced. Computed once and shared, because both the catalogue and the
// supplemented reference are measured against the same evidence.
process RNA_JUNCTIONS {
  tag "splice junctions from the run's own BAMs"
  publishDir "${params.outdir}/rnaseq", mode: "copy"
  input:
    path bams
  output:
    path "rnaseq_junctions.tsv", emit: tsv
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  export SAMTOOLS=${params.env_annot}/bin/samtools
  bash ${projectDir}/bin/junctions_from_bam.sh rnaseq_junctions.tsv ${task.cpus} ${bams}
  """
}

// Per-TRANSCRIPT structural support. The tier annotation is locus-level and
// says which KINDS of evidence a gene has; this says whether THIS transcript's
// structure is the one the reads support, which is the question asked of a
// revised model. Applied to whatever GFF it is given, so the supplemented
// reference gets it on BOTH sides -- annotating only the added transcripts
// would tilt every comparison toward the side carrying numbers.
process RNA_SUPPORT {
  tag { label }
  publishDir "${params.outdir}", mode: "copy"
  input:
    // One invocation over a channel of (gff, label) pairs: DSL2 forbids calling
    // a process twice, and the catalogue and the supplemented reference both
    // need the same measure from the same junctions.
    tuple path(gff), val(label)
    path junctions
  output:
    path "${label}.rna_support.gff3", emit: gff
    path "${label}.rna_support_report.tsv"
    path "${label}.rna_support_comparison.tsv", optional: true
  script:
  """
  python3 ${projectDir}/bin/annotate_rna_support.py \
    --gff ${gff} --junctions ${junctions} \
    --out ${label}.rna_support.gff3 \
    --report ${label}.rna_support_report.tsv \
    --compare ${label}.rna_support_comparison.tsv \
    --min-reads ${params.rna_min_junction_reads}
  """
}

process PROVENANCE {
  publishDir "${params.outdir}", mode: "copy"
  input:
    path final_gff
    val  aug_model
    // BUSCO's own output directory, so the completeness score lands IN the
    // provenance file. The first real run of this process produced a
    // provenance with an empty BUSCO line, which is the one number a reader
    // looks for first. `assets/NO_BUSCO` stands in when --run_qc is off.
    path busco_dir
  output:
    path "provenance.txt"
  script:
  """
  export GM_PV_PROJECTDIR='${projectDir}'
  export GM_PV_ANNOT='${params.env_annot}'       GM_PV_AUGUSTUS='${params.env_augustus}'
  export GM_PV_RMOD='${params.masked_genome || params.repeat_gff ? "" : params.env_rmod}'
  export GM_PV_BUSCO='${params.env_busco}'
  export GM_PV_EGGNOG='${asBool(params.run_qc) ? params.env_eggnog : ""}'
  export GM_PV_EGGNOG_DB='${asBool(params.run_qc) ? params.eggnog_db : ""}'
  export GM_PV_GENEMARK_ETP='${asBool(params.run_genemark) ? params.genemark_etp_dir : ""}'
  export GM_PV_LINEAGE='${params.busco_lineage}'
  export GM_PV_AUG_MODEL='${aug_model}'
  export GM_PV_GENOME='${params.genome}'         GM_PV_PROTEOME='${params.proteome}'
  export GM_PV_REPEAT_LIB='${params.repeat_lib ?: ""}'
  export GM_PV_MASKED_GENOME='${params.masked_genome ?: ""}'
  export GM_PV_REPEAT_GFF='${params.repeat_gff ?: ""}'
  export GM_PV_READS='${params.reads}'
  export PATH=/usr/bin:/bin:\$PATH
  BUSCO_ARG=${busco_dir}
  [ -d "\$BUSCO_ARG" ] || BUSCO_ARG=""
  bash ${projectDir}/bin/write_provenance.sh provenance.txt ${final_gff} \$BUSCO_ARG
  """
}

process VALIDATE_MASKED_INPUT {
  // Guards the --masked_genome/--masked_rmout injection path. When MASK_GENOME
  // derives the masked genome from `genome` itself, coordinate agreement is
  // structural; here the caller supplies both files from elsewhere, and a
  // seqid/assembly mismatch or a hard- instead of soft-masked file would
  // otherwise fail silently downstream (plausible gene counts and BUSCO
  // scores from two coordinate systems merged by BUILD_UNION; an empty
  // TE_FILTER overlap that reads as "no TE-derived genes"). See
  // bin/validate_masked_input.py for exactly what it checks.
  tag "coordinate check"
  input:
    path genome
    path masked_genome
    path masked_rmout
  output:
    path "masked_genome.checked.fa", emit: fasta
    path "masked_rmout.checked.out", emit: rmout
  script:
  """
  python3 ${projectDir}/bin/validate_masked_input.py ${genome} ${masked_genome} ${masked_rmout}
  cp ${masked_genome} masked_genome.checked.fa
  cp ${masked_rmout} masked_rmout.checked.out
  """
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
  REPEAT_ENV=${params.env_rmod} PA=\$(( ${task.cpus} / 4 )) \
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
  # perl (+CPAN) and python3 for gmetp.pl come from the `genemark` conda env;
  # GeneMark-ETP itself (bundled tools + ProtHint) lives in genemark_etp_dir.
  export PATH=${params.env_genemark}/bin:\$PATH
  export GENEMARK_ETP_DIR=${params.genemark_etp_dir}
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
  python3 ${projectDir}/bin/build_union.py --prefix ${params.gene_prefix} ${gm_arg} --min-score ${params.aug_score} --min-aa ${params.aug_min_aa} > union_summary.txt 2>&1
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

process SPROT_HITS {
  publishDir "${params.outdir}/qc", mode: "copy"
  input:
    path pep
    path proteome
  output:
    path "sprot_hits.outfmt6", emit: hits
  script:
  """
  export PATH=${params.env_augustus}/bin:\$PATH
  diamond makedb --in ${proteome} -d prot -p ${task.cpus} > makedb.log 2>&1
  diamond blastp -q ${pep} -d prot -p ${task.cpus} -e 1e-5 -k1 --outfmt 6 \
    -o sprot_hits.outfmt6 > blastp.log 2>&1
  """
}

process ANNOTATE_TIERS {
  publishDir "${params.outdir}", mode: "copy"
  input:
    path gff
    path td
    path ann
    path hits
  output:
    path "union.final.tiers.gff3", emit: gff
    path "confidence_tiers.tsv",   emit: tsv
  script:
  """
  export PATH=${params.env_annot}/bin:\$PATH
  python3 ${projectDir}/bin/annotate_tiers.py ${gff} union.final.tiers.gff3 \
    --td ${td} --eggnog ${ann} --sprot ${hits} --summary confidence_tiers.tsv > tiers.log 2>&1
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
  run_busco_proteins.sh ${pep} ${params.busco_lineage} ${task.cpus} busco_union
  """
}

// Run bin/preflight.sh before any task is submitted. Everything it checks is
// cheap, and every one of the failures it covers otherwise arrives hours in,
// after the expensive stages have been paid for: a wrong GM_CONDA_BASE surfaces
// at BUSCO_TRAIN, a missing eggNOG DB and an emapper that cannot start surface
// after BUILD_UNION, an unreachable ezlab surfaces after 17 h of alignment and
// assembly. Skippable with --skip_preflight for a host where a check itself is
// the problem; the skip is printed, so it cannot be mistaken for a pass.
def preflight() {
  if( asBool(params.skip_preflight) ) {
    log.warn "preflight SKIPPED (--skip_preflight). Nothing below was verified."
    return
  }
  def script = "${projectDir}/bin/preflight.sh"
  if( !file(script).exists() ) {
    log.warn "preflight SKIPPED: ${script} not found (older checkout?). Nothing was verified."
    return
  }
  def args = ['bash', script,
    '--genome',            params.genome?.toString() ?: '',
    '--reads',             params.reads?.toString() ?: '',
    '--proteome',          params.proteome?.toString() ?: '',
    '--busco-lineage',     params.busco_lineage?.toString() ?: '',
    '--env-annot',         params.env_annot?.toString() ?: '',
    '--env-augustus',      params.env_augustus?.toString() ?: '',
    '--env-rmod',          params.env_rmod?.toString() ?: '',
    '--env-busco',         params.env_busco?.toString() ?: '',
    '--env-eggnog',        params.env_eggnog?.toString() ?: '',
    '--env-genemark',      params.env_genemark?.toString() ?: '',
    '--genemark-etp-dir',  params.genemark_etp_dir?.toString() ?: '',
    '--eggnog-db',         params.eggnog_db?.toString() ?: '',
    '--run-genemark',      asBool(params.run_genemark) ? '1' : '0',
    '--run-qc',            asBool(params.run_qc) ? '1' : '0']
  if( params.augustus_species ) args += ['--augustus-species', params.augustus_species.toString()]
  if( params.masked_genome )    args += ['--masked-genome',    params.masked_genome.toString()]
  if( params.masked_rmout )     args += ['--masked-rmout',     params.masked_rmout.toString()]
  if( params.repeat_lib )       args += ['--repeat-lib',       params.repeat_lib.toString()]
  if( params.repeat_gff )       args += ['--repeat-gff',       params.repeat_gff.toString()]
  if( params.ref_gff )          args += ['--ref-gff',          params.ref_gff.toString()]

  def out = new StringBuilder()
  def err = new StringBuilder()
  def proc = args.execute()
  proc.consumeProcessOutput(out, err)
  proc.waitForOrKill(180_000)
  if( out ) print out.toString()
  if( err ) System.err.print err.toString()
  if( proc.exitValue() != 0 )
    exit 1, "preflight failed -- nothing was submitted. Fix the FAIL lines above, or pass --skip_preflight to run anyway."
}

workflow {
  preflight()
  genome   = file(reqd(params.genome,"genome"))
  proteome = file(reqd(params.proteome,"proteome"))
  reads_ch = Channel.fromFilePairs(reqd(params.reads,"reads"))

  // A caller who already has a soft-masked genome (and its RepeatMasker .out,
  // needed by TE_FILTER) can supply both and skip RepeatMasker/RepeatModeler
  // entirely. --repeat_lib only skips de-novo RepeatModeler library building
  // inside run_repeatmasker.sh; RepeatMasker itself still runs unless the
  // genome is injected here.
  if( params.masked_genome ) {
    VALIDATE_MASKED_INPUT(genome, file(params.masked_genome), file(reqd(params.masked_rmout,"masked_rmout")))
    masked_fa    = VALIDATE_MASKED_INPUT.out.fasta
    masked_rmout = VALIDATE_MASKED_INPUT.out.rmout
  } else if( params.repeat_gff ) {
    // One published file in; the converted pair goes through the same gate.
    CONV = REPEAT_GFF_TO_RMOUT(genome, file(params.repeat_gff))
    VALIDATE_MASKED_INPUT(genome, CONV.fasta, CONV.rmout)
    masked_fa    = VALIDATE_MASKED_INPUT.out.fasta
    masked_rmout = VALIDATE_MASKED_INPUT.out.rmout
  } else {
    MASK_GENOME(genome)
    masked_fa    = MASK_GENOME.out.fasta
    masked_rmout = MASK_GENOME.out.rmout
  }
  HISAT2_BUILD(genome)
  HISAT2_ALIGN(reads_ch, HISAT2_BUILD.out.idx.collect())
  STRINGTIE(HISAT2_ALIGN.out.bam)
  STRINGTIE_MERGE(STRINGTIE.out.gtf.collect(), genome)
  TRANSDECODER(STRINGTIE_MERGE.out.transcripts, STRINGTIE_MERGE.out.aligngff, proteome)

  // AUGUSTUS model: native species model if given (e.g. rice), else BUSCO-trained
  def use_busco = !params.augustus_species || params.augustus_species == 'auto'
  aug_cfg = use_busco ? BUSCO_TRAIN(genome).config : file("${params.env_augustus}/config")
  aug_sp  = use_busco ? 'auto' : params.augustus_species
  AUGUSTUS_PREDICT(masked_fa, aug_cfg, aug_sp)

  // optional 3rd ab-initio stream: GeneMark-ETP (precomputed gtf, or run in-pipeline)
  no_gm = file("${projectDir}/assets/NO_GENEMARK")
  no_busco = file("${projectDir}/assets/NO_BUSCO")
  if( params.genemark_gtf ) {
    gm_ch = Channel.value(file(params.genemark_gtf))
  } else if( asBool(params.run_genemark) ) {
    GENEMARK_ETP(masked_fa, HISAT2_ALIGN.out.bam.collect(), proteome)
    gm_ch = GENEMARK_ETP.out.gtf
  } else {
    gm_ch = Channel.value(no_gm)
  }
  BUILD_UNION(genome, AUGUSTUS_PREDICT.out.gff, TRANSDECODER.out.gff, gm_ch)

  if( asBool(params.run_qc) ) {
    EGGNOG(BUILD_UNION.out.pep)
    FILTER_TAXONOMY(EGGNOG.out.ann, BUILD_UNION.out.gff, BUILD_UNION.out.pep)
    TE_FILTER(masked_rmout, FILTER_TAXONOMY.out.gff, FILTER_TAXONOMY.out.pep)
    SPROT_HITS(TE_FILTER.out.pep, proteome)
    ANNOTATE_TIERS(TE_FILTER.out.gff, TRANSDECODER.out.gff, EGGNOG.out.ann, SPROT_HITS.out.hits)
    BUSCO(TE_FILTER.out.pep)
    PROVENANCE(TE_FILTER.out.gff, aug_sp, BUSCO.out.out)
    if( params.ref_gff )
      MERGE_REFERENCE(TE_FILTER.out.gff, ANNOTATE_TIERS.out.gff, file(params.ref_gff))
    if( asBool(params.rna_support) ) {
      RNA_JUNCTIONS(HISAT2_ALIGN.out.bam.collect())
      rna_targets = ANNOTATE_TIERS.out.gff.map{ [it, 'union.final.tiers'] }
      if( params.ref_gff )
        rna_targets = rna_targets.mix( MERGE_REFERENCE.out[0].map{ [it, 'reference_updated'] } )
      // RNA_JUNCTIONS is fed a collected (value) channel, so its output is a
      // value channel too and both targets can read it. No .first() needed --
      // nextflow warns that it is useless here, and it is.
      RNA_SUPPORT(rna_targets, RNA_JUNCTIONS.out.tsv)
    }
  } else {
    // Without the QC stage the union is the final product, so provenance is
    // written against that instead of being skipped: a run with --run_qc off
    // still needs to say which versions produced it. There is no BUSCO result
    // to report, hence the placeholder.
    PROVENANCE(BUILD_UNION.out.gff, aug_sp, no_busco)
    if( params.ref_gff )
      MERGE_REFERENCE(BUILD_UNION.out.gff, BUILD_UNION.out.gff, file(params.ref_gff))
  }
}
