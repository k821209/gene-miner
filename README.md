# gene-miner

Mine a **clean protein-coding gene catalog** from a non-model genome assembly —
even a messy, whole-body, gut-microbiome-laden, modest-coverage one. The
pipeline takes a genome (+ RNA-seq + a protein DB), maximises true gene recovery
by **unioning** independent evidence streams, then **QC-filters out the junk**
(residual bacterial genes, TE-derived genes) that a raw annotation silently
keeps.

```
RNA-seq  ─HISAT2─StringTie─TransDecoder─┐
                                        ├─ UNION ─→ raw gene set ─QC filters─→ CLEAN catalog
masked   ─AUGUSTUS (BUSCO-trained)──────┘            (over-counts)   │           (union.final.*)
contigs                                                              ├ eggNOG taxonomy  → −bacterial
                                                                     └ RepeatMasker     → −TE
```

The scripts are simple; the **Lessons** below are the point of the repo — they
are the gotchas that cost real time on real data. Worked example: a Korean
*Zeugodacus scutellatus* (pumpkin-flower fruit fly; 호박꽃과실파리) assembly,
mined to **11,463** genes (raw union 14,311 − 2,010 bacterial − 838 TE ≈ 86% of
a published congener's 13,327, matching its BUSCO-protein 88%).

---

## Why "mining" and not just "annotation"

A single annotation tool gives a biased answer, and a raw merge over-counts.
This pipeline is built around two facts learned the hard way:

1. **No single evidence stream is complete.** AUGUSTUS alone over-predicts
   (23,505 here); EVidenceModeler consensus alone *under*-counts (5,919–7,282,
   because it drops ab-initio-only genes). RNA-seq contributed 6,350 loci absent
   from high-confidence AUGUSTUS (1,439 predicted by AUGUSTUS nowhere). The right
   answer is the **union** of high-confidence AUGUSTUS and RNA-only loci.
2. **A union still contains non-target genes.** Whole-body insect DNA is largely
   gut bacteria; contig-level decontamination removes bacterial *contigs* but
   low-GC ones leak through as gene models. Reference pipelines also exclude
   TE-derived genes. So the union gets a **functional QC pass** before it is a
   real catalog: eggNOG taxonomy drops bacterial genes, RepeatMasker overlap
   drops TE genes. Here that was 15% of the raw union.

---

## Pipeline stages

| # | Stage | Script | Output |
|---|---|---|---|
| 1 | RNA-seq evidence (HISAT2 → StringTie → TransDecoder + DIAMOND/Pfam) | `bin/run_rnaseq_transdecoder.sh` | `genome.transdecoder.gff3` |
| 2 | AUGUSTUS ab initio (BUSCO-trained model, per-contig parallel) | `bin/run_augustus.sh`, `bin/lift_agp.py` | `augustus_scaffold.gff3` |
| 3 | **UNION** of usable AUGUSTUS (score≥0.8, ≥100 aa) ∪ RNA-only loci | `bin/build_union.py` | `union.gff3` |
| 4 | Translate to proteins (dependency-free fallback for `gffread -y` segfaults) | `bin/extract_pep.py` | `union.pep.fa` |
| 5 | **QC-a** eggNOG taxonomy → drop residual bacterial genes | `bin/run_eggnog.sh`, `bin/filter_taxonomy.py` | `union.clean.*` |
| 6 | **QC-b** RepeatMasker overlap → drop TE-derived genes | `bin/run_repeatmasker.sh`, `bin/te_filter_genes.py` | `union.final.*` |
| 7 | BUSCO-protein validation (run where the lineage DB lives) | — | completeness check |

Optional helpers: `bin/run_miniprot.sh` (protein-spliced alignment evidence),
`bin/run_evm3.sh` (EVidenceModeler consensus — kept for comparison; **not** the
final set, see Lesson 1), `bin/filter_models.py` (tiered high-confidence subset
of the RNA-seq models), `bin/union_genes.py` (standalone union utility).

---

## Setup

Create the conda environments the pipeline uses, in one command (needs conda or
mamba on PATH):

```bash
bash setup_envs.sh     # builds: annot, augustus, rmod, eggnog, busco
```

This covers the default **two-stream** run (AUGUSTUS + RNA-seq + QC). The
optional **3rd stream, GeneMark-ETP** (`--run_genemark true`) is *not*
conda-installable on a clean machine — the bioconda `braker3` recipe requires
`genomethreader`, which bioconda no longer ships. For it, use the official
BRAKER container (`docker://teambraker/braker3`, which bundles GeneMark-ETP) and
point `run_genemark_etp.sh` at it, or run GeneMark-ETP standalone.

Databases (eggNOG ~50 GB, BUSCO lineage, optional Pfam) are fetched on first use
— see the notes `setup_envs.sh` prints. The pipeline finds the envs under
`$HOME/miniconda3/envs` by default; point elsewhere with
`export GM_CONDA_BASE=/path/to/miniconda`.

## Run

**One command, end to end (Nextflow).** `main.nf` runs every stage — masking,
RNA-seq assembly, ab-initio prediction, GeneMark-ETP, the union, the two QC
filters and BUSCO — and writes the clean catalogue to `<outdir>/union.final.*`:

```bash
nextflow run main.nf -c nextflow.config \
  --genome    genome.fa \
  --reads     'rnaseq/*_{1,2}.fastq.gz' \
  --proteome  uniprot_sprot.fasta \
  --augustus_species rice \      # native AUGUSTUS model; omit / 'auto' => BUSCO-train one
  --repeat_lib repeats.fa \      # supplied library;       omit          => de-novo RepeatModeler
  --run_genemark true \          # add GeneMark-ETP as a 3rd ab-initio stream
  --busco_lineage poales_odb10 \
  --outdir gm_out
```

The workflow's branch points are parameters (all defaults in `nextflow.config`):

| Choice | Parameter | Options |
|---|---|---|
| AUGUSTUS model | `--augustus_species` | a native model name (e.g. `rice`) — or omit / `auto` to **BUSCO-train** one for a new clade |
| Repeat library | `--repeat_lib` | a supplied `.fa` (fast) — or omit to build one **de novo** with RepeatModeler (slow) |
| GeneMark-ETP 3rd stream | `--run_genemark` | `true` / `false` (or pass a precomputed `--genemark_gtf`) |
| QC filters | `--run_qc` | `true` / `false` (stop at the raw union if false) |

**Alternative — the shell driver** (`run_gene_mining.sh`) is the original,
battle-tested implementation aimed at *messy, contaminated* assemblies (per-contig
AUGUSTUS + RagTag AGP lift; see Lesson 3). Edit its CONFIG block, then:

```bash
bash run_gene_mining.sh          # stages 1-7 -> union.final.gff3 / union.final.pep.fa
# RUN_QC=0 stops at the raw union; REPEAT_LIB=<lib.fa> skips de-novo RepeatModeler
```

### Requirements (conda envs)
- core: `HISAT2, StringTie, TransDecoder, DIAMOND, HMMER(hmmscan), AUGUSTUS (+BUSCO-retrained model), samtools, gffread, miniprot, EvidenceModeler`
- QC:   `eggNOG-mapper (+DB), RepeatMasker/RepeatModeler, BUSCO (+lineage DB)`

The QC runners (`run_eggnog.sh`, `run_repeatmasker.sh`) create their conda env /
fetch their DB on first run.

### Repeat library: de-novo vs supplied

The TE QC step (stage 6) is driven by RepeatMasker. The `REPEAT_LIB` config
variable controls where the repeat library comes from:

- **`REPEAT_LIB=` (empty, default)** — build one *de novo* with RepeatModeler
  (slow, hours). `bin/run_repeatmasker.sh` does this automatically when no
  library is passed.
- **`REPEAT_LIB=<lib.fa>`** — reuse an existing/curated library and skip
  RepeatModeler (fast; also the right choice when comparing TE content
  *between* genomes, which needs the same `-lib` on both).

TE-derived genes are then flagged purely by **RepeatMasker CDS overlap**
(`TE_THRESH`, default 0.5 of the CDS in interspersed repeats) — there is no
Pfam-domain-based TE branch. Pfam (`hmmscan`) is used only in stage 1 to
*retain* coding ORFs (`TransDecoder --retain_pfam_hits`), not to identify
repeats.

## Tests

Lightweight smoke tests exercise the union and QC-filter logic on tiny
synthetic fixtures — **`python3` only**, no external predictors or test-data
download, under a second to run:

```bash
python3 -m unittest discover -s test -v     # or: pytest test/
```

See `test/README.md` for what each test covers.

---

## Lessons (read before reusing)

1. **The final gene set is a UNION, not an evidence-only consensus.** EVM
   (transcript+protein consensus) dropped ab-initio-only genes → under-count;
   AUGUSTUS alone over-predicts. Keep **usable AUGUSTUS** (gene score ≥ 0.8,
   ≥100 aa) and **union** it with RNA-seq/TransDecoder loci AUGUSTUS missed. Each
   evidence stream finds genes the other misses. See `bin/build_union.py`. **But
   that raw count is _provisional_ — clean it (Lesson 4) before comparing.**

2. **Reuse the BUSCO-retrained AUGUSTUS model.** BUSCO genome mode with
   `--augustus`/`--long` retrains AUGUSTUS on the found single-copy orthologs and
   leaves a species model in `run_*/augustus_output/retraining_parameters/`. Copy
   it into `AUGUSTUS_CONFIG_PATH/species/<name>/`. Old (3.2.2) models work with
   AUGUSTUS 3.5.0. No retraining = bad ab-initio predictions on a new clade.

3. **Keep coordinate systems consistent.** RNA-seq/TransDecoder ran on the RagTag
   scaffold; AUGUSTUS ran on masked contigs. Lift AUGUSTUS contig→scaffold via the
   RagTag AGP (`bin/lift_agp.py`) **before** unioning — otherwise the union
   double-counts the same locus under two coordinate systems.

4. **A union gene set needs a functional QC pass — eggNOG caught 15% junk.**
   eggNOG-mapper of the union exposed two non-target classes the union had kept:
   **2,010 residual bacterial genes** (gut microbiome surviving *contig*
   decontam — `groL`/`tuf`/ribosomal; a clean reference had 0) and **838
   TE-derived genes** (CDS ≥50% in interspersed repeats; reference pipelines mask
   these out). Dropping both: 14,311 → 12,301 → **11,463** (≈86% of the
   reference's 13,327, now consistent with BUSCO-protein 88% — the raw
   14,311≈13,327 was coincidental). Filters: `bin/filter_taxonomy.py`,
   `bin/te_filter_genes.py`.

5. **Contamination inflates the gene count upstream too.** Gene-dense bacterial
   *contigs* (35.8 Mb here) once inflated an earlier annotation to 42,728 genes.
   Decontaminate the assembly (drop non-Arthropoda contigs) and validate with
   BUSCO *before* annotating — removing 188 bacterial contigs left BUSCO
   unchanged (97.9 vs 98.1%) with *lower* duplication (1.1→0.7%), proving only
   non-fly sequence went. The eggNOG taxonomy filter (Lesson 4) is the
   gene-level mop-up for what survives contig-level cleaning.

6. **Comparing gene/TE content between genomes? Use the SAME library/DB on both.**
   An apparent "lineage A has expanded family X" was pure annotation-filter
   confound. Run the same RepeatModeler library and same eggNOG DB on **both**
   genomes; here TE content was 29.1% vs 32.9% (near-identical → no expansion).
   `run_repeatmasker.sh` takes a shared `-lib` for exactly this.

7. **Infrastructure gotchas.** `gffread -y` can segfault on some valid gff3 —
   `bin/extract_pep.py` is a dependency-free fallback translator. `stringtie` has
   no `--dta` (that's a HISAT2 flag). `download_eggnog_data.py` may silently fail
   (0-byte DB files) — `wget -c` the DB from `eggnog6.embl.de` directly
   (`run_eggnog.sh` does this). Conda `RepeatMasker` needs its env **activated**
   (`REPEATMASKER_DIR`), not just the binary on `$PATH`. Long remote jobs need
   `setsid nohup … </dev/null &` (a plain background SSH dies on SIGHUP if the
   connection drops).
