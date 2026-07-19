# Reproduce the rice catalogue (one command, end to end)

This reproduces the paper's rice result with a single `nextflow run`. It is a
full pipeline run (masking → RNA-seq → AUGUSTUS → GeneMark-ETP → union → QC →
BUSCO), not a quick unit test — budget ~6 h on a 24-core node (the fast
per-script checks are `python3 -m unittest discover -s test`).

Validated on deevo-219 (2026-07-18, `nextflow` exit 0, ~6 h): **41,513 genes /
78,373 transcripts**, all-isoform BUSCO `poales_odb10` **96.5 %** — matching the
paper's step-wise catalogue (41,509 / 78,328; the ±4 genes are RepeatMasker
re-run boundary noise).

## Public inputs

| Input | Source |
|---|---|
| Genome | *Oryza sativa* Nipponbare **IRGSP-1.0** (Ensembl Plants / RAP-DB) |
| RNA-seq (4 tissues) | `SRR26136061` flag leaf, `SRR26136051` panicle, `SRR26136064` seedling, `SRR26136038` endosperm (BioProject PRJNA1019705) |
| Proteins | UniProtKB/Swiss-Prot (`uniprot_sprot.fasta`) |

```bash
# 1. genome
#    download IRGSP-1.0 (e.g. from Ensembl Plants) -> genome/IRGSP.fa
# 2. reads: fetch the four libraries and name them <tissue>_{1,2}.fastq.gz
for a in SRR26136061 SRR26136051 SRR26136064 SRR26136038; do
  prefetch $a && fasterq-dump --split-files -O rnaseq $a       # sra-tools
done
# rename to leaf_/panicle_/seedling_/endosperm_{1,2}.fastq.gz (any names work)
# 3. proteins: download uniprot_sprot.fasta -> db/uniprot_sprot.fasta
```

## Run

First point the `env_*` prefixes in `nextflow.config` at **your** conda envs
(annot, augustus, rmod, eggnog, busco) and `eggnog_db` at your eggNOG DB — the
committed values are the deevo-219 host paths.

```bash
nextflow run main.nf -c nextflow.config -profile rice \
  --genome    genome/IRGSP.fa \
  --reads     'rnaseq/*_{1,2}.fastq.gz' \
  --proteome  db/uniprot_sprot.fasta \
  --repeat_lib repeats/rice.repeatlib.fa \   # optional: skip de-novo RepeatModeler
  --outdir    gm_out
```

`-profile rice` sets `--augustus_species rice`, `--busco_lineage poales_odb10`,
`--gene_prefix RICE`, `--run_genemark true`, `--run_qc true` (see
`nextflow.config`). Omit `--repeat_lib` to build the repeat library de novo (adds
hours). Omit `--augustus_species` (or pass `auto`) on a clade with no native
AUGUSTUS model to BUSCO-train one instead.

## Expected output

`gm_out/union.final.gff3` + `gm_out/union.final.pep.fa`, ~41.5 k genes; the
`gm_out/busco/` summary reports ~96.5 % all-isoform completeness on
`poales_odb10` (~93.6 % on representative transcripts). Small deviations from the
exact counts are expected (RepeatMasker/RNA-seq assembly are not bit-deterministic).
