# Re-annotation analysis scripts

Small, dependency-light scripts that reproduce the supplementary tables and the
per-locus figure comparing a Gene-Miner catalogue against a reference
annotation. They all read the pipeline's own outputs (`union.gff3` before QC,
`union.final.gff3`, `genome.transdecoder.gff3`, `union.emapper.annotations`) plus
the coordinate-matched reference GFF3, and use only `python3` (with `pysam` /
`numpy` / `matplotlib` for the figure). Matching everywhere is same-strand
CDS-footprint overlap, reciprocal ≥ 0.5 of the shorter locus — the same rule as
`bin/gm_compare.py`.

Supplementary numbers below are the manuscript's rendered numbers (Supplementary
Table 5–8, Supplementary Figure 1); update them if the supplement is renumbered.

| Script | Produces | Reads |
|---|---|---|
| `cds_to_exon_gtf.py` | inputs for the coding-level gffcompare of **STable 6** (recall/precision/F1) | a GFF3 |
| `classify_missed.py` | **STable 5** — why reference genes were not recovered | raw union.gff3, union.final.gff3, reference |
| `revised_junctions.py` | **STable 7** — splice-junction support of revised models | union.final.gff3, reference, junctions.tsv |
| `added_intron_streams.py` | **STable 7** decomposition — GM-added introns split RNA-derived vs *ab initio*, with junction-confirmation % | union.final.gff3, reference, junctions.tsv, transdecoder.gff3, genemark.gtf, augustus.gff3 |
| `novel_tiers.py` | **STable 8** novel-locus confidence tiers + the per-locus `STable_novel_confidence_tiers.tsv` | union.final.gff3, reference, transdecoder.gff3, eggNOG, blastp |
| `pick_example.py` | picks an example revised locus for **SFigure 1a** (GM better-supported) | as `novel_tiers.py` inputs |
| `pick_worse_example.py` | picks an example revised locus for **SFigure 1b** (GM worse: added *ab initio* introns unconfirmed) | as `added_intron_streams.py` inputs |
| `../bin/plot_splice_support.py` | one panel of **SFigure 1** — spliced-read view of one locus | GFF3s, junctions.tsv, BAMs |
| `compose_sfig.py` | stacks the two panels into **SFigure 1** (a)/(b) | two panel PNGs |

## STable 6 — structural recall/precision/F1 (gffcompare, coding level)

Restrict both annotations to protein-coding CDS (drops ncRNA, ignores UTRs) then
run gffcompare, which reports sensitivity (recall) and precision at the
nucleotide/exon/transcript/locus levels; F1 is their harmonic mean.

```bash
python3 cds_to_exon_gtf.py reference.gff3       ref.cds.gtf
python3 cds_to_exon_gtf.py union.final.gff3     gm.cds.gtf
gffcompare -r ref.cds.gtf -o cmp gm.cds.gtf     # read Sn|Sp from cmp.stats
```

## STable 5 — reasons reference genes were not recovered

Classifies each not-recovered reference gene as prediction-failure / structural-
mismatch / wrong-strand / QC-removed (the last needs the pre-QC `union.gff3`).

```bash
python3 classify_missed.py union.gff3 union.final.gff3 reference.gff3
```

## STable 7 + SFigure 1 — splice-junction support of revised models

First extract splice junctions from the pipeline's BAMs and merge to a
`chrom<TAB>intron_start<TAB>intron_end<TAB>read_count` TSV (1-based introns):

```bash
for b in sorted_*.bam; do regtools junctions extract -s XS -o $b.junc.bed $b; done
# convert each BED12 row: istart = chromStart+blockSizes[0]+1, iend = chromEnd-blockSizes[1],
# summing BED col 5 (read count) across libraries -> junctions.tsv
```

Intron-level shared / added / dropped confirmation (STable 7):

```bash
python3 revised_junctions.py union.final.gff3 reference.gff3 junctions.tsv 3
```

Split the **GM-added** introns by originating stream and report each class's
junction-confirmation rate — RNA-seq-derived (in a TransDecoder transcript) is
self-confirming, whereas *ab initio*-predicted (AUGUSTUS / GeneMark-ETP, not in
any assembled transcript) is mostly unconfirmed. This is what makes a locus score
"GM worse" per-locus: added but expression-unconfirmed *ab initio* structure, not
lost support (STable 7 caption; the rice split is 38.8% RNA-derived @ 97.2%
confirmed vs 61.2% *ab initio* @ 7.8%):

```bash
python3 added_intron_streams.py union.final.gff3 reference.gff3 junctions.tsv \
        genome.transdecoder.gff3 genemark.gtf augustus_raw.gff3 3
```

### SFigure 1 — per-locus spliced-read view, (a) confirmed and (b) unconfirmed added structure

`pick_example.py` finds a clean "GM better" locus (added exons bridged by spliced
reads) for panel (a); `pick_worse_example.py` finds a "GM worse" locus (GM keeps
the reference's confirmed introns but appends *ab initio* introns with no spliced
reads) for panel (b). Render each with `../bin/plot_splice_support.py`, then stack
them with `compose_sfig.py`:

```bash
python3 pick_worse_example.py union.final.gff3 reference.gff3 junctions.tsv \
        genome.transdecoder.gff3 genemark.gtf augustus_raw.gff3 3     # candidate loci for panel (b)

# panel (a): GM better (e.g. rice Os_A001065 vs gene:Os10g0117000)
python3 ../bin/plot_splice_support.py --chrom 10 --start 1093271 --end 1096798 --strand - \
        --gm-gene Os_A001065 --ref-gene gene:Os10g0117000 \
        --gm-gff union.final.gff3 --ref-gff reference.gff3 --ref-label "RAP-DB reference" \
        --junctions junctions.tsv --out panelA.png \
        --bam sorted_1.bam --bam sorted_2.bam --bam sorted_3.bam --bam sorted_4.bam

# panel (b): GM worse (e.g. rice Os_A012705 vs gene:Os02g0686700)
python3 ../bin/plot_splice_support.py --chrom 2 --start 28132995 --end 28134455 --strand - \
        --gm-gene Os_A012705 --ref-gene gene:Os02g0686700 \
        --gm-gff union.final.gff3 --ref-gff reference.gff3 --ref-label "RAP-DB reference" \
        --junctions junctions.tsv --out panelB.png \
        --bam sorted_1.bam --bam sorted_2.bam --bam sorted_3.bam --bam sorted_4.bam

python3 compose_sfig.py panelA.png panelB.png sfigure_1.png       # needs Pillow
```

## STable 8 — novel-locus confidence tiers

Scores each novel locus for three independent evidence types (RNA-seq expression
= CDS overlaps a TransDecoder transcript; Swiss-Prot `blastp` hit E<1e-5; eggNOG
ortholog) and bins into High (≥2), Medium (1), Low (0). Prints the tier counts;
the trailing two args write a per-locus TSV.

```bash
# blastp hits: diamond blastp -q union.final.pep.fa -d sprot_db.dmnd -e 1e-5 -k1 --outfmt 6 -o sprot_hits.tsv
python3 novel_tiers.py union.final.gff3 reference.gff3 genome.transdecoder.gff3 \
        union.emapper.annotations sprot_hits.tsv  Rice rice_tiers.tsv
```
Note: the eggNOG file must cover **all** streams. If the GeneMark (`_E`) stream was
run through a separate eggNOG job, concatenate the base-union and `_E` annotation
files before passing them in.

## Main-text figures 2 and 4

Standalone matplotlib scripts that redraw the two data figures from their
(hard-coded, from the tables) values — kept so the layouts can be regenerated.

```bash
python3 plot_figure2_venn.py          # Figure 2 (revised-model Venn, rice/soybean)
python3 plot_figure4_completeness.py  # Figure 4 (completeness-vs-evidence + isoforms)
```
`plot_figure2_venn.py` needs `matplotlib-venn`.
