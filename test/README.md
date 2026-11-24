# Tests

Lightweight smoke tests for the pure-Python catalogue-building logic. They run
the real `bin/` scripts on tiny synthetic fixtures (a few genes on one
contig), so they need **only `python3`** — no AUGUSTUS, GeneMark-ETP, eggNOG,
RepeatMasker or BUSCO, and no test data download. The whole suite finishes in
under a second.

## Run

```bash
python3 -m unittest discover -s test -v        # stdlib only
# or, if you use pytest:
pytest test/
```

## What is covered

| Test file | Script under test | Checks |
|---|---|---|
| `test_build_union.py` | `bin/build_union.py` | usable-AUGUSTUS score/length gates; RNA-only and GeneMark-only loci become their own `_R` / `_E` genes; overlapping RNA/GeneMark transcripts are folded onto the AUGUSTUS locus as extra isoforms; GeneMark `status != complete` transcripts are skipped |
| `test_te_filter.py` | `bin/te_filter_genes.py` | a gene whose CDS lies inside an interspersed repeat is dropped; a clean gene and a gene overlapping only a *simple* repeat are kept; the CDS-fraction threshold is respected |
| `test_filter_taxonomy.py` | `bin/filter_taxonomy.py` | a gene with only bacterial-rooted eggNOG orthogroups is dropped; eukaryotic and unassigned genes are kept |
| `test_merge_placement.py` | `bin/merge_into_reference.py` | which reference gene an added transcript goes under (its own CDS overlap, not the locus's pairing); a read-through model gets its own record instead of being filed under one of the two genes it spans; a widening that would reach into a neighbouring reference gene is refused and the transcript written separately; a widening into empty sequence still happens; no reference gene ends up holding coding bases of another |
| `test_merge_placement.py` (second half) | `bin/check_merge_placement.py` | the independent auditor: it passes a correct merge, fails a hand-built file with a transcript under the wrong gene, and returns **COULD NOT ASSESS (exit 2)** rather than "clean" for a file with no additions or a wrong `--source` |
| `test_lift_agp.py` | `bin/lift_agp.py` | contig→scaffold coordinate lift used by the shell driver: forward contig unchanged, reverse-oriented contig mirrored with a strand flip, contig absent from the AGP dropped, header preserved |

These exercise the decisions that define the final catalogue — the union and
the two QC filters — not the heavy external predictors.

## Full pipeline reproduction (rice)

For an end-to-end check on real data — a one-command Nextflow run that
reproduces the paper's rice catalogue (~41.5 k genes, ~96.5 % all-isoform BUSCO,
~6 h) — see [`reproduce_rice.md`](reproduce_rice.md). It lists the public inputs
(IRGSP-1.0 genome + four SRA libraries + Swiss-Prot) and the exact `-profile rice`
invocation.
