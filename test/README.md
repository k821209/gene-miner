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

These exercise the decisions that define the final catalogue — the union and
the two QC filters — not the heavy external predictors, which are integration
concerns covered by a full pipeline run (see the top-level `README.md`).
