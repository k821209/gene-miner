# Revision analyses (round 1, BMC Bioinformatics)

Scripts behind the analyses added during the first revision. They read the
pipeline outputs of a finished run and write the tables of the manuscript;
paths at the top of each script point at the run directories used here and are
the only thing to change for another dataset. Supplementary table numbers below
are those of the revised manuscript.

| Script | Produces |
|---|---|
| `regen.sh`, `rice_rebuild.sh` | rebuilt catalogues for the four genomes (Tables 1-3 and most supplementary tables) |
| `sensitivity.sh` | AUGUSTUS posterior-probability sweep (STable 1) |
| `split_fusion.py`, `fusion_origin.py`, `structural_benchmark.sh` | gene split/fusion counts and their source stream (STable 2) |
| `te_threshold_audit.py` | TE-overlap threshold sweep (STable 4) |
| `extract_classes.py`, `kaks_pipeline.py`, `kaks_by_tier.py` | dN/dS against wild Oryza by category and tier (STable 5) |
| `heldout_junctions.sh`, `heldout_test.py` | held-out-library intron check (STable 6) |
| `subsample_rice.sh`, `busco_primary.sh`, `collect_sub.py`, `nextflow.226.config` | RNA-seq subsampling series (STable 8) |
| `runtime_run.sh`, `collect_runtime.py` | per-stage runtime, CPU-hours, peak memory (STable 10) |
| `run_evm_rice.sh`, `run_evm_h100.sh`, `evm_te_filter_h100.sh`, `evm_busco_h100.sh` | EVidenceModeler comparison (STable 14) |
| `revised_features.py`, `isoform_origin.py`, `novel_by_stream.py` | revision types and the origin of added models (STables 16, 21) |
| `contamination.sh`, `organellar.sh` | non-target taxa and NUMT/NUPT screen (STable 22) |
| `heldout_novel_introns.py` | held-out confirmation of the introns the reference lacks (Results) |
| `wild_gene_overlap.py`, `wild_overlap_run.sh` | whether the sampled loci fall on genes the wild relatives' own annotations record (Results) |
| `codon_usage.py` | codon usage by comparison class and tier against an intergenic-ORF null (Results) |
| `paired_codon.py`, `paired_spread.py` | check that the revised class' weaker codon bias belongs to the reference genes it revises (not reported in the paper) |
| `transfer_eggnog.py` | reuse eggNOG rows for unchanged proteins when a catalogue is rebuilt |

`te_threshold_audit.py`, `extract_classes.py`, `heldout_test.py` and
`paired_codon.py` import `gm_compare.py` from `bin/` and `gm_overlap.py` from
`analysis/`; `paired_codon.py` and `paired_spread.py` also import
`codon_usage.py` from this directory.
