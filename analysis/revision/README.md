# Revision analyses (round 1, BMC Bioinformatics)

Scripts behind the analyses added during the first revision. They read the
pipeline outputs of a finished run and write the tables of the manuscript;
paths at the top of each script point at the run directories used here and are
the only thing to change for another dataset.

| Script | Produces |
|---|---|
| `regen.sh`, `rice_rebuild.sh` | rebuilt catalogues for the four genomes (Tables 1, 2, STables 1, 5-13) |
| `sensitivity.sh` | AUGUSTUS score sweep (STable 14) |
| `split_fusion.py`, `fusion_origin.py` | gene split/fusion counts and their source stream (STable 15) |
| `te_threshold_audit.py` | TE-overlap threshold sweep (STable 16) |
| `runtime_run.sh`, `collect_runtime.py` | per-stage runtime, CPU-hours, peak memory (STable 17) |
| `revised_features.py`, `isoform_origin.py`, `novel_by_stream.py` | revision types and the origin of added models (STables 12, 18) |
| `extract_classes.py`, `kaks_pipeline.py`, `kaks_by_tier.py` | dN/dS against wild Oryza by category and tier (STable 19) |
| `contamination.sh`, `organellar.sh` | non-target taxa and NUMT/NUPT screen (STable 20) |
| `run_evm_rice.sh`, `run_evm_h100.sh`, `evm_te_filter_h100.sh`, `evm_busco_h100.sh` | EVidenceModeler comparison (STable 21) |
| `heldout_junctions.sh`, `heldout_test.py` | held-out-library intron check (STable 22) |
| `subsample_rice.sh`, `busco_primary.sh`, `collect_sub.py` | RNA-seq subsampling series (STable 23) |
| `transfer_eggnog.py` | reuse eggNOG rows for unchanged proteins when a catalogue is rebuilt |

`te_threshold_audit.py`, `extract_classes.py` and `heldout_test.py` import
`gm_compare.py` from `bin/` and `gm_overlap.py` from `analysis/`.
