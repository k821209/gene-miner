#!/usr/bin/env bash
# After run_evm_h100.sh: rename EVM IDs to the GM-style pattern te_filter_genes.py
# recognises (<prefix>_[AER]<digits>), apply the same TE filter as Gene-Miner, re-score.
# Run in /data/k821209/gm_evm/run on h100. The python renaming block is the one used
# on 2026-09-18 (gene EVM_E%06d, mRNA <gene>.1).
set -uo pipefail
B=/data/k821209/gm_evm; REF=$B/rice/genome/IRGSP.ref.gff3
export PATH=$B/env/bin:$PATH
python3 $B/gene-miner/bin/te_filter_genes.py $B/rm.out evm.ren.gff3 evm.ren.pep.fa evm.te.gff3 evm.te.pep.fa > te.log 2>&1
PYTHONPATH=$B/gene-miner/bin python3 $B/gene-miner/bin/gm_compare.py --gm evm.te.gff3 --ref $REF --min_ro 0.5 --out gm_compare.te.tsv
python3 $B/revision/split_fusion.py $REF evm.te.gff3 0.5 > split_fusion.te.txt
python3 $B/gene-miner/analysis/cds_to_exon_gtf.py evm.te.gff3 evm.te.cds.gtf
gffcompare -r $B/revision/structural/rice.ref.cds.gtf -o cmpte evm.te.cds.gtf
