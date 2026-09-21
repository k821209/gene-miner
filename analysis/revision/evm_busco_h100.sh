export PATH=/data/k821209/gm_evm/busco_env/bin:$PATH
busco --version
for n in evm.te evm.ren gm.primary; do
  busco -i $n.pep.fa -l poales_odb10 -m proteins -o busco_$n -f --offline --download_path /data/k821209/gm_evm/busco_dl -c 16 > busco_$n.log 2>&1
  echo "$n $(grep -h "C:" busco_$n/short_summary*.txt)"
done
echo BUSCO_DONE
