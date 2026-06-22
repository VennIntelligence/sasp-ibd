#!/usr/bin/env bash
set -u
cd /Users/ujs/Downloads/lzy
export PYTHONPATH=src

mkdir -p data/raw/gwas outputs/mr

url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002398/harmonised/32888494-GCST90002398-EFO_0004833.h.tsv.gz"
out="data/raw/gwas/LTL_Codd2021_GCST90002398.h.tsv.gz"
log="outputs/mr/ltl_download_run.log"

{
  echo "[$(date)] LTL download started"
  echo "URL: $url"
  echo "OUT: $out"
} | tee -a "$log"

curl -L --fail -C - --progress-bar -o "$out" "$url"
rc=$?
echo "[$(date)] curl exit code: $rc" | tee -a "$log"

if [ "$rc" -eq 0 ]; then
  echo "[$(date)] gzip test" | tee -a "$log"
  gzip -t "$out"
  gzrc=$?
  echo "[$(date)] gzip exit code: $gzrc" | tee -a "$log"
  if [ "$gzrc" -eq 0 ]; then
    echo "[$(date)] running src/20_ltl_mr.py" | tee -a "$log"
    .venv/bin/python src/20_ltl_mr.py 2>&1 | tee -a "$log"
  fi
fi

echo "[$(date)] LTL window finished." | tee -a "$log"
