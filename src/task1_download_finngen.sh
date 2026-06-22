#!/usr/bin/env bash
set -u
cd /Users/ujs/Downloads/lzy
export PYTHONPATH=src

mkdir -p data/raw/gwas outputs/mr
log="outputs/mr/finngen_download_run.log"

echo "[$(date)] FinnGen R12 download + replication started" | tee -a "$log"
.venv/bin/python src/21_replicate_finngen.py --download 2>&1 | tee -a "$log"
rc=${PIPESTATUS[0]}
echo "[$(date)] FinnGen script exit code: $rc" | tee -a "$log"
echo "[$(date)] FinnGen window finished." | tee -a "$log"
