#!/usr/bin/env bash
# 把工作区(代码/脚本/results)推到东京 devbox 的 ~/lzy 镜像，用于"提交前先在远端试跑"。
# 用 .gitignore 过滤 → 永远不会把 data/、outputs/、.venv/ 推过去（不带垃圾过去）。
#
# 用法:
#   tools/sync_push.sh            # 增量推送（不删远端多余文件，安全）
#   tools/sync_push.sh -n         # dry-run，只看会传什么
#   tools/sync_push.sh --delete   # 镜像式：删掉远端多余的(已被 ignore 的 data/、outputs/ 受保护，不会删)
# 远端仓库位置：~/mycode/sasp-ibd（单细胞数据在其 data/scrna/，被 ignore，不会被同步覆盖）
#
# 注意：这是"快速迭代"通道，会让远端 git 树变 dirty。正式定稿请走 git commit/push + 远端 git pull。
set -euo pipefail
cd "$(dirname "$0")/.."
rsync -avz --human-readable \
  --filter=':- .gitignore' \
  --exclude='.git/' \
  "$@" \
  ./ tokyo:~/mycode/sasp-ibd/
