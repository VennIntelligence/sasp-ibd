#!/usr/bin/env bash
# 从东京 devbox 把结果(或指定子路径)拉回本地 Mac，用于查看/挑选/晋升进 results/ 再提交。
# 不套 .gitignore 过滤 —— 这样你可以按需抓被 ignore 的 outputs/ 里的某张图/某张表。
# 只挡缓存类垃圾(__pycache__、.DS_Store)。
#
# 用法:
#   tools/sync_pull.sh                    # 默认拉 results/
#   tools/sync_pull.sh outputs/14_mr/     # 拉某个具体 output 子目录
#   tools/sync_pull.sh outputs/Fig_MAIN_senescence.png   # 拉单个文件
set -euo pipefail
cd "$(dirname "$0")/.."
SUB="${1:-results/}"
rsync -avz --human-readable --mkpath \
  --exclude='__pycache__/' --exclude='*.pyc' \
  --exclude='.DS_Store' --exclude='.ipynb_checkpoints/' \
  "tokyo:~/lzy/${SUB}" "./${SUB}"
