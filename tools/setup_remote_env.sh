#!/usr/bin/env bash
# task2 (单细胞 + 基础模型) 的 GPU 环境，在远端 devbox 上构建。
# 远端：Ubuntu24 / 2×RTX3090 / 驱动580 / CUDA12+ / uv（无 conda）。
# 用法：在 ~/mycode/sasp-ibd 下  bash tools/setup_remote_env.sh
set -uo pipefail
cd "$(dirname "$0")/.."

uv venv --python 3.11 .venv
source .venv/bin/activate

# torch：PyPI 默认即 CUDA 构建；若拿到 CPU 版则回退 cu121 索引
uv pip install torch
python -c "import torch,sys;sys.exit(0 if torch.cuda.is_available() else 3)" \
  || uv pip install torch --index-url https://download.pytorch.org/whl/cu121
python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available(),'ngpu',torch.cuda.device_count())"

# 单细胞栈
uv pip install scanpy anndata celltypist huggingface_hub
uv pip install scrublet || echo "scrublet 跳过(scanpy 内置可替代)"

# Geneformer：免 git-lfs，用 hf 只拉包+最小权重(V1-10M)，跳过 V2 大模型
python - <<'PY'
from huggingface_hub import snapshot_download
snapshot_download("ctheodoris/Geneformer", local_dir="models/Geneformer",
    ignore_patterns=["Geneformer-V2-*","model.safetensors","fine_tuned_models/*","*.ipynb"])
PY
uv pip install -e models/Geneformer
# 关键：Geneformer 要 transformers==4.46；不钉死会被装成 v5 而报 SpecialTokensMixin 导入错
uv pip install "transformers==4.46.3"
python -c "import geneformer;print('geneformer OK')"

echo "[ENV READY] 冒烟测试： python tools/smoke_test_task2.py"
# 注：cellphonedb(第4步细胞通讯)依赖常与 scanpy 冲突，单独建 venv，勿装进此环境。
