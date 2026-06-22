import time, torch, pandas as pd
print("="*60)
print("[1] torch / GPU")
print("  torch", torch.__version__, "| cuda", torch.cuda.is_available(),
      "| ngpu", torch.cuda.device_count())
for i in range(torch.cuda.device_count()):
    print("   ", i, torch.cuda.get_device_name(i),
          f"{torch.cuda.get_device_properties(i).total_memory/1e9:.0f}GB")
x = torch.randn(4096, 4096, device="cuda")
t = time.time(); y = (x @ x).sum().item(); dt = time.time()-t
print(f"  cuda matmul 4096^2 ok ({dt*1000:.0f} ms), checksum finite: {y==y}")

print("[2] scanpy 读下载的数据")
import scanpy as sc, anndata
print("  scanpy", sc.__version__, "anndata", anndata.__version__)
df = pd.read_csv("data/scrna/GSE116222_Expression_matrix.txt.gz", sep="\t", nrows=4)
print(f"  Smillie(UC) 矩阵: 读到 {df.shape[0]} 基因 × {df.shape[1]-1} 细胞列, 头列={list(df.columns[:2])}")

print("[3] Geneformer 载模型 + 3090 前向")
from transformers import AutoModel, AutoConfig
mp = "models/Geneformer/Geneformer-V1-10M"
cfg = AutoConfig.from_pretrained(mp)
model = AutoModel.from_pretrained(mp).to("cuda").half().eval()
nparam = sum(p.numel() for p in model.parameters())/1e6
print(f"  模型载入: {cfg.model_type}, {nparam:.1f}M 参数, vocab={cfg.vocab_size}, hidden={cfg.hidden_size}")
ids = torch.randint(0, cfg.vocab_size, (4, 512), device="cuda")
with torch.no_grad():
    out = model(ids)
emb = out.last_hidden_state.mean(1)
print(f"  GPU 前向 OK: last_hidden={tuple(out.last_hidden_state.shape)} -> 细胞嵌入={tuple(emb.shape)}")
print("="*60)
print("[SMOKE_PASS] 全部通过：GPU 计算 + 数据读取 + 基础模型前向")
