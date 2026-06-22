# IBD 黏膜「细胞衰老 / SASP」与生物制剂应答

公开数据 + 纯 Python 流水线（本机未装 R）。核心问题：**能否把"细胞衰老"从"泛炎症"中
因果地剥离出来——IBD 黏膜里到底有没有真·衰老细胞？**

- **转录组层**：SenMayo SASP 负荷在活动期 IBD 升高，且**治疗前越高越难治**
  （预测生物制剂应答，AUC 0.85 / 0.74），独立于内镜严重度。
- **因果遗传层**：MR（eQTLGen→de Lange 2017）+ coloc 把关后，**CCL8（风险）**、
  **CXCR2（保护）** 两个基因稳健过关；多数 MR 命中是 LD 假象被正确剔除。

> 头条结果、诚实局限、下一步并行任务见 `research_journal/README.md`。

## 目录地图

| 目录 | 内容 |
|---|---|
| `data/raw/` | 只读源数据（GWAS / eQTLGen / GEO / GTEx）——不可变 |
| `data/external/` | 参考资源（gene sets） |
| `data/interim/` | 可重生派生数据（scored.tsv、parquet 矩阵） |
| `src/` | `01–17` 编号流水线 + `paths.py` 中枢路径 |
| `notebooks/` | 仅探索 / 出图原型（不进流水线） |
| `outputs/` | **机器自动吐**、可弃、按脚本名归档 |
| `results/` | **晋升而来**的论文最终 `figures/` + `tables/` |
| `paper/` | 正文 / 参考文献 / 补充材料 |
| `research_journal/` | 实验室笔记：方法叙事(`docs/`)、任务书(`planning/`) |

**`outputs/` vs `results/`**：脚本只写 `outputs/`（可随便重跑）；你**主动挑中**的终版
才晋升进 `results/`。两层模型与工程标准见 **`AGENTS.md`**。

## 怎么跑

```bash
source .venv/bin/activate            # Python 3.14
python src/02_build_clock.py         # 例：构建衰老时钟
python src/paths.py                  # 打印中枢路径表，自检布局
```

脚本按编号顺序（`01_gtex_samples` → … → `17_integrated_figure`）构成流水线；
新代码用 `from paths import P` 取路径，不要硬编码绝对路径。
