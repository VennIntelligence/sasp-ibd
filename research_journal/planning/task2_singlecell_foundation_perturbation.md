# 任务2 · 单细胞分辨率确认"真·衰老细胞" + 用基础模型对 CCL8/CXCR2 做计算扰动
# Task 2 · Single-cell: identify bona-fide senescent cells in IBD gut + in-silico perturb causal genes with foundation models

> 本文件是**自包含**任务书。执行者无需任何聊天记录。所有背景、数据 accession、工具、验证标准写在下面。
> 本任务是**方法学创新点**，需 **GPU（未来工作站 2×RTX 3090 + 多 TB 磁盘）**。可与 task1、task3 完全并行、互不依赖。

---

## 项目背景（必读·已压缩）

主题：肠黏膜**细胞衰老 / SASP** 与炎症性肠病(IBD)，仅用公开数据，Python 栈（机器未装 R），venv 在 `/Users/ujs/Downloads/lzy/.venv`。执行者：临床检验医师、生信新手，求一篇 2026-09 前**不烂大街**的中上论文，看重新颖性与算力效率。

**已有结果：**
- 转录组层：黏膜 SASP(SenMayo 124 基因) 在活动期 IBD 升高，**基线 SASP 预测生物制剂无应答**（AUC 0.85 GSE16879/英夫利西，0.74 GSE73661/维得利珠），独立于内镜严重度，应答者治疗后消退。
- **关键诚实点**：能预测应答的是衰老的**分泌(SASP)臂** = 与泛炎症重叠；**经典停滞臂(p16/CDKN2A、p21/CDKN1A) 并不升高**，黏膜反而在增殖(MKI67↑)。所以"是不是真有衰老细胞"在 bulk 层面**分不清**——这正是本任务要在单细胞层面解决的。
- 因果遗传层：MR+coloc 锁定 **CCL8（OR≈2.35 风险，PP4=0.955）** 与 **CXCR2（OR≈0.75 保护，PP4=0.951）**；CXCR2 遗传保护却在发炎黏膜升高（疑似代偿）。CCL8 是单核/巨噬趋化因子，CXCR2 是中性粒细胞受体。

**本任务两个目标：**
1. **(a) 确认 IBD 肠道是否存在 BONA FIDE（真·停滞型）衰老细胞**（不只是 SASP 分泌），以及**在哪类细胞**——直接回答 bulk 层面的解耦疑问。
2. **(b) 对因果基因 CCL8/CXCR2 做 IN-SILICO 扰动**（敲除/过表达），预测细胞状态向"炎症态/缓解态"的移动，**在机制层验证遗传方向**（CCL8 风险、CXCR2 保护）。

---

## 背景与动机（Why）

1. bulk 转录组无法区分"少量真衰老细胞"与"大量增殖细胞背景里的 SASP 分泌"——**只有单细胞能定位是哪类细胞、是真停滞还是纯分泌**。这把项目从"信号"推进到"细胞与机制"。
2. **方法学新颖、不撞车**：用**单细胞基础模型(foundation model, Geneformer / scGPT)** 做迁移学习 + **计算扰动(in-silico perturbation)**，在 IBD×衰老领域罕见。**关键好处：基础模型在数千万细胞上预训练，小样本本地数据不致命**——正好契合"公开数据+样本量有限"。
3. 把遗传层的因果方向(CCL8 risk / CXCR2 protective) 在**细胞机制层做独立验证**：扰动预测应与"应答者 vs 无应答者"差异、与遗传方向一致 → 三角闭环。

---

## 具体做法（How：数据 / 工具 / 步骤）

### 0. 取数据（公开 IBD 单细胞图谱）
下载到 `/Users/ujs/Downloads/lzy/data/scrna/`：
- **UC：Smillie 2019**，GEO **GSE116222** / Broad Single-Cell Portal **SCP259**（结肠上皮+基质+免疫，含炎症 vs 非炎症 vs 健康）。
- **CD：Martin 2019**，GEO **GSE134809**（回肠，含炎症/非炎症）。
- 可补：**Kong 2023**（多区段 IBD 图谱）、儿童 IBD 图谱等（在 GEO/SCP 搜索确认 accession 与元数据：病变状态、治疗应答标签若有）。
- 优先选**带炎症/非炎症标注**、最好**带治疗应答信息**的数据集。

### 1. QC 与预处理（scanpy，CPU 即可）
- `scanpy` 标准流程：过滤低质量细胞(min genes/counts、线粒体% 阈值)、双胞检测(scrublet)、归一化(log1p / 或 SCT 思路)、HVG、PCA、邻域、UMAP、Leiden 聚类。
- **细胞类型注释**：marker 基因 + 参考映射（如 CellTypist）。区分上皮(EPCAM)、成纤维/基质(COL1A1/PDGFRB)、髓系(LYZ/CD68)、T/B/浆细胞、内皮等。

### 2. 衰老评分：分离"停滞臂" vs "SASP 臂"（每细胞类型）
- 用**专用衰老工具**而非单一 signature：
  - **SenCID**（机器学习衰老分类器）、**hUSI**（human universal senescence index）——给每个细胞"是否衰老"的概率/指数。
  - **SenMayo / CellAge** 评分（scanpy `score_genes`）——SASP/分泌臂。
  - **停滞臂**单独评分：CDKN2A(p16)、CDKN1A(p21)、CDKN2B、GLB1、SERPINE1、LMNB1↓、MKI67↓。
- **关键分析**：在每个细胞类型里**分别**看"停滞评分"与"SASP 评分"——**是否存在某细胞群同时停滞高 + 增殖低 + SASP 高 = 真·衰老细胞**？还是 SASP 只来自增殖型炎症细胞？这直接回答 bulk 的解耦疑问。
- 衰老细胞分数 vs **疾病活动度/炎症状态** 的关联（按细胞类型）。

### 3. 基础模型迁移学习 + 计算扰动（GPU，方法学核心）
- **模型**：**Geneformer**（rank-value 编码、in-silico perturbation 原生支持）和/或 **scGPT**（细胞状态嵌入 + 扰动预测）。用**预训练权重**(tens of millions of cells)，对本地 IBD 数据做 zero-shot 嵌入或轻量微调。
- **细胞状态嵌入**：得到炎症态 vs 缓解态/健康态的 latent 表征；定位 CCL8/CXCR2 高表达细胞所处状态。
- **IN-SILICO 扰动**：对目标细胞**敲除/过表达 CCL8、CXCR2**，用模型预测嵌入移动：
  - 预期 **敲除 CCL8 → 远离炎症态**（与"CCL8 风险"一致）；
  - **敲除/抑制 CXCR2 → 是否移向炎症态？**（若是，则与"CXCR2 遗传保护、盲目抑制可能有害"一致——机制层验证该反向洞见）。
  - 也可做 SenMayo 命中基因(MMP3/CXCL8/MMP9/CXCL10 等)的扰动作对照。

### 4. 细胞-细胞通讯（SASP 发送→接收轴）
- **CellPhoneDB**（Python，优先，无需 R）或 CellChat（若装 R 则用），刻画 **SASP 发送细胞 → 接收细胞** 的配受体轴：重点 CCL8→CCR(单核/巨噬募集)、CXCL8/CXCL1/2/3 → **CXCR2**(中性粒细胞) 轴。
- 比较炎症 vs 非炎症、应答 vs 无应答（若有标签）的通讯强度差异。

---

## 预期结果（Expected）

- 一张**细胞类型 × (停滞评分, SASP 评分)** 图：明确"真·衰老细胞"是否存在、在哪类细胞（预期可能在基质/成纤维或上皮的少数亚群，而髓系主要贡献 SASP/炎症）。
- 衰老细胞分数随疾病活动升高（按细胞类型）。
- **扰动预测**：敲除 CCL8 使细胞远离炎症态；抑制 CXCR2 的方向与"保护性、勿盲目抑制"一致。
- **通讯图**：清晰的 SASP 发送→接收轴（CCL8/CXCR2 居中）。

---

## 如何验证（Validation）

- **衰老细胞比例 vs 疾病活动度**正相关（多数据集一致）。
- **跨工具一致性**：SenCID / hUSI / SenMayo 对"哪类细胞衰老"结论应大体一致；不一致则记录。
- **扰动预测 ↔ 观测**：模型预测的 CCL8/CXCR2 扰动方向，应与**应答者 vs 无应答者**的观测差异、以及**遗传方向(CCL8 risk, CXCR2 protective)** 一致——三方吻合 = 强证据。
- **跨数据集复现**：在 UC(Smillie) 与 CD(Martin) 至少两套数据上结论一致。
- **诚实底线**：若发现 IBD 肠道**几乎没有真停滞型衰老细胞**、信号纯由炎症细胞 SASP 贡献——**如实报告**，这恰恰强化了 bulk 的"解耦"结论，是有价值的阴性发现。

---

## 交付物（Deliverables）

写到 `/Users/ujs/Downloads/lzy/outputs/scrna/`（表/嵌入）与 `results/figures/`（图）：
1. 处理后的 `.h5ad`（每数据集）+ 细胞类型注释。
2. `senescence_per_celltype.tsv` —— 各细胞类型停滞/SASP 评分（多工具）。
3. `insilico_perturbation.tsv` + 图 —— CCL8/CXCR2 敲除/过表达的嵌入移动量与方向。
4. `cellcomm_SASP_axis.tsv/png` —— 发送→接收通讯（CCL8/CXCR2 轴）。
5. `Fig_task2_singlecell.png` —— 整合主图。
6. 诚实小结 md：IBD 是否存在真衰老细胞、在哪类细胞、扰动是否支持遗传方向。
新脚本建议放 `src/sc_*.py`（如 `sc_01_qc.py`、`sc_02_senescence.py`、`sc_03_foundation_perturb.py`、`sc_04_cellcomm.py`）。

---

## 资源与注意（Compute / pitfalls）

### 数据量估算（不靠"大数据"——基础模型已预训练于数千万细胞）
| 档位 | 细胞数 | 数据集 | 磁盘 |
|---|---|---|---|
| 最小可行 | 20–50 万 | Smillie UC(~36万) + Martin CD(~8万) + 对照 | ~5–15 GB |
| 舒适 | 50万–100万 | 再加 Kong 2023(~72万)、儿童 IBD、治疗前后队列 | ~20–50 GB |
| 模型权重 | — | Geneformer / scGPT 预训练权重(下载) | ~0.3–2 GB |

→ **总计 ~30–80 GB 即可**；**现有 200 GB 机器已能跑"最小可行版"**（QC+评分+零样本扰动），TB 工作站用于多模型/多实验的宽裕。
→ 真正吃**系统内存 RAM**（百万细胞 h5ad 读入），建议 **≥64 GB RAM**；不是显存。

### 显存估算（in-silico 扰动=推理，最省显存；微调才耗显存）
| 模型 | 参数 | 推理显存 | 微调显存 | 单张 3090(24G) |
|---|---|---|---|---|
| Geneformer(base) | ~10M | ~2–4 GB | ~8 GB | 轻松 ✅ |
| Geneformer-V2(large) | ~316M | ~8–12 GB | ~20–24 GB(需梯度技巧) | 1–2 张 |
| scGPT | ~50M | ~4–8 GB | ~12–20 GB | 单张可 ✅ |
| scFoundation | ~100M | ~8–12 GB | 宜双卡 | 2 张更稳 |
| UCE | 650M | ~16–22 GB | 多卡 | 双卡推理 |

→ **单张 3090(24G) 足以起步**（Geneformer/scGPT 嵌入+扰动+轻量微调）；**2×3090(48G)=舒适**（更大 batch、跑大模型、数据并行、多模型对照）。**无需 A100/80G。**
→ **时间(2×3090)**：嵌入~100万细胞 几十分钟–2小时；轻量微调 几小时；CCL8/CXCR2 扰动 分钟–小时级 → **几天出第一版**。
→ **执行策略**：先在**现有机器/单卡**跑最小可行版（零样本嵌入+扰动，验证流程通），再上工作站双卡扩规模/微调。

### 安装与坑
- **安装**：`torch`(CUDA 对应版本)、`scanpy`、`anndata`、`scrublet`、`celltypist`、`cellphonedb`，以及 **geneformer**(Hugging Face 权重) 和/或 **scgpt**(权重)。提前下载基础模型权重(大文件)。3090 为 Ampere，启用 **fp16/bf16** 省显存。
- **坑**：(1) 基础模型对基因 ID 体系(Ensembl/symbol、token 词表)敏感，输入需按模型词表对齐；(2) zero-shot 扰动的绝对数值不可尽信，**看方向与相对移动**；(3) 批次效应(多数据集)用 Harmony/scVI 或基础模型自身嵌入缓解；(4) CellPhoneDB 优先(避免 R 依赖)；(5) 评分用 control gene set 做背景校正(scanpy score_genes 自带)。
- **不要跑 task1/task3 的遗传/bulk-meta 内容**；本任务专注单细胞 + 基础模型。
