# 任务3 · 多队列严格验证"基线黏膜 SASP 预测生物制剂应答"，并定量剥离泛炎症
# Task 3 · Rigorous multi-cohort validation of "baseline mucosal SASP predicts biologic response" + isolate senescence-specific signal from inflammation

> 本文件是**自包含**任务书。执行者无需任何聊天记录。所有背景、accession、工具、验证标准写在下面。
> 全程 **CPU、Python（机器未装 R）**。可与 task1、task2 完全并行、互不依赖。

---

## 项目背景（必读·已压缩）

主题：肠黏膜**细胞衰老 / SASP** 与炎症性肠病(IBD)，仅用公开数据，Python 栈，venv 在 `/Users/ujs/Downloads/lzy/.venv`。执行者：临床检验医师、生信新手，求 2026-09 前一篇**不烂大街**的中上论文，看重算力效率。硬件：现 ~200GB 磁盘。

**已有核心结果（本任务要把它做硬、做广）：**
- 黏膜 **SASP(SenMayo, 124 基因)** 在活动期 IBD 升高(p~1e-8)，**基线 SASP 越高越倾向无应答**：
  - **GSE16879**（133 样本，英夫利西 IFX 治前/后，平台 U133+2.0/GPL570）：基线预测无应答 **AUC 0.85**。
  - **GSE73661**（178 样本，维得利珠 VDZ/IFX，Gene1.0ST/GPL6244）：**AUC 0.74**。
  - 应答者治疗后 SASP 消退；预测力**独立于内镜严重度**（与基线 Mayo 不共线 rho=0.12，增量 ΔAUC≈+0.40）。
- **诚实点（本任务要定量回应）**：SASP = 分泌臂，与**泛炎症重叠**；经典停滞臂(p16/CDKN2A、p21/CDKN1A) 并不升高。必须证明 SASP 在**调整炎症严重度后仍有预测增量**。

**本任务目标：**
1. **跨多个公开 IBD 治疗应答转录组队列**，严格检验"基线 SASP 预测应答"的**泛化性**（留一队列/Meta）。
2. **定量把"衰老特异信号"从"泛炎症"中分离**（在转录组层呼应 task1 的遗传层检验）。

---

## 背景与动机（Why）

1. 现有结论仅基于两队列两药，审稿人会问泛化性与"是不是只是炎症"。**多队列留一法 + 炎症调整**能把结论从"两队列观察"升级为"跨队列稳健、且对炎症有增量"的硬证据。
2. 与 task1（遗传 MVMR 调整炎症）形成**双层呼应**：遗传层与转录组层都做"调整炎症后衰老是否仍有效"——叙事闭环、抗审稿。
3. CPU 即可、不烫手，**性价比最高的"加固"任务**。

---

## 具体做法（How：数据 / 工具 / 步骤）

### 1. 收集更多公开 IBD 治疗应答转录组队列（GEO）
- 在 **GEO** 搜索并**确认 accession 与应答标签**（抗 TNF 与抗整合素均要）。候选（执行者须逐一核实是否含黏膜表达谱 + 应答标签 + 治前/后）：
  - 抗 TNF（英夫利西/阿达木）：**GSE12251、GSE23597、GSE52746、GSE92415（戈利木单抗/PURSUIT）、GSE16879**（已有）、以及 **RISK / RISC** 队列(儿童 CD，注意其 GEO 编号)。
  - 抗整合素（维得利珠）：**GSE73661**（已有）等。
  - 其它可用的 IBD 黏膜应答队列(如 GSE14580、GSE3629 等)按检索结果纳入。
- 下载 series matrix + 平台注释到 `/Users/ujs/Downloads/lzy/data/geo/`，命名 `GSExxxxx_series_matrix.txt.gz`、平台 `GPLxxxxx.annot.gz`。
- **务必核对**：应答定义(内镜/临床缓解)、时间点(基线必须有)、组织(黏膜活检)、药物。不一致的队列单独标注或剔除。

### 2. 统一 probe→gene 映射与解析（复用现有流程）
- **复用 `src/03_parse_geo.py` 的模式**：用 GEO `.annot.gz`(列 `ID`, `Gene symbol`) 做 probe→symbol，多探针**按平均表达最大塌缩到基因**；线性值(max>100)自动 `log2(x+1)`；输出 `data/interim/GSExxxxx_expr.parquet` + `_meta.tsv`。把 `COHORTS` 字典扩成新队列的 `{GSE: GPL}`。
- 元数据里抽出 **基线/治后、应答/无应答、活动度(Mayo/CRP 若有)** 字段。

### 3. 评分：SenMayo + 停滞臂 + 炎症代理
- **复用 `data/genesets/senescence_sets.json`**（含 `SenMayo` 124 基因 + 停滞 core 等）做单样本评分（z-score/ssGSEA 风格，与现有 `src/07/10` 一致）。
- 同时算**停滞臂**(p16/CDKN2A、p21/CDKN1A、CDKN2B、GLB1、SERPINE1、MKI67 反向) 评分。
- **炎症代理**：CRP 样基因签名 + **免疫/中性粒细胞去卷积**（**CIBERSORTx** signature matrix，或 Python 的去卷积如 `nnls`/`scaden` 思路、xCell-like）估各队列炎症/中性粒细胞负荷，作为"泛炎症"协变量。

### 4. 跨队列预测与炎症调整
- **留一队列法(LOCO, leave-one-cohort-out)**：在 N-1 队列训练、留出队列测 AUC；或更简单：每队列内基线 SASP 单变量 LOO-AUC，再 **Meta 合并**(随机效应)得 pooled AUC + 95%CI + I²(异质性)。
- **炎症调整(核心)**：logistic 回归 `应答 ~ SASP + 炎症代理(CRP签名/中性粒去卷积/Mayo)`，看 **SASP 系数是否仍显著、ΔAUC 增量是否>0**。这在转录组层定量回答"衰老特异 vs 纯炎症"，**镜像 task1 的遗传检验**。
- 报告每队列 + 合并的 AUC、森林图、异质性。

---

## 大文件 / 效率规则（务必遵守）

- 若需对大表(如全基因表达矩阵、外部去卷积参考)做 ID/样本连接，**用 awk 哈希连接(单遍 O(1))，不要 `grep -f`**：
  ```bash
  gunzip -c BIG.tsv.gz | awk -F'\t' -v OFS='\t' \
    'NR==FNR{a[$1];next} ($1 in a){print $0}' keys.txt -
  ```
- 纯 **Python**(pandas/numpy/scipy/sklearn/statsmodels)，**不引入 R**（CIBERSORTx 可用其网页/容器或 Python 等价去卷积，避免 R 依赖；若必须 R 则单独隔离并注明）。
- 表达矩阵存 parquet（已是项目惯例）。

---

## 复用的现有文件（直接拿来用）

```
/Users/ujs/Downloads/lzy/src/03_parse_geo.py          GEO 解析 + probe→gene 模板
/Users/ujs/Downloads/lzy/src/07_senescence.py / 10_senescence_main.py  评分/预测模式
/Users/ujs/Downloads/lzy/src/11_incremental.py        增量(调整严重度)检验模式
/Users/ujs/Downloads/lzy/data/genesets/senescence_sets.json  SenMayo + 停滞 gene sets
/Users/ujs/Downloads/lzy/data/interim/GSE16879_scored.tsv, GSE73661_scored.tsv  既有评分(对照)
/Users/ujs/Downloads/lzy/data/interim/GSE{16879,73661}_expr.parquet + _meta.tsv  既有队列
```

---

## 预期结果（Expected）

- **pooled AUC + 95%CI**（多队列合并），预期与 0.74–0.85 同量级；给出**跨队列异质性(I²)**。
- **炎症调整后**：SASP 是否保留独立预测增量（ΔAUC>0、系数显著）。预期保留一部分增量(呼应已有 ΔAUC≈+0.40)，但须诚实报告剥离后剩多少。
- 停滞臂预期仍**弱**(对照)，强化解耦论点。

---

## 如何验证（Validation）

- **方向一致性**：各队列基线 SASP 高 → 无应答，方向应一致；个别相反需查队列质量/标签。
- **留一稳健**：LOCO-AUC 不应被单一队列主导(留一后仍 >0.65 较稳)。
- **炎症调整**：SASP 在加入 CRP 签名/中性粒去卷积/Mayo 后**仍显著**= 衰老特异；若被吸收为零 **如实报告**（与 task1 遗传层结论对照）。
- **平台稳健**：结论应跨 Affy U133/Gene-ST 等平台成立(批次/平台作敏感性)。

---

## 交付物（Deliverables）

写到 `/Users/ujs/Downloads/lzy/outputs/`（表/统计）与 `results/figures/`（图）、`results/tables/`：
1. `data/interim/GSExxxxx_expr.parquet` + `_meta.tsv` + `_scored.tsv`（每个新队列）。
2. `multicohort_auc.tsv` —— 每队列 + pooled AUC/CI/I²。
3. `inflammation_adjusted.tsv` —— SASP 调整炎症后的系数/ΔAUC。
4. `Fig_task3_multicohort.png` —— 森林图 + 调整前后 AUC。
5. 诚实小结 md：泛化性结论 + "衰老特异 vs 纯炎症"在转录组层的最终判读。
新脚本建议命名 `src/22_harvest_cohorts.py`、`23_score_all.py`、`24_loco_meta.py`、`25_inflammation_adjust.py`。

---

## 资源与注意（Compute / pitfalls）

- **纯 CPU**，现 ~200GB 磁盘足够（GEO series matrix 每个数十~数百 MB）。
- **坑**：(1) 应答定义不统一(内镜 vs 临床缓解 vs CRP) → 统一口径或分层报告；(2) 跨平台基因覆盖不同 → SenMayo 缺失基因要记录、按可用基因子集评分并注明；(3) 批次/平台效应 → 队列内相对评分 + 留一法缓解，必要时 ComBat(Python 版 `pycombat`)；(4) 去卷积参考矩阵需与组织匹配；(5) 小队列 AUC 不稳 → 报 CI、Meta 随机效应。
- **不要跑 task1/task2 的遗传/单细胞内容**；本任务专注 bulk 多队列。
