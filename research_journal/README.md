# 研究日志索引 · Research Journal Index

> IBD 黏膜「细胞衰老 / SASP」与生物制剂应答 —— 公开数据 + 纯 Python 流水线（电脑未装 R）。
> 本目录是整个项目的**实验室笔记（过程档案）**：决策、死胡同、方法叙事、任务书。
> 先看本页，再按需要钻进 `docs/`、`planning/`。
> 论文最终图表已不在本目录，统一晋升到仓库根 `results/`（见根 `README.md` / `AGENTS.md`）。

---

## 🚩 START HERE · 新窗口接手先读这里

**如果你是新会话/新执行者，按此顺序读，3 分钟即可接上进度：**
1. 本 README（项目全貌 + 头条结果 + 诚实局限，见下文）。
2. `docs/RESULTS_SUMMARY.md`（结果细节）、`docs/讲解_给你听.md`（因果/MR/coloc 大白话）。
3. `planning/task1_*.md`、`planning/task2_*.md`、`planning/task3_*.md`（三个**并行、互不依赖**的下一步任务书，自包含可直接开干）。

**🚩 当前方向（最新决定，2026-06-23，已重大转型）**：**退出「衰老/SASP」框架**，转向 **「难治性 IBD 炎症-分泌模块的因果解剖」**。
完整理由与新论点见 **[`docs/PIVOT_2026-06-23_causal_refractory_module.md`](docs/PIVOT_2026-06-23_causal_refractory_module.md)**（**新北极星，先读它**）。
一句话：第一性原理探索证明「SASP 预测应答」只是泛炎症的别名（不新颖）、衰老停滞臂为零；项目真正稀缺的资产是 **MR+coloc 因果引擎**，现把它从「82 个衰老基因」**重新对准数据/文献真正指认的难治炎症模块**（OSM/OSMR/TREM1/IL13RA2/CXCR2/CCL8/中性粒基因…），并改用**肠道 eQTL**——问「这些公认炎症标志物里哪些是因果、哪些是路人」。
- 旧 task1/2/3 任务书与下方衰老/SASP 头条**保留为历史证据链**，定位已按 PIVOT 文档改写。

**注意**：本项目目录已配持久记忆（`MEMORY.md`），新会话会**自动加载**用户画像、项目状态与"awk 优于 grep"等约定——不必从头交代。数据/代码在仓库根：`src/01-17`、`data/raw|external|interim`（eQTLGen 4.6G、GWAS、单细胞待下）、机器产物 `outputs/`、论文终版 `results/`。

---

## 一句话项目简介（大白话）

我们想搞清楚：**肠子里的"细胞衰老"到底是炎症性肠病(IBD)的"真凶"，还是跟着炎症一起涨的"路人"。**
分两层做：(1) **转录组层**——用公开肠黏膜表达谱算"衰老分泌信号(SASP)"，发现它在活动期 IBD 升高、且**治疗前 SASP 越高越难治（能预测生物制剂应答）**；(2) **因果遗传层**——用孟德尔随机化(MR)+共定位(coloc) 把"相关"升级为"因果"，锁定 **CCL8（风险）** 与 **CXCR2（保护）** 两个真·因果基因。全部公开数据、Python 自实现。

---

## 目录结构

```
research_journal/                 ← 实验室笔记（过程，不是产品）
├── README.md          ← 本页（索引）
├── docs/              ← 文字记录：计划 / 方法 / 结果叙事 / 小白讲解
└── planning/          ← 下一步并行任务书（task1/2/3，自包含可直接开干）

# 产物不在这里，按"出身"分在仓库根：
#   outputs/  机器自动吐、可弃（按脚本名）
#   results/  晋升而来的论文最终 figures/ + tables/
```

### docs/ —— 文字记录

| 文件 | 一句话说明 |
|---|---|
| `IBD_aging_project_plan.md` | 最初的项目计划与立意（衰老×IBD 选题、阶段划分）。 |
| `METHODS_AND_NARRATIVE.md` | 方法学与叙事主线：数据来源、时钟构建、SASP 评分、MR/coloc 流程的方法描述。 |
| `RESULTS_SUMMARY.md` | **结果汇总（诚实版）**：支持/不支持的结果都如实记录，含 SASP 预测应答 + MR/coloc 因果发现 + 局限。 |
| `讲解_给你听.md` | **小白友好讲解**：用比喻把"相关 vs 因果""MR""共定位"讲清楚，可直接拿来汇报。 |

### 论文最终成图（已晋升到根 `results/figures/`）

| 文件 | 一句话说明 |
|---|---|
| `OVERVIEW.png` | 全项目总览图。 |
| `Fig1_clock_accuracy.png` | GTEx 肠道转录组衰老时钟交叉验证精度（r=0.70，对照用）。 |
| `Fig_MAIN_senescence.png` | **主图**：SASP 在活动期升高、随应答消退、预测应答(AUC)。 |
| `Fig_CAUSAL_integrated.png` | **因果整合图**：MR + coloc + 三角验证（CCL8/CXCR2）。 |
| `MR_forest_IBD.png` | IBD 各衰老基因 MR 效应森林图。 |

### 论文最终结果表（已晋升到根 `results/tables/`）

| 文件 | 一句话说明 |
|---|---|
| `instruments.tsv` | 82 个衰老基因的 MR 工具变量（rsid、assessed/other 等位、eaf、beta/se，由 eQTLGen Z 经 Zhu 2016 公式换算）。 |
| `mr_IBD.tsv` / `mr_CD.tsv` / `mr_UC.tsv` | 三个结局(IBD/CD/UC) 的 Wald-ratio MR 结果（OR、p、FDR）。 |
| `coloc_IBD.tsv` | 共定位(coloc ABF) 后验概率 PP.H0–H4，用于剔除 LD 假象。 |
| `triangulation.tsv` | 遗传(MR/coloc) × 转录组(黏膜表达) 三角验证汇总表。 |

> 产物按"出身"分两层（见根 `AGENTS.md`）：机器自动吐的图/`*_stats.json`/MR scratch 在
> `outputs/`（含 `outputs/mr/`），上表这批**挑中的终版**已晋升到 `results/figures|tables/`；
> 评分表（可重生派生数据）在 `data/interim/GSE*_scored.tsv`。

---

## 头条结果（Headline）

1. **SASP 预测生物制剂应答**：治疗前黏膜 SenMayo(124 基因) 负荷越高越倾向**无应答**，**AUC 0.85（GSE16879/英夫利西单抗）、0.74（GSE73661/维得利珠单抗）**；应答者治疗后 SASP 消退。预测力**独立于内镜严重度**（与基线 Mayo 不共线 rho=0.12，增量 ΔAUC≈+0.40）。
2. **因果遗传锁定 2 个基因**：MR(eQTLGen→de Lange 2017 GWAS) + coloc 把关后，仅 **CCL8（OR≈2.35，风险，PP4=0.955）** 与 **CXCR2（OR≈0.75，保护，PP4=0.951）** 稳健过关；多数 MR 命中是 PP.H3≈1 的 LD 假象，被正确剔除。**TNFRSF1A（抗 TNF 靶点）作为阳性对照在 MR 中被揪出**。
3. **方向有正有反**：CXCR2 遗传上保护、却在发炎黏膜中升高 → 其升高更像代偿/修复；提示"盲目抑制 CXCR2 可能有害"——单纯 signature 得不出的临床反向洞见。

---

## 诚实局限 / 待解问题（Limitations & Open Questions）

- **"衰老→IBD"框架部分是被"硬套"的**：能预测应答的是衰老的**分泌(SASP)臂**，与经典**停滞臂**(p16/CDKN2A、p21/CDKN1A) **解耦**——后者在 IBD 中并不升高（黏膜在增殖，MKI67↑）。所以"衰老驱动 IBD"部分是给"炎症"换了标签。
- **稳健因果命中(CCL8/CXCR2) 都是趋化因子**，可能"只是炎症"。→ 需在遗传层(MVMR 调整炎症/CRP) 与转录组层(去卷积调整炎症) 同时把"衰老特异信号"从"泛炎症"中分离出来。
- **转录组时钟是阴性结果**（IBD 黏膜并非全局转录组"变老"），已如实保留为对照，非失败。
- **eQTLGen 是血液 eQTL**，肠道特异性有限；诱导型 SASP 因子(IL6/IL8) 在健康组织不表达、找不到工具变量。
- **样本量/跨平台有限**（两队列两药）；需多队列留一法验证 + 独立 GWAS(FinnGen R12) 复制。
- **机制未在单细胞分辨率定位**：是否存在真·衰老细胞、哪类细胞、CCL8/CXCR2 的发送→接收轴尚待单细胞 + 计算扰动(foundation model) 验证。

> 上述待解问题对应 `planning/` 下三个并行任务：MVMR 因果分离 / 单细胞+基础模型扰动 / 多队列验证。
