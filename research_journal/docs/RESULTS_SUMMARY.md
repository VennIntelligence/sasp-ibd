# 结果汇总（诚实版）— IBD 黏膜衰老/SASP 与生物制剂应答

> ⚠️ **2026-06-23 已重大转型**：项目已退出「衰老/SASP」框架，转向**「难治性 IBD 炎症-分泌模块的因果解剖」**。
> 本文档记录的是**转型前**的衰老/SASP 分析——它们没作废，但**定位已改**（SASP「预测应答」被证明只是泛炎症的别名、并非新颖；衰老停滞臂为零）。
> **先读新北极星：[`PIVOT_2026-06-23_causal_refractory_module.md`](PIVOT_2026-06-23_causal_refractory_module.md)**，它解释了为什么转、新论点、各结果的新定位。

> 自动化分析于 2026-06-21/22 完成。本文档如实记录**支持**与**不支持**的结果，未做美化。
> 共两大层：**因果遗传层（第3阶段，新增，最有冲击力）** + **转录组层（第1-2阶段）**。
> 小白友好的完整讲解见根目录 `讲解_给你听.md`。

---

## 第3阶段（旗舰·新增）：因果遗传层 —— MR + 共定位

**问题**：衰老/SASP 基因是 IBD 的"因果驱动者"还是"伴随路人"？
**方法**：用 eQTLGen（3.1万人血液 eQTL）为 82 个衰老基因找遗传工具变量 → 孟德尔随机化(MR) 检验对 de Lange 2017 IBD/CD/UC GWAS 的因果效应 → 共定位(coloc)排除连锁假象 → 与转录组三角验证。Python 自实现（无需 R）。

**关键结果（FDR 校正）**：
- **IBD 9 个因果基因**：`TNFRSF1A`(OR1.16,风险，**抗TNF靶点=阳性对照**)、`CCL8`(2.35)、`PTBP1`(2.13)、`HGF`(1.83)、`PLAUR`、`ICAM1`(风险)；`CXCR2`(0.75)、`MMP9`(0.79)、`IL1B`(0.59)(保护)。
- **CD 7 个 / UC 5 个**；`CXCR2`、`MMP9`、`TNFRSF1A` 跨亚型反复出现 → 稳健。
- **方向有正有反** → "衰老信号≠一律有害"，比单纯 signature 深刻。
**共定位(coloc)把关后——只有 2 个基因稳健过关（PP.H4>0.95）**：
- **CCL8**（OR≈2.35 风险，PP4=0.955）、**CXCR2**（OR≈0.75 保护，PP4=0.951）。
- 多数 MR 命中（ICAM1/FGF2/ETS2/PTBP1/GDF15…）PP.H3≈1 = **连锁假象，被 coloc 正确剔除**（体现严谨，不是抓到就信）。
- TNFRSF1A coloc 不确定（PP4=0.18，PP1/PP3 各 0.4）——blood eQTL + 该位点复杂所致，如实标注。

**三角验证（遗传 × 转录组）**：CXCR2 因果保护、coloc 过关、但黏膜活动期↑（FC +3.3/+1.9）——提示其升高是**代偿/修复**而非致病；CCL8 因果风险 + coloc 过关 + 黏膜↑，方向一致。

**可成药性 + 反向洞见**：CXCR2 有拮抗剂，但因其**遗传上保护**，盲目抑制可能有害——这是单纯 signature 得不出的临床提醒；CCL8/CCR 轴是更合理的风险靶点。

**为何新颖不撞车**：系统性地对整个衰老基因集做 MR+coloc 筛 IBD 因果基因（含跨 UC/CD/IBD），IBD×衰老领域尚无；把"相关"升级为"因果"，并得出方向性/可成药的具体结论。
图：`outputs/mr/Fig_CAUSAL_integrated.png`；表：`coloc_IBD.tsv`、`triangulation.tsv`、`mr_{IBD,CD,UC}.tsv`。

### 第3阶段补强（2026-06-22）：FinnGen R12 独立复制 + 炎症调整 MR

**FinnGen R12 独立复制（task1）**：
- 已下载并校验 FinnGen R12 三个 strict endpoint：`K11_IBD_STRICT`、`K11_CD_STRICT2`、`K11_UC_STRICT2`。输出表：`outputs/mr/replication_finngen.tsv`、`outputs/mr/finngen_replication_summary.tsv`。
- **CCL8 在 IBD 中强复现**：de Lange discovery 中 OR≈2.35、FDR≈0.030、PP4=0.955；FinnGen R12 IBD 中 OR≈4.37、p≈2.9e-6、PP4=0.950。结论：CCL8-IBD 风险信号有独立队列方向一致 + 共定位支持，是目前最硬的遗传因果点。
- FinnGen UC 也显示 CCL8 风险/共享因果信号（OR≈6.28、p≈1.9e-6、PP4=0.952），但 de Lange UC discovery 未达 FDR 显著，因此应写作 **FinnGen-positive extension**，不是正式 UC 复现。
- CCL8 在 FinnGen CD 不复现（OR≈0.87、p≈0.84、PP4=0.016），需如实报告。
- **CXCR2 只做方向性复现，不做 coloc 复现表述**：FinnGen IBD/UC Wald 方向仍为保护（IBD OR≈0.80，UC OR≈0.74），但 PP4 不过阈值（IBD PP4=0.619，UC PP4=0.438，CD PP4=0.010）。结论：CXCR2 可写成方向性支持，但不能写成独立共定位复制。

**CRP 调整 MVMR（task1，探索性）**：
- CRP GWAS 已下载并用于两暴露 MVMR。CXCR2 调整 CRP 后在 IBD/CD/UC 中仍保持保护方向（IBD OR≈0.749，p≈1.1e-10）。
- 但 CRP 条件 F 很弱（约 2.09），且 CCL8 可用工具不足，故该结果只能作为"调整泛炎症后 CXCR2 方向仍稳"的探索性支持，不能当作强 MVMR 结论。

**LTL→IBD 机体衰老 MR（task1，Codd 2021）**：
- LTL harmonised GWAS `GCST90002398` 已下载并校验；`src/20_ltl_mr.py` 可 headless 重跑。输出：`outputs/mr/mr_LTL.tsv`、`outputs/mr/mr_LTL_instruments.tsv`、`outputs/mr/mr_LTL_leave_one_out.tsv`，另有 MHC 排除敏感性 `outputs/mr/mr_LTL_sensitivity.tsv`。
- 593 个 genome-wide significant LTL SNP 经 500kb distance-clump 后进入候选；与 de Lange IBD/CD/UC 可 harmonise 的 SNP 数分别为 439/436/439。
- 主结果方向**反直觉**：遗传预测的**更长 LTL**与更高 IBD/CD/UC 风险相关（IVW：IBD OR≈1.26, p≈2.3e-15；CD OR≈1.43, p≈1.4e-21；UC OR≈1.12, p≈0.0028）。换成"更短 LTL/机体衰老"方向则是保护性估计（IBD OR≈0.79；CD OR≈0.70；UC OR≈0.90），不支持原先"短端粒增加 IBD 风险"假设。
- 稳健性：MR-Egger 斜率支持 IBD/CD，Egger 截距不显著（无明显方向性水平多效性）；leave-one-out 不翻转主方向。排除 chr6:25-34Mb MHC 后，IBD/CD 仍显著（IBD IVW OR≈1.21, p≈6.0e-10；CD OR≈1.30, p≈2.7e-11），UC 变弱（IVW OR≈1.12, p≈0.0045；weighted median≈1.00, p≈0.96）。
- 判读：LTL 不是支持"整体机体衰老驱动 IBD"的正交证据；更像是端粒生物学/免疫增殖/遗传多效性与 IBD 风险的复杂关系。可作为诚实的阴性/反向证据，与 bulk 里的"不是全局衰老时钟，而是 SASP/趋化分泌臂"相互呼应。

---

## TL;DR（一句话）

**肠黏膜「衰老相关分泌表型(SASP)」负荷在活动期 IBD 显著升高，且治疗前 SASP 越高越倾向"无应答"（refractory），可预测生物制剂应答（AUC 0.74–0.85，两队列两药验证），并在应答者治疗后消退；该预测力独立于内镜严重度（与基线 Mayo 不共线 rho=0.12，增量 ΔAUC=+0.40）。但信号来自衰老的"分泌(SASP)臂"而非经典"细胞周期停滞臂"(p16/p21)——二者在 IBD 中解耦。** 原计划的"转录组衰老时钟"为**阴性结果**，已如实保留为对照。

---

## 1. 衰老时钟：阴性结果（不回避）

- 时钟本身在 GTEx 正常肠道有效：交叉验证 r=0.70，MAE=7.8 岁，525 基因（`Fig1_clock_accuracy.png`）。
- **但套到 IBD 后不成立**：
  - IBD 基线 vs 对照 age acceleration 无差异（GSE16879 p=0.71，GSE73661 p=0.32），方向不一致。
  - 与疾病活动度(Mayo)无关（rho=0.007，p=0.96）。
  - 预测应答 ≈ 随机（AUC 0.54 / 0.56）。
- **解读**：IBD 黏膜并非"全局转录组年龄变老"。既往甲基化时钟显示的"加速衰老"未在转录组层面重现——这是一个有价值的阴性/对照结论，而非失败。

## 2. 衰老/SASP（SenMayo）：强且一致 ✅

| 指标 | GSE16879(IFX) | GSE73661(VDZ/IFX) |
|---|---|---|
| 基线 IBD vs 对照（升高） | p=4.8e-8 | p=3.3e-8 |
| 应答者配对前后下降 | p=4.5e-8 | p=9.5e-5 |
| 无应答者配对前后下降 | p=2.9e-4 | p=6.1e-3 |
| ΔR vs ΔNR（应答者降更多） | **p=0.003** | p=0.091（趋势） |
| **基线预测无应答 AUC** | **0.85** | **0.74** |

- 应答者治疗后下调最显著的 SASP 基因：**MMP3, CXCL8(IL8), MMP1, MMP10, MMP12, CXCR2, CXCL10, CXCL1, MMP9, CCL4**（经典 SASP/趋化-基质重塑，`Fig4_senescence.png`）。
- 图：`Fig_MAIN_senescence.png`（主图）。

## 3. 关键诚实点：SASP 臂 vs 停滞臂解耦 ⚠️➡️亮点

- 经典衰老停滞标志（p16/CDKN2A、p21/CDKN1A、CDKN2B、GLB1、SERPINE1 上调 − 增殖标志 MKI67/LMNB1 下调）**在 IBD 中并不升高**（基线 IBD 反而低于对照：GSE16879 −0.11 vs +1.02）——因为活动期黏膜在**增殖**（MKI67↑，隐窝再生）。
- 这些停滞标志预测应答很弱：core AUC 0.52–0.58；p16 单基因 AUC 0.53–0.64；p21 ~0.52。
- **含义**：能预测应答的"衰老信号"其实是**SASP=分泌/炎症臂**，不是真正的细胞周期停滞性衰老。
- **为什么这是创新而非缺陷**：既往 IBD"衰老基因 signature"文章把两臂混为一谈；我们用数据**把两臂解耦**——这正是与现有文献的差异化论点，也回应了"是不是只是炎症"的质疑（部分是，且我们诚实地证明了它来自分泌臂）。

---

## 4. 这对发表意味着什么（务实判断）

**能立住的论文（建议主线）：**
> *"在 IBD 黏膜中，衰老相关分泌表型(SASP)与细胞周期停滞程序解耦：SASP 负荷预测并随生物制剂应答消退，而经典衰老停滞标志不然。"*
- 卖点：① 临床有用的**应答预测**（AUC 0.85/0.74，两队列两药）；② 机制上的**SASP–停滞解耦**新观点；③ 诚实的时钟阴性对照增加严谨性。
- 现实定位：IF 3–6（*J Transl Med / Inflamm Bowel Dis / Frontiers Immunol / Clin Transl Med*）。**不撞车**，因为既往是"诊断 signature"，本文是"应答预测 + 两臂解耦 + 可逆性"。

**需要警惕的审稿质疑（已部分回应）：**
1. "SenMayo 预测应答 = 基线炎症/严重度预测应答？" —— **已做增量检验并基本排除**：在 GSE73661 中，SASP 与基线内镜 Mayo 评分**不共线**（rho=0.12, p=0.31）；SASP 单独 LOO-AUC=0.70，而内镜严重度单独仅 0.28，SASP 提供 **ΔAUC=+0.40** 的独立预测力。即 SASP 反映的是一个**独立于内镜严重度**的"难治性分泌表型"。（仍建议补组织学严重度/CRP 等做更全调整。）
2. 方向：两队列一致——**基线 SASP 越高越倾向无应答**（R 中位 ≈ −0.01 vs NR ≈ +0.30~0.38），提示"SASP-high"为难治亚型，临床可解释。
3. 样本量与跨平台——已用两队列/两药 + 队列内相对比较缓解。

---

## 5. 下一步（让结论更硬，1–2 周内可做）

1. **增量价值检验**：logistic 中加入基线 Mayo/炎症评分，看 SASP 是否仍显著（剥离"纯炎症"解释）。
2. **单细胞定位**：哪类细胞(上皮/成纤维/髓系)表达预测性 SASP——把"炎症"细化到细胞来源，提升机制深度。
3. **可进血标志物**：从预测性 SASP 基因里挑分泌型(MMP3、CXCL8、MMP9、CXCL10)，作为你**血样验证**靶点（ELISA）。

---

## 6. 产出文件

- 图：`Fig1_clock_accuracy.png`、`Fig2_accelerated_aging.png`、`Fig3_reversal.png`、`Fig4_senescence.png`、`Fig5_predictor.png`、`Fig_MAIN_senescence.png`、`ALL_FIGURES.png`
- 统计：`outputs/*_stats.json`、`outputs/clock_model.json`
- 评分表：`data/interim/GSE16879_scored.tsv`、`GSE73661_scored.tsv`
- 脚本：`src/01..10`，方法学：`METHODS_AND_NARRATIVE.md`
