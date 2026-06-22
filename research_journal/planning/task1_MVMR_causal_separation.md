# 任务1 · 因果上把"衰老/SASP"与"泛炎症"分离，并检验"机体衰老"对 IBD 的因果性
# Task 1 · Causally separate senescence/SASP from generic inflammation + test organismal-aging causality

> 本文件是**自包含**任务书。执行者无需任何聊天记录。所有需要的背景、数据路径、accession、公式、效率规则都写在下面。
> 全程 **summary-stat only，纯 CPU，Python（机器未装 R）**。可与 task2、task3 完全并行、互不依赖。

---

## 项目背景（必读·已压缩）

主题：肠黏膜**细胞衰老 / 衰老相关分泌表型(SASP)** 与炎症性肠病(IBD)，仅用公开数据。
执行者画像：临床检验医师、前 IBD 研究者、生信新手，目标 2026-09 前出一篇中等偏上、**不烂大街**的论文。硬件：现 ~200GB 磁盘；未来工作站 多 TB + 2×RTX 3090。venv 在 `/Users/ujs/Downloads/lzy/.venv`。

**已有结果（两层）：**
- 转录组层：黏膜 SASP(SenMayo, 124 基因) 在活动期 IBD 升高，**基线 SASP 预测生物制剂无应答**（AUC 0.85 GSE16879/英夫利西，0.74 GSE73661/维得利珠），独立于内镜严重度。**诚实点**：信号来自 SASP/分泌臂 = 与泛炎症重叠；经典停滞标志(p16/CDKN2A、p21/CDKN1A) 并不升高（黏膜在增殖，MKI67↑）。
- 因果遗传层（本任务的基础）：用 **eQTLGen 血液 cis-eQTL** 为 82/124 个 SenMayo 基因建工具变量 → **de Lange 2017 IBD/CD/UC GWAS**(GWAS Catalog GCST004131/2/3，harmonised GRCh38) 做 **Wald-ratio MR**，FDR 校正。命中：IBD 9、CD 7、UC 5。**TNFRSF1A（抗 TNF 靶点）= 阳性对照**。**coloc(自实现 ABF) 后仅 CCL8（OR≈2.35 风险，PP4=0.955）与 CXCR2（OR≈0.75 保护，PP4=0.951）过关**，其余多为 PP.H3≈1 的 LD 假象。

**本任务要回答的核心问题（项目最大软肋）：**
> 稳健因果命中 CCL8/CXCR2 都是**趋化因子**，很可能"只是炎症"。**这到底是真·衰老，还是泛炎症？** 本任务用遗传方法直接把二者分开，并把"机体衰老"本身也拉来做因果检验。

---

## 背景与动机（Why）

1. **审稿人会问**："SASP 预测应答会不会只是炎症严重度预测应答？基因层面 CCL8/CXCR2 会不会只是 IBD 风险位点顺带富集在趋化因子上？"——必须有因果层面的"调整后仍显著"证据。
2. **纯统计/MR 方法本身在变得烂大街**——靠"多做几条 MR"不够新；本任务的差异化在于**多变量 MR(MVMR) 做条件因果分离** + **衰老硬指标(端粒) 的器官层因果** + **独立 GWAS 复制**，三管齐下回答"是不是真衰老"。
3. 把项目从"又一个 signature/又一条 MR"升级为"**能在调整泛炎症后仍站得住的、方向明确的因果结论**"。

---

## 具体做法（How：数据 / 工具 / 步骤）

### (i) 多变量 MR(MVMR)：衰老 vs 泛炎症，谁还因果？
- **目标**：把 SASP 基因暴露与"泛炎症"暴露放进同一个模型，看衰老在调整炎症后是否仍因果。
- **暴露 A（衰老/SASP）**：复用 `/Users/ujs/Downloads/lzy/outputs/mr/instruments.tsv`（含 rsid、assessed/other 等位、eaf、beta_eqtl、se_eqtl），重点 CCL8、CXCR2，必要时纳入其它 SenMayo MR 命中基因。
- **暴露 B（泛炎症）**：从公开 GWAS 取 **CRP**（如 Said 2022 / GCST90029070 系列，或 CHARGE CRP）作为泛炎症代表；可补 **IL6 / TNF**（如 Folkersen 2020 / Ahola-Olli 2017 循环蛋白 pQTL，eQTL Catalogue blood QTD000356，或 deCODE/UKB-PPP TNF/IL6 pQTL）。在 GWAS Catalog 确认最终 accession 与 harmonised 文件。
- **方法**：MVMR-IVW（多暴露线性，残差相关需用 GMM 或 MVMR 包思路自实现；新手可先做"两步法/条件 F + MVMR-IVW"）。条件 F 统计量检验工具强度，Q 检验异质性。输出每个暴露在**互相调整后**的条件因果效应与 p。
- **判读**：若 CCL8/CXCR2 在调整 CRP/IL6/TNF 后**仍显著且方向不变** → 支持"衰老特异、非纯炎症"；若被吸收为零 → 诚实记录"以炎症为中介/不可分"。

### (ii) 双向 & 反向 MR：钉死方向
- **正向**：衰老基因表达 → IBD（已有，复用）。
- **反向**：IBD → 衰老基因表达。用 de Lange IBD GWAS 的**独立显著 SNP（distance-clump，见下）** 作工具，以 eQTLGen 中 CCL8/CXCR2（及其它命中基因）为结局，IVW。
- **判读**：若反向不显著、正向显著 → 方向稳固"衰老在前、病在后"；若反向也显著 → 标注双向/反向因果可能（炎症推高 SASP）。

### (iii) 衰老硬指标 MR：白细胞端粒长度(LTL) → IBD
- **暴露**：**Codd 2021 LTL GWAS，GWAS Catalog GCST90002398**（~47 万 UKB）。下载 harmonised summary stats。
- **工具**：取全基因组显著（p<5e-8）SNP，**distance-clump**（按 rsid 的染色体/位置，窗口 ±500kb 取每窗口最显著 SNP；因无 R/无 LD 面板，用距离剪枝近似 LD-clump，并在局限里注明）。
- **结局**：de Lange IBD/CD/UC（同下文路径）。
- **方法**：**IVW + MR-Egger（截距测多效性）+ 加权中位数**三法一致性；留一法敏感性。
- **判读**：LTL 短(衰老) 是否因果增加 IBD 风险 → 这是"机体层面衰老"对 IBD 的直接因果检验，独立于 SASP 基因，**正交证据**。

### (iv) 独立 GWAS 复制 CCL8/CXCR2 的 cis-MR + coloc
- **独立结局**：**FinnGen R12** 端点 `K11_IBD` / `K11_CD` / `K11_UC`（或 `CHRONSMALL`/`ULCERCOLI` 视版本命名，从 FinnGen 官网 manifest 确认精确端点名与 summary-stat 下载 URL）；如可，再加 UKB IBD。
- **流程**：对 CCL8、CXCR2 重跑 cis-MR(Wald) + coloc(下面 ABF 公式)，看方向与 PP.H4 是否在独立队列复现。
- **判读**：方向一致 + PP.H4 仍高 → 复制成功，结论硬。

---

## 关键公式（直接照抄，已在本项目验证过）

**Zhu 2016：Z → beta/se**（eQTLGen 给的是 Z 和样本量 N，需配 assessed 等位频率 p）：
```
se   = 1 / sqrt( 2 * p * (1-p) * (N + Z^2) )
beta = Z * se
```
（p = assessed 等位频率；由 eQTLGen AF 文件给 AlleleB + AlleleB_all freq，按 assessed/other 对齐取 p 或 1-p。）

**coloc ABF（Giambartolomei 2014，self-implemented，已跑通）**：每个 SNP 的近似贝叶斯因子
```
labf(z, V, W) = 0.5 * ( log(1 - r) + r * z * z ),  其中 r = W / (V + W)
```
eQTL 端 V_e = 1 / (2 * N_e * f * (1-f))，W_eqtl=0.15^2；GWAS 端 V_g = se^2，z_g=beta/se，W_gwas=0.2^2。
先验 P1=P2=1e-4, P12=1e-5。
按 H1=logsumexp(l1)、H2=logsumexp(l2)、H4=logsumexp(l1+l2)，H3 用 log-diff-exp 稳定计算；
最后 logs=[0, logP1+H1, logP2+H2, logP1+logP2+H3, logP12+H4]，softmax 得 PP.H0..H4。
**PP.H4>0.9 = 共享同一因果变异（真）；PP.H3≈1 = LD 假象（剔除）。**

---

## 大文件效率规则（务必遵守）

- **大 GWAS / eQTLGen 全量文件用 awk 哈希连接（单遍 O(1) 查找），不要用 `grep -f`**（后者对大文件极慢）。模板（本项目已用过）：
  ```bash
  # 先把需要的 rsid 写进 _rsids.txt，再：
  gunzip -c BIG.tsv.gz | awk -F'\t' -v OFS='\t' \
    'NR==FNR{a[$1];next} ($RSCOL in a){print $RSCOL,$BETACOL,$SECOL}' _rsids.txt -
  ```
- **eQTLGen 是 hg19，GWAS 是 GRCh38 → 一律用 rsid 作连接键**（坐标不同，rsid 安全）。
- 先 `gunzip -c FILE | head -1` 定位列号（awk 用 1-based），再传 `-v` 变量给 awk。
- 纯 Python：numpy / pandas / scipy；**不要引入 R 依赖**。

---

## 复用的现有文件（直接拿来用）

```
/Users/ujs/Downloads/lzy/outputs/mr/instruments.tsv   工具变量(rsid/等位/eaf/beta/se)
/Users/ujs/Downloads/lzy/outputs/mr/mr_{IBD,CD,UC}.tsv  既有 MR 结果(OR/p/FDR)
/Users/ujs/Downloads/lzy/outputs/mr/coloc_IBD.tsv      既有 coloc PP
/Users/ujs/Downloads/lzy/outputs/mr/triangulation.tsv  遗传×转录组三角验证
/Users/ujs/Downloads/lzy/data/eqtlgen/cis_sig.txt.gz   显著 cis-eQTL(列含 SNP,GeneSymbol,Zscore,AssessedAllele,OtherAllele,NrSamples,Pvalue)
/Users/ujs/Downloads/lzy/data/eqtlgen/cis_full.txt.gz  全 cis(4.6GB，区域 coloc 用)
/Users/ujs/Downloads/lzy/data/eqtlgen/snp_af.txt.gz    AF(列: SNP,...,AlleleA,AlleleB(第5列),...,AlleleB_all(第9列))
/Users/ujs/Downloads/lzy/data/gwas/{IBD,CD,UC}.h.tsv.gz de Lange 2017 harmonised(列含 hm_rsid,hm_beta,beta,standard_error)
/Users/ujs/Downloads/lzy/data/genesets/senescence_sets.json  含 SenMayo 基因列表
```
参考脚本风格：`src/13_build_instruments.py`（Z→beta/se）、`src/14_mr.py`（Wald MR/FDR）、`src/15_coloc.py`（ABF + awk 哈希连接）。

新增 GWAS（CRP/IL6/TNF/LTL/FinnGen）请下载到 `data/gwas/`，建议命名 `CRP.h.tsv.gz`、`LTL_Codd2021.tsv.gz`、`finngen_R12_K11_IBD.tsv.gz` 等。

---

## 预期结果（Expected）

- **MVMR 表**：CCL8/CXCR2 在调整 CRP(±IL6/TNF) 后的条件 OR/p。预期 CXCR2 的保护方向较可能稳健（机制独立于泛炎症方向）；CCL8 部分可能被炎症吸收——无论哪种，**都是对核心问题的实质回答**。
- **反向 MR**：预期正向显著、反向弱/不显著（方向稳固）；若反向亦显著则记录双向。
- **LTL→IBD**：给出 IVW/Egger/加权中位数三法 OR 与一致性；这是"机体衰老"对 IBD 的正交因果证据。
- **复制**：CCL8/CXCR2 在 FinnGen R12 的方向与 PP.H4 复现情况。

---

## 如何验证（Validation）

- **阳性对照**：在新管线里 TNFRSF1A 仍应为 IBD 风险（OR>1）——验证流程正确。
- **MVMR 工具强度**：条件 F 统计量 >10；报告 Q/异质性。
- **LTL**：三法方向一致、Egger 截距不显著（无明显多效性）才采信。
- **复制一致性**：独立 GWAS 方向与发现队列一致 = 通过；不一致则诚实标注。
- **诚实底线**：若衰老在调整炎症后被吸收为零，**如实写"信号以炎症为中介/不可分"**——这本身就是对项目软肋的有价值回答，不得粉饰。

---

## 交付物（Deliverables）

写到 `/Users/ujs/Downloads/lzy/outputs/mr/`（表）与 `/Users/ujs/Downloads/lzy/results/figures/`（图）：
1. `mvmr_results.tsv` —— CCL8/CXCR2 ± 炎症暴露的条件因果效应表 + 条件 F/Q。
2. `reverse_mr.tsv` —— IBD→衰老基因 反向 MR 结果。
3. `mr_LTL.tsv` —— LTL→IBD/CD/UC（IVW/Egger/加权中位数 + 留一）。
4. `replication_finngen.tsv` —— CCL8/CXCR2 在 FinnGen R12 的 cis-MR + coloc。
5. `Fig_task1_causal_separation.png` —— 一张整合图（MVMR 条件效应 + LTL 森林 + 复制）。
6. 一段诚实小结(写入 `research_journal/docs/RESULTS_SUMMARY.md` 的新小节或单独 md)：核心问题"是真衰老还是炎症"的最终判读。
新脚本建议命名 `src/18_mvmr.py`、`19_reverse_mr.py`、`20_ltl_mr.py`、`21_replicate_finngen.py`。

---

## 资源与注意（Compute / pitfalls）

- 纯 CPU、数小时内可完成；磁盘主要被新下载 GWAS 占用（CRP/LTL/FinnGen 各数百 MB~数 GB），现有 ~200GB 足够。
- **坑**：(1) 无 R/无 LD 面板 → distance-clump 近似，务必在局限注明；(2) MVMR 暴露间工具重叠/弱工具会偏倚，先查条件 F；(3) CRP/pQTL GWAS 的等位与坐标需与结局**按 rsid 对齐并校正等位方向**（assessed 等位翻转时 beta 取负）；(4) FinnGen 端点命名各版本不同，务必核对 R12 manifest；(5) pQTL 若用 cis 工具更可信(减少多效性)。
- **不要跑 task2/task3 的内容**；本任务只做遗传/summary-stat 层。
