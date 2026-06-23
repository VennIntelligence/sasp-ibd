# 写作交接文档 · FINDINGS FOR WRITE-UP（2026-06-23 定稿）

> 给写作智能体：**本文件自包含**，写论文只需读它 + 它引用的 `results/tables/*` 与 `results/figures/*`，
> 不需要聊天记录。所有数字均来自已落盘的结果表（括号内给了文件名）。
> **务必照"诚实框架/不可夸大"一节执行**——本项目的价值正在于诚实的因果分诊，夸大会被审稿当场打掉。

---

## 0. 一句话主线（HEADLINE）

> **难治性 IBD 由一个中性粒/髓系炎症-分泌程序标记；对该程序系统做 MR+coloc（血液/肠/免疫细胞多语境 eQTL）后，
> 众多"明星 marker"中只有 CXCR2 是稳健因果——且反直觉地保护——据此预测 CXCR2 拮抗剂在 IBD 注定碰壁
> （与真实失败试验吻合），而 OSM/TREM1/IL13RA2 等被坐实为非因果路人。**

**论文不是**：衰老论文、新预测器论文、或多基因药靶图谱（见第 5 节）。

---

## 1. 四个头条贡献（按重要性）

1. **因果分诊（去伪存真）**：用 MR+coloc 把难治炎症模块的"相关 marker"分成因果 vs 路人。**只有 CXCR2 稳健因果**；
   **OSM/OSMR/TREM1/IL13RA2 是路人**（TREM1 是范例：中性粒 eQTL p=2.4e-27 极强，但 MR OR=1.002, p=0.92，coloc PP4≈0）。
2. **CXCR2 可成药悖论（translational 卖点）**：CXCR2 是药靶（拮抗剂 navarixin/AZD5069 等），但**遗传上保护**
   → 药靶 MR 预测抑制它"碰壁/有害" → **与真实 CXCR2 拮抗剂 IBD/炎症试验失败吻合**；anti-TNF(TNFRSF1A) 作阳性对照方向吻合。
3. **机制（深化悖论）**：CXCR2↑ → 循环中性粒↑；中性粒计数本身是 IBD 因果风险因子；**可 CXCR2 仍保护** →
   保护**不是靠减少中性粒**，而是中性粒的**定位/稳态**层面（朴素"抑制 CXCR2 减少肠道中性粒"的逻辑因此是错的）。
4. **诚实阴性弧（增严谨）**：衰老停滞臂为零、转录组衰老时钟阴性、LTL MR 反向/阴性；SASP/"预测器"只是泛炎症的别名（见第4节）。

---

## 2. 因果核心（精确数字）

### CXCR2 —— 稳健因果保护（论文主角）
来源表：`results/tables/module_causal_map_multicontext.tsv`、`finngen_cxcr2.tsv`、`mr_sensitivity.tsv`、`mvmr_crp.tsv`、`reverse_mr_proper.tsv`、`steiger.tsv`、`coloc_IBD.tsv`
- 血液 eQTLGen：MR **OR 0.753**（p 1.6e-09, FDR 1.24e-07），**coloc PP4 0.951**。
- 中性粒 BLUEPRINT(QTD000026)：MR **OR 0.851**（FDR 1.05e-07），**coloc PP4 0.935**（细胞类型匹配的独立确认）。
- 独立队列 **FinnGen R12 复制**：IBD OR 0.796（p 1.3e-05）、UC OR 0.736（p 1.6e-06）；CD 不显著(0.968, p 0.77)。
- 多工具敏感性（中性粒, 5 SNP）：IVW 0.842 / MR-Egger 0.841（截距无方向性多效性）/ weighted-median 0.846 / LOO ~0.84。
- MVMR 调 CRP：方向稳（IBD OR 0.749, p 1.1e-10）；注：CRP 条件 F 弱(~2.1)。
- 反向因果**已排除**：排 cis 区后 IBD→CXCR2 反向无工具/无信号（先前 p=3e-200 是共定位 cis 变异伪反向）。Steiger 支持 exposure→outcome。

### CCL8 —— 提示性风险（**降级，需谨慎措辞**）
来源表：`triangulation.tsv`、`ccl8_mvmr.tsv`、`ccl8_mr_sensitivity.tsv`、`pqtl_ccl8_v2.tsv`
- 单 lead eQTL：MR OR 2.349（p 0.0023），coloc PP4 0.955；FinnGen IBD 复制 OR≈4.37（p 2.9e-6, coloc 0.95)。
- **但多工具削弱**：放宽取 4–5 工具后 IVW 对 IBD p=0.10、UC p=0.50（不显著），异质 + LOO 不稳。
- **MVMR 调 CRP 后不显著**（与泛炎症纠缠）；**SCALLOP-INF 血浆 pQTL(MCP-2) coloc PP4=0.007**（无蛋白层支持）。
- → **写作措辞**：CCL8 = "共定位 + 单队列复制支持的**提示性**风险信号，但在多工具/调炎症/蛋白层面不稳健"。**不要**与 CXCR2 平起平坐当"两个确证因果基因"。

### 路人（bystanders）
来源表：`module_causal_map_multicontext.tsv`、`concordance_map.tsv`
- OSM/OSMR/IL13RA2：无可用工具或无信号；**TREM1**：强 eQTL(p 2.4e-27) 但 MR OR 1.002, p 0.92, PP4≈0 → 明确路人。

---

## 3. 药靶 MR + 试验对账（B1，第 2 头条）
来源表：`drugtarget_mr_predictions.tsv`、`trial_outcomes.tsv`、`concordance_map.tsv`；图：`Fig_drugtarget_concordance.png`
- 框架：cis-MR = "基因型模拟药物抑制靶点"；升高表达→风险 ⇒ 抑制有效；升高表达→保护 ⇒ 抑制碰壁。
- **可严格评判的靶 = 2 个，一致 2/2**（**务必写明 N=2，不可吹"100% 一致图谱"**）：
  - **TNFRSF1A**（anti-TNF 阳性对照）：升高→风险(OR 1.16) ⇒ 抑制有效 ⇒ 对账"已批有效"。**方向性对照**（coloc PP4 0.18 不过严格关，写作如实标注）。
  - **CXCR2**（关键预警）：升高→保护(OR 0.753, PP4 0.951) ⇒ 抑制碰壁 ⇒ 对账 CXCR2 拮抗剂炎症试验 failed/no_efficacy。**这是干净、过 coloc 的核心案例。**
- **多数靶评不了**：已批生物药靶(IL12/23、α4β7 整合素、IL23R)**无 cis-eQTL 工具**；IL1B/IL6/ICAM1/CCL2/MMP9 **MR 有信号但 coloc 全不过**（写作要诚实，这是"遗传学只能做分诊层，做不成完整药理图谱"）。

---

## 4. 机制（B2，第 3 头条）
来源表：`cxcr2_to_bloodtraits_mr.tsv`、`bloodtrait_to_ibd_mr.tsv`、`mediation_summary.tsv`；图：`Fig_cxcr2_mechanism.png`
- CXCR2↑ → 中性粒计数↑：theta +0.052（p 1.2e-10, blood, **coloc PP4 0.998**）；+0.041（p 8e-19, neutrophil）。
- 中性粒计数 → IBD：**OR 1.23 (IBD) / 1.39 (CD) / 1.10 (UC)**（439 工具，p 5.5e-13）——中性粒计数是 IBD **因果风险因子**（独立可引用的干净结果）。
- 中介乘积为正、**与 CXCR2 保护方向相反** → **简单计数中介解释不了悖论**（codex 诚实标 "opposes_protection"）。
- **写作解读**：排除了"CXCR2 靠减少中性粒来保护"；保护应在中性粒**定位/稳态/功能**层面（如留在循环、不浸润肠道；count-MR 量不到）。这与 B1 自洽——朴素"抑制 CXCR2 减少肠道中性粒"逻辑因此失败。**如实写"机制部分解释、未完全解开"，不要硬说已解。**

---

## 5. 描述层 + 诚实阴性（支撑/对照，**不当头条**）
- **多队列应答预测**（`multicohort_auc.tsv`、`inflammation_adjusted.tsv`）：5 队列(GSE12251/16879/23597/73661/92415, n=237) SenMayo/SASP 池化 **AUC 0.78**（单队列 0.67–0.93），方向=高→无应答。**但**：与泛炎症 r≈0.88、**零增量**(LRT p≈0.10)；OSM/IL13RA2/TREM1/中性粒 全在 0.80–0.82 饱和家族。→ **不能当"新预测器"卖**。
- **单细胞反卷积**（`celltype_fraction_vs_response.tsv`、`deconv_incremental_vs_senmayo.tsv`；`Fig_deconv_response.png`）：**髓系丰度**预测无应答(OR/SD 2.65, p 4.8e-07)，但**不比 SASP 多增量**(ΔAUC<0.011, n.s.)。→ 把信号定位到髓系，是机制支撑，不是独立预测力。
- **衰老阴性对照**：停滞臂(p16/p21/CDKN)为零；转录组衰老时钟阴性；LTL→IBD MR 反向/阴性。→ 写作放"诚实阴性"小节增严谨，并解释为何退出衰老框架（见 `PIVOT_2026-06-23_causal_refractory_module.md`）。
- **单细胞衰老分诊**（`senescence_per_celltype.tsv`）：真衰老候选细胞少数(8–11%)、主要髓系 SASP，非全局。

---

## 6. 诚实框架 / 不可夸大（写作红线）
- ❌ 不写"衰老驱动 IBD"——停滞臂为零，信号是炎症-分泌。
- ❌ 不写"新颖应答预测器"——饱和、与炎症冗余、无增量。
- ❌ 不写"多基因因果/药靶图谱"——只有 CXCR2 干净；阳性对照 coloc 还弱。
- ❌ 不把 CCL8 与 CXCR2 并列为"两个确证因果基因"——CCL8 是提示性。
- ❌ "药靶一致 2/2 = 100%" 必须带 **N=2** 限定。
- ❌ 机制不说"已解开"——说"排除了数量解释、指向定位/稳态、部分解释"。
- ✅ 卖点是：**因果分诊（debunk marker 因果）+ CXCR2 可成药悖论（trial-concordant 的临床警示）+ 中性粒计数→IBD 因果 + 诚实阴性弧**。

---

## 7. 建议论文结构（写作智能体可据此搭）
- **定位**：整合型 translational，IF 3–6（如 *J Crohns Colitis / Gut(短) / J Transl Med / Clin Transl Gastroenterol / Cell Rep Med(冲)*）。
- **Title 方向**：如 *"Genetic causal triage of the refractory-IBD inflammatory module: CXCR2 is protective and its blockade is genetically contraindicated"*。
- **Figures**（已就绪，在 `results/figures/`）：
  - Fig1 = 难治模块描述 + 应答预测饱和家族（`Fig_deconv_response.png` / `multicohort_auc.tsv` 作图）。
  - Fig2 = 多语境因果地图 + 分诊（`Fig_module_causal_multicontext.png`）。
  - Fig3 = CXCR2 加固（FinnGen 复制 + 敏感性 + 反向澄清，`Fig_causal_hardening.png`）。
  - Fig4 = 药靶 MR vs 试验一致性（`Fig_drugtarget_concordance.png`）。
  - Fig5 = CXCR2 机制（`Fig_cxcr2_mechanism.png`）。
- **章节**：背景(难治 IBD 的炎症-分泌 marker 乱象)→ 模块定义+应答预测(描述)→ 多语境因果分诊→ CXCR2 加固→ 药靶 MR+试验对账→ 机制→ 诚实局限(衰老阴性/CCL8 提示/coloc 缺口)。
- **方法学溯源**：脚本 `src/12–31`（MR/coloc/多语境/加固/药靶/机制）；纯 Python 无 R；数据=de Lange 2017 + FinnGen R12 + eQTLGen + GTEx colon + BLUEPRINT 免疫 eQTL + Vuckovic/Astle 血细胞 + SCALLOP pQTL + 5 个 GEO 应答队列 + Smillie/Martin 单细胞。

---

## 8. 待补/局限（写作时如实列入 Limitations）
- 只有 CXCR2 一个稳健因果；CCL8 提示性。anti-TNF 阳性对照 coloc 弱(0.18)。
- 已批生物药靶无 eQTL 工具，药靶对账 N 小。
- 机制悖论未完全解开（缺肠道 vs 循环中性粒的分室数据）。
- 肠 allpairs phenotype-id 映射技术债（肠 coloc 仍近似；不影响结论，但严格版未跑通）。
- CXCR2 跨族裔复制未做（讨论后判定"加 robustness 不加维度"，可列为 future work）。
