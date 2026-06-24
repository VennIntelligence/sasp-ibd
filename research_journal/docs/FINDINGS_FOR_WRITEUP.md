# 写作交接文档 · FINDINGS FOR WRITE-UP（2026-06-23 定稿；同日重绘/分层更新）

> **（2026-06-23 三次更新：locus 特异性已结论——CXCR2 不可分辨于 CXCR1，主角改为 CXCL8/IL-8 受体 CXCR1/CXCR2 轴）**
> 因果信号本身不变（coloc ~0.95、MR OR 0.75、FinnGen IBD/UC 复制、中性粒 coloc 0.93）——只是从基因级降为**轴级**。
> 药靶预警因此**更干净**：navarixin/reparixin 本就是 **CXCR1/CXCR2 双拮抗剂**，与"遗传不看好阻断 IL-8 受体轴"直接对账。

> 给写作智能体：**本文件自包含**，写论文只需读它 + 它引用的 `results/tables/*` 与 `results/figures/*`，
> 不需要聊天记录。所有数字均来自已落盘的结果表（括号内给了文件名）。
> **务必照"诚实框架/不可夸大"一节执行**——本项目的价值正在于诚实的因果分诊，夸大会被审稿当场打掉。

---

## 0. 一句话主线（HEADLINE）

> **难治性 IBD 由一个中性粒/髓系炎症-分泌程序标记；对该程序系统做 MR+coloc（血液/肠/免疫细胞多语境 eQTL）后，
> 众多"明星 marker"大多被过滤为遗传学沉默、LD 可疑或非因果路人；唯一可作为正文主角的严格 actionable 锚点是
> 反直觉保护性的 **CXCL8/IL-8 受体（CXCR1/CXCR2）轴**，由此给出 IL-8 受体轴拮抗剂在 IBD 中不宜乐观的遗传预警。**
>
> 该轴**稳健因果 + 保护**（血液 coloc 0.951、MR OR 0.753、FinnGen IBD/UC 复制、中性粒 coloc 0.935）；
> 但**不可解析到单一基因**——血液、中性粒（细胞类型匹配）、条件分析三关一致显示 **CXCR1≈CXCR2**，
> 共读同一 IBD 主变异 **rs62183956（chr2:219,046,122, hg19）**，CXCR2 与 CXCR1 互为共线、遗传上不可分辨（见第 2 节）。
> 边界照旧：无 IBD 三期试验，预警措辞为"遗传不看好"而非"已证实失败"。

**论文不是**：衰老论文、新预测器论文、或多基因药靶图谱（见第 5 节）。

---

## 1. 头条贡献（按重要性；**药靶对账已从第 2 头条降为讨论段**，见第 3 节）

1. **因果分诊（去伪存真）**：用 MR+coloc 把难治炎症模块的"相关 marker"分成严格 actionable、方向性对照、
   提示性信号、LD 可疑、遗传学沉默/无工具和路人。**正文主角是 IL-8 受体轴（CXCR1/CXCR2）**（locus 测定不可解析到单一基因，见第 2 节）。
   negative triage 的两类**务必分开**：**TREM1 = tested-null**（中性粒 eQTL p=2.4e-27 极强，但 MR OR=1.002, p=0.92, coloc PP4≈0，干净证伪范例）；
   **OSM/OSMR/IL13RA2 = 遗传学沉默/无 cis 工具**（**不可评判，≠ 非因果**；OSM 尚有 anti-OSM Crohn 试验在跑，措辞别写"坐实非因果"）。
2. **机制（中性粒计数悖论，*真·第 2 头条*）**：IL-8 受体轴↑ → 循环中性粒↑（coloc PP4 0.998）；中性粒计数本身是 IBD 因果风险因子(OR 1.23)；**可该轴仍保护** →
   保护**不是靠减少中性粒**，而是中性粒的**定位/稳态**层面（朴素"抑制 IL-8 受体减少肠道中性粒"的逻辑因此是错的）。详见第 4 节。
3. **IL-8 受体轴可成药推论（translational coda，*非正文头条*）**：该轴是药靶（navarixin/reparixin 等**双拮抗剂**），遗传上保护 ⇒ 抑制"遗传上不被看好"；
   与**非 IBD** 炎症适应症拮抗剂"结局不鼓舞"对账（**无 IBD 三期，不可写"试验已证实失败"**）；anti-TNF(TNFRSF1A) 仅作*方向性*阳性对照(coloc 弱)。详见第 3 节。
4. **诚实阴性弧（增严谨）**：衰老停滞臂为零、转录组衰老时钟阴性、LTL MR 反向/阴性；SASP/"预测器"只是泛炎症的别名（见第 5 节）。

---

## 2. 因果核心（精确数字）

### IL-8 受体轴（CXCR1/CXCR2）—— 稳健因果保护（论文主角；轴级，非单基因）
来源表：`results/tables/module_causal_map_multicontext.tsv`、`finngen_cxcr2.tsv`、`mr_sensitivity.tsv`、`mvmr_crp.tsv`、`reverse_mr_proper.tsv`、`steiger.tsv`、`coloc_IBD.tsv`；locus 测定见 `src/45_cxcr2_locus_specificity.py`、`src/46_cxcr2_neutrophil_conditional.py`、`outputs/cxcr2_locus_specificity/SUMMARY.md` + `SUMMARY_conditional.md`
> 下列 MR/coloc/FinnGen/中性粒数字现读作**轴级**信号（信号经 CXCR2 探针测得，但 locus 测定证明不能归到单一基因，见末条）。
- 血液 eQTLGen：MR **OR 0.753**（p 1.6e-09, FDR 1.24e-07），**coloc PP4 0.951**。
- 中性粒 BLUEPRINT(QTD000026)：MR **OR 0.851**（FDR 1.05e-07），**coloc PP4 0.935**（细胞类型匹配的独立确认）。
- 独立队列 **FinnGen R12 复制**：IBD OR 0.796（p 1.3e-05）、UC OR 0.736（p 1.6e-06）；CD 不显著(0.968, p 0.77)。
- 多工具敏感性（中性粒, 5 SNP）：IVW 0.842 / MR-Egger 0.841（截距无方向性多效性）/ weighted-median 0.846 / LOO ~0.84。
- MVMR 调 CRP：方向稳（IBD OR 0.749, p 1.1e-10）；注：CRP 条件 F 弱(~2.1)。
- 反向因果**已排除**：排 cis 区后 IBD→轴反向无工具/无信号（先前 p=3e-200 是共定位 cis 变异伪反向）。Steiger 支持 exposure→outcome。
- **✅ locus 特异性已结论 — 信号是 IL-8 受体轴、不可解析到单一基因**（`src/45`、`src/46`；曾是首要硬化项，现已测定）：
  - **血液 eQTLGen coloc 扫描**：CXCR2 PP4_IBD **0.951** / CXCR1(IL8RA) **0.889** / GPBAR1 **0.677**，三者读**同一** IBD 主变异 **rs62183956 @ chr2:219,046,122 (hg19)** → 非 CXCR2 专属。
  - **中性粒（细胞类型匹配, BLUEPRINT QTD000026）**：CXCR1 PP4 **0.937** ≈ CXCR2 **0.935**（阳性对照复现），**同一 eQTL lead rs6737563**；GPBAR1 0 SNP 出局。匹配细胞类型**仍不能区分**二者（CXCR1 若有差别反略胜）。
  - **条件分析（eQTLGen 血液 + 1000G EUR LD, 单 SNP COJO）**：在 IBD 变异处，CXCR2 调 CXCR1 后其 eQTL 坍塌（p 3.2e-200 → 0.38），CXCR1 调 CXCR2 后亦坍塌（p 0 → 0.98）→ **互为共线，遗传上不可解析。**
  - **结论**：报告为 **CXCR1/CXCR2 轴**；GPBAR1 在中性粒出局，把单元干净收窄到两个 IL-8 受体。基因级归因需功能/等位特异表达(ASE)/MPRA（湿实验，超本计算论文范围）。**（"单因果变量假设"的待定限定已删除——已解析，结论即"不是单一基因"。）**

### CCL8 —— 提示性风险（**降级，需谨慎措辞**）
来源表：`triangulation.tsv`、`ccl8_mvmr.tsv`、`ccl8_mr_sensitivity.tsv`、`pqtl_ccl8_v2.tsv`
- 单 lead eQTL：MR OR 2.349（p 0.0023），coloc PP4 0.955；FinnGen IBD 复制 OR≈4.37（p 2.9e-6, coloc 0.95)。
- **但多工具削弱**：放宽取 4–5 工具后 IVW 对 IBD p=0.10、UC p=0.50（不显著），异质 + LOO 不稳。
- **MVMR 调 CRP 后不显著**（与泛炎症纠缠）；**SCALLOP-INF 血浆 pQTL(MCP-2) coloc PP4=0.007**（无蛋白层支持）。
- → **写作措辞**：CCL8 = "共定位 + 单队列复制支持的**提示性**风险信号，但在多工具/调炎症/蛋白层面不稳健"。
  正文只作降级案例或补充结果，**不要**与 IL-8 受体轴平起平坐当"确证因果基因"。

### 路人 / 不可评判（triage 两类，**务必分开写**）
来源表：`module_causal_map_multicontext.tsv`、`concordance_map.tsv`、`drugtarget_evidence_tiers.tsv`
- **证伪（tested-null）**：**TREM1**——强 eQTL(p 2.4e-27) 但 MR OR 1.002, p 0.92, PP4≈0 → 明确非因果路人（干净 debunk 范例）。
- **遗传学沉默（no instrument，不可评判）**：OSM/OSMR/IL13RA2——无可用 cis 工具/无信号 → **只能写"genetics is silent"，绝不能写"非因果"**（OSM 有 anti-OSM Crohn 试验在跑，硬说非因果会被掀桌）。

---

## 3. 药靶 MR + 试验对账（**讨论段 / translational coda —— 非正文第 2 头条**）
来源表：`drugtarget_evidence_tiers.tsv`、`drugtarget_mr_predictions.tsv`、`trial_outcomes.tsv`、`concordance_map.tsv`；图：`Fig4_drugtarget_concordance.png`
- 框架：cis-MR = "基因型模拟药物抑制靶点"；升高表达→风险 ⇒ 抑制有效；升高表达→保护 ⇒ 抑制碰壁。
- **可严格评判的靶 = 2 个，一致 2/2**（**务必写明 N=2，不可吹"100% 一致图谱"**）：
  - **TNFRSF1A**（anti-TNF 阳性对照）：升高→风险(OR 1.16) ⇒ 抑制有效 ⇒ 对账"已批有效"。**方向性对照**（coloc PP4 0.18 不过严格关，写作如实标注）。
  - **IL-8 受体轴（CXCR1/CXCR2）**（关键预警）：升高→保护(OR 0.753, PP4 0.951) ⇒ 抑制遗传上不被看好 ⇒ 对账 IL-8 受体拮抗剂**非 IBD** 炎症适应症（AZD5069 支气管扩张、danirixin COPD）结局不鼓舞。**临床在研的 navarixin / reparixin 本就是双 CXCR1/CXCR2 拮抗剂——预警因此为轴级、且与药物机制直接对账，更干净。**过 coloc；务必注明无 IBD 三期试验，不可写"已证实失败"。
- **多数靶评不了**：已批生物药靶(IL12/23、α4β7 整合素、IL23R)**无 cis-eQTL 工具**；IL1B/IL6/ICAM1/CCL2/MMP9 **MR 有信号但 coloc 全不过**（写作要诚实，这是"遗传学只能做分诊层，做不成完整药理图谱"）。
- `drugtarget_evidence_tiers.tsv` 是写作主表（**2026-06-24 扩展至 27 个靶点**）：
  - **strict_actionable_warning**: CXCR2（正文读作 CXCR1/CXCR2 轴级预警）
  - **directional_positive_control**: TNFRSF1A
  - **downgraded_hint_supplementary**: CCL8
  - **genetics_silent_approved_pathway (14 个)**: JAK1/JAK3/TYK2/S1PR1/MADCAM1/PDE4B/PDE4D/IL6R/TNFSF15（新增，血液 eQTLGen 已测，coloc PP4≈0）+ IL12B/IL23A/IL23R/ITGA4/ITGB7（原有，无工具）
  - **ld_suspect_not_actionable**: CCL2/ICAM1/IL1B/IL6/MMP9 + **JAK2**（新增，MR OR 0.665 p≈0，F=638，但 coloc PP4=0 → LD 存疑）
  - **marker_bystander_no_tool**: OSM/OSMR/IL13RA2
  - **tested_bystander_null**: TREM1
- **新增 10 靶（eQTLGen 血液 cis-MR benchmark #1）**（`outputs/48_drugtarget_expand/expanded_targets_blood.tsv`，脚本 `src/48`）：全部 coloc PP4≈0；JAK2 有强 MR 但 coloc 失败（LD-suspect）；TNFSF15（TL1A）在 eQTLGen 中缺失（genetics-silent）。**写作结论**：遗传学对 27 个靶做了系统筛查，只有 CXCR1/CXCR2 轴通过严格 PP4+MR 门槛；JAK 类/S1P 调节剂/MAdCAM1/PDE4/TL1A 均为遗传学沉默——这反而强化"系统筛查结果就是 CXCR1/CXCR2 是唯一遗传学锚点"的叙事（不要说"遗传学反对 JAK 疗法"，而是"遗传学对这些靶沉默，临床证据才是这些靶的依据"）。

---

## 3b. PheWAS 安全性地图（CXCR1/CXCR2 轴两个 lead 变异）
来源：`outputs/49_phewas_axis/`；脚本 `src/49_phewas_axis.py`；图：`results/figures/Fig_phewas_blood_lead_rs62183956.png`；表：`results/tables/phewas_axis_significant.tsv`
- **FinnGen R12 phenome-wide 扫描**（2470 个疾病终点）：血液 lead rs62183956（2-218181399-C-T）+ 中性粒 lead rs6737563（2-218080951-T-C）两个变异独立确认（两者 IBD 锚定方向正确，p 8.3e-6 / 1.4e-5）。
- **Bonferroni 显著（两变异一致）= 6 个终点**：
  - IBD/UC 保护（4 个，oriented beta < 0）：K11_IBD_STRICT、K11_UC_STRICT2、ULCERNAS、MUCOPROCT——**预期内，确认方向**
  - 胆囊结石/胆囊切除（2 个，oriented beta > 0，**p 5.8e-10 / 5.1e-10**，最强信号！）：K11_CHOLELITH + K11_CHOLECYSTECTOMY——**意外的肝胆信号，但生物学合理**（胆道上皮富含 CXCL8；中性粒介导的胆道炎症是胆石症机制；IBD 保护性等位基因(更高受体活性)→ 略升胆石风险 OR~1.04，不是大风险，但显著）
- **无感染/呼吸道/血液免疫疾病终点通过 Bonferroni**：FinnGen 疾病表型谱中，该 locus 的多效性很窄（非广泛免疫枢纽）；但这**不排除**中性粒抗菌功能的生物学作用——FinnGen 疾病终点捕捉不到中性粒定量特征（后者在 GWAS Catalog 血细胞定量 GWAS 中，超出本 PheWAS 范围）。
- **写作策略**：
  - ✅ 可写：该 locus 在 FinnGen 中多效性窄，主要信号是 IBD 保护 + 肝胆风险；无感染/呼吸枢纽（在疾病表型层面）
  - ✅ 可写：胆囊结石信号（beta +0.043，Bonferroni p<1e-9）表明 CXCL8 轴在胆道亦有生物活性，与已知中性粒-胆道炎症生物学一致
  - ❌ 不写"axis is not an immune hub"——中性粒功能枢纽有丰富机制证据；PheWAS 说的是**这个 locus 在 FinnGen 疾病表型中多效性窄**，不是"axis 无免疫功能"
  - ❌ 不说"PheWAS 表明拮抗剂安全"——胆石信号恰好是潜在安全顾虑

## 4. 机制（B2，第 2 头条 —— 中性粒计数悖论）
来源表：`cxcr2_to_bloodtraits_mr.tsv`、`bloodtrait_to_ibd_mr.tsv`、`mediation_summary.tsv`；图：`Fig_cxcr2_mechanism.png`
- IL-8 受体轴↑ → 中性粒计数↑：theta +0.052（p 1.2e-10, blood, **coloc PP4 0.998**）；+0.041（p 8e-19, neutrophil）。
- 中性粒计数 → IBD：**OR 1.23 (IBD) / 1.39 (CD) / 1.10 (UC)**（439 工具，p 5.5e-13）——中性粒计数是 IBD **因果风险因子**（独立可引用的干净结果）。
- 中介乘积为正、**与该轴保护方向相反** → **简单计数中介解释不了悖论**（codex 诚实标 "opposes_protection"）。
- **写作解读**：排除了"IL-8 受体轴靠减少中性粒来保护"；保护应在中性粒**定位/稳态/功能**层面（如留在循环、不浸润肠道；count-MR 量不到）。这与 B1 自洽——朴素"抑制 IL-8 受体减少肠道中性粒"逻辑因此失败。**如实写"机制部分解释、未完全解开"，不要硬说已解。**

### 单细胞机制定位（Smillie UC atlas SCP259 Imm compartment，2026-06-24）
来源：`outputs/51_scrnaseq_cxcr2/`；脚本 `src/51_scrnaseq_cxcr2_mechanism.py`；图：`outputs/Fig_cxcr2_dotplot.png`、`outputs/Fig_cxcr2_health_bar.png`（待促推至 results/）
- **数据**：365,493 个细胞（Imm 共 210,614），三种 Health 状态（Healthy/Non-inflamed/Inflamed），使用 SCP259 基因矩阵（IL8 = CXCL8，旧 HGNC 名）。
- **CXCL8/IL8 配体来源（最强信号）**：
  - 炎症单核细胞(Inflammatory Monocytes)：mean 0.892，**30.4%** 细胞表达（n=906）
  - 巡回单核细胞(Cycling Monocytes)：mean 1.027，**42.3%**（n=78）
  - 巨噬细胞(Macrophages)：mean 0.855，**29.0%**（n=7,162）
  - DC2：mean 0.807，**32.7%**；DC1：mean 0.585，**29.9%**
  - **结论：CXCL8 由黏膜髓系细胞（单核/巨噬/DC）产生，这是 IBD 炎症最主要的 IL-8 来源**
- **CXCR1/CXCR2 受体**：在所有捕获的免疫细胞中表达率均 <2%（接近零）——**符合预期**：受体主要在中性粒上，中性粒因脆性在 scRNA-seq 中缺失；NKs 和 DC2 有极微量表达（CXCR2 pct 2%/1.4%）可能是噪声。
- **分层（Health 状态）**：myeloid 单核/巨噬细胞在各 Health 状态均有相近水平的 IL8 表达（单核 42–62% pct across states）——**CXCL8 产生是 constitutive（组成性），非仅在炎症期上调**；与稳态黏膜免疫监控功能（homeostatic surveillance）一致。
- **写作整合**：单细胞数据提供 ligand-side 证据：myeloid → CXCL8 → CXCR2 on neutrophils → neutrophil recruitment。遗传学提供 receptor-side 证据（血液 eQTL）。二者拼出完整回路，且 CXCL8 的 constitutive myeloid 产生支持"阻断受体会干扰正常黏膜免疫稳态"的预警逻辑。
- **注意**：中性粒在 scRNA-seq 缺失，因此无法直接展示 CXCR2+ 细胞；图表展示的是 ligand-source 证据，receptor-side 来自 eQTL 遗传学（已 coloc）。写作中如实说明这一互补关系。
- ⚠️ **审稿风险（主动拆弹）**：MR 方向与其自身中介相反，正是审稿人识别"水平多效性/工具无效"的标准签名。务必先用第 2 节的 locus 特异性（已结论为轴级）+ Steiger + MR-Egger 截距堵死"伪信号"解释，再讲"反直觉生物学"。

---

## 5. 描述层 + 诚实阴性（支撑/对照，**不当头条**）
- **多队列应答预测**（`multicohort_auc.tsv`、`inflammation_adjusted.tsv`）：5 队列(GSE12251/16879/23597/73661/92415, n=237) SenMayo/SASP 池化 **AUC 0.78**（单队列 0.67–0.93），方向=高→无应答。**但**：与泛炎症 r≈0.88、**零增量**(LRT p≈0.10)；OSM/IL13RA2/TREM1/中性粒 全在 0.80–0.82 饱和家族。→ **不能当"新预测器"卖**。
- **单细胞反卷积**（`celltype_fraction_vs_response.tsv`、`deconv_incremental_vs_senmayo.tsv`；`Fig_deconv_response.png`）：**髓系丰度**预测无应答(OR/SD 2.65, p 4.8e-07)，但**不比 SASP 多增量**(ΔAUC<0.011, n.s.)。→ 把信号定位到髓系，是机制支撑，不是独立预测力。
- **衰老阴性对照**：停滞臂(p16/p21/CDKN)为零；转录组衰老时钟阴性；LTL→IBD MR 反向/阴性。→ 写作放"诚实阴性"小节增严谨，并解释为何退出衰老框架（见 `PIVOT_2026-06-23_causal_refractory_module.md`）。
- **单细胞衰老分诊**（`senescence_per_celltype.tsv`、`FigS6`）：严格 bona fide 衰老比例**上限仅 8–11%**，集中在上皮分泌(M 细胞/杯状/未成熟肠)、内皮、成纤维，髓系约 7%，淋巴/增殖最低 → **衰老是局灶 ~10% 现象、非全局驱动**（"诚实退出衰老框架"的实证地基）。注：SASP **程序**表达由炎症成纤维/基质主导(score 最高 ~0.29)，与"严格衰老比例"是两件事，别混。
- **肠 GTEx eQTL 无严格锚点 → 反向写成腔室特异性**（`FigS3_gut_eqtl_null`）：IL-8 受体轴因果信号在**血液/免疫腔室**(中性粒 BLUEPRINT coloc 0.935)、不在 bulk 肠上皮——与中性粒/髓系机制自洽，写成"信号在免疫腔室而非上皮"的特异性，不是纯缺陷。
- **单细胞 foundation-model 扰动（Geneformer）= 正交验证、*已干净重跑、未复现 → 诚实弃用***（`perturb_consensus.tsv`、`insilico_perturbation.tsv`）：V1-10M 下 CXCR2 过表达把细胞推**离**炎症中心（martin −0.0157、smillie ≈0、方向与 MR 一致）；但 **V2-104M（104M 大模型）干净重跑后方向反转**（martin +0.0273、smillie +0.0448、**朝向**炎症），`perturb_consensus.tsv` 中 (CXCR2, overexpress) 的 `projection_consistent_significant`=**False**。→ 这是**真·跨模型不一致**（已排除前次 GPU 故障/空面板伪阴性；四分类器 AUC 0.54–0.71 均 >0.5 排除"炎症向量无意义"解释）；按预注册裁决规则(`tmp/gpu_rerun.md`)**诚实弃用、不作正文证据**，最多 Limitations 一句。主线 CXCR1/CXCR2 保护依赖 MR+coloc+中性粒机制、**不依赖此项**。

---

## 6. 诚实框架 / 不可夸大（写作红线）
- ❌ 不写"衰老驱动 IBD"——停滞臂为零，信号是炎症-分泌。
- ❌ 不写"新颖应答预测器"——饱和、与炎症冗余、无增量。
- ❌ 不写"多基因因果/药靶图谱"——只有 IL-8 受体轴（CXCR1/CXCR2）干净；阳性对照 coloc 还弱。
- ❌ 不把 CCL8 与 IL-8 受体轴并列为"确证因果基因"——CCL8 是提示性。
- ❌ "药靶一致 2/2 = 100%" 必须带 **N=2** 限定，且药靶对账是**讨论段**不是头条。
- ❌ **不把 OSM/OSMR/IL13RA2 写成"非因果"**——它们是"无工具/遗传学沉默"；只有 TREM1 是真·tested-null。
- ❌ **不写"IL-8 受体拮抗剂 IBD 试验已失败"**——无 IBD 试验；写"遗传上不被看好 + 非 IBD 拮抗结局不鼓舞"。
- ❌ 不写 CXCR2 是单一因果基因——locus 已测定不可分辨于 CXCR1，写 **CXCL8/IL-8 受体（CXCR1/CXCR2）轴**。
- ❌ 机制不说"已解开"——说"排除了数量解释、指向定位/稳态、部分解释"。
- ✅ 卖点是：**negative causal triage（相关 marker 不等于因果靶点）+ IL-8 受体轴（CXCR1/CXCR2）可成药悖论（遗传预警）+ 中性粒计数→IBD 因果 + 诚实阴性弧**。

---

## 7. 建议论文结构（写作智能体可据此搭）
- **定位**：整合型 translational，IF 3–6（如 *J Crohns Colitis / Gut(短) / J Transl Med / Clin Transl Gastroenterol / Cell Rep Med(冲)*）。
- **Title 方向**：如 *"Genetic causal triage of the refractory-IBD inflammatory module: a protective CXCL8/IL-8-receptor (CXCR1/CXCR2) axis whose blockade is genetically cautioned"*（**用 "cautioned / genetically disfavored"，不用 "contraindicated"**——无 IBD 试验，措辞不可过强；主角为轴级、非单基因 CXCR2）。
- **Figures**（已就绪并对齐到轴级，在 `results/figures/`；**6 主图 Fig1–6 + 7 补充图 FigS1–S7**，每张含 `.pdf/.png/.svg`，CNS 顶刊风格走 `src/figstyle.py`）：
  - Fig1 = 总览图：预测信号降级（ΔAUC 标 n.s.）+ 因果分诊 + IL-8 受体轴保护 + blockade warning（`Fig1_causal_overview.*`）。
  - Fig2 = 详细多语境因果地图 + 分诊审计（**CXCR2=轴级实线框、CCL8=降级虚线框**，不平起平坐）（`Fig2_causal_triage.*`）。
  - Fig3 = IL-8 受体轴加固：FinnGen 复制 + 多工具敏感性 + MVMR×CRP + 反向澄清（标题已改轴级、不写"causal gene"）（`Fig3_cxcr2_hardening.*`）。
  - **Fig4 = locus 特异性（新·keystone，堵"单基因归因"）：coloc 扫描三基因共读 rs62183956 + 条件互坍塌(3.2e-200→0.38 / 0→0.98) + 中性粒匹配 trio(CXCR1 0.937≈CXCR2 0.935) → CXCR1≈CXCR2 遗传不可解析、报告为轴级**（`Fig4_locus_specificity.*`；脚本 `src/47`、数据 `outputs/cxcr2_locus_specificity/*.tsv`）。
  - Fig5 = 药靶 MR vs 证据分层/试验对账（轴级拮抗剂 navarixin/reparixin 预警、**N=2**、措辞"genetically disfavoured"非"trials failed"）（`Fig5_drugtarget_concordance.*` + `drugtarget_evidence_tiers.tsv`）。
  - Fig6 = IL-8 受体轴机制：中性粒计数路径与保护方向相反（`Fig6_mechanism.*`）。
  - 补充层 = `FigS1` coloc PP3/PP4、`FigS2` CCL8 downgrade、`FigS3` gut eQTL null（免疫腔室特异性）、`FigS4` deconvolution response、`FigS5` neutrophil IBD MR、`FigS6` celltype senescence、`FigS7` biomarker deflation。
- **章节**：背景(难治 IBD 的炎症-分泌 marker 乱象)→ 模块定义+应答预测(描述)→ 多语境因果分诊→ IL-8 受体轴加固（含 locus 特异性=轴级）→ 药靶 MR+试验对账→ 机制→ 诚实局限(衰老阴性/CCL8 提示/coloc 缺口/基因级不可解析)。
- **方法学溯源**：脚本 `src/12–31` 主流水线 + `src/45–46`（locus 特异性/中性粒条件分析）（MR/coloc/多语境/加固/药靶/机制）；纯 Python 无 R；数据=de Lange 2017 + FinnGen R12 + eQTLGen + GTEx colon + BLUEPRINT 免疫/中性粒 eQTL(QTD000026) + 1000G EUR LD + Vuckovic/Astle 血细胞 + SCALLOP pQTL + 5 个 GEO 应答队列 + Smillie/Martin 单细胞。

---

## 8. 待补/局限（写作时如实列入 Limitations）
- **locus 特异性已结论（首要硬化项，已测定）**：2q35 簇内 CXCR2 与高同源 CXCR1 **遗传上不可解析**——血液(CXCR2 0.951/CXCR1 0.889) + 中性粒(CXCR1 0.937≈CXCR2 0.935，同 lead) + 条件分析(互相条件后双双坍塌) 一致；因此主结论报告为 **CXCL8/IL-8 受体（CXCR1/CXCR2）轴**，**基因级归因从当前遗传学无法实现**，需功能/等位特异表达(ASE)/MPRA 等湿实验（future work，超本计算论文范围）。
- 稳健因果是 IL-8 受体轴（轴级、非单基因）；CCL8 提示性。anti-TNF 阳性对照 coloc 弱(0.18) 且 TNF/TNFRSF1A 方向自相矛盾。
- 药靶-试验对账 N=2，且"碰壁"证据来自**非 IBD** 适应症——已降为讨论段；已批生物药靶无 eQTL 工具。
- 机制悖论未完全解开（缺肠道 vs 循环中性粒的分室数据）。
- 肠 GTEx 无严格因果锚点（已正面重构为免疫腔室特异性，见第 5 节）；allpairs phenotype-id 映射技术债（肠 coloc 仍近似；不影响结论）。
- **Geneformer 扰动加固已完成、未复现 → 诚实弃用**：V1-10M 方向支持 CXCR2 保护，但 V2-104M 干净重跑后方向反转（朝向炎症）、跨模型不一致(`projection_consistent_significant`=False)；in-silico 扰动证据**模型版本依赖、不稳健**，故不入正文，仅此一句列入 Limitations。主线结论不依赖此项。
- IL-8 受体轴跨族裔复制未做（讨论后判定"加 robustness 不加维度"，可列为 future work）。
