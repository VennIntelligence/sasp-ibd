"""Assemble all results into results/RESULTS_SUMMARY.md and a figure contact sheet."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

RES = "/Users/ujs/Downloads/lzy/outputs"


def load(name):
    p = f"{RES}/{name}"
    return json.load(open(p)) if os.path.exists(p) else {}

clock = load("clock_model.json")
act = load("activity_stats.json")
rev = load("reversal_stats.json")
sen = load("senescence_stats.json")
pred = load("predictor_stats.json")

md = []
md.append("# IBD 肠黏膜转录组衰老时钟 — 结果汇总\n")
md.append("> 主题：肠黏膜转录组『生物学年龄』是 IBD 治疗应答的**可逆**生物标志。\n")

md.append("## 1. 衰老时钟构建（GTEx 正常肠道）")
md.append(f"- 时钟基因数：**{clock.get('n_clock_genes','?')}**")
md.append(f"- 交叉验证（受试者分组，无泄漏）：Pearson r = **{clock.get('cv_r',float('nan')):.3f}**, "
          f"MAE = **{clock.get('cv_mae',float('nan')):.2f} 岁**")
md.append(f"- 超参数：alpha={clock.get('alpha')}, l1_ratio={clock.get('l1_ratio')}")
md.append("- 图：`Fig1_clock_accuracy.png`\n")

md.append("## 2. IBD 加速衰老 + 活动度（Fig2）")
if act:
    md.append(f"- GSE16879 基线：UC vs 对照 p={act.get('16879_UC_vs_ctrl_p',float('nan')):.2e}; "
              f"CD vs 对照 p={act.get('16879_CD_vs_ctrl_p',float('nan')):.2e}")
    md.append(f"- GSE73661 基线：UC vs 对照 p={act.get('73661_UC_vs_ctrl_p',float('nan')):.2e}")
    md.append(f"- 活动度梯度：年龄加速 vs Mayo 内镜评分 Spearman rho="
              f"{act.get('73661_ageaccel_vs_mayo_rho',float('nan')):.2f} "
              f"(p={act.get('73661_ageaccel_vs_mayo_p',float('nan')):.2e})")
md.append("- 图：`Fig2_accelerated_aging.png`\n")

md.append("## 3. 可逆性（核心，Fig3）")
if rev:
    md.append(f"- GSE16879 H1 基线 IBD>对照：p={rev.get('16879_H1_baseline_vs_control_p',float('nan')):.2e}")
    md.append(f"- 应答者(R) 配对前后：p={rev.get('16879_H2_R_paired_p',float('nan')):.3f}, "
              f"Δ中位={rev.get('16879_H2_R_delta_median',float('nan')):.2f}")
    md.append(f"- 无应答者(NR) 配对前后：p={rev.get('16879_H2_NR_paired_p',float('nan')):.3f}, "
              f"Δ中位={rev.get('16879_H2_NR_delta_median',float('nan')):.2f}")
    md.append(f"- ΔR vs ΔNR：p={rev.get('16879_H2_deltaRvsNR_p',float('nan')):.3f}")
    md.append(f"- 验证队列 GSE73661 R 配对：p={rev.get('73661_H2_R_paired_p',float('nan')):.3f}")
md.append("- 图：`Fig3_reversal.png`\n")

md.append("## 4. 衰老/SASP 机制（Fig4）")
if sen:
    md.append(f"- GTEx 正常肠道 SenMayo vs 年龄：rho={sen.get('GTEx_SenMayo_vs_age_rho',float('nan')):.2f} "
              f"(p={sen.get('GTEx_SenMayo_vs_age_p',float('nan')):.2e})")
    md.append(f"- 时钟年龄加速 vs SenMayo：GSE16879 rho={sen.get('GSE16879_clock_vs_senmayo_rho',float('nan')):.2f}, "
              f"GSE73661 rho={sen.get('GSE73661_clock_vs_senmayo_rho',float('nan')):.2f}")
    if "top_down_SASP_in_responders" in sen:
        genes = list(sen["top_down_SASP_in_responders"].keys())[:8]
        md.append(f"- 应答者治疗后下调的主要 SASP 基因：{', '.join(genes)}")
md.append("- 图：`Fig4_senescence.png`\n")

md.append("## 5. 预测应答（Fig5）")
if pred:
    for gse in ["GSE16879", "GSE73661"]:
        a = pred.get(f"{gse}_AUC_ageaccel"); s = pred.get(f"{gse}_AUC_senmayo")
        c = pred.get(f"{gse}_AUC_combined_LOOCV")
        if a is not None:
            md.append(f"- {gse}: 基线年龄加速 AUC={a:.2f}, SenMayo AUC={s:.2f}, 组合(LOO) AUC={c:.2f}")
md.append("- 图：`Fig5_predictor.png`\n")

with open(f"{RES}/RESULTS_SUMMARY.md", "w") as f:
    f.write("\n".join(md))
print("wrote results/RESULTS_SUMMARY.md")

# contact sheet
figs = ["Fig1_clock_accuracy.png", "Fig2_accelerated_aging.png", "Fig3_reversal.png",
        "Fig4_senescence.png", "Fig5_predictor.png"]
figs = [f for f in figs if os.path.exists(f"{RES}/{f}")]
if figs:
    fig, axes = plt.subplots(len(figs), 1, figsize=(13, 4.3 * len(figs)))
    if len(figs) == 1:
        axes = [axes]
    for ax, f in zip(axes, figs):
        ax.imshow(mpimg.imread(f"{RES}/{f}")); ax.axis("off"); ax.set_title(f, fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{RES}/ALL_FIGURES.png", dpi=120)
    print(f"wrote results/ALL_FIGURES.png ({len(figs)} panels)")
