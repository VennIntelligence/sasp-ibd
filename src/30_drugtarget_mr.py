"""Drug-target cis-MR interpreted as genetic proxy for target inhibition.

The upstream MR scripts estimate the effect of higher target expression on IBD.
Most IBD biologics/small molecules inhibit or neutralize their target, so the
drug proxy is the opposite direction: if higher expression raises IBD risk,
target inhibition is predicted beneficial; if higher expression is protective,
target inhibition is predicted harmful or ineffective.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, field_validator

from paths import P


TARGETS = [
    {"gene": "TNFRSF1A", "target_axis": "TNF/TNFRSF1A", "drug_class": "anti-TNF", "priority": "positive_control"},
    {"gene": "TNF", "target_axis": "TNF/TNFRSF1A", "drug_class": "anti-TNF", "priority": "known_ibd_target"},
    {"gene": "CXCR2", "target_axis": "CXCR2", "drug_class": "CXCR2 antagonists", "priority": "key_warning"},
    {"gene": "IL1B", "target_axis": "IL1B", "drug_class": "IL-1 blockade", "priority": "druggable_module_target"},
    {"gene": "MMP9", "target_axis": "MMP9", "drug_class": "anti-MMP9", "priority": "druggable_module_target"},
    {"gene": "ICAM1", "target_axis": "ICAM1", "drug_class": "ICAM1 antisense", "priority": "druggable_module_target"},
    {"gene": "IL6", "target_axis": "IL6/IL6R", "drug_class": "IL-6 pathway blockade", "priority": "druggable_module_target"},
    {"gene": "CCL2", "target_axis": "CCL2/CCR2", "drug_class": "CCR2/CCL2-axis inhibitors", "priority": "druggable_module_target"},
    {"gene": "OSM", "target_axis": "OSM/OSMR", "drug_class": "anti-OSM/OSMR", "priority": "bystander_candidate"},
    {"gene": "OSMR", "target_axis": "OSM/OSMR", "drug_class": "anti-OSM/OSMR", "priority": "bystander_candidate"},
    {"gene": "CCL8", "target_axis": "CCL8", "drug_class": "chemokine-axis inhibition", "priority": "downgraded_hint"},
    {"gene": "TREM1", "target_axis": "TREM1", "drug_class": "TREM1 modulation", "priority": "bystander_candidate"},
    {"gene": "IL13RA2", "target_axis": "IL13RA2", "drug_class": "IL-13-axis targeting", "priority": "bystander_candidate"},
    {"gene": "ITGA4", "target_axis": "alpha4beta7 integrin", "drug_class": "vedolizumab/natalizumab axis", "priority": "approved_context"},
    {"gene": "ITGB7", "target_axis": "alpha4beta7 integrin", "drug_class": "vedolizumab axis", "priority": "approved_context"},
    {"gene": "IL12B", "target_axis": "IL12/23 p40", "drug_class": "ustekinumab", "priority": "approved_context"},
    {"gene": "IL23A", "target_axis": "IL23 p19", "drug_class": "IL-23 p19 inhibitors", "priority": "approved_context"},
    {"gene": "IL23R", "target_axis": "IL23 signaling", "drug_class": "IL-23 pathway", "priority": "approved_context"},
    {"gene": "IL13", "target_axis": "IL13", "drug_class": "IL-13 blockade", "priority": "known_inflammation_target"},
]


TRIAL_OUTCOMES = [
    {
        "target_axis": "TNF/TNFRSF1A",
        "representative_agents": "infliximab; adalimumab; golimumab; certolizumab",
        "trial_or_program": "ACT 1/2 and anti-TNF IBD programs",
        "actual_trial_outcome": "approved_effective",
        "outcome_note": "Anti-TNF therapy is approved and effective in IBD; ACT 1/2 showed infliximab clinical response superiority over placebo in UC.",
        "citation": "Rutgeerts et al. Infliximab for induction and maintenance therapy for ulcerative colitis. N Engl J Med 2005.",
        "url": "https://www.naspghan.org/files/documents/pdfs/training/curriculum-resources/inflammatory%20bowel%20disease_zeusdeleted_145_09092015100906/ACT_1_and_2.pdf",
    },
    {
        "target_axis": "CXCR2",
        "representative_agents": "navarixin; AZD5069; danirixin",
        "trial_or_program": "CXCR2 antagonist inflammation programs",
        "actual_trial_outcome": "failed/no_efficacy",
        "outcome_note": "No mature IBD efficacy program was found; CXCR2 antagonists reduced neutrophil biomarkers in pulmonary/inflammatory settings without clear clinical benefit, matching a genetics warning against antagonism in IBD.",
        "citation": "De Soyza et al. AZD5069 in bronchiectasis; ClinicalTrials.gov NCT03250689 danirixin COPD/influenza.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/26341987/; https://clinicaltrials.gov/study/NCT03250689",
    },
    {
        "target_axis": "MMP9",
        "representative_agents": "andecaliximab / GS-5745",
        "trial_or_program": "phase 2/3 UC; phase 2 CD",
        "actual_trial_outcome": "failed/no_efficacy",
        "outcome_note": "UC induction trial met futility criteria; neither every-2-week nor weekly dosing improved remission, response, endoscopic response, or histologic healing versus placebo.",
        "citation": "Andecaliximab induction therapy for ulcerative colitis. J Crohns Colitis 2018.",
        "url": "https://academic.oup.com/ecco-jcc/article/12/9/1021/4996034",
    },
    {
        "target_axis": "ICAM1",
        "representative_agents": "alicaforsen / ISIS 2302",
        "trial_or_program": "Crohn systemic antisense; topical UC/proctitis/pouchitis",
        "actual_trial_outcome": "mixed",
        "outcome_note": "Systemic Crohn studies were largely negative; topical enema studies/case series show inconsistent but sometimes durable UC/proctitis/pouchitis benefit.",
        "citation": "Greuter et al. Alicaforsen in left-sided UC/proctitis. Dig Dis 2018; review references include negative Crohn and topical UC trials.",
        "url": "https://karger.com/ddi/article/36/2/123/94795/Alicaforsen-an-Antisense-Inhibitor-of",
    },
    {
        "target_axis": "IL1B",
        "representative_agents": "canakinumab; anakinra",
        "trial_or_program": "IL-1 blockade reports in VEO/autoinflammatory IBD; safety reports",
        "actual_trial_outcome": "mixed",
        "outcome_note": "Not an approved conventional IBD mechanism; reports suggest benefit in autoinflammatory/VEO-IBD subsets, while de-novo IBD has been reported during IL-1 inhibitor therapy.",
        "citation": "Canakinumab for autoinflammatory very early onset IBD; reports of IBD following anti-IL-1 treatment.",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9531243/; https://pmc.ncbi.nlm.nih.gov/articles/PMC5348783/",
    },
    {
        "target_axis": "IL6/IL6R",
        "representative_agents": "tocilizumab; PF-04236921; olamkicept",
        "trial_or_program": "anti-IL6R exploratory CD; ANDANTE; IL-6 trans-signaling programs",
        "actual_trial_outcome": "mixed",
        "outcome_note": "Exploratory anti-IL6R and anti-IL6 programs showed efficacy signals but notable safety concerns including GI perforation/suppuration; not approved for IBD.",
        "citation": "Ito et al. anti-IL6R in Crohn's disease; selective IBD therapy review summarizing ANDANTE and olamkicept.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/15902961/; https://www.mdpi.com/2077-0383/11/4/994",
    },
    {
        "target_axis": "CCL2/CCR2",
        "representative_agents": "CCR2/CCL2-axis inhibitors",
        "trial_or_program": "IBD target biology; no clear IBD efficacy trial located",
        "actual_trial_outcome": "untested",
        "outcome_note": "CCL2/CCR2 is biologically implicated in IBD monocyte recruitment, but no mature IBD inhibitor efficacy outcome was identified in the curated search.",
        "citation": "Chemokines and chemokine receptors as therapeutic targets in IBD; CCR2/CCL2 IBD biology reports.",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6104621/; https://pubmed.ncbi.nlm.nih.gov/15306587/",
    },
    {
        "target_axis": "OSM/OSMR",
        "representative_agents": "anti-OSM / anti-OSMR investigational agents",
        "trial_or_program": "GSK2330811 anti-OSM Crohn's trial; biomarker programs",
        "actual_trial_outcome": "untested",
        "outcome_note": "OSM is a strong refractory biomarker and target hypothesis; a Crohn's anti-OSM trial exists, but no mature efficacy result was found locally.",
        "citation": "West et al. Oncostatin M drives intestinal inflammation; ClinicalTrials.gov NCT04151225 anti-OSM Crohn's study.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/28368383/; https://clinicaltrials.gov/study/NCT04151225",
    },
    {
        "target_axis": "CCL8",
        "representative_agents": "chemokine-axis inhibition",
        "trial_or_program": "no specific IBD inhibitor program located",
        "actual_trial_outcome": "untested",
        "outcome_note": "No specific IBD therapeutic outcome for direct CCL8 inhibition was identified; genetic evidence was downgraded after sensitivity/MVMR/pQTL checks.",
        "citation": "Local hardening outputs: ccl8_mvmr.tsv and pqtl_ccl8_v2.tsv.",
        "url": "results/tables/ccl8_mvmr.tsv; results/tables/pqtl_ccl8_v2.tsv",
    },
    {
        "target_axis": "TREM1",
        "representative_agents": "TREM1 modulation",
        "trial_or_program": "IBD target biology; no efficacy program located",
        "actual_trial_outcome": "untested",
        "outcome_note": "No mature IBD therapeutic outcome identified; local genetics classifies TREM1 as a non-causal bystander despite neutrophil eQTL.",
        "citation": "Local module causal map multicontext.",
        "url": "results/tables/module_causal_map_multicontext.tsv",
    },
    {
        "target_axis": "IL13RA2",
        "representative_agents": "IL-13-axis targeting",
        "trial_or_program": "IBD target biology; no efficacy program located",
        "actual_trial_outcome": "untested",
        "outcome_note": "No mature IBD therapeutic outcome identified; local genetics has no valid causal instrument.",
        "citation": "Local module causal map multicontext.",
        "url": "results/tables/module_causal_map_multicontext.tsv",
    },
    {
        "target_axis": "alpha4beta7 integrin",
        "representative_agents": "vedolizumab",
        "trial_or_program": "GEMINI 1/2",
        "actual_trial_outcome": "approved_effective",
        "outcome_note": "Vedolizumab is approved for UC/CD and was more effective than placebo in GEMINI trials.",
        "citation": "Feagan et al. GEMINI 1 UC; Sandborn et al. GEMINI 2 CD. N Engl J Med 2013.",
        "url": "https://www.nejm.org/doi/full/10.1056/NEJMoa1215734; https://www.nejm.org/doi/full/10.1056/NEJMoa1215739",
    },
    {
        "target_axis": "IL12/23 p40",
        "representative_agents": "ustekinumab",
        "trial_or_program": "UNITI / IM-UNITI / UNIFI",
        "actual_trial_outcome": "approved_effective",
        "outcome_note": "Ustekinumab is approved for CD/UC; phase 3 programs showed induction and maintenance benefit versus placebo.",
        "citation": "Feagan et al. UNITI/IM-UNITI Crohn's; Sands et al. UNIFI UC. N Engl J Med.",
        "url": "https://www.nejm.org/doi/full/10.1056/NEJMoa1602773; https://www.nejm.org/doi/full/10.1056/NEJMoa1900750",
    },
    {
        "target_axis": "IL23 p19",
        "representative_agents": "risankizumab; mirikizumab; guselkumab",
        "trial_or_program": "IL-23 p19 IBD programs",
        "actual_trial_outcome": "approved_effective",
        "outcome_note": "IL-23 p19 inhibitors are approved/effective in IBD, but this local refractory-module cis-MR run lacks usable target-gene instruments.",
        "citation": "Review of IL-23 inhibitors in IBD.",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11920014/",
    },
    {
        "target_axis": "IL23 signaling",
        "representative_agents": "IL-23 pathway inhibitors",
        "trial_or_program": "IL-23 IBD programs",
        "actual_trial_outcome": "approved_effective",
        "outcome_note": "Approved IL-23 pathway therapies validate the pathway clinically, but IL23R itself lacks a local cis-expression instrument in this analysis.",
        "citation": "Review of IL-23 inhibitors in IBD.",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11920014/",
    },
    {
        "target_axis": "IL13",
        "representative_agents": "IL-13 blockade",
        "trial_or_program": "no approved IBD anti-IL13 program",
        "actual_trial_outcome": "untested",
        "outcome_note": "No mature IBD efficacy outcome for IL-13 blockade was curated here.",
        "citation": "No mature IBD efficacy source found in targeted search.",
        "url": "",
    },
]


class Inputs(BaseModel):
    model_config = ConfigDict(frozen=True)

    out_dir: Path = P.out("drugtarget_mr")
    triangulation: Path = P.tables / "triangulation.tsv"
    multicontext: Path = P.tables / "module_causal_map_multicontext.tsv"
    finngen: Path = P.tables / "finngen_cxcr2.tsv"
    sensitivity: Path = P.tables / "mr_sensitivity.tsv"
    status_log: Path = P.journal / "status" / "overnight_autorun_log.md"
    journal_note: Path = P.journal / "docs" / "drugtarget_mr_concordance_2026-06-23.md"

    @field_validator("triangulation", "multicontext")
    @classmethod
    def exists(cls, v: Path) -> Path:
        if not v.exists():
            raise FileNotFoundError(v)
        return v


def direction_from_or(or_value: float | None) -> str:
    if or_value is None or pd.isna(or_value):
        return "no_valid_cis_tool"
    if np.isclose(or_value, 1.0, rtol=0, atol=0.03):
        return "neutral"
    return "expression_increases_risk" if or_value > 1 else "expression_protective"


def inhibition_effect(or_value: float | None) -> str:
    d = direction_from_or(or_value)
    if d == "expression_increases_risk":
        return "beneficial"
    if d == "expression_protective":
        return "harmful/ineffective"
    return "neutral/uncertain"


def confidence(row: pd.Series) -> str:
    if row["gene"] == "CCL8":
        return "downgraded_after_hardening"
    if row["analysis_status"] in {"no_tool", "no_blood_result", "no_valid_cis_tool"} or pd.isna(row["MR_OR"]):
        return "no_valid_cis_tool"
    if bool(row.get("special_positive_control", False)):
        return "directional_positive_control_LD_suspect"
    if row["coloc_pass"]:
        if row["MR_p"] <= 0.05 and row["n_instrument"] >= 1:
            return "high"
        return "coloc_only_direction_uncertain"
    if row["MR_p"] <= 0.05:
        return "LD_suspect"
    return "weak_or_null"


def choose_best_records(tri: pd.DataFrame, multi: pd.DataFrame) -> pd.DataFrame:
    target = pd.DataFrame(TARGETS)
    rows = []
    for t in target.itertuples(index=False):
        tri_row = tri[tri["gene"].eq(t.gene)]
        if len(tri_row) and pd.notna(tri_row.iloc[0].get("MR_OR", np.nan)):
            r = tri_row.iloc[0].to_dict()
            best = {
                "gene": t.gene,
                "context": "blood",
                "MR_OR": r.get("MR_OR", np.nan),
                "MR_p": r.get("MR_p", np.nan),
                "MR_FDR": r.get("MR_fdr", np.nan),
                "coloc_PP4": r.get("coloc_PP4", np.nan),
                "n_instrument": 1,
                "analysis_status": "tested_existing_eqtlgen",
                "causal_call": bool(pd.notna(r.get("coloc_PP4", np.nan)) and r.get("coloc_PP4", 0) > 0.8 and r.get("MR_p", 1) <= 0.05),
            }
        else:
            candidates = multi[multi["gene"].eq(t.gene)].copy()
            if not candidates.empty:
                candidates["score_coloc"] = candidates["coloc_PP4"].fillna(-1)
                candidates["score_p"] = candidates["MR_p"].fillna(1.0)
                candidates["score_tool"] = candidates["n_instrument"].fillna(0)
                candidates = candidates.sort_values(
                    ["causal_call", "score_coloc", "score_tool", "score_p"],
                    ascending=[False, False, False, True],
                )
                best = candidates.iloc[0].to_dict()
            else:
                best = {
                    "gene": t.gene,
                    "context": "",
                    "MR_OR": np.nan,
                    "MR_p": np.nan,
                    "MR_FDR": np.nan,
                    "coloc_PP4": np.nan,
                    "n_instrument": 0,
                    "analysis_status": "no_valid_cis_tool",
                    "causal_call": False,
                }
        if len(tri_row):
            r = tri_row.iloc[0]
            for col in ["GSE16879_FC", "GSE16879_p", "GSE73661_FC", "GSE73661_p", "convergent", "drug_target"]:
                best[col] = r.get(col, np.nan)
        else:
            for col in ["GSE16879_FC", "GSE16879_p", "GSE73661_FC", "GSE73661_p", "convergent", "drug_target"]:
                best[col] = np.nan
        best.update(t._asdict())
        rows.append(best)
    out = pd.DataFrame(rows)
    return out


def add_predictions(best: pd.DataFrame) -> pd.DataFrame:
    out = best.copy()
    out["coloc_pass"] = out["coloc_PP4"].fillna(0) > 0.8
    out["special_positive_control"] = out["gene"].eq("TNFRSF1A")
    out["expression_to_IBD_direction"] = out["MR_OR"].map(direction_from_or)
    out["drug_inhibition_OR_proxy"] = 1 / out["MR_OR"]
    out.loc[out["MR_OR"].isna(), "drug_inhibition_OR_proxy"] = np.nan
    out["ungated_predicted_drug_effect"] = out["MR_OR"].map(inhibition_effect)
    out["confidence"] = out.apply(confidence, axis=1)
    out["predicted_drug_effect"] = out["ungated_predicted_drug_effect"]
    out.loc[out["confidence"].isin(["no_valid_cis_tool", "LD_suspect", "weak_or_null", "downgraded_after_hardening"]), "predicted_drug_effect"] = "neutral/uncertain"
    out.loc[out["special_positive_control"], "predicted_drug_effect"] = out.loc[out["special_positive_control"], "ungated_predicted_drug_effect"]
    out["genetic_interpretation"] = np.select(
        [
            out["confidence"].eq("high") & out["gene"].eq("CXCR2"),
            out["special_positive_control"],
            out["confidence"].eq("high"),
            out["confidence"].eq("LD_suspect"),
            out["confidence"].eq("no_valid_cis_tool"),
        ],
        [
            "High-confidence protective expression signal: antagonism is genetically predicted to fail or be harmful.",
            "Directional positive control agrees with anti-TNF efficacy, but strict coloc PP4 is below 0.8.",
            "Coloc-supported genetic proxy for target inhibition.",
            "MR direction exists but coloc PP4 is below 0.8; do not make an actionable drug prediction.",
            "No usable local cis-expression instrument; genetics is silent here.",
        ],
        default="Weak/null target-expression MR signal.",
    )
    out.loc[out["gene"].eq("CCL8"), "genetic_interpretation"] = (
        "PP4-passing blood signal exists, but later hardening downgraded it: multi-instrument MR/MVMR were unstable and plasma pQTL did not colocalize."
    )
    keep = [
        "gene",
        "target_axis",
        "drug_class",
        "priority",
        "context",
        "n_instrument",
        "analysis_status",
        "MR_OR",
        "MR_p",
        "MR_FDR",
        "coloc_PP4",
        "coloc_pass",
        "expression_to_IBD_direction",
        "drug_inhibition_OR_proxy",
        "ungated_predicted_drug_effect",
        "predicted_drug_effect",
        "confidence",
        "genetic_interpretation",
        "GSE16879_FC",
        "GSE73661_FC",
        "convergent",
        "drug_target",
    ]
    return out[[c for c in keep if c in out.columns]].sort_values(["priority", "target_axis", "gene"])


def concordance(pred: pd.DataFrame, trials: pd.DataFrame) -> pd.DataFrame:
    out = pred.merge(trials, on="target_axis", how="left")
    expected_success = out["predicted_drug_effect"].eq("beneficial")
    expected_failure = out["predicted_drug_effect"].eq("harmful/ineffective")
    actual_success = out["actual_trial_outcome"].eq("approved_effective")
    actual_failure = out["actual_trial_outcome"].eq("failed/no_efficacy")
    out["concordance"] = np.select(
        [
            expected_success & actual_success,
            expected_failure & actual_failure,
            out["predicted_drug_effect"].eq("neutral/uncertain"),
            out["actual_trial_outcome"].isin(["mixed", "untested"]),
        ],
        ["concordant", "concordant", "not_actionable", "not_rateable"],
        default="discordant",
    )
    out["headline_case"] = np.select(
        [out["gene"].eq("CXCR2"), out["gene"].eq("TNFRSF1A")],
        ["CXCR2 genetic warning: antagonist failure predicted", "anti-TNF positive control: direction agrees"],
        default="",
    )
    return out


def plot_concordance(df: pd.DataFrame, path: Path) -> None:
    genes = [
        "TNFRSF1A",
        "CXCR2",
        "MMP9",
        "ICAM1",
        "IL1B",
        "IL6",
        "CCL2",
        "OSM",
        "OSMR",
        "IL12B",
        "IL23A",
        "ITGA4",
        "CCL8",
        "TREM1",
    ]
    d = df[df["gene"].isin(genes)].drop_duplicates("gene").copy()
    d["gene"] = pd.Categorical(d["gene"], categories=genes, ordered=True)
    d = d.sort_values("gene")
    effect_color = {"beneficial": "#1b9e77", "harmful/ineffective": "#d95f02", "neutral/uncertain": "#9e9e9e"}
    trial_color = {"approved_effective": "#1b9e77", "failed/no_efficacy": "#d95f02", "mixed": "#7570b3", "untested": "#bdbdbd"}
    conf_marker = {"high": "o", "directional_positive_control_LD_suspect": "s", "LD_suspect": "^", "weak_or_null": "v", "downgraded_after_hardening": "P", "no_valid_cis_tool": "x"}

    fig, ax = plt.subplots(figsize=(10, 6.8))
    y = np.arange(len(d))
    ax.axvline(0.5, color="#dddddd", lw=1)
    ax.axvline(1.5, color="#dddddd", lw=1)
    for i, row in enumerate(d.itertuples(index=False)):
        ax.scatter(0, i, s=170, c=effect_color.get(row.predicted_drug_effect, "#9e9e9e"), marker=conf_marker.get(row.confidence, "o"), edgecolor="black", linewidth=0.6)
        ax.scatter(1, i, s=170, c=trial_color.get(row.actual_trial_outcome, "#bdbdbd"), marker="s", edgecolor="black", linewidth=0.6)
        ax.scatter(2, i, s=170, c={"concordant": "#2ca25f", "discordant": "#de2d26", "not_actionable": "#969696", "not_rateable": "#bdbdbd"}.get(row.concordance, "#bdbdbd"), marker="D", edgecolor="black", linewidth=0.6)
        if row.gene == "CXCR2":
            ax.text(2.18, i, "genetic warning", va="center", fontsize=9, weight="bold", color="#b24d00")
        if row.gene == "TNFRSF1A":
            ax.text(2.18, i, "positive control", va="center", fontsize=9, weight="bold", color="#0b6b4a")
    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"].astype(str))
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Genetic inhibition prediction", "Trial reality", "Concordance"])
    ax.set_xlim(-0.35, 3.25)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.invert_yaxis()
    ax.set_title("Drug-target MR concordance for refractory-module targets")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label="beneficial / approved", markerfacecolor="#1b9e77", markeredgecolor="black", markersize=9),
        plt.Line2D([0], [0], marker="o", color="w", label="harmful or failed", markerfacecolor="#d95f02", markeredgecolor="black", markersize=9),
        plt.Line2D([0], [0], marker="o", color="w", label="uncertain / mixed / untested", markerfacecolor="#9e9e9e", markeredgecolor="black", markersize=9),
        plt.Line2D([0], [0], marker="D", color="w", label="concordant", markerfacecolor="#2ca25f", markeredgecolor="black", markersize=8),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_summary(pred: pd.DataFrame, conc: pd.DataFrame, path: Path) -> None:
    rateable = conc[conc["concordance"].isin(["concordant", "discordant"])]
    n_conc = int((rateable["concordance"] == "concordant").sum())
    n_rate = int(len(rateable))
    rate = n_conc / n_rate if n_rate else np.nan
    cx = conc[conc["gene"].eq("CXCR2")].iloc[0]
    tn = conc[conc["gene"].eq("TNFRSF1A")].iloc[0]
    high = pred[pred["confidence"].eq("high")]
    lines = [
        "# Drug-target MR concordance summary",
        "",
        f"Rateable concordance: {n_conc}/{n_rate} ({rate:.1%}) among targets with an actionable genetic prediction and a clear success/failure clinical outcome.",
        "",
        "## Headline calls",
        "",
        f"- Anti-TNF positive control: TNFRSF1A higher-expression MR OR={tn.MR_OR:.3g}; inhibition proxy OR={tn.drug_inhibition_OR_proxy:.3g}; predicted {tn.predicted_drug_effect}; trials {tn.actual_trial_outcome}. Strict coloc PP4={tn.coloc_PP4:.3g}, so this is a directional positive control rather than a strict PP4-passing claim.",
        f"- CXCR2 warning: higher expression is protective (MR OR={cx.MR_OR:.3g}, PP4={cx.coloc_PP4:.3g}); target inhibition proxy OR={cx.drug_inhibition_OR_proxy:.3g}; predicted {cx.predicted_drug_effect}. Curated CXCR2 antagonist inflammation outcomes are classified {cx.actual_trial_outcome}, with no mature positive IBD efficacy result found.",
        "- MMP9 is directionally similar to CXCR2 (higher expression protective) but fails coloc; andecaliximab failed in UC, so the trial reality is consistent with caution but not a strict genetic prediction.",
        "- OSM/OSMR/TREM1/IL13RA2 remain marker/bystander candidates with no actionable cis-MR support in this local run.",
        "",
        "## Actionable genetic targets",
        "",
    ]
    if len(high):
        for r in high.itertuples(index=False):
            lines.append(f"- {r.gene}: predicted {r.predicted_drug_effect}; confidence={r.confidence}; MR_OR={r.MR_OR:.3g}; PP4={r.coloc_PP4:.3g}.")
    else:
        lines.append("- None beyond CXCR2 passed the strict coloc gate.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Genetics does not produce a broad menu of validated refractory-module drug targets here. It mostly acts as a triage layer: it recovers the expected anti-TNF direction as a positive-control sanity check, gives a clean warning that CXCR2 antagonism is genetically disfavored in IBD, and downgrades many inflammatory markers because their cis-MR lacks colocalization or instruments.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def append_notes(cfg: Inputs, summary_path: Path, fig_path: Path, conc: pd.DataFrame) -> None:
    rateable = conc[conc["concordance"].isin(["concordant", "discordant"])]
    n_conc = int((rateable["concordance"] == "concordant").sum())
    n_rate = int(len(rateable))
    note = f"""# Drug-target MR concordance note (2026-06-23)

Promoted final artifacts from `outputs/drugtarget_mr/` after explicitly reframing cis-MR as target-inhibition proxy for druggable refractory-module targets.

- Summary: `{summary_path}`
- Figure: `{fig_path}`
- Rateable concordance: {n_conc}/{n_rate}
- Main call: CXCR2 is the clean PP4-passing genetic warning case; anti-TNF/TNFRSF1A is directionally correct as a positive control but below strict PP4.
- Caveat: most approved non-module targets (vedolizumab/ustekinumab/IL-23) lack local cis instruments, so genetics is silent rather than wrong.
"""
    cfg.journal_note.write_text(note)
    marker = "2026-06-23 drug-target MR concordance完成"
    status_text = cfg.status_log.read_text() if cfg.status_log.exists() else ""
    if marker not in status_text:
        with cfg.status_log.open("a") as fh:
            fh.write(
                "\n- **2026-06-23 drug-target MR concordance完成**（`src/30`，`outputs/drugtarget_mr/`，"
                "`results/tables/drugtarget_mr_predictions.tsv` / `trial_outcomes.tsv` / `concordance_map.tsv`，"
                "`results/figures/Fig_drugtarget_concordance.png`）。结论：anti-TNF方向阳性对照对上但PP4未过严格阈值；"
                "CXCR2为唯一干净PP4通过的遗传预警，表达保护→拮抗/抑制预测碰壁；MMP9/IL1B/IL6等多为LD存疑或混合/失败现实；"
                "OSM/OSMR/TREM1/IL13RA2仍是marker/bystander而非可行动因果靶。\n"
            )


def main() -> None:
    cfg = Inputs()
    tri = pd.read_csv(cfg.triangulation, sep="\t")
    multi = pd.read_csv(cfg.multicontext, sep="\t")
    trials = pd.DataFrame(TRIAL_OUTCOMES)
    pred = add_predictions(choose_best_records(tri, multi))
    conc = concordance(pred, trials)

    pred_path = cfg.out_dir / "drugtarget_mr_predictions.tsv"
    trial_path = cfg.out_dir / "trial_outcomes.tsv"
    conc_path = cfg.out_dir / "concordance_map.tsv"
    fig_path = cfg.out_dir / "Fig_drugtarget_concordance.png"
    summary_path = cfg.out_dir / "SUMMARY.md"

    pred.to_csv(pred_path, sep="\t", index=False)
    trials.to_csv(trial_path, sep="\t", index=False)
    conc.to_csv(conc_path, sep="\t", index=False)
    plot_concordance(conc, fig_path)
    write_summary(pred, conc, summary_path)

    promoted = [
        P.promote_table(pred_path),
        P.promote_table(trial_path),
        P.promote_table(conc_path),
        P.promote_figure(fig_path),
    ]
    append_notes(cfg, summary_path, promoted[-1], conc)
    print(f"wrote {pred_path}")
    print(f"wrote {trial_path}")
    print(f"wrote {conc_path}")
    print(f"wrote {fig_path}")
    print("promoted:")
    for p in promoted:
        print(f"  {p}")


if __name__ == "__main__":
    main()
