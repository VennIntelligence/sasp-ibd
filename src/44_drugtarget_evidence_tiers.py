#!/usr/bin/env python3
"""Build the strict drug-target evidence tier table for the manuscript.

The table is deliberately conservative. It separates:
  - strict actionable genetics (CXCR2 only);
  - directional controls (TNFRSF1A/anti-TNF);
  - genetics-silent approved pathways (IL12/23, integrin targets);
  - LD-suspect MR signals that should not be interpreted as drug predictions;
  - bystander or downgraded marker signals.

Run from repo root:
    .venv/bin/python src/44_drugtarget_evidence_tiers.py
"""
from __future__ import annotations

from typing import Optional

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from paths import P


class EvidenceRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    tier_order: int = Field(ge=1)
    evidence_tier: str
    gene: str
    target_axis: str
    representative_agents: str
    main_text_role: str
    n_instrument: int
    context: str
    mr_or: Optional[float]
    mr_fdr: Optional[float]
    coloc_pp4: Optional[float]
    predicted_drug_effect: str
    clinical_outcome: str
    concordance: str
    manuscript_interpretation: str
    caveat: str
    citation: str
    url: str


def _str(v) -> str:
    if pd.isna(v):
        return ""
    return str(v)


def _float(v) -> Optional[float]:
    x = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
    return None if pd.isna(x) else float(x)


def _int(v) -> int:
    x = pd.to_numeric(pd.Series([v]), errors="coerce").iloc[0]
    return 0 if pd.isna(x) else int(x)


def classify(row: pd.Series) -> tuple[int, str, str, str, str]:
    gene = _str(row["gene"])
    priority = _str(row["priority"])
    confidence = _str(row["confidence"])
    outcome = _str(row["actual_trial_outcome"])
    n_inst = _int(row["n_instrument"])

    if gene == "CXCR2":
        return (
            1,
            "strict_actionable_warning",
            "Main text key warning",
            "High-confidence protective expression signal predicts that CXCR2 antagonism is harmful or ineffective.",
            "Clinical no-efficacy evidence is not a mature positive IBD endpoint; state as a genetic warning, not a completed IBD trial validation.",
        )
    if gene == "TNFRSF1A":
        return (
            2,
            "directional_positive_control",
            "Main text positive control",
            "Anti-TNF direction agrees with known IBD efficacy, but the coloc result is below the strict threshold.",
            "Use as a directional sanity check only; PP4 is 0.18.",
        )
    if gene == "CCL8":
        return (
            3,
            "downgraded_hint_supplementary",
            "Supplementary cautionary signal",
            "CCL8 has a coloc-positive single-lead signal and FinnGen support, but hardening weakens it.",
            "Do not present as a co-equal causal target with CXCR2; multi-instrument, MVMR, and pQTL checks are not robust.",
        )
    if priority == "approved_context" and (n_inst == 0 or confidence in {"weak_or_null", "no_valid_cis_tool"}):
        return (
            4,
            "genetics_silent_approved_pathway",
            "Limitations / context",
            "The pathway is clinically validated; this local cis-expression analysis either has no usable target-gene instrument or the instrument fails to produce a colocalized causal signal.",
            "Genetics is silent or uninformative, not discordant with clinical efficacy.",
        )
    if gene in {"OSM", "OSMR", "IL13RA2"}:
        return (
            5,
            "marker_bystander_no_tool",
            "Negative triage",
            "Strong refractory marker biology, but no usable causal expression instrument in tested contexts.",
            "Absence of an instrument is not proof of no biology; it blocks an actionable MR claim.",
        )
    if gene == "TREM1":
        return (
            6,
            "tested_bystander_null",
            "Negative triage",
            "TREM1 has a strong neutrophil eQTL but near-null MR and no coloc support.",
            "Useful as a concrete example that marker strength is not causality.",
        )
    if confidence == "LD_suspect":
        return (
            7,
            "ld_suspect_not_actionable",
            "Supplementary filtered signal",
            "MR direction exists but coloc fails, so the signal is treated as LD-suspect rather than actionable.",
            "Do not use for drug prediction.",
        )
    if confidence in {"weak_or_null", "no_valid_cis_tool"} or outcome == "untested":
        return (
            8,
            "unrateable_or_weak",
            "Supplementary context",
            "No actionable local cis-MR prediction is available.",
            "Mention only if needed for completeness.",
        )
    return (
        9,
        "unclassified_context",
        "Supplementary context",
        "Not part of the strict rateable set.",
        "Manual review required before any manuscript claim.",
    )


def main() -> None:
    df = pd.read_csv(P.tables / "concordance_map.tsv", sep="\t")
    rows: list[EvidenceRow] = []
    for _, r in df.iterrows():
        tier_order, tier, role, interpretation, caveat = classify(r)
        rows.append(EvidenceRow(
            tier_order=tier_order,
            evidence_tier=tier,
            gene=_str(r["gene"]),
            target_axis=_str(r["target_axis"]),
            representative_agents=_str(r["representative_agents"]),
            main_text_role=role,
            n_instrument=_int(r["n_instrument"]),
            context=_str(r["context"]),
            mr_or=_float(r["MR_OR"]),
            mr_fdr=_float(r["MR_FDR"]),
            coloc_pp4=_float(r["coloc_PP4"]),
            predicted_drug_effect=_str(r["predicted_drug_effect"]),
            clinical_outcome=_str(r["actual_trial_outcome"]),
            concordance=_str(r["concordance"]),
            manuscript_interpretation=interpretation,
            caveat=caveat,
            citation=_str(r["citation"]),
            url=_str(r["url"]),
        ))

    out_df = pd.DataFrame([row.model_dump() for row in rows])
    out_df = out_df.sort_values(["tier_order", "gene"]).reset_index(drop=True)

    outdir = P.out("44_drugtarget_evidence_tiers")
    out_tsv = outdir / "drugtarget_evidence_tiers.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    promoted = P.promote_table(out_tsv)

    summary = out_df.groupby("evidence_tier", sort=False)["gene"].apply(lambda s: ", ".join(s)).reset_index()
    summary_path = outdir / "SUMMARY.md"
    summary_lines = [
        "# Drug-target evidence tier summary",
        "",
        "Strict rateable/actionable genetics is intentionally narrow.",
        "",
    ]
    for _, s in summary.iterrows():
        summary_lines.append(f"- **{s['evidence_tier']}**: {s['gene']}")
    summary_lines.extend([
        "",
        "Writing rule: only CXCR2 is a strict coloc-passing actionable warning; TNFRSF1A is a directional positive control; CCL8 is supplementary and downgraded.",
    ])
    summary_path.write_text("\n".join(summary_lines) + "\n")

    print(f"wrote {out_tsv}")
    print(f"promoted {promoted}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
