"""Identify GTEx colon/ileum samples and attach decade-bracket ages.
Outputs data/gtex/gut_samples.tsv  (SAMPID, SUBJID, tissue, age_bracket, age_mid, sex)
"""
import pandas as pd

BASE = "/Users/ujs/Downloads/lzy/data/raw/gtex"

GUT_TISSUES = {
    "Colon - Sigmoid": "colon_sigmoid",
    "Colon - Transverse": "colon_transverse",
    "Small Intestine - Terminal Ileum": "ileum",
}

attr = pd.read_csv(f"{BASE}/SampleAttributes.txt", sep="\t", low_memory=False)
attr = attr[attr["SMTSD"].isin(GUT_TISSUES)][["SAMPID", "SMTSD"]].copy()
attr["tissue"] = attr["SMTSD"].map(GUT_TISSUES)
# subject id = first two dash-fields, e.g. GTEX-1117F
attr["SUBJID"] = attr["SAMPID"].str.split("-").str[:2].str.join("-")

pheno = pd.read_csv(f"{BASE}/SubjectPhenotypes.txt", sep="\t")
# AGE like "60-69" -> midpoint 64.5
pheno["age_mid"] = pheno["AGE"].str.split("-").apply(lambda x: (int(x[0]) + int(x[1])) / 2)

out = attr.merge(pheno[["SUBJID", "AGE", "age_mid", "SEX"]], on="SUBJID", how="left")
out = out.rename(columns={"AGE": "age_bracket", "SEX": "sex"})[
    ["SAMPID", "SUBJID", "tissue", "age_bracket", "age_mid", "sex"]
]
out.to_csv(f"{BASE}/gut_samples.tsv", sep="\t", index=False)

print("gut samples:", len(out))
print(out["tissue"].value_counts().to_string())
print("\nage bracket distribution:")
print(out["age_bracket"].value_counts().sort_index().to_string())
