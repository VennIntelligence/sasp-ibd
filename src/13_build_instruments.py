"""Build MR instruments for SenMayo genes from eQTLGen significant cis-eQTL.
For each gene -> lead cis-eQTL SNP (max |Z|). Convert Z to beta/se (Zhu 2016)
using assessed-allele frequency from the eQTLGen AF file. Match later to GWAS
by rsid (eQTLGen=hg19, GWAS=GRCh38, so rsid is the safe key).
"""
import gzip, json
import numpy as np
import pandas as pd

BASE = "/Users/ujs/Downloads/lzy"
EG = f"{BASE}/data/raw/eqtlgen"
sets = json.load(open(f"{BASE}/data/external/genesets/senescence_sets.json"))
SENMAYO = set(sets["SenMayo"])

# ---- 1. lead cis-eQTL per SenMayo gene ----
print("scanning eQTLGen significant cis-eQTL...")
keep = []
with gzip.open(f"{EG}/cis_sig.txt.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    print("eQTLGen columns:", header)
    iz, isym, isnp = idx["Zscore"], idx["GeneSymbol"], idx["SNP"]
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[isym] in SENMAYO:
            keep.append(p)
cols = header
df = pd.DataFrame(keep, columns=cols)
for c in ["Pvalue", "Zscore", "SNPChr", "SNPPos", "NrSamples", "FDR"]:
    if c in df:
        df[c] = pd.to_numeric(df[c], errors="coerce")
print(f"SenMayo cis-eQTL rows: {len(df)} across {df['GeneSymbol'].nunique()} genes")
# lead per gene = max |Z|
df["absZ"] = df["Zscore"].abs()
lead = df.sort_values("absZ", ascending=False).drop_duplicates("GeneSymbol")
print(f"lead instruments: {len(lead)} genes")

# ---- 2. assessed-allele frequency from AF file ----
need = set(lead["SNP"])
af = {}
print("scanning AF file for instrument SNPs...")
with gzip.open(f"{EG}/snp_af.txt.gz", "rt") as fh:
    afhead = fh.readline().rstrip("\n").split("\t")
    aidx = {c: i for i, c in enumerate(afhead)}
    print("AF columns:", afhead)
    isnp_af = aidx.get("SNP", 0)
    # find allele + freq columns
    a_alleleB = aidx.get("AlleleB")
    a_freqB = aidx.get("AlleleB_all", aidx.get("AlleleB_freq"))
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[isnp_af] in need:
            af[p[isnp_af]] = {"AlleleB": p[a_alleleB] if a_alleleB is not None else None,
                              "freqB": float(p[a_freqB]) if a_freqB is not None and p[a_freqB] not in ("", "NA") else np.nan}
            if len(af) == len(need):
                break
print(f"AF found for {len(af)}/{len(need)} instrument SNPs")

# ---- 3. Z -> beta/se (Zhu 2016): need freq of ASSESSED allele ----
rows = []
for _, r in lead.iterrows():
    a = af.get(r["SNP"])
    if a is None or np.isnan(a["freqB"]):
        continue
    # eQTLGen assessed allele freq: AF file gives AlleleB freq; align
    if a["AlleleB"] == r["AssessedAllele"]:
        p_assessed = a["freqB"]
    elif a["AlleleB"] == r["OtherAllele"]:
        p_assessed = 1 - a["freqB"]
    else:
        continue
    z, n = r["Zscore"], r["NrSamples"]
    denom = 2 * p_assessed * (1 - p_assessed) * (n + z**2)
    if denom <= 0:
        continue
    se = 1 / np.sqrt(denom)
    beta = z * se
    rows.append({"gene": r["GeneSymbol"], "rsid": r["SNP"],
                 "assessed": r["AssessedAllele"], "other": r["OtherAllele"],
                 "eaf": p_assessed, "z": z, "n": n, "p_eqtl": r["Pvalue"],
                 "beta_eqtl": beta, "se_eqtl": se})
inst = pd.DataFrame(rows)
inst.to_csv(f"{BASE}/outputs/mr/instruments.tsv", sep="\t", index=False)
print(f"\nfinal instruments with beta/se: {len(inst)} genes")
print(inst[["gene", "rsid", "eaf", "z", "beta_eqtl", "se_eqtl"]].head(15).to_string(index=False))
