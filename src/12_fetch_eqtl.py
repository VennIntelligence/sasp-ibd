"""Fetch cis-eQTL associations for SenMayo genes from the eQTL Catalogue API.
Caches one TSV per (dataset, gene). Colon (sigmoid) primary, blood secondary.
Usage: python 12_fetch_eqtl.py [test|full]
"""
import json, sys, time, os
import urllib.request, urllib.error
import pandas as pd

BASE = "/Users/ujs/Downloads/lzy"
OUT = f"{BASE}/data/raw/eqtl"
os.makedirs(OUT, exist_ok=True)
ensg = json.load(open(f"{BASE}/data/external/genesets/senmayo_ensg.json"))
DATASETS = {"colon": "QTD000226", "blood": "QTD000356"}
FIELDS = ["rsid", "chromosome", "position", "ref", "alt", "beta", "se", "pvalue", "maf", "an"]
API = "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{ds}/associations"


def fetch_gene(ds_id, gene_id, max_pages=40, size=1000):
    recs, start = [], 0
    for _ in range(max_pages):
        url = f"{API.format(ds=ds_id)}?gene_id={gene_id}&size={size}&start={start}"
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    data = json.load(r)
                break
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
                time.sleep(2 * (attempt + 1))
        else:
            break
        if not isinstance(data, list) or not data:
            break
        recs.extend(data)
        if len(data) < size:
            break
        start += size
        time.sleep(0.05)
    return recs


def run(genes):
    summary = []
    for i, (sym, gid) in enumerate(genes.items(), 1):
        for tag, ds in DATASETS.items():
            path = f"{OUT}/{tag}_{sym}.tsv"
            if os.path.exists(path):
                df = pd.read_csv(path, sep="\t") if os.path.getsize(path) > 5 else pd.DataFrame()
            else:
                recs = fetch_gene(ds, gid)
                df = pd.DataFrame(recs)
                if len(df):
                    df = df[[c for c in FIELDS if c in df.columns]]
                df.to_csv(path, sep="\t", index=False)
            minp = df["pvalue"].min() if len(df) and "pvalue" in df else None
            n = len(df)
            summary.append((sym, tag, n, minp))
        if i % 10 == 0:
            print(f"  ...{i}/{len(genes)} genes done")
    s = pd.DataFrame(summary, columns=["gene", "tissue", "n_snps", "min_p"])
    s.to_csv(f"{OUT}/_fetch_summary.tsv", sep="\t", index=False)
    sig = s[(s.min_p < 5e-8)]
    print(f"\ngene-tissue pairs with genome-wide cis-eQTL (p<5e-8): {len(sig)}")
    print(sig.sort_values('min_p').head(25).to_string(index=False))
    return s


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    if mode == "test":
        test = {k: ensg[k] for k in ["ERAP2", "IL6", "CXCL8", "MMP9", "ICAM1"] if k in ensg}
        print("TEST genes:", list(test))
        run(test)
    else:
        run(ensg)
