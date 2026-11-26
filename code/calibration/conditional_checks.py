#!/usr/bin/env python3
"""Conditional-dependency supplement for the drug-target calibration run.

The whole-cancer pipeline cannot see dependencies that exist only in a
mutation-defined subset. These checks test the three conditional cases directly
from the raw DepMap 24Q4 files (CRISPRGeneEffect.csv + OmicsSomaticMutations.csv,
LikelyLoF=True, deduplicated by ModelID):
  1. PARP1 in BRCA1/2-LoF lines vs the rest (pan-cancer; and within ovarian/PDAC/TNBC).
  2. MDM2 in TP53-WT lines vs TP53-LoF lines (pan-cancer; and within melanoma/RCC).
  3. EPAS1 in VHL-LoF vs VHL-WT RCC lines.
Also: identity of the single RET- and NTRK1-dependent lines genome-wide.
Output: conditional_checks.csv (one row per comparison).
"""
import csv, os
import pandas as pd
from scipy import stats

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
CAL = os.path.join(BASE, "drug-target-calibration")
VRK = os.path.expanduser("~/workspace/research_notes/vrk1-codependency-20261108")
MUT = os.path.expanduser("~/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv")
TNBC_CODEP = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-codependency-20261120")

def col_for(df, symbol):
    hits = [c for c in df.columns if c.split(" (")[0] == symbol]
    assert len(hits) == 1, (symbol, hits)
    return hits[0]

def summarise(s):
    s = s.dropna()
    return dict(n_scored=len(s), mean=round(float(s.mean()), 3),
                dep_n=int((s < -0.5).sum()),
                dep_frac=round(float((s < -0.5).mean()), 3))

def main():
    print("loading ...")
    crispr = pd.read_csv(os.path.join(VRK, "CRISPRGeneEffect.csv"), index_col=0)
    model = pd.read_csv(os.path.join(VRK, "Model.csv"), index_col="ModelID")
    disease = model["OncotreePrimaryDisease"].to_dict()
    cellname = model["CellLineName"].to_dict()
    mut = pd.read_csv(MUT, usecols=["ModelID", "HugoSymbol", "LikelyLoF"], low_memory=True)
    lof = mut[mut["LikelyLoF"] == True]
    lof_sets = {g: set(lof[lof["HugoSymbol"] == g]["ModelID"]) for g in ["BRCA1", "BRCA2", "TP53", "VHL"]}
    for g, s in lof_sets.items():
        print(f"  {g} LikelyLoF lines: {len(s)}")
    brca_lof = lof_sets["BRCA1"] | lof_sets["BRCA2"]
    crispr_ids = set(crispr.index)

    def group_by_disease(label):
        return sorted(m for m in crispr_ids if disease.get(m) == label)
    tnbc_rows = list(csv.DictReader(open(os.path.join(TNBC_CODEP, "tnbc_lines.csv"))))
    tnbc = sorted({r["ModelID"] for r in tnbc_rows if r.get("Core_TNBC") == "True"} & crispr_ids)
    groups = {"Pan-cancer": sorted(crispr_ids),
              "Ovarian": group_by_disease("Ovarian Epithelial Tumor"),
              "PDAC": group_by_disease("Pancreatic Adenocarcinoma"),
              "TNBC": tnbc,
              "Melanoma": group_by_disease("Melanoma"),
              "RCC": group_by_disease("Renal Cell Carcinoma")}

    rows = []
    def compare(label, gene, in_ids, out_ids):
        col = col_for(crispr, gene)
        a = crispr.loc[sorted(set(in_ids) & crispr_ids), col]
        b = crispr.loc[sorted(set(out_ids) & crispr_ids), col]
        sa, sb = summarise(a), summarise(b)
        p = stats.ttest_ind(a.dropna(), b.dropna(), equal_var=False).pvalue
        rows.append(dict(Comparison=label, Gene=gene,
                         In_n=sa["n_scored"], In_mean=sa["mean"], In_dep_n=sa["dep_n"], In_dep_frac=sa["dep_frac"],
                         Out_n=sb["n_scored"], Out_mean=sb["mean"], Out_dep_n=sb["dep_n"], Out_dep_frac=sb["dep_frac"],
                         Welch_p=f"{p:.3g}"))
        print(f"{label}: {gene} in={sa['mean']} ({sa['dep_n']}/{sa['n_scored']}) vs out={sb['mean']} "
              f"({sb['dep_n']}/{sb['n_scored']}) p={p:.3g}")

    # 1. PARP1 x BRCA1/2 LoF
    print("\n1. PARP1 dependency by BRCA1/2 loss-of-function")
    for gname in ["Pan-cancer", "Ovarian", "PDAC", "TNBC"]:
        g = groups[gname]
        compare(f"PARP1 | BRCA1/2-LoF vs rest | {gname}", "PARP1",
                [m for m in g if m in brca_lof], [m for m in g if m not in brca_lof])
    # 2. MDM2 x TP53
    print("\n2. MDM2 dependency by TP53 loss-of-function (WT = no TP53 LikelyLoF call)")
    tp53_lof = lof_sets["TP53"]
    for gname in ["Pan-cancer", "Melanoma", "RCC"]:
        g = groups[gname]
        compare(f"MDM2 | TP53-WT vs TP53-LoF | {gname}", "MDM2",
                [m for m in g if m not in tp53_lof], [m for m in g if m in tp53_lof])
    # 3. EPAS1 x VHL within RCC
    print("\n3. EPAS1 dependency by VHL loss-of-function within RCC")
    rcc = groups["RCC"]
    compare("EPAS1 | VHL-LoF vs VHL-WT | RCC", "EPAS1",
            [m for m in rcc if m in lof_sets["VHL"]], [m for m in rcc if m not in lof_sets["VHL"]])
    # also: EPAS1 in VHL-LoF lines pan-cancer
    compare("EPAS1 | VHL-LoF vs rest | Pan-cancer", "EPAS1",
            [m for m in crispr_ids if m in lof_sets["VHL"]],
            [m for m in crispr_ids if m not in lof_sets["VHL"]])

    # 4. single-line dependencies
    print("\n4. The single RET- and NTRK1-dependent lines")
    for gene in ["RET", "NTRK1"]:
        col = col_for(crispr, gene)
        s = crispr[col].dropna()
        for mid in s[s < -0.5].index:
            print(f"  {gene}: {cellname.get(mid, mid)} [{disease.get(mid)}] Chronos={s[mid]:.3f}")

    with open(os.path.join(CAL, "conditional_checks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        w.writerows(rows)
    print(f"\nwrote conditional_checks.csv ({len(rows)} rows)")

if __name__ == "__main__":
    main()
