"""Independent verification for the KLF5/TAZ partner programme.
Recomputes headline numbers from the RAW files with separate code from the
pipeline scripts. Run: python3 verify_headlines.py — exits non-zero on any fail.
"""
import pandas as pd, numpy as np, json, sys

BASE = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120"
KLF = f"{BASE}/klf5-taz-partners"
fails = []
def check(name, got, want, tol=0.005):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: got {got}, want {want}")
    if not ok:
        fails.append(name)

# --- 1. DepMap 24Q4 codependency, recomputed from CRISPRGeneEffect.csv ---
eff = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", index_col=0)
eff.columns = [c.split(" (")[0] for c in eff.columns]
for anchor, tops in [("WWTR1", {"PKN2": 0.428, "FERMT2": 0.426, "TEAD1": 0.424, "PTK2": 0.413,
                                "ITGAV": 0.411, "AMOTL2": -0.424, "LATS2": -0.418, "NF2": -0.376}),
                    ("KLF5", {"STX4": 0.372, "EGFR": 0.365, "GRHL2": 0.346, "TCF7L2": 0.337,
                              "ITGA3": 0.327, "CTNNB1": 0.291, "YAP1": 0.258, "ERBB2": 0.247})]:
    for g, want in tops.items():
        check(f"codepend {anchor}-{g}", float(eff[anchor].corr(eff[g])), want)
check("codepend WWTR1-KLF5", float(eff["WWTR1"].corr(eff["KLF5"])), -0.051)
check("codepend KLF5-FERMT2 (negative)", float(eff["KLF5"].corr(eff["FERMT2"])), -0.227)

# --- 2. Cell-line expression co-expression, recomputed from the expression matrix ---
expr = pd.read_csv("/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
expr.columns = [c.split(" (")[0] for c in expr.columns]
check("expr WWTR1-KLF5", float(expr["WWTR1"].corr(expr["KLF5"])), 0.213)
check("expr WWTR1-YAP1", float(expr["WWTR1"].corr(expr["YAP1"])), 0.710)
check("expr KLF5-EGFR", float(expr["KLF5"].corr(expr["EGFR"])), 0.565)

# --- 3. TCGA co-expression, recomputed from the fetched cBioPortal matrices ---
for study, f, pairs in [
    ("BRCA", "tcga_brca_mrna.csv", {("KLF5","WWTR1"): 0.301, ("WWTR1","CDK6"): 0.529, ("WWTR1","FERMT2"): 0.445}),
    ("COAD", "tcga_coadread_mrna.csv", {("KLF5","WWTR1"): -0.129, ("WWTR1","ITGB1"): 0.617, ("WWTR1","FERMT2"): 0.607}),
    ("SKCM", "tcga_skcm_mrna.csv", {("KLF5","TP63"): 0.750, ("KLF5","GRHL2"): 0.707, ("WWTR1","CDK6"): 0.458}),
]:
    piv = pd.read_csv(f"{KLF}/{f}", index_col=0)
    for (a, b), want in pairs.items():
        check(f"TCGA {study} {a}-{b}", float(piv[a].corr(piv[b])), want)

# --- 4. Relaxed-bar passes, recomputed from the per-cancer scan files ---
scans = {
 "TNBC": "/home/hatch/workspace/research_notes/depmap-tnbc-selective-targets-20261120/dependency_scan.csv",
 "HCC": f"{BASE}/hcc/dependency_scan.csv", "Melanoma": f"{BASE}/melanoma/dependency_scan.csv",
 "Colon": f"{BASE}/colon/dependency_scan.csv", "Glioma": f"{BASE}/glioma/dependency_scan.csv",
 "PDAC": f"{BASE}/pancreatic/dependency_scan.csv", "HNSC": f"{BASE}/head-and-neck/dependency_scan.csv",
 "Ovarian": f"{BASE}/ovarian/dependency_scan.csv",
}
def relaxed(row):
    return (row["Enrichment"] <= -0.05 and row["_gfrac"] >= 0.20 and row["FracDiff"] >= 0.10
            and row["_gw"] <= 0.60 and not bool(row["_pan"]))
expected_pass = {
 ("PTK2","TNBC"): True, ("PTK2","Melanoma"): True, ("PTK2","Glioma"): True, ("PTK2","Ovarian"): True,
 ("PTK2","PDAC"): False, ("PTK2","HNSC"): False,
 ("TEAD1","Melanoma"): True, ("TEAD1","Glioma"): False, ("TEAD1","Ovarian"): False,
 ("EGFR","HNSC"): True, ("EGFR","Colon"): True, ("EGFR","PDAC"): True,
 ("CTNNB1","Colon"): True, ("CTNNB1","PDAC"): True,
 ("CDK6","TNBC"): True, ("CDK6","HNSC"): True,
 ("ERBB2","HNSC"): True, ("ERBB2","Ovarian"): True,
 ("CCNE1","Ovarian"): True, ("MAPK1","Colon"): True, ("MAPK1","PDAC"): True,
 ("PIK3CA","Colon"): True, ("PIK3CA","PDAC"): True, ("PIK3CA","HNSC"): True,
 ("WWTR1","TNBC"): True, ("WWTR1","PDAC"): True, ("WWTR1","Ovarian"): True,
 ("YAP1","HNSC"): True, ("YAP1","Colon"): True,
 ("CDK2","Ovarian"): False, ("CCND1","Melanoma"): False, ("MYC","Colon"): False,
}
for (gene, cancer), want in expected_pass.items():
    df = pd.read_csv(scans[cancer])
    df.columns = [c.strip() for c in df.columns]
    gfrac = [c for c in df.columns if c.lower().endswith("dep_frac") or c.lower().endswith("depfrac")][0]
    gwc = [c for c in df.columns if "genomewide" in c.lower()][0]
    pan = [c for c in df.columns if c.lower() == "panessential"][0]
    fd = [c for c in df.columns if c.lower() == "fracdiff"][0]
    row = df[df.Gene == gene].iloc[0]
    got = relaxed({"Enrichment": row["Enrichment"], "_gfrac": row[gfrac], "FracDiff": row[fd],
                   "_gw": row[gwc], "_pan": row[pan]})
    check(f"relaxed-bar {gene} in {cancer}", got, want)

# --- 5. Genome-wide dependent fractions (>0.60 = excluded by the programme bar) ---
tnbc = pd.read_csv(scans["TNBC"])
gwcol = [c for c in tnbc.columns if "genomewide" in c.lower()][0]
gw = dict(zip(tnbc.Gene, tnbc[gwcol]))
for g, want in {"CDK2": 0.688, "CCND1": 0.649, "MYC": 0.948, "CRKL": 0.749,
                "RAC1": 0.798, "PRMT5": 0.941, "BRD4": 0.952, "PTK2": 0.495,
                "CDK6": 0.505, "TEAD1": 0.167, "EGFR": 0.190}.items():
    check(f"genome-wide dep frac {g}", float(gw[g]), want)

# --- 6. Partner CSVs agree with the computed tables ---
for f, anchor in [("taz_partners.csv", "WWTR1"), ("klf5_partners.csv", "KLF5")]:
    df = pd.read_csv(f"{KLF}/{f}")
    print(f"[INFO] {f}: {len(df)} rows; Tier 1 = {', '.join(df[df.tier=='1'].partner)}")
    for _, r in df.iterrows():
        if pd.notna(r["dep_r_vs_anchor"]):
            check(f"{f} dep_r {r['partner']}", float(r["dep_r_vs_anchor"]),
                  float(eff[anchor].corr(eff[r["partner"]])))

print()
if fails:
    print(f"FAILED: {len(fails)} checks — {fails}")
    sys.exit(1)
print("ALL CHECKS PASSED")
