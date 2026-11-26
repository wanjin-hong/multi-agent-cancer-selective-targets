"""Independent verification: recompute headline numbers from the RAW files (not the CSVs),
then compare against dependency_scan.csv / overexpression.csv / candidates_both.csv.
"""
import json, numpy as np, pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/colon"
BASE = "/home/hatch/workspace/research_notes"
model = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/Model.csv", low_memory=False)
colon = set(model.loc[model["OncotreePrimaryDisease"] == "Colorectal Adenocarcinoma", "ModelID"])
eff = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/CRISPRGeneEffect.csv")
expr = pd.read_csv(f"{BASE}/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
scan = pd.read_csv(f"{OUT}/dependency_scan.csv").set_index("Gene")
ov = pd.read_csv(f"{OUT}/overexpression.csv").set_index("Gene")
cand = pd.read_csv(f"{OUT}/candidates_both.csv")

eids = eff["Unnamed: 0"]; is_c = eids.isin(colon).values
xids = expr["Unnamed: 0"]; is_cx = xids.isin(colon).values

def check(name, got, want, tol=1e-6):
    ok = abs(got - want) <= tol * max(1, abs(want))
    print(f"{'OK ' if ok else 'MISMATCH'} {name}: recomputed={got:.6g} csv={want:.6g}")
    return ok

allok = True
for g in ["CTNNB1", "KRAS", "KLF5", "TCF7L2", "YAP1", "WRN", "WWTR1", "PTK2", "TEAD1"]:
    col = [c for c in eff.columns if c.startswith(g + " (")][0]
    v = eff[col].values.astype(float)
    gc, cc = v[is_c], v[~is_c]
    gc, cc = gc[~np.isnan(gc)], cc[~np.isnan(cc)]
    r = scan.loc[g]
    allok &= check(f"{g} colon mean", gc.mean(), r["Colon_mean"])
    allok &= check(f"{g} enrichment", gc.mean() - cc.mean(), r["Enrichment"])
    allok &= check(f"{g} colon dep frac", (gc < -0.5).mean(), r["Colon_dep_frac"])
    allok &= check(f"{g} noncolon dep frac", (cc < -0.5).mean(), r["NonColon_dep_frac"])
    print(f"    {g}: dep {int((gc<-0.5).sum())}/{len(gc)} vs {int((cc<-0.5).sum())}/{len(cc)}; colon mean {gc.mean():.3f}")

# cell-line RNA for CDX2, KLF5, FERMT1, YAP1
for g in ["CDX2", "KLF5", "FERMT1", "YAP1", "SOX9"]:
    col = [c for c in expr.columns if c.startswith(g + " (")][0]
    v = expr[col].values.astype(float)
    gc, cc = v[is_cx], v[~is_cx]
    gc, cc = gc[~np.isnan(gc)], cc[~np.isnan(cc)]
    fc = (2 ** gc.mean() - 1) / (2 ** cc.mean() - 1)
    r = ov.loc[g]
    allok &= check(f"{g} RNA delta log2", gc.mean() - cc.mean(), r["RNA_delta_log2"])
    allok &= check(f"{g} RNA TPM FC", fc, r["RNA_TPM_FC"])
    print(f"    {g}: colon RNA n={len(gc)} mean={gc.mean():.3f}; comp n={len(cc)} mean={cc.mean():.3f}")

# TCGA RNA for CDX2, KLF5, FERMT1, CTNNB1 from raw JSONs
coad = json.load(open(f"{OUT}/downloads/cbioportal_pca_coadread_rna.json"))
agg = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/pca_comparator_agg.json")).items()}
m = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/entrez2symbol.json")).items()}
inv = {v: k for k, v in m.items()}
def stype(s):
    try: return s.split("-")[3][:2]
    except Exception: return "?"
for g in ["CDX2", "KLF5", "FERMT1", "CTNNB1", "YAP1"]:
    eid = inv[g]
    x = np.array([np.log2(d["value"] + 1) for d in coad if d["entrezGeneId"] == eid and stype(d["sampleId"]) == "01"])
    s, ss, n2 = agg[eid]
    fc = x.mean() - s / n2
    allok &= check(f"{g} TCGA log2FC", fc, ov.loc[g, "TCGA_RNA_log2FC"])
    print(f"    {g}: TCGA coad n={len(x)} mean={x.mean():.3f}; other n={n2} mean={s/n2:.3f}")

# counts
strict = set(open(f"{OUT}/downloads/strict_genes.txt").read().split())
relaxed = set(open(f"{OUT}/downloads/relaxed_genes.txt").read().split())
print("strict count:", len(strict), "| relaxed count:", len(relaxed))
print("candidates (relaxed+overexp):", len(cand), "| Tier A in top20:", int((cand.head(20)["Tier"].str.startswith("A")).sum()))
print("top20 genes:", cand.head(20)["Gene"].tolist())
print("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED")
