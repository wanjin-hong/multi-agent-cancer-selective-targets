"""EGFR-mutant lung selective-targets scan — DepMap 24Q4, TNBC pipeline replica.

Primary group: NSCLC lines with ACTIVATING EGFR mutation (CRISPR n=6, RNA n=9).
Secondary (sensitivity): all EGFR-mutant NSCLC lines (CRISPR n=21, RNA n=28).
Comparator: EGFR-wild-type NSCLC lines (CRISPR n=77, RNA n=106). SCLC excluded
entirely; 9 expression-panel NSCLC lines with no mutation-file record excluded
(unknown status).

Activating rule (ProteinChange, from OmicsSomaticMutations.csv):
  INCLUDE any of: exon-19 in-frame deletion/delins p.E746_A750del,
  p.L747_P753delinsS, p.L747_E749del; p.L858R; p.T790M.
  In lung models there is NO exon-20 insertion, NO G719S/S768I/L861Q
  (those categories are empty, not merely uncounted).
  EXCLUDE: p.R309G and other VUS (C307R, C311F, V592F, Q1159H, D1127N,
  S1130G, Q1174R); p.G857_L858delinsDM (non-canonical delins, LK-2);
  LikelyLoF frameshifts (ETCC-016, HCC2935); SCLC lines.
  NOTE: AMClass/Hotspot alone cannot define this — R309G is AMClass
  likely_pathogenic, and LOU-NH91's L747_P753delinsS is Hotspot=False.
"""
import numpy as np
import pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/egfr-mutant-lung"
GE = "/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv"
MODEL = "/home/hatch/workspace/research_notes/vrk1-codependency-20261108/Model.csv"
EXPR = "/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MUT = "/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv"

ACT_PC = {"p.E746_A750del", "p.L747_P753delinsS", "p.L747_E749del", "p.L858R", "p.T790M"}

def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p); order = np.argsort(p); ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1.0)
    return out

def welch_stats(X, grp):
    """Vectorised Welch t-test, group (bool) vs complement. NaN-aware."""
    a = X[grp]; b = X[~grp]
    m1 = np.nanmean(a, axis=0); m2 = np.nanmean(b, axis=0)
    s1 = np.nanstd(a, axis=0, ddof=1); s2 = np.nanstd(b, axis=0, ddof=1)
    n1 = np.sum(~np.isnan(a), axis=0); n2 = np.sum(~np.isnan(b), axis=0)
    se = np.sqrt(s1**2 / n1 + s2**2 / n2)
    t = (m1 - m2) / se
    dfree = (s1**2 / n1 + s2**2 / n2) ** 2 / ((s1**2 / n1) ** 2 / (n1 - 1) + (s2**2 / n2) ** 2 / (n2 - 1))
    p = 2 * stats.t.sf(np.abs(t), dfree)
    return m1, m2, n1, n2, p

# ---------- load ----------
model = pd.read_csv(MODEL)
md = model.set_index("ModelID")
ge = pd.read_csv(GE, index_col=0)
assert ge.shape == (1178, 17916), ge.shape
genes = [c.split(" (")[0] for c in ge.columns]
assert len(set(genes)) == len(genes)
ge.columns = genes
X = ge.values.astype(float)
model_ids = list(ge.index)
print("CRISPR lines:", len(model_ids), "genes:", len(genes), "missing entries:", int(np.isnan(X).sum()))

mut = pd.read_csv(MUT, usecols=["HugoSymbol", "ModelID", "ProteinChange"])
eg = mut[mut["HugoSymbol"] == "EGFR"]
eg_models = set(eg["ModelID"])
act_models = set(eg[eg["ProteinChange"].isin(ACT_PC)]["ModelID"])
seq_models = set(mut["ModelID"].unique())

def in_panel(ids, panel): return [i for i in ids if i in set(panel)]

nsclc_ids = set(model[(model["OncotreeLineage"] == "Lung") &
                      (model["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer")]["ModelID"])
# --- CRISPR groups ---
crispr_nsclc = [m for m in model_ids if m in nsclc_ids]
act_c = [m for m in crispr_nsclc if m in act_models]
broad_c = [m for m in crispr_nsclc if m in eg_models]
wt_c = [m for m in crispr_nsclc if m not in eg_models]
assert len(crispr_nsclc) == 98 and len(act_c) == 6 and len(broad_c) == 21 and len(wt_c) == 77, \
    (len(crispr_nsclc), len(act_c), len(broad_c), len(wt_c))
print("CRISPR NSCLC:", len(crispr_nsclc), "| activating:", len(act_c), "| broad mutant:", len(broad_c), "| WT:", len(wt_c))
print("activating CRISPR:", [(m, md.loc[m, "CellLineName"]) for m in act_c])
print("broad CRISPR:", [(m, md.loc[m, "CellLineName"]) for m in broad_c])

def scan(grp_ids, comp_ids, label):
    grp = np.array([m in set(grp_ids) for m in model_ids])
    comp = np.array([m in set(comp_ids) for m in model_ids])
    pool = grp | comp
    Xp = X[pool]; g = grp[pool]
    m1, m2, n1, n2, p = welch_stats(Xp, g)
    dep1 = np.nansum((Xp[g] < -0.5), axis=0); dep2 = np.nansum((Xp[~g] < -0.5), axis=0)
    f1 = dep1 / n1; f2 = dep2 / n2
    gw_dep = np.nansum((X < -0.5), axis=0) / np.sum(~np.isnan(X), axis=0)
    enrich = m1 - m2
    pan = (f1 >= 0.70) & (f2 >= 0.70)
    q = bh(np.where(np.isnan(p), 1.0, p))
    q = np.where(np.isnan(p), np.nan, q)
    df = pd.DataFrame({"Gene": genes, "GroupMean": m1, "GroupN": n1, "GroupDep": dep1,
                       "GroupDepFrac": f1, "CompMean": m2, "CompN": n2, "CompDep": dep2,
                       "CompDepFrac": f2, "Enrichment": enrich, "FracDiff": f1 - f2,
                       "GenomeWideDepFrac": gw_dep, "WelchP": p, "BH_Q": q,
                       "PanEssential": pan})
    strict = (df["Enrichment"] <= -0.10) & (df["GroupDepFrac"] >= 0.25) & \
             (df["FracDiff"] >= 0.10) & (~df["PanEssential"]) & (df["GenomeWideDepFrac"] <= 0.60)
    relaxed = (df["Enrichment"] <= -0.05) & (df["GroupDepFrac"] >= 0.20) & \
              (df["FracDiff"] >= 0.10) & (~df["PanEssential"]) & (df["GenomeWideDepFrac"] <= 0.60)
    df["StrictBar"] = strict; df["RelaxedBar"] = relaxed
    print(f"[{label}] strict={strict.sum()} relaxed={relaxed.sum()} q<0.05={(df['BH_Q']<0.05).sum()} "
          f"min q={np.nanmin(q):.3g} pan-essential={pan.sum()}")
    return df

prim = scan(act_c, wt_c, "PRIMARY 6 vs 77")
sec = scan(broad_c, wt_c, "SECONDARY 21 vs 77")
prim.to_csv(f"{OUT}/dependency_scan.csv", index=False)
sec.to_csv(f"{OUT}/dependency_scan_secondary.csv", index=False)

# ---------- cell-line RNA ----------
ex = pd.read_csv(EXPR, index_col=0)
print("expr shape:", ex.shape, "missing:", int(np.isnan(ex.values).sum()))
egenes = [c.split(" (")[0] for c in ex.columns]
assert len(set(egenes)) == len(egenes)
ex.columns = egenes
E = ex.values.astype(float)
expr_ids = list(ex.index)
expr_nsclc = [m for m in expr_ids if m in nsclc_ids]
act_e = [m for m in expr_nsclc if m in act_models]
broad_e = [m for m in expr_nsclc if m in eg_models]
wt_e = [m for m in expr_nsclc if m in seq_models and m not in eg_models]
unseq_e = [m for m in expr_nsclc if m not in seq_models]
assert len(expr_nsclc) == 143 and len(act_e) == 9 and len(broad_e) == 28 and len(wt_e) == 106, \
    (len(expr_nsclc), len(act_e), len(broad_e), len(wt_e))
print("RNA NSCLC:", len(expr_nsclc), "| activating:", len(act_e), "| broad:", len(broad_e),
      "| WT:", len(wt_e), "| unsequenced excluded:", len(unseq_e))
print("activating RNA:", [(m, md.loc[m, "CellLineName"]) for m in act_e])

def rna_stats(grp_ids, comp_ids, label):
    grp = np.array([m in set(grp_ids) for m in expr_ids])
    comp = np.array([m in set(comp_ids) for m in expr_ids])
    pool = grp | comp
    Ep = E[pool]; g = grp[pool]
    m1, m2, n1, n2, p = welch_stats(Ep, g)
    # TPM fold change from back-transform of group means of log2(TPM+1)
    fc = (2.0 ** m1 - 1) / (2.0 ** m2 - 1)
    q = bh(np.where(np.isnan(p), 1.0, p)); q = np.where(np.isnan(p), np.nan, q)
    df = pd.DataFrame({"Gene": egenes, "GroupMeanLog2": m1, "GroupN": n1,
                       "CompMeanLog2": m2, "CompN": n2, "DeltaLog2": m1 - m2,
                       "TPM_FC": fc, "WelchP": p, "BH_Q": q})
    print(f"[{label}] RNA genes passing bar: "
          f"{((df['TPM_FC']>=1.5)&(df['DeltaLog2']>=0.5)&(df['WelchP']<0.05)).sum()} of {len(df)}")
    return df

rna_prim = rna_stats(act_e, wt_e, "RNA PRIMARY 9 vs 106")
rna_sec = rna_stats(broad_e, wt_e, "RNA SECONDARY 28 vs 106")
rna_prim.to_csv(f"{OUT}/celline_rna.csv", index=False)
rna_sec.to_csv(f"{OUT}/celline_rna_secondary.csv", index=False)

# spot checks: TEAD programme genes, EGFR itself
for gname in ["EGFR", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "PTK2", "PKN2", "GRB2", "ERBB2", "MET"]:
    r = prim[prim["Gene"] == gname].iloc[0]
    rr = rna_prim[rna_prim["Gene"] == gname].iloc[0]
    print(f"{gname}: dep mean {r['GroupMean']:.3f} ({int(r['GroupDep'])}/{int(r['GroupN'])}) vs "
          f"{r['CompMean']:.3f} ({int(r['CompDep'])}/{int(r['CompN'])}), enrich {r['Enrichment']:.3f}, "
          f"q={r['BH_Q']:.3g} | RNA FC {rr['TPM_FC']:.2f}x dlog2 {rr['DeltaLog2']:.2f} p={rr['WelchP']:.3g}")
# group lists for the report
pd.DataFrame({"ModelID": act_c, "CellLine": [md.loc[m, "CellLineName"] for m in act_c],
              "EGFR_ProteinChange": ["; ".join(sorted(set(eg[eg['ModelID']==m]['ProteinChange'].dropna()))) for m in act_c],
              "Panel": "CRISPR+RNA"}).to_csv(f"{OUT}/activating_lines.csv", index=False)
pd.DataFrame({"ModelID": broad_c, "CellLine": [md.loc[m, "CellLineName"] for m in broad_c],
              "EGFR_ProteinChange": ["; ".join(sorted(set(eg[eg['ModelID']==m]['ProteinChange'].dropna()))) for m in broad_c],
              "Activating": [m in act_models for m in broad_c]}).to_csv(f"{OUT}/broad_mutant_lines.csv", index=False)
print("done")
