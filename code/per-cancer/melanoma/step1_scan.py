"""Melanoma selective-targets pipeline, step 1: genome-wide dependency scan + cell-line RNA.
DepMap 24Q4. Group = OncotreePrimaryDisease == 'Melanoma' (Model.csv).
Comparator = all non-melanoma lines in the same matrix.
Dependency = Chronos < -0.5. NaN-aware denominators, Welch t, BH.
"""
import numpy as np, pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/melanoma"
BASE = "/home/hatch/workspace/research_notes"

model = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/Model.csv", low_memory=False)
mel_ids = set(model.loc[model["OncotreePrimaryDisease"] == "Melanoma", "ModelID"])
print("melanoma Model lines:", len(mel_ids))

# ---------- CRISPR ----------
eff = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/CRISPRGeneEffect.csv")
print("effect shape:", eff.shape, "missing:", int(eff.iloc[:, 1:].isna().sum().sum()))
# first col is ModelID (check)
print("first col:", eff.columns[0])
eff = eff.set_index(eff.columns[0])
genes = [c.split(" (")[0] for c in eff.columns]
print("n genes:", len(set(genes)), "n cols:", len(eff.columns))
mel_cols = [c for c in eff.index if c in mel_ids]
non_cols = [c for c in eff.index if c not in mel_ids]
print("melanoma CRISPR lines:", len(mel_cols), "non-melanoma CRISPR lines:", len(non_cols))

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    ok = ~np.isnan(p)
    q = np.empty(ok.sum()); o = np.argsort(p[ok]); ranked = p[ok][o]
    q = ranked * ok.sum() / (np.arange(ok.sum()) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[np.where(ok)[0][o]] = np.clip(q, 0, 1)
    return out

rows = []
for col in eff.columns:
    g = col.split(" (")[0]
    a = eff.loc[mel_cols, col].to_numpy(dtype=float)
    b = eff.loc[non_cols, col].to_numpy(dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    gw = eff[col].to_numpy(dtype=float); gw = gw[~np.isnan(gw)]
    ma = a.mean() if len(a) else np.nan; mb = b.mean() if len(b) else np.nan
    da = (a < -0.5).sum(); db = (b < -0.5).sum()
    fa = da / len(a) if len(a) else np.nan; fb = db / len(b) if len(b) else np.nan
    gwf = (gw < -0.5).mean() if len(gw) else np.nan
    if len(a) >= 3 and len(b) >= 3:
        p = stats.ttest_ind(a, b, equal_var=False).pvalue
    else:
        p = np.nan
    pan = (fa >= 0.70 and fb >= 0.70) if not (np.isnan(fa) or np.isnan(fb)) else False
    rows.append(dict(Gene=g, Column=col, Mel_mean=ma, Mel_dep_n=da, Mel_N_scored=len(a),
                     Mel_dep_frac=fa, NonMel_mean=mb, NonMel_dep_n=db, NonMel_N_scored=len(b),
                     NonMel_dep_frac=fb, Enrichment=ma - mb, FracDiff=fa - fb,
                     GenomeWide_dep_frac=gwf, p_welch=p, PanEssential=pan))
scan = pd.DataFrame(rows)
scan["q_BH"] = bh(scan["p_welch"].values)
scan["EnrichRank"] = scan["Enrichment"].rank(method="min").astype(int)
scan = scan.sort_values("Enrichment")
scan.to_csv(f"{OUT}/dependency_scan.csv", index=False)
print("scan done", scan.shape)
print("q<0.05:", int((scan.q_BH < 0.05).sum()), "q<0.1:", int((scan.q_BH < 0.1).sum()),
      "q<0.2:", int((scan.q_BH < 0.2).sum()))
print("pan-essential:", int(scan.PanEssential.sum()))
strict = scan[(scan.Enrichment <= -0.10) & (scan.Mel_dep_frac >= 0.25) & (scan.FracDiff >= 0.10)
              & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relax = scan[(scan.Enrichment <= -0.05) & (scan.Mel_dep_frac >= 0.20) & (scan.FracDiff >= 0.10)
             & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print("strict:", len(strict), "relaxed:", len(relax))
print(scan.head(15)[["Gene","Mel_mean","Mel_dep_n","Mel_N_scored","NonMel_mean","NonMel_dep_n","NonMel_N_scored","Enrichment","FracDiff","q_BH","PanEssential"]].to_string())
for g in ["TEAD1","YAP1","WWTR1","PTK2","PKN2"]:
    r = scan[scan.Gene == g]
    print(g, r[["EnrichRank","Mel_mean","Mel_dep_frac","NonMel_mean","NonMel_dep_frac","Enrichment","FracDiff","GenomeWide_dep_frac","q_BH"]].to_dict("records"))
scan.to_pickle(f"{OUT}/_scan.pkl")
strict[["Gene"]].to_csv(f"{OUT}/_strict_genes.csv", index=False)
relax[["Gene"]].to_csv(f"{OUT}/_relax_genes.csv", index=False)

# ---------- Cell-line RNA ----------
expr = pd.read_csv(f"{BASE}/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
print("expr shape:", expr.shape)
expr = expr.set_index(expr.columns[0])
mel_e = [c for c in expr.index if c in mel_ids]
non_e = [c for c in expr.index if c not in mel_ids]
print("melanoma RNA lines:", len(mel_e), "non-melanoma RNA lines:", len(non_e))
erows = []
for col in expr.columns:
    g = col.split(" (")[0]
    a = expr.loc[mel_e, col].to_numpy(dtype=float)
    b = expr.loc[non_e, col].to_numpy(dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    ma = a.mean(); mb = b.mean(); d = ma - mb
    fc = (2**ma - 1) / (2**mb - 1) if (2**mb - 1) > 0 else np.nan
    p = stats.ttest_ind(a, b, equal_var=False).pvalue if len(a) >= 3 and len(b) >= 3 else np.nan
    erows.append(dict(Gene=g, Column_cl=col, Mel_RNA_mean_log2TPM1=ma, NonMel_RNA_mean_log2TPM1=mb,
                      RNA_delta_log2=d, RNA_TPM_FC=fc, RNA_p_welch=p,
                      Mel_RNA_n=len(a), NonMel_RNA_n=len(b)))
rna = pd.DataFrame(erows)
rna["RNA_q_BH"] = bh(rna["RNA_p_welch"].values)
rna.to_csv(f"{OUT}/celline_rna.csv", index=False)
print("rna done", rna.shape, "pass RNA bar:",
      int(((rna.RNA_TPM_FC >= 1.5) & (rna.RNA_delta_log2 >= 0.5) & (rna.RNA_p_welch < 0.05)).sum()))
for g in ["TEAD1","YAP1","WWTR1","PTK2","PKN2"]:
    print(g, rna[rna.Gene == g][["RNA_TPM_FC","RNA_delta_log2","RNA_p_welch"]].to_dict("records"))
rna.to_pickle(f"{OUT}/_rna.pkl")
print("STEP1 OK")
