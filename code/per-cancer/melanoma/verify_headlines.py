"""Independent verification: recompute headline numbers straight from the raw files,
not from the pipeline's intermediate tables."""
import json, numpy as np, pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/melanoma"
BASE = "/home/hatch/workspace/research_notes"
model = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/Model.csv", low_memory=False)
mel_ids = set(model.loc[model["OncotreePrimaryDisease"] == "Melanoma", "ModelID"])
eff = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/CRISPRGeneEffect.csv").set_index("Unnamed: 0")
mel = [c for c in eff.index if c in mel_ids]; non = [c for c in eff.index if c not in mel_ids]
print("VERIFY lines:", len(mel), "mel /", len(non), "non-mel; total", len(eff))

def check(gene):
    col = [c for c in eff.columns if c.split(" (")[0] == gene][0]
    a = eff.loc[mel, col].dropna(); b = eff.loc[non, col].dropna()
    gw = eff[col].dropna()
    print(f"{gene}: mel mean {a.mean():.4f} dep {(a < -0.5).sum()}/{len(a)} | "
          f"non mean {b.mean():.4f} dep {(b < -0.5).sum()}/{len(b)} | "
          f"enrich {a.mean() - b.mean():.4f} | gw dep frac {(gw < -0.5).mean():.4f}")
for g in ["SOX10", "BRAF", "DUSP4", "MITF", "TEAD1", "WWTR1", "PTK2", "YAP1", "PKN2", "MAPK1"]:
    check(g)

# strict / relaxed counts, recomputed from scratch
def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = ~np.isnan(p)
    o = np.argsort(p[ok]); r = p[ok][o]
    q = np.minimum.accumulate((r * ok.sum() / (np.arange(ok.sum()) + 1))[::-1])[::-1]
    out[np.where(ok)[0][o]] = np.clip(q, 0, 1); return out

enr, frac, depf, gwf, pan = [], [], [], [], []
for col in eff.columns:
    a = eff.loc[mel, col].dropna(); b = eff.loc[non, col].dropna(); gw = eff[col].dropna()
    fa = (a < -0.5).mean(); fb = (b < -0.5).mean()
    enr.append(a.mean() - b.mean()); frac.append(fa - fb); depf.append(fa)
    gwf.append((gw < -0.5).mean()); pan.append(fa >= 0.70 and fb >= 0.70)
enr = np.array(enr); frac = np.array(frac); depf = np.array(depf); gwf = np.array(gwf); pan = np.array(pan)
strict = (enr <= -0.10) & (depf >= 0.25) & (frac >= 0.10) & (~pan) & (gwf <= 0.60)
relax = (enr <= -0.05) & (depf >= 0.20) & (frac >= 0.10) & (~pan) & (gwf <= 0.60)
print("VERIFY strict:", int(strict.sum()), "relaxed:", int(relax.sum()), "pan-essential:", int(pan.sum()))

# TCGA SOX10 recompute from raw JSON
e2s = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/entrez2symbol.json")).items()}
sk = json.load(open(f"{OUT}/downloads/cbioportal_skcm_rna_seq_v2.json"))
oth = []
for s in ["brca_tcga", "coadread_tcga", "gbm_tcga", "ov_tcga", "luad_tcga",
          "ucec_tcga", "kirc_tcga", "hnsc_tcga", "lgg_tcga", "thca_tcga"]:
    oth += json.load(open(f"{OUT}/downloads/cbioportal_{s}_rna_seq_v2.json"))
def st(sid):
    try: return sid.split("-")[3][:2]
    except Exception: return "?"
for g in ["SOX10", "MITF", "WWTR1"]:
    x = np.array([np.log2(d["value"] + 1) for d in sk if e2s.get(d["entrezGeneId"]) == g and st(d["sampleId"]) in ("01", "06")])
    y = np.array([np.log2(d["value"] + 1) for d in oth if e2s.get(d["entrezGeneId"]) == g and st(d["sampleId"]) in ("01", "06")])
    print(f"VERIFY TCGA {g}: SKCM n={len(x)} mean {x.mean():.4f} | other n={len(y)} mean {y.mean():.4f} | log2FC {x.mean()-y.mean():.4f} | p {stats.ttest_ind(x, y, equal_var=False).pvalue:.2e}")

# top-20 score recompute from overexpression.csv components
m = pd.read_csv(f"{OUT}/overexpression.csv")
rel = m[m.RelaxedBar].copy()
def z(s):
    s = s.fillna(0.0); return (s - s.mean()) / s.std(ddof=1)
score = z(-rel.Enrichment) + z(rel.FracDiff) + z(rel.TCGA_RNA_log2FC.fillna(0)) + z(rel.RNA_delta_log2.fillna(0))
top = rel[rel.Pass_OE].assign(Score2=score).sort_values("Score2", ascending=False).head(20)
saved = pd.read_csv(f"{OUT}/_top20_raw.csv")
print("VERIFY top20 genes match saved:", list(top.Gene) == list(saved.Gene))
print("VERIFY top20 scores match (max abs diff):", float(np.abs(top.Score2.values - saved.Score.values).max()))
print("VERIFY top gene:", top.iloc[0].Gene, "score", round(top.iloc[0].Score2, 2))
print("VERIFY strict&OE:", int(((m.StrictBar) & (m.Pass_OE)).sum()), "relaxed&OE:", int(((m.RelaxedBar) & (m.Pass_OE)).sum()))
print("VERIFY missing values in effect matrix:", int(eff.isna().sum().sum()))
print("ALL VERIFY DONE")
