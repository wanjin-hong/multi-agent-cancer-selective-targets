"""Step 3: TCGA stats + merge + bars + combined score + top 20 + CSVs.
TCGA RNA: log2(RSEM+1); log2FC = mean SKCM tumours - mean pooled other-cancer tumours.
Tumour samples only (TCGA barcode type 01 primary, 06 metastatic; normals 11 excluded).
RPPA: values already log2-scale; delta = SKCM mean - other-cancer mean.
BH across the fetched genes. Protein bar: RPPA delta >= 0.5 and p < 0.05
(CPTAC unavailable for SKCM — no CPTAC columns, recorded in report).
"""
import json, numpy as np, pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/melanoma"
OTHERS = ["brca_tcga", "coadread_tcga", "gbm_tcga", "ov_tcga", "luad_tcga",
          "ucec_tcga", "kirc_tcga", "hnsc_tcga", "lgg_tcga", "thca_tcga"]
entrez2sym = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/entrez2symbol.json")).items()}

def bh(p):
    p = np.asarray(p, dtype=float); out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    q = p[ok] * ok.sum() / (np.arange(ok.sum()) + 1)
    q = np.minimum.accumulate(q[np.argsort(p[ok])][::-1])[::-1] if False else None
    o = np.argsort(p[ok]); ranked = p[ok][o]
    qq = ranked * ok.sum() / (np.arange(ok.sum()) + 1)
    qq = np.minimum.accumulate(qq[::-1])[::-1]
    out[np.where(ok)[0][o]] = np.clip(qq, 0, 1)
    return out

def stype(sid):
    try: return sid.split("-")[3][:2]
    except Exception: return "?"

def load_json(tag):
    return json.load(open(f"{OUT}/downloads/cbioportal_{tag}.json"))

def to_frame(records, value_transform):
    rows = []
    for d in records:
        g = entrez2sym.get(d["entrezGeneId"])
        if g is None: continue
        rows.append((g, d["sampleId"], value_transform(d["value"])))
    return pd.DataFrame(rows, columns=["Gene", "Sample", "Value"])

# ---- RNA ----
sk = to_frame(load_json("skcm_rna_seq_v2"), lambda v: np.log2(v + 1))
print("SKCM RNA sample types:", pd.Series([stype(s) for s in sk.Sample.unique()]).value_counts().to_dict())
sk_t = sk[sk.Sample.map(stype).isin(["01", "06"])]
print("SKCM tumour RNA samples:", sk_t.Sample.nunique(), "primary:", sk_t[sk_t.Sample.map(stype) == "01"].Sample.nunique(),
      "metastatic:", sk_t[sk_t.Sample.map(stype) == "06"].Sample.nunique())
oth_frames = []
for s in OTHERS:
    f = to_frame(load_json(f"{s}_rna_seq_v2"), lambda v: np.log2(v + 1))
    f = f[f.Sample.map(stype).isin(["01", "06"])]
    oth_frames.append(f)
    print(s, "tumour RNA samples:", f.Sample.nunique())
oth = pd.concat(oth_frames)
print("pooled other-cancer tumour RNA samples:", oth.Sample.nunique())

def per_gene_stats(a, b):
    out = []
    for g in sorted(set(a.Gene) | set(b.Gene)):
        x = a.loc[a.Gene == g, "Value"].to_numpy(float)
        y = b.loc[b.Gene == g, "Value"].to_numpy(float)
        if len(x) == 0 or len(y) == 0:
            out.append(dict(Gene=g, mean_a=np.nan, mean_b=np.nan, delta=np.nan, p=np.nan, n_a=len(x), n_b=len(y)))
            continue
        p = stats.ttest_ind(x, y, equal_var=False).pvalue
        out.append(dict(Gene=g, mean_a=x.mean(), mean_b=y.mean(), delta=x.mean() - y.mean(),
                        p=p, n_a=len(x), n_b=len(y)))
    return pd.DataFrame(out)

rna = per_gene_stats(sk_t, oth)
rna["q_BH"] = bh(rna["p"].values)
rna = rna.rename(columns={"mean_a": "TCGA_SKCM_mean_log2RSEM1", "mean_b": "TCGA_Other_mean_log2RSEM1",
                          "delta": "TCGA_RNA_log2FC", "p": "TCGA_RNA_p", "q_BH": "TCGA_RNA_q",
                          "n_a": "TCGA_RNA_n_SKCM", "n_b": "TCGA_RNA_n_Other"})
rna.to_csv(f"{OUT}/tcga_rna.csv", index=False)
print("TCGA RNA genes:", len(rna), "pass bar (log2FC>=0.5 & p<0.05):",
      int(((rna.TCGA_RNA_log2FC >= 0.5) & (rna.TCGA_RNA_p < 0.05)).sum()))

# secondary: SKCM primary vs metastatic
prim = sk_t[sk_t.Sample.map(stype) == "01"]; met = sk_t[sk_t.Sample.map(stype) == "06"]
pm = per_gene_stats(met, prim).rename(columns={"delta": "Met_vs_Prim_log2FC", "p": "Met_vs_Prim_p",
                                               "n_a": "n_met", "n_b": "n_prim"})
pm["Met_vs_Prim_q"] = bh(pm["Met_vs_Prim_p"].values)
pm.to_csv(f"{OUT}/tcga_skcm_primary_vs_metastatic.csv", index=False)
print("primary-vs-met done; genes with |log2FC|>=0.5 & q<0.05:",
      int(((pm.Met_vs_Prim_log2FC.abs() >= 0.5) & (pm.Met_vs_Prim_q < 0.05)).sum()))

# ---- RPPA ----
skp = to_frame(load_json("skcm_rppa"), lambda v: float(v))
print("SKCM RPPA genes:", sorted(skp.Gene.unique()), "samples:", skp.Sample.nunique())
othp = pd.concat([to_frame(load_json(f"{s}_rppa"), lambda v: float(v)) for s in OTHERS])
othp = othp[othp.Sample.map(stype).isin(["01", "06"])] if othp.Sample.map(stype).isin(["01", "06"]).any() else othp
# RPPA sampleIds may lack barcode type info in some studies; keep all if barcode parse fails
print("pooled other RPPA samples:", othp.Sample.nunique(), "genes:", sorted(othp.Gene.unique()))
rp = per_gene_stats(skp, othp)
rp["q_BH"] = bh(rp["p"].values)
rp = rp.rename(columns={"mean_a": "RPPA_SKCM_mean", "mean_b": "RPPA_Other_mean",
                        "delta": "RPPA_delta", "p": "RPPA_p", "q_BH": "RPPA_q",
                        "n_a": "RPPA_n_SKCM", "n_b": "RPPA_n_Other"})
rp.to_csv(f"{OUT}/tcga_rppa.csv", index=False)
print("RPPA pass bar (delta>=0.5 & p<0.05):", int(((rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)).sum()))

# ---- Merge ----
scan = pd.read_pickle(f"{OUT}/_scan.pkl")
crna = pd.read_pickle(f"{OUT}/_rna.pkl")
m = scan.merge(crna, on=["Gene"], how="left", suffixes=("", "_rna"))
m = m.merge(rna[["Gene", "TCGA_RNA_log2FC", "TCGA_RNA_p", "TCGA_RNA_q", "TCGA_RNA_n_SKCM", "TCGA_RNA_n_Other"]],
            on="Gene", how="left")
m = m.merge(rp[["Gene", "RPPA_delta", "RPPA_p", "RPPA_q", "RPPA_n_SKCM", "RPPA_n_Other"]],
            on="Gene", how="left")
# over-expression bars
m["Pass_CL_RNA"] = (m.RNA_TPM_FC >= 1.5) & (m.RNA_delta_log2 >= 0.5) & (m.RNA_p_welch < 0.05)
m["Pass_TCGA_RNA"] = (m.TCGA_RNA_log2FC >= 0.5) & (m.TCGA_RNA_p < 0.05)
m["Pass_Protein"] = (m.RPPA_delta >= 0.5) & (m.RPPA_p < 0.05)  # CPTAC unavailable
m["Pass_OE"] = m.Pass_CL_RNA | m.Pass_TCGA_RNA | m.Pass_Protein
m["StrictBar"] = (m.Enrichment <= -0.10) & (m.Mel_dep_frac >= 0.25) & (m.FracDiff >= 0.10) & (~m.PanEssential) & (m.GenomeWide_dep_frac <= 0.60)
m["RelaxedBar"] = (m.Enrichment <= -0.05) & (m.Mel_dep_frac >= 0.20) & (m.FracDiff >= 0.10) & (~m.PanEssential) & (m.GenomeWide_dep_frac <= 0.60)
print("strict:", int(m.StrictBar.sum()), "relaxed:", int(m.RelaxedBar.sum()))
print("strict & OE:", int((m.StrictBar & m.Pass_OE).sum()), "relaxed & OE:", int((m.RelaxedBar & m.Pass_OE).sum()))
print("strict & CL-RNA:", int((m.StrictBar & m.Pass_CL_RNA).sum()),
      "strict & TCGA:", int((m.StrictBar & m.Pass_TCGA_RNA).sum()),
      "strict & protein:", int((m.StrictBar & m.Pass_Protein).sum()))

# combined score across relaxed-bar genes
rel = m[m.RelaxedBar].copy()
def z(s):
    s = s.fillna(0.0)
    sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd and sd > 0 else s * 0
rel["Score"] = (z(-rel.Enrichment) + z(rel.FracDiff) + z(rel.TCGA_RNA_log2FC.fillna(0)) +
                z(rel.RNA_delta_log2.fillna(0)) + rel.Pass_Protein.fillna(False).astype(float) * 1.0)
rel = rel.sort_values("Score", ascending=False)
rel["Tier"] = np.where(rel.StrictBar, "A (strict dep bar)", "B (relaxed dep bar only)")
top = rel[rel.Pass_OE].head(20).copy()
print("top20 (relaxed & OE) n =", len(top))
print(top[["Gene", "Tier", "Mel_mean", "Mel_dep_frac", "NonMel_dep_frac", "Enrichment",
           "RNA_TPM_FC", "TCGA_RNA_log2FC", "RPPA_delta", "Score"]].to_string())
print("ranks 21-25 of relaxed&OE:")
print(rel[rel.Pass_OE].iloc[20:25][["Gene", "Tier", "Enrichment", "TCGA_RNA_log2FC", "Score"]].to_string())
# TEAD axis
for g in ["TEAD1", "YAP1", "WWTR1", "PTK2", "PKN2", "SOX10", "BRAF"]:
    print(g, m[m.Gene == g][["EnrichRank", "Enrichment", "Mel_dep_frac", "NonMel_dep_frac",
          "GenomeWide_dep_frac", "RNA_TPM_FC", "RNA_delta_log2", "TCGA_RNA_log2FC", "TCGA_RNA_q",
          "RPPA_delta", "Pass_OE", "StrictBar", "RelaxedBar"]].to_dict("records"))
m.to_csv(f"{OUT}/overexpression.csv", index=False)
rel.to_csv(f"{OUT}/_relaxed_scored.csv", index=False)
top.to_csv(f"{OUT}/_top20_raw.csv", index=False)
print("STEP3 OK")
