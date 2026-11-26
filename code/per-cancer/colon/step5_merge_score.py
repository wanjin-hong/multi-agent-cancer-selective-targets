"""Step 5: merge DepMap + TCGA PanCan RNA + CPTAC + RPPA; bars; combined score; top 20; Q4.
Welch t-tests; BH correction (dependency: across all 17,916 genes — already in
dependency_scan.csv; TCGA RNA: across fetched genes with data; CPTAC/RPPA: across genes
with data in that assay). NaN p-values excluded from BH input and left NaN.
"""
import json, numpy as np, pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/colon"
entrez2sym = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/entrez2symbol.json")).items()}

scan = pd.read_csv(f"{OUT}/dependency_scan.csv")
rna_cl = pd.read_csv(f"{OUT}/celline_rna.csv")
relaxed = [g.strip() for g in open(f"{OUT}/downloads/relaxed_genes.txt") if g.strip()]
strict = [g.strip() for g in open(f"{OUT}/downloads/strict_genes.txt") if g.strip()]
extras = ["TEAD1", "WWTR1", "PTK2", "PKN2", "WRN"]
fetch_genes = relaxed + [g for g in extras if g not in relaxed]

def welch_from_stats(m1, v1, n1, m2, v2, n2):
    if n1 < 2 or n2 < 2 or v1 is None or v2 is None or np.isnan(v1) or np.isnan(v2):
        return np.nan, np.nan
    se = v1 / n1 + v2 / n2
    if se <= 0: return np.nan, np.nan
    t = (m1 - m2) / np.sqrt(se)
    df = se**2 / ((v1 / n1)**2 / (n1 - 1) + (v2 / n2)**2 / (n2 - 1))
    p = 2 * stats.t.sf(abs(t), df)
    return t, p

def bh(pvals):
    p = np.array(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    mask = ~np.isnan(p)
    if mask.sum():
        out[mask] = stats.false_discovery_control(p[mask], method="bh")
    return out

def stype(sid):
    try: return sid.split("-")[3][:2]
    except Exception: return "?"

# ---------- TCGA PanCan RNA ----------
coad = json.load(open(f"{OUT}/downloads/cbioportal_pca_coadread_rna.json"))
agg = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/pca_comparator_agg.json")).items()}
coad_by_gene = {}
for d in coad:
    if stype(d["sampleId"]) != "01": continue
    coad_by_gene.setdefault(d["entrezGeneId"], []).append(np.log2(d["value"] + 1))
print("COADREAD PanCan RNA genes with data:", len(coad_by_gene))
print("Comparator genes aggregated:", len(agg))
missing_rna = [g for g in fetch_genes if entrez2sym and g not in
               [entrez2sym.get(e, "") for e in set(coad_by_gene) | set(agg)]]
print("Fetched genes with NO TCGA RNA data anywhere:", missing_rna)

tcga_rows = []
for g in fetch_genes:
    eid = [k for k, v in entrez2sym.items() if v == g]
    if not eid: continue
    eid = eid[0]
    x = np.array(coad_by_gene.get(eid, []), dtype=float)
    a = agg.get(eid)
    if len(x) == 0 or a is None:
        tcga_rows.append(dict(Gene=g, TCGA_RNA_log2FC=np.nan, TCGA_RNA_p=np.nan,
                              TCGA_RNA_n_coad=0, TCGA_RNA_n_other=0,
                              TCGA_RNA_mean_coad=np.nan, TCGA_RNA_mean_other=np.nan))
        continue
    s, ss, n2 = a
    m1, m2 = x.mean(), s / n2
    v1 = x.var(ddof=1) if len(x) > 1 else np.nan
    v2 = (ss - s * s / n2) / (n2 - 1) if n2 > 1 else np.nan
    _, p = welch_from_stats(m1, v1, len(x), m2, v2, n2)
    tcga_rows.append(dict(Gene=g, TCGA_RNA_log2FC=m1 - m2, TCGA_RNA_p=p,
                          TCGA_RNA_n_coad=len(x), TCGA_RNA_n_other=n2,
                          TCGA_RNA_mean_coad=m1, TCGA_RNA_mean_other=m2))
tcga = pd.DataFrame(tcga_rows)
tcga["TCGA_RNA_q"] = bh(tcga["TCGA_RNA_p"].values)
tcga.to_csv(f"{OUT}/tcga_rna.csv", index=False)
print("TCGA RNA bar pass (log2FC>=0.5 & p<0.05):",
      int(((tcga["TCGA_RNA_log2FC"] >= 0.5) & (tcga["TCGA_RNA_p"] < 0.05)).sum()), "of", len(tcga))
print(tcga.sort_values("TCGA_RNA_log2FC", ascending=False).head(12).to_string(index=False))

# ---------- Protein: CPTAC + RPPA (coadread vs pooled BRCA+OV) ----------
def protein_stats(tag):
    data = json.load(open(f"{OUT}/downloads/cbioportal_pca_{tag}.json"))
    coad_vals, comp_vals = {}, {}
    for st, recs in data.items():
        for d in recs:
            if stype(d["sampleId"]) != "01": continue
            (coad_vals if st.startswith("coadread") else comp_vals).setdefault(d["entrezGeneId"], []).append(d["value"])
    rows = []
    for g in fetch_genes:
        eid = [k for k, v in entrez2sym.items() if v == g]
        if not eid: continue
        eid = eid[0]
        x = np.array(coad_vals.get(eid, []), dtype=float)
        y = np.array(comp_vals.get(eid, []), dtype=float)
        if len(x) < 2 or len(y) < 2:
            rows.append(dict(Gene=g, delta=np.nan, p=np.nan, n_coad=len(x), n_other=len(y)))
            continue
        _, p = stats.ttest_ind(x, y, equal_var=False)
        rows.append(dict(Gene=g, delta=x.mean() - y.mean(), p=p, n_coad=len(x), n_other=len(y)))
    df = pd.DataFrame(rows)
    df["q"] = bh(df["p"].values)
    df.columns = ["Gene", f"{tag}_delta", f"{tag}_p", f"{tag}_n_coad", f"{tag}_n_other", f"{tag}_q"]
    return df

cptac = protein_stats("cptac")
# --- CPTAC SCALE CHECK (added after first run) ---
# PanCancer Atlas coadread protein_quantification values are on an uncentred intensity
# scale (mean ~19.2, range 0–34.5) while BRCA/OV are log-ratio scale (mean ~0, range
# -5 to +4). Cross-study CPTAC deltas (e.g. +28) are therefore scale artefacts, NOT
# biology. CPTAC is EXCLUDED from the protein bar and from the combined score.
# tcga_cptac.csv is rewritten below to per-group means only, with no delta/p, as the
# audit record of the attempted (invalid) comparison.
_raw = json.load(open(f"{OUT}/downloads/cbioportal_pca_cptac.json"))
_cptac_means = []
for g in fetch_genes:
    eid = [k for k, v in entrez2sym.items() if v == g]
    if not eid: continue
    eid = eid[0]
    _row = {"Gene": g, "Comparable_across_studies": "No - scale mismatch (coadread intensity scale vs BRCA/OV log-ratio scale; see report)"}
    for _st, _key in [("coadread_tcga_pan_can_atlas_2018", "CPTAC_mean_coadread"),
                      ("brca_tcga_pan_can_atlas_2018", "CPTAC_mean_BRCA"),
                      ("ov_tcga_pan_can_atlas_2018", "CPTAC_mean_OV")]:
        _vals = [x["value"] for x in _raw[_st] if x["entrezGeneId"] == eid and stype(x["sampleId"]) == "01"]
        _row[_key] = float(np.mean(_vals)) if _vals else np.nan
        _row[_key.replace("_mean_", "_n_")] = len(_vals)
    _cptac_means.append(_row)
pd.DataFrame(_cptac_means).to_csv(f"{OUT}/tcga_cptac.csv", index=False)
print("CPTAC EXCLUDED from protein bar (scale mismatch). Per-group means saved to tcga_cptac.csv.")
print("CPTAC genes with coadread data:", int(pd.DataFrame(_cptac_means)["CPTAC_mean_coadread"].notna().sum()))
cptac = pd.DataFrame({"Gene": fetch_genes,
                      "cptac_delta": np.nan, "cptac_p": np.nan,
                      "cptac_n_coad": np.nan, "cptac_n_other": np.nan, "cptac_q": np.nan})
rppa = protein_stats("rppa"); rppa.to_csv(f"{OUT}/tcga_rppa.csv", index=False)
print("RPPA genes with data:", int(rppa["rppa_p"].notna().sum()),
      "bar pass:", int(((rppa["rppa_delta"] >= 0.5) & (rppa["rppa_p"] < 0.05)).sum()))

# ---------- Merge: overexpression.csv (all 17,916 genes) ----------
ov = scan.merge(rna_cl[["Gene", "Colon_RNA_mean_log2TPM1", "NonColon_RNA_mean_log2TPM1",
                        "RNA_delta_log2", "RNA_TPM_FC", "RNA_p_welch", "RNA_q_BH"]],
                on="Gene", how="left")
ov = ov.merge(tcga[["Gene", "TCGA_RNA_log2FC", "TCGA_RNA_p", "TCGA_RNA_q",
                    "TCGA_RNA_n_coad", "TCGA_RNA_n_other"]], on="Gene", how="left")
ov = ov.merge(cptac[["Gene", "cptac_delta", "cptac_p", "cptac_q", "cptac_n_coad", "cptac_n_other"]],
              on="Gene", how="left")
ov = ov.merge(rppa[["Gene", "rppa_delta", "rppa_p", "rppa_q", "rppa_n_coad", "rppa_n_other"]],
              on="Gene", how="left")
ov.columns = [c.replace("cptac_", "CPTAC_").replace("rppa_", "RPPA_") for c in ov.columns]
ov.to_csv(f"{OUT}/overexpression.csv", index=False)
print("overexpression.csv rows:", len(ov), "| genes with no cell-line RNA:",
      int(ov["RNA_delta_log2"].isna().sum()))

# ---------- Bars & intersection ----------
def passes_celline(df):
    return (df["RNA_TPM_FC"] >= 1.5) & (df["RNA_delta_log2"] >= 0.5) & (df["RNA_p_welch"] < 0.05)
def passes_tcga(df):
    return (df["TCGA_RNA_log2FC"] >= 0.5) & (df["TCGA_RNA_p"] < 0.05)
def passes_protein(df):
    # RPPA only: CPTAC cross-study comparison is invalid (scale mismatch) — see above.
    return ((df["RPPA_delta"] >= 0.5) & (df["RPPA_p"] < 0.05))
ov["Pass_celline_RNA"] = passes_celline(ov).fillna(False)
ov["Pass_TCGA_RNA"] = passes_tcga(ov).fillna(False)
ov["Pass_protein"] = passes_protein(ov).fillna(False)
ov["Pass_overexpression_any"] = ov["Pass_celline_RNA"] | ov["Pass_TCGA_RNA"] | ov["Pass_protein"]
ov.to_csv(f"{OUT}/overexpression.csv", index=False)

rel = ov[ov["Gene"].isin(relaxed)].copy()
strc = ov[ov["Gene"].isin(strict)].copy()
print("Strict bar:", len(strc), "of which pass over-expression bar:", int(strc["Pass_overexpression_any"].sum()))
print("Relaxed bar:", len(rel), "of which pass over-expression bar:", int(rel["Pass_overexpression_any"].sum()))
print("Strict∩overexp genes:", sorted(strc.loc[strc["Pass_overexpression_any"], "Gene"].tolist()))
print("Relaxed∩overexp genes:", sorted(rel.loc[rel["Pass_overexpression_any"], "Gene"].tolist()))

# ---------- Combined score (across relaxed-bar genes) ----------
def z(s):
    s = s.astype(float)
    m, sd = np.nanmean(s), np.nanstd(s, ddof=1)
    return (s - m) / sd
rel["z_neg_enrich"] = z(-rel["Enrichment"])
rel["z_fracdiff"] = z(rel["FracDiff"])
rel["z_tcga"] = z(rel["TCGA_RNA_log2FC"].fillna(0))  # missing TCGA -> z of 0 contribution handled below
# TNBC rule: missing TCGA or cell-line values are z-scored as 0. Implement exactly:
rel.loc[rel["TCGA_RNA_log2FC"].isna(), "z_tcga"] = 0.0
rel["z_celline"] = z(rel["RNA_delta_log2"])
rel.loc[rel["RNA_delta_log2"].isna(), "z_celline"] = 0.0
rel["CombinedScore"] = (rel["z_neg_enrich"] + rel["z_fracdiff"] + rel["z_tcga"] + rel["z_celline"]
                        + np.where(rel["Pass_protein"], 1.0, 0.0))
rel["Tier"] = np.where(rel["Gene"].isin(strict), "A (strict dep bar)", "B (relaxed dep bar only)")
cand = rel[rel["Pass_overexpression_any"]].sort_values("CombinedScore", ascending=False).copy()
cand.to_csv(f"{OUT}/candidates_both.csv", index=False)
print("\nCandidates (relaxed + over-expression), ranked:")
print(cand[["Gene", "Tier", "Colon_mean", "Colon_dep_frac", "NonColon_dep_frac", "Enrichment",
            "RNA_TPM_FC", "TCGA_RNA_log2FC", "CPTAC_delta", "CombinedScore"]].to_string(index=False))
print("\nRanks 21+ of candidates:")
print("total candidates:", len(cand))
