"""Colon selective-targets pipeline, step 1: DepMap 24Q4 dependency scan + cell-line RNA.
Group = Colorectal Adenocarcinoma lines (OncotreePrimaryDisease) in Model.csv.
Comparator = all other lines. Dependency = Chronos < -0.5.
"""
import numpy as np, pandas as pd
from scipy import stats
import os

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/colon"
BASE = "/home/hatch/workspace/research_notes"

model = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/Model.csv", low_memory=False)
colon_ids = set(model.loc[model["OncotreePrimaryDisease"] == "Colorectal Adenocarcinoma", "ModelID"])
print("Model.csv rows:", len(model))
print("Colorectal Adenocarcinoma lines in Model.csv:", len(colon_ids))
print("  colon subtype counts:", model.loc[model["OncotreePrimaryDisease"]=="Colorectal Adenocarcinoma","OncotreeSubtype"].value_counts().to_dict())
print("  bowel non-colorectal excluded:", model.loc[(model["OncotreeLineage"]=="Bowel") & (model["OncotreePrimaryDisease"]!="Colorectal Adenocarcinoma"), ["CellLineName","OncotreePrimaryDisease"]].values.tolist())

# ---- Dependency scan ----
eff = pd.read_csv(f"{BASE}/vrk1-codependency-20261108/CRISPRGeneEffect.csv")
ids = eff["Unnamed: 0"]
print("Effect matrix shape:", eff.shape, "missing entries:", int(eff.iloc[:,1:].isna().sum().sum()))
is_colon = ids.isin(colon_ids).values
print("CRISPR lines:", len(ids), "colon with CRISPR:", int(is_colon.sum()), "comparator:", int((~is_colon).sum()))
# sanity: model ids not in effect matrix
print("Colon ModelIDs missing from CRISPR matrix:", len(colon_ids - set(ids)))

def gene_sym(col):
    return col.split(" (")[0]

rows = []
cols = eff.columns[1:]
for col in cols:
    v = eff[col].values.astype(float)
    g = v[is_colon]; c = v[~is_colon]
    g = g[~np.isnan(g)]; c = c[~np.isnan(c)]
    gm, cm = g.mean(), c.mean()
    gdep = (g < -0.5).sum(); cdep = (c < -0.5).sum()
    gw = eff[col].values.astype(float); gw = gw[~np.isnan(gw)]
    gwfrac = (gw < -0.5).mean()
    if len(g) >= 3 and len(c) >= 3:
        _, p = stats.ttest_ind(g, c, equal_var=False)
    else:
        p = np.nan
    rows.append(dict(Gene=gene_sym(col), Column=col,
        Colon_mean=gm, Colon_dep_n=gdep, Colon_N_scored=len(g), Colon_dep_frac=gdep/len(g),
        NonColon_mean=cm, NonColon_dep_n=cdep, NonColon_N_scored=len(c), NonColon_dep_frac=cdep/len(c),
        Enrichment=gm-cm, FracDiff=gdep/len(g)-cdep/len(c),
        GenomeWide_dep_frac=gwfrac, p_welch=p,
        PanEssential=(gdep/len(g) >= 0.70 and cdep/len(c) >= 0.70)))
scan = pd.DataFrame(rows)
scan["q_BH"] = bh = stats.false_discovery_control(scan["p_welch"].values, method="bh")
scan["EnrichRank"] = scan["Enrichment"].rank(method="min").astype(int)
scan = scan.sort_values("Enrichment")
scan.to_csv(f"{OUT}/dependency_scan.csv", index=False)
print("Genes scanned:", len(scan))
print("Pan-essential (dep>=70% both):", int(scan["PanEssential"].sum()))
print("q<0.05:", int((scan["q_BH"]<0.05).sum()), " q<0.1:", int((scan["q_BH"]<0.1).sum()))

def bar(df, enr, frac, diff):
    return df[(df["Enrichment"]<=enr) & (df["Colon_dep_frac"]>=frac) & (df["FracDiff"]>=diff)
              & (~df["PanEssential"]) & (df["GenomeWide_dep_frac"]<=0.60)]
strict = bar(scan, -0.10, 0.25, 0.10)
relaxed = bar(scan, -0.05, 0.20, 0.10)
print("Strict bar genes:", len(strict), " Relaxed bar genes:", len(relaxed))
strict["Gene"].to_csv(f"{OUT}/downloads/strict_genes.txt", index=False, header=False)
relaxed["Gene"].to_csv(f"{OUT}/downloads/relaxed_genes.txt", index=False, header=False)
print("\nTop 20 by enrichment (non-pan-essential):")
print(scan[~scan["PanEssential"]].head(20)[["Gene","Colon_mean","Colon_dep_n","Colon_N_scored","Colon_dep_frac","NonColon_mean","NonColon_dep_frac","Enrichment","FracDiff","GenomeWide_dep_frac","q_BH"]].to_string(index=False))
print("\nTop 20 by enrichment (all):")
print(scan.head(20)[["Gene","Colon_mean","Colon_dep_frac","NonColon_mean","NonColon_dep_frac","Enrichment","FracDiff","GenomeWide_dep_frac","q_BH","PanEssential"]].to_string(index=False))

# ---- Cell-line RNA genome-wide ----
expr = pd.read_csv(f"{BASE}/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
eids = expr["Unnamed: 0"]
print("\nExpression matrix shape:", expr.shape)
e_is_colon = eids.isin(colon_ids).values
print("Expression lines:", len(eids), "colon with RNA:", int(e_is_colon.sum()), "comparator:", int((~e_is_colon).sum()))
print("Colon ModelIDs missing from RNA matrix:", len(colon_ids - set(eids)))
rrows = []
for col in expr.columns[1:]:
    v = expr[col].values.astype(float)
    g = v[e_is_colon]; c = v[~e_is_colon]
    g = g[~np.isnan(g)]; c = c[~np.isnan(c)]
    gm, cm = g.mean(), c.mean()
    d = gm - cm
    fc = (2**gm - 1) / (2**cm - 1) if (2**cm - 1) != 0 else np.nan
    _, p = stats.ttest_ind(g, c, equal_var=False)
    rrows.append(dict(Gene=gene_sym(col), Column=col,
        Colon_RNA_mean_log2TPM1=gm, NonColon_RNA_mean_log2TPM1=cm,
        RNA_delta_log2=d, RNA_TPM_FC=fc, RNA_p_welch=p,
        Colon_RNA_n=len(g), NonColon_RNA_n=len(c)))
rna = pd.DataFrame(rrows)
rna["RNA_q_BH"] = stats.false_discovery_control(rna["RNA_p_welch"].values, method="bh")
rna.to_csv(f"{OUT}/celline_rna.csv", index=False)
print("RNA genes:", len(rna))
print("Cell-line RNA bar pass (FC>=1.5 & d>=0.5 & p<0.05):",
      int(((rna["RNA_TPM_FC"]>=1.5)&(rna["RNA_delta_log2"]>=0.5)&(rna["RNA_p_welch"]<0.05)).sum()))
print("\nTop 15 cell-line RNA over-expressed (by delta log2):")
print(rna.sort_values("RNA_delta_log2", ascending=False).head(15)[["Gene","Colon_RNA_mean_log2TPM1","NonColon_RNA_mean_log2TPM1","RNA_delta_log2","RNA_TPM_FC","RNA_p_welch"]].to_string(index=False))
