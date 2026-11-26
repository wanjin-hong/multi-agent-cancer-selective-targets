"""Ovarian selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: Ovarian Epithelial Tumor lines
(OncotreePrimaryDisease == 'Ovarian Epithelial Tumor', 71 lines in
Model.csv, 58 with CRISPR data). Excludes the 1 germ-cell line (PA-1, Mixed
Germ Cell Tumor) and 1 'Ovarian Cancer, Other' line (HSKTC), mirroring the
PDAC/HNSC exclusions of non-epithelial lines.
Comparator: all 1,120 non-ovarian lines in the effect matrix (the matrix
contains no non-cancerous lines at all, so the comparator is entirely
cancerous).
Same-patient pairs in the group: Hey-A8 / Hey, OC 314 / OC 316, PEO1 / PEO4,
NIH:OVCAR-3 / ONCO-DG-1 (all four pairs in the CRISPR panel — effective
CRISPR n ~54, not 58).
Subtypes in the CRISPR panel: High-Grade Serous 18, Serous 16, Clear Cell 8,
Endometrioid 6, Mucinous 5, Small Cell 3, Brenner 1, Mixed 1. No HGSOC-only
split as primary (the 'Serous' annotation is ambiguous vs HGSOC); an
HGSOC-only split is a follow-up, not a second primary analysis.
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/ovarian/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
ov_ids = model.loc[model['OncotreePrimaryDisease'] == 'Ovarian Epithelial Tumor', 'ModelID'].tolist()
ov = [x for x in ov_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(ov_ids)]
print('OV CRISPR n=', len(ov), 'comparator n=', len(comp))
print('OV lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in ov])
print('OV missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in ov_ids if i not in eff.index])
print('OV subtype counts (CRISPR panel):')
print(model.loc[model.ModelID.isin(ov), 'OncotreeSubtype'].value_counts(dropna=False).to_dict())

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p)
    pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

rows = []
for g in eff.columns:
    a = eff.loc[ov, g].dropna().values
    b = eff.loc[comp, g].dropna().values
    gw = eff[g].dropna().values
    hm = a.mean() if len(a) else np.nan
    cm = b.mean() if len(b) else np.nan
    hdep = (a < -0.5).sum(); cdep = (b < -0.5).sum()
    hfrac = hdep / len(a) if len(a) else np.nan
    cfrac = cdep / len(b) if len(b) else np.nan
    gwfrac = (gw < -0.5).mean() if len(gw) else np.nan
    if len(a) >= 2 and len(b) >= 2:
        t, p = stats.ttest_ind(a, b, equal_var=False)
    else:
        p = np.nan
    rows.append(dict(Gene=g, OV_mean=hm, OV_n=len(a), OV_dep=hdep, OV_dep_frac=hfrac,
                     Comp_mean=cm, Comp_n=len(b), Comp_dep=cdep, Comp_dep_frac=cfrac,
                     Enrichment=hm - cm, FracDiff=hfrac - cfrac, GenomeWide_dep_frac=gwfrac,
                     PanEssential=bool(hfrac >= 0.70 and cfrac >= 0.70), p=p))
scan = pd.DataFrame(rows)
scan['q'] = bh(scan['p'].values)
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
print('genes', len(scan))
print('q<0.05:', int((scan.q < 0.05).sum()), ' q<0.1:', int((scan.q < 0.1).sum()))
print('q<0.05 selective-direction (Enrichment<0):', int(((scan.q < 0.05) & (scan.Enrichment < 0)).sum()))
strict = scan[(scan.Enrichment <= -0.10) & (scan.OV_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.OV_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any OV line:', int(scan.OV_n.eq(0).sum()))
print('genes scored in <58 OV lines:', int((scan.OV_n < 58).sum()))
print(scan.head(25)[['Gene', 'OV_mean', 'OV_dep', 'OV_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'CCNE1', 'FOXM1', 'PAX8', 'WT1']:
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'OV_mean', 'OV_dep', 'OV_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
