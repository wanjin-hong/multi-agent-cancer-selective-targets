"""HNSC selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: Head and Neck Squamous Cell Carcinoma lines
(OncotreePrimaryDisease == 'Head and Neck Squamous Cell Carcinoma', 90 lines in
Model.csv, 72 with CRISPR data). Excludes 3 'Head and Neck Carcinoma, Other'
(RPMI 2650, Ty-82 — NUT midline carcinoma; DOK) and 2 salivary carcinoma lines,
mirroring the PDAC exclusion of the single neuroendocrine line QGP-1.
Comparator: all 1,106 non-HNSC lines in the effect matrix (the matrix contains
no non-cancerous lines at all, so the comparator is entirely cancerous).
Same-patient pairs in the group: UPCI-SCC-090 / UPCI-SCC-152 (both in CRISPR)
=> effective CRISPR n ~71, not 72.
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/head-and-neck/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
hnsc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Head and Neck Squamous Cell Carcinoma', 'ModelID'].tolist()
hnsc = [x for x in hnsc_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(hnsc_ids)]
print('HNSC CRISPR n=', len(hnsc), 'comparator n=', len(comp))
print('HNSC lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in hnsc])
print('HNSC missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in hnsc_ids if i not in eff.index])

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
    a = eff.loc[hnsc, g].dropna().values
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
    rows.append(dict(Gene=g, HNSC_mean=hm, HNSC_n=len(a), HNSC_dep=hdep, HNSC_dep_frac=hfrac,
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
strict = scan[(scan.Enrichment <= -0.10) & (scan.HNSC_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.HNSC_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any HNSC line:', int(scan.HNSC_n.eq(0).sum()))
print('genes scored in <72 HNSC lines:', int((scan.HNSC_n < 72).sum()))
print(scan.head(25)[['Gene', 'HNSC_mean', 'HNSC_dep', 'HNSC_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA']:
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'HNSC_mean', 'HNSC_dep', 'HNSC_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
