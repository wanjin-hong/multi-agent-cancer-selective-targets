"""Small-cell lung cancer selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: Lung Neuroendocrine Tumor lines
(OncotreePrimaryDisease == 'Lung Neuroendocrine Tumor', 26 lines with CRISPR
data). This is the SCLC group; subtype annotation (ASCL1/NEUROD1) is checked
in the report from OncotreeSubtype.
Comparator: all non-LNET lines in the effect matrix (the matrix contains no
non-cancerous lines at all, so the comparator is entirely cancerous).
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/small-cell-lung/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
sclc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Lung Neuroendocrine Tumor', 'ModelID'].tolist()
sclc = [x for x in sclc_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(sclc_ids)]
print('SCLC CRISPR n=', len(sclc), 'comparator n=', len(comp))
print('SCLC lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in sclc])
print('SCLC missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in sclc_ids if i not in eff.index])
print('SCLC subtype counts (CRISPR panel):')
print(model.loc[model.ModelID.isin(sclc), 'OncotreeSubtype'].value_counts(dropna=False).to_dict())
print('SCLC subtype counts (all Model lines):')
print(model.loc[model.ModelID.isin(sclc_ids), 'OncotreeSubtype'].value_counts(dropna=False).to_dict())
# same-patient / duplicate check: CellLineName duplicates and PatientID
sub = model.loc[model.ModelID.isin(sclc)].copy()
print('PatientID counts:', sub['PatientID'].value_counts().to_dict())
print('duplicate CellLineName:', sub[sub.duplicated('CellLineName', keep=False)][['CellLineName', 'ModelID']].to_string() if sub.duplicated('CellLineName').any() else 'none')

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
    a = eff.loc[sclc, g].dropna().values
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
    rows.append(dict(Gene=g, SCLC_mean=hm, SCLC_n=len(a), SCLC_dep=hdep, SCLC_dep_frac=hfrac,
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
strict = scan[(scan.Enrichment <= -0.10) & (scan.SCLC_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.SCLC_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any SCLC line:', int(scan.SCLC_n.eq(0).sum()))
print('genes scored in <26 SCLC lines:', int((scan.SCLC_n < 26).sum()))
print('low-n breakdown:', scan[scan.SCLC_n < 26]['SCLC_n'].value_counts().to_dict())
print(scan.head(25)[['Gene', 'SCLC_mean', 'SCLC_dep', 'SCLC_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'MYCL', 'MYCN', 'ASCL1', 'NEUROD1', 'POU2F3', 'INSM1', 'DLL3', 'BCL2', 'EZH2', 'LSD1', 'KDM1A', 'RB1', 'TP53', 'SOX10', 'KLF5', 'FOSL1']:
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'SCLC_mean', 'SCLC_dep', 'SCLC_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
