"""Esophagogastric adenocarcinoma (ESGA) selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: OncotreePrimaryDisease == 'Esophagogastric Adenocarcinoma', 74 lines in
Model.csv, 45 with CRISPR data. Comparator: all 1,133 non-ESGA lines in the effect
matrix (the matrix contains no non-cancerous lines at all, so the comparator is
entirely cancerous).
Site split inside the group (CRISPR panel): gastric adenocarcinoma 35 lines,
esophageal adenocarcinoma 9 lines, GEJ adenocarcinoma 1 line. The all-group scan
is primary; the site split is reported descriptively (esophageal n=9 is too small
for its own scan).
Exclusions: none needed for the CRISPR panel — the two small-cell carcinoma of the
stomach lines (ECC12, ECC10) and the adenosquamous line (MKN1) are all missing
CRISPR data anyway.
Same-patient pairs in the CRISPR panel: IS076P / IS076A (PT-qcon1I) — effective
CRISPR n ~44, not 45. (PT-i9nBPr, PT-T2KFTw and PT-L7CBMa pairs have at most one
line in the CRISPR panel — harmless.)
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/esophagogastric/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
esga_ids = model.loc[model['OncotreePrimaryDisease'] == 'Esophagogastric Adenocarcinoma', 'ModelID'].tolist()
esga = [x for x in esga_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(esga_ids)]
print('ESGA CRISPR n=', len(esga), 'comparator n=', len(comp))
print('ESGA lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in esga])
print('ESGA missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in esga_ids if i not in eff.index])
# site split within the CRISPR panel
sub = model.set_index('ModelID').loc[esga, ['CellLineName', 'OncotreeSubtype']]
print(sub['OncotreeSubtype'].value_counts().to_string())
gastric = sub[sub['OncotreeSubtype'].str.contains('Stomach', na=False)].index.tolist()
eso = sub[sub['OncotreeSubtype'].str.contains('Esophageal|Gastroesophageal', na=False)].index.tolist()
print('gastric CRISPR n=', len(gastric), 'esophageal+GEJ CRISPR n=', len(eso))
print('esophageal+GEJ lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in eso])

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
    a = eff.loc[esga, g].dropna().values
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
    # site-split descriptive counts
    ag = eff.loc[gastric, g].dropna().values
    ae = eff.loc[eso, g].dropna().values
    rows.append(dict(Gene=g, ESGA_mean=hm, ESGA_n=len(a), ESGA_dep=hdep, ESGA_dep_frac=hfrac,
                     Comp_mean=cm, Comp_n=len(b), Comp_dep=cdep, Comp_dep_frac=cfrac,
                     Enrichment=hm - cm, FracDiff=hfrac - cfrac, GenomeWide_dep_frac=gwfrac,
                     PanEssential=bool(hfrac >= 0.70 and cfrac >= 0.70), p=p,
                     Gastric_dep=int((ag < -0.5).sum()), Gastric_n=len(ag),
                     EsoGEJ_dep=int((ae < -0.5).sum()), EsoGEJ_n=len(ae)))
scan = pd.DataFrame(rows)
scan['q'] = bh(scan['p'].values)
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
print('genes', len(scan))
print('q<0.05:', int((scan.q < 0.05).sum()), ' q<0.1:', int((scan.q < 0.1).sum()))
print('q<0.05 selective-direction (Enrichment<0):', int(((scan.q < 0.05) & (scan.Enrichment < 0)).sum()))
strict = scan[(scan.Enrichment <= -0.10) & (scan.ESGA_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.ESGA_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any ESGA line:', int(scan.ESGA_n.eq(0).sum()))
print('genes scored in <45 ESGA lines:', int((scan.ESGA_n < 45).sum()))
print('genes scored in <=24 ESGA lines:', int((scan.ESGA_n <= 24).sum()))
print(scan.head(25)[['Gene', 'ESGA_mean', 'ESGA_dep', 'ESGA_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'KLF5', 'SOX9', 'HNF1B', 'FOSL1', 'CDX2', 'ERBB2', 'CDH17', 'GATA6']:
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'ESGA_mean', 'ESGA_dep', 'ESGA_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q', 'Gastric_dep', 'Gastric_n', 'EsoGEJ_dep', 'EsoGEJ_n']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
