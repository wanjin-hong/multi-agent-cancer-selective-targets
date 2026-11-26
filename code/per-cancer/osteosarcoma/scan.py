"""Osteosarcoma selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: Osteosarcoma lines (OncotreePrimaryDisease == 'Osteosarcoma',
28 lines in Model.csv, 17 with CRISPR data). Smallest group in the
programme (n = 17) — read as a hypothesis generator only.
Comparator: all 1,161 non-osteosarcoma lines in the effect matrix,
including the 14 Non-Cancerous lines (verified 21 Nov 2026 from Model.csv
+ the effect matrix; 14 of 1,161, negligible effect). The HNSC report's
claim that the matrix contains no non-cancerous lines is wrong.
Parent/derivative pair in the group: 143B was derived from HOS (KRAS
transformation of the same patient's tumour) — effective CRISPR n ~16.
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/osteosarcoma/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
os_ids = model.loc[model['OncotreePrimaryDisease'] == 'Osteosarcoma', 'ModelID'].tolist()
os_ = [x for x in os_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(os_ids)]
print('OS CRISPR n=', len(os_), 'comparator n=', len(comp))
nc = model.loc[model['OncotreePrimaryDisease'] == 'Non-Cancerous', 'ModelID'].tolist()
print('non-cancerous lines in effect matrix:', len([x for x in nc if x in eff.index]),
      'of which in comparator:', len([x for x in nc if x in eff.index and x not in set(os_ids)]))
print('OS lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in os_])
print('OS missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in os_ids if i not in eff.index])
print('OS age/sex (CRISPR panel):')
print(model.loc[model.ModelID.isin(os_), ['CellLineName', 'Age', 'Sex']].to_string(index=False))
print('OS subtype counts (CRISPR panel):')
print(model.loc[model.ModelID.isin(os_), 'OncotreeSubtype'].value_counts(dropna=False).to_dict())

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
    a = eff.loc[os_, g].dropna().values
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
    sym = g.split(' (')[0]
    rows.append(dict(Gene=g, Symbol=sym, OS_mean=hm, OS_n=len(a), OS_dep=hdep, OS_dep_frac=hfrac,
                     Comp_mean=cm, Comp_n=len(b), Comp_dep=cdep, Comp_dep_frac=cfrac,
                     Enrichment=hm - cm, FracDiff=hfrac - cfrac, GenomeWide_dep_frac=gwfrac,
                     PanEssential=bool(hfrac >= 0.70 and cfrac >= 0.70), p=p))
scan = pd.DataFrame(rows)
scan['q'] = bh(scan['p'].values)
# EnrichRank: rank by enrichment ascending (1 = most selective), ties broken by order
scan = scan.sort_values('Enrichment').reset_index(drop=True)
scan['EnrichRank'] = np.arange(1, len(scan) + 1)
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
print('genes', len(scan))
print('q<0.05:', int((scan.q < 0.05).sum()), ' q<0.1:', int((scan.q < 0.1).sum()))
print('q<0.05 selective-direction (Enrichment<0):', int(((scan.q < 0.05) & (scan.Enrichment < 0)).sum()))
strict = scan[(scan.Enrichment <= -0.10) & (scan.OS_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.OS_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any OS line:', int(scan.OS_n.eq(0).sum()))
print('genes scored in <17 OS lines:', int((scan.OS_n < 17).sum()))
print('genes scored in <10 OS lines:', int((scan.OS_n < 10).sum()))
print(scan.head(25)[['Symbol', 'OS_mean', 'OS_dep', 'OS_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'CCNE1', 'TP53', 'RB1', 'MDM2', 'CDK4', 'CDK6', 'FGFR1', 'PDGFRA', 'IGF1R', 'MET', 'SRC', 'RUNX2', 'SATB2', 'SP7', 'DLX5']:
    r = scan[scan.Symbol == g]
    if len(r):
        print(g, r[['Symbol', 'OS_mean', 'OS_n', 'OS_dep', 'OS_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int(r.iloc[0].EnrichRank))
    else:
        print(g, 'NOT FOUND')
