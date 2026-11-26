"""Renal cell carcinoma selective-targets dependency scan — DepMap 24Q4,
TNBC pipeline replica.

Group: Renal Cell Carcinoma lines (OncotreePrimaryDisease == 'Renal Cell
Carcinoma', 57 lines in Model.csv, 30 with CRISPR data, 37 with RNA data).
Subtype split inside the CRISPR group: 16 clear cell, 12 unspecified
'Renal Cell Carcinoma', 1 papillary (Caki-2), 1 medullary
(CCLF_PEDS_0001_T). Unusual members kept per the pre-registered group rule
and flagged in the report: G-402 (age 0, rhabdoid tumour of kidney) and
CCLF_PEDS_0001_T (age 15, renal medullary carcinoma).
Comparator: all 1,148 non-RCC lines in the effect matrix (the matrix
contains no non-cancerous lines at all, so the comparator is entirely
cancerous).
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/renal-cell/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
rcc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Renal Cell Carcinoma', 'ModelID'].tolist()
rcc = [x for x in rcc_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(rcc_ids)]
print('RCC CRISPR n=', len(rcc), 'comparator n=', len(comp))
print('RCC lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in rcc])
print('RCC missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in rcc_ids if i not in eff.index])
# same-patient check via PatientID
sub = model.set_index('ModelID')
print('PatientID dupes in CRISPR group:', sub.loc[rcc, 'PatientID'].duplicated().sum())
print('subtypes in CRISPR group:')
print(sub.loc[rcc, 'OncotreeSubtype'].value_counts(dropna=False))

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
    a = eff.loc[rcc, g].dropna().values
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
    rows.append(dict(Gene=g, RCC_mean=hm, RCC_n=len(a), RCC_dep=hdep, RCC_dep_frac=hfrac,
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
strict = scan[(scan.Enrichment <= -0.10) & (scan.RCC_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.RCC_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any RCC line:', int(scan.RCC_n.eq(0).sum()))
print('genes scored in <30 RCC lines:', int((scan.RCC_n < 30).sum()))
print('scored-line-count distribution (RCC_n):')
print(scan['RCC_n'].value_counts().sort_index().to_string())
print(scan.head(25)[['Gene', 'RCC_mean', 'RCC_dep', 'RCC_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'VHL', 'HIF1A', 'EPAS1', 'EGLN1', 'KLF5', 'FERMT2', 'EGFR', 'HNF1B',
         'FOSL1', 'SOX9', 'MT1E', 'ITGAV', 'TCF7L2', 'PARD3', 'PARD6B',
         'ARHGAP29', 'ITGB5', 'CDK6', 'TIPARP', 'KRAS', 'BRAF', 'MAPK1',
         'CCND1', 'MET', 'FLT1', 'KDR', 'PDGFRA', 'FGFR1', 'CA9', 'NDUFA4L2',
         'SLC2A1', 'VEGFA', 'PAX8', 'WT1', 'TP63', 'SOX10', 'CDX2', 'GRB2']
for g in CHECK:
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'RCC_mean', 'RCC_dep', 'RCC_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
