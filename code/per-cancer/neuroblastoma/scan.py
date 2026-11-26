"""Neuroblastoma selective-targets dependency scan — DepMap 24Q4,
TNBC pipeline replica.

Group: Neuroblastoma lines (OncotreePrimaryDisease == 'Neuroblastoma',
51 lines in Model.csv, 39 with CRISPR data, 36 with RNA data).
Comparator: all 1,139 non-neuroblastoma lines in the effect matrix —
including the 14 Non-Cancerous lines in CRISPRGeneEffect.csv (verified
21 Nov 2026; 14 of 1,139, negligible effect, stated explicitly).
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/neuroblastoma/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
nb_ids = model.loc[model['OncotreePrimaryDisease'] == 'Neuroblastoma', 'ModelID'].tolist()
nb = [x for x in nb_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(nb_ids)]
print('NB CRISPR n=', len(nb), 'comparator n=', len(comp))
print('NB lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in nb])
print('NB missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in nb_ids if i not in eff.index])
sub = model.set_index('ModelID')
crispr_sub = sub.loc[nb]
print('PatientID dupes in CRISPR group:', int(crispr_sub['PatientID'].duplicated().sum()))
dups = crispr_sub[crispr_sub['PatientID'].duplicated(keep=False)].sort_values('PatientID')
print(dups[['CellLineName', 'PatientID', 'Age']].to_string())
print('subtypes in CRISPR group:')
print(crispr_sub['OncotreeSubtype'].value_counts(dropna=False))
print('ages in CRISPR group:', crispr_sub['Age'].describe().to_string())
print('Non-Cancerous lines in comparator:', int(sum(1 for i in comp if model.loc[model.ModelID == i, 'OncotreePrimaryDisease'].iloc[0] == 'Non-Cancerous')))

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
    a = eff.loc[nb, g].dropna().values
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
    rows.append(dict(Gene=g, NB_mean=hm, NB_n=len(a), NB_dep=hdep, NB_dep_frac=hfrac,
                     Comp_mean=cm, Comp_n=len(b), Comp_dep=cdep, Comp_dep_frac=cfrac,
                     Enrichment=hm - cm, FracDiff=hfrac - cfrac, GenomeWide_dep_frac=gwfrac,
                     PanEssential=bool(hfrac >= 0.70 and cfrac >= 0.70), p=p))
scan = pd.DataFrame(rows)
scan['q'] = bh(scan['p'].values)
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan_raw.csv', index=False)
print('genes', len(scan))
print('q<0.05:', int((scan.q < 0.05).sum()), ' q<0.05 selective-direction (Enrichment<0):', int(((scan.q < 0.05) & (scan.Enrichment < 0)).sum()))
strict = scan[(scan.Enrichment <= -0.10) & (scan.NB_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.NB_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any NB line:', int(scan.NB_n.eq(0).sum()))
print('scored-line-count distribution (NB_n):')
print(scan['NB_n'].value_counts().sort_index().to_string())
print(scan.head(25)[['Gene', 'NB_mean', 'NB_dep', 'NB_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'MYCN', 'MYC', 'MYCL', 'ALK', 'PHOX2B', 'ATRX', 'KLF5', 'FERMT2', 'EGFR', 'HNF1B',
         'FOSL1', 'SOX9', 'MT1E', 'ITGAV', 'TCF7L2', 'PARD3', 'PARD6B',
         'ARHGAP29', 'ITGB5', 'CDK6', 'TIPARP', 'KRAS', 'BRAF', 'MAPK1',
         'CCND1', 'TP63', 'SOX10', 'CDX2', 'GRB2', 'PHOX2A', 'HAND2', 'GATA3',
         'ISL1', 'ASCL1', 'NEUROD1', 'BCL2', 'MDM2', 'CDK4', 'AURKA', 'CHEK1',
         'EZH2', 'HDAC2', 'LIN28B', 'TERT', 'SOX11', 'DLK1', 'DBH', 'TH',
         'NCAM1', 'CHGA', 'SYP', 'RET', 'NTRK1', 'NTRK2', 'NGFR', 'SOX9']
for g in dict.fromkeys(CHECK):
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'NB_mean', 'NB_dep', 'NB_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
