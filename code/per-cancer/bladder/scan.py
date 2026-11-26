"""Bladder urothelial carcinoma selective-targets dependency scan — DepMap 24Q4,
TNBC pipeline replica.

Group: Bladder Urothelial Carcinoma lines (OncotreePrimaryDisease ==
'Bladder Urothelial Carcinoma', 35 lines in Model.csv, 30 with CRISPR data,
33 with RNA data). No luminal/basal annotation exists in Model.csv (all 35
carry subtype 'Bladder Urothelial Carcinoma'), so no subtype split is
possible from the metadata; FGFR3 mutation status is checked separately
as a luminal-papillary proxy (see fgfr3_check.py).
Comparator: all 1,148 non-bladder lines in the effect matrix — 1,134
cancerous + 14 non-cancerous (Hs 729, HA1E, WPE1-NA22, BPH-1, P4E6,
Shmac 4, Shmac 5, ARH-77, IM-9, RPE1-ss48, RPE1-ss77, RPE1-ss6,
RPE1-ss119, RPE1-ss51). The 14 non-cancerous lines are 1.2% of the
comparator; they sit in it per the pre-registered pan-cancer comparator
rule and their effect on rankings is negligible, but the composition is
recorded exactly here (correcting the old claim, repeated in earlier
reports, that the matrix contains no non-cancerous lines).
Same-patient pairs in the CRISPR group (3): 253J / 253J-BV, RT-112 /
RT112/84, T24 / SLR 20 => effective CRISPR n ~27, not 30.
Excluded, mirroring the PDAC/HNSC exclusions: 2 'Bladder Squamous Cell
Carcinoma' lines (not urothelial).
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/bladder/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
blad_ids = model.loc[model['OncotreePrimaryDisease'] == 'Bladder Urothelial Carcinoma', 'ModelID'].tolist()
blad = [x for x in blad_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(blad_ids)]
print('Bladder CRISPR n=', len(blad), 'comparator n=', len(comp))
nc = model[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(eff.index))]
print('non-cancerous lines in comparator:', len(nc), nc['CellLineName'].tolist())
print('Bladder lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in blad])
print('Bladder missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in blad_ids if i not in eff.index])
sub = model.set_index('ModelID')
print('PatientID dupes in CRISPR group:', sub.loc[blad, 'PatientID'].duplicated().sum())

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
    a = eff.loc[blad, g].dropna().values
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
    rows.append(dict(Gene=g, Bladder_mean=hm, Bladder_n=len(a), Bladder_dep=hdep, Bladder_dep_frac=hfrac,
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
strict = scan[(scan.Enrichment <= -0.10) & (scan.Bladder_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.Bladder_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any bladder line:', int(scan.Bladder_n.eq(0).sum()))
print('genes scored in <30 bladder lines:', int((scan.Bladder_n < 30).sum()))
print('scored-line-count distribution (Bladder_n):')
print(scan['Bladder_n'].value_counts().sort_index().to_string())
print(scan.head(25)[['Gene', 'Bladder_mean', 'Bladder_dep', 'Bladder_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q']].to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'FGFR3', 'FGFR1', 'ERBB2', 'ERBB3', 'EGFR', 'KLF5', 'FERMT2', 'HNF1B',
         'FOSL1', 'SOX9', 'MT1E', 'ITGAV', 'TCF7L2', 'PARD3', 'PARD6B',
         'ARHGAP29', 'ITGB5', 'CDK6', 'TIPARP', 'KRAS', 'BRAF', 'MAPK1',
         'CCND1', 'TP63', 'GATA3', 'FOXA1', 'PPARG', 'KRT5', 'KRT20',
         'CDKN2A', 'RB1', 'TP53', 'PIK3CA', 'GRB2', 'MYBL2', 'CXADR']
for g in CHECK:
    r = scan[scan.Gene.str.startswith(g + ' (')]
    if len(r):
        print(g, r[['Gene', 'Bladder_mean', 'Bladder_dep', 'Bladder_dep_frac', 'Comp_dep_frac', 'Enrichment', 'GenomeWide_dep_frac', 'q']].to_string(index=False),
              ' enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum()) + 1)
    else:
        print(g, 'NOT FOUND')
