"""Ewing sarcoma selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: Ewing Sarcoma lines (OncotreePrimaryDisease == 'Ewing Sarcoma',
51 lines in Model.csv, 27 with CRISPR data). Comparator: all 1,151
non-Ewing lines in the effect matrix — 1,137 cancerous + 14 non-cancerous
lines (they sit in it per the pre-registered pan-cancer comparator rule;
1.2% of the comparator, negligible effect, composition recorded exactly —
correcting the old claim, repeated in earlier reports, that the matrix
contains no non-cancerous lines).
Same-patient pair in the CRISPR group (1): SK-N-MC / MC-IXC share
PatientID PT-wGDMGe => effective CRISPR n ~26, not 27. The four other
PT-IAac7b lines (CCLF A-673 derivatives) are not in the CRISPR matrix.
Excluded from the group by construction: the four engineered A673 STAG2
knockout / non-targeting lines (A673STAG2KO16, A673STAG2KO45,
A673STAG2NT14, A673STAG2NT23) — EngineeredModel lines, and none of them
is in the CRISPR matrix anyway.
Fusion annotation: Model.csv carries no fusion annotation, and the
DepMap 24Q4 somatic-mutation file records point mutations/indels, not
fusions — no EWSR1-FLI1 vs EWSR1-ERG split is possible from the data on
disk (see fusion_check.py). The EWSR1-FLI1 programme (FLI1, EWSR1 and
downstream targets) is checked explicitly in Q4.
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/ewing/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
ew_ids = model.loc[model['OncotreePrimaryDisease'] == 'Ewing Sarcoma', 'ModelID'].tolist()
ew = [x for x in ew_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(ew_ids)]
print('Ewing CRISPR n=', len(ew), 'comparator n=', len(comp))
nc = model[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(eff.index))]
print('non-cancerous lines in comparator:', len(nc), nc['CellLineName'].tolist())
sub = model.set_index('ModelID')
print('Ewing lines:', [sub.loc[i, 'CellLineName'] for i in ew])
print('PatientID dupes in CRISPR group:', int(sub.loc[ew, 'PatientID'].duplicated().sum()))

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
    a = eff.loc[ew, g].dropna().values
    b = eff.loc[comp, g].dropna().values
    gw = eff[g].dropna().values
    hm = a.mean() if len(a) else np.nan
    cm = b.mean() if len(b) else np.nan
    hdep = int((a < -0.5).sum()); cdep = int((b < -0.5).sum())
    hfrac = hdep / len(a) if len(a) else np.nan
    cfrac = cdep / len(b) if len(b) else np.nan
    gwfrac = (gw < -0.5).mean() if len(gw) else np.nan
    if len(a) >= 2 and len(b) >= 2:
        t, p = stats.ttest_ind(a, b, equal_var=False)
    else:
        p = np.nan
    rows.append(dict(Gene=g, Ewing_mean=hm, Ewing_n=len(a), Ewing_dep=hdep, Ewing_dep_frac=hfrac,
                     Comp_mean=cm, Comp_n=len(b), Comp_dep=cdep, Comp_dep_frac=cfrac,
                     Enrichment=hm - cm, FracDiff=hfrac - cfrac, GenomeWide_dep_frac=gwfrac,
                     PanEssential=bool(hfrac >= 0.70 and cfrac >= 0.70), p=p))
scan = pd.DataFrame(rows)
scan['q'] = bh(scan['p'].values)
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
print('genes', len(scan))
print('q<0.05:', int((scan.q < 0.05).sum()), ' q<0.05 selective-direction (Enrichment<0):',
      int(((scan.q < 0.05) & (scan.Enrichment < 0)).sum()))
strict = scan[(scan.Enrichment <= -0.10) & (scan.Ewing_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.Ewing_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('low-n genes: scored in <27 Ewing lines:', int((scan.Ewing_n < 27).sum()))
print(scan[scan.Ewing_n < 27][['Gene','Ewing_n']].head(20).to_string())
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)
