"""AML selective-targets dependency scan — DepMap 24Q4, TNBC pipeline replica.

Group: Acute Myeloid Leukemia lines
(OncotreePrimaryDisease == 'Acute Myeloid Leukemia', 61 lines in Model.csv,
28 with CRISPR data).

TWO comparators (the haematology adaptation, stated plainly):
- PRIMARY (programme-comparable): all 1,150 non-AML lines in the effect
  matrix. This comparator is mostly solid-tumour lines and includes the 14
  Non-Cancerous lines present in CRISPRGeneEffect.csv (verified 21 Nov 2026;
  14 of 1,150 — negligible effect, but recorded exactly).
- SENSITIVITY (honest control): the 104 other blood-cancer lines with CRISPR
  data (Mature B-Cell Neoplasms 57, B-Lymphoblastic Leukemia/Lymphoma 13,
  Myeloproliferative Neoplasms 11, Mature T and NK Neoplasms 10,
  T-Lymphoblastic Leukemia/Lymphoma 6, Non-Hodgkin Lymphoma 4, Hodgkin
  Lymphoma 3; Acute Leukemias of Ambiguous Lineage 0 with CRISPR data).
  The solid-tumour comparator inflates lineage/haematopoietic genes; the
  blood-only comparator is the honest control for those.

AML subtypes in the CRISPR panel: Acute Myeloid Leukemia 21, Acute
Monoblastic/Monocytic Leukemia 3, Acute Megakaryoblastic Leukemia 3,
Myeloid Leukemia Associated with Down Syndrome 1. No same-patient pairs
found in the CRISPR panel (checked below).
"""
import pandas as pd, numpy as np
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/aml/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
print('genes:', eff.shape[1])

aml_ids = model.loc[model['OncotreePrimaryDisease'] == 'Acute Myeloid Leukemia', 'ModelID'].tolist()
aml = [x for x in aml_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(aml_ids)]
print('AML CRISPR n=', len(aml), 'pan-cancer comparator n=', len(comp))
print('AML lines:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in aml])
print('AML missing CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in aml_ids if i not in eff.index])
print('AML subtype counts (CRISPR panel):')
print(model.loc[model.ModelID.isin(aml), 'OncotreeSubtype'].value_counts(dropna=False).to_dict())
# same-patient pairs inside the group?
sub = model.loc[model.ModelID.isin(aml)]
print('distinct PatientIDs in CRISPR panel:', sub['PatientID'].nunique(), 'of', len(sub))
print('duplicated PatientIDs:', sub[sub.duplicated('PatientID', keep=False)][['ModelID', 'CellLineName', 'PatientID']].to_string() if sub.duplicated('PatientID').any() else 'none')

BLOOD = ['Mature B-Cell Neoplasms', 'B-Lymphoblastic Leukemia/Lymphoma',
         'Myeloproliferative Neoplasms', 'Mature T and NK Neoplasms',
         'T-Lymphoblastic Leukemia/Lymphoma', 'Non-Hodgkin Lymphoma',
         'Hodgkin Lymphoma', 'Acute Leukemias of Ambiguous Lineage']
blood_ids = model.loc[model['OncotreePrimaryDisease'].isin(BLOOD), 'ModelID'].tolist()
haem = [x for x in blood_ids if x in eff.index and x not in set(aml_ids)]
print('haematology comparator n=', len(haem))
print('haem composition:', model.loc[model.ModelID.isin(haem), 'OncotreePrimaryDisease'].value_counts().to_dict())
nc = model.loc[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(eff.index)), 'ModelID'].tolist()
print('Non-Cancerous lines in effect matrix:', len(nc), '(in pan-cancer comparator)')

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p)
    pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

def scan_against(comp_ids, tag):
    rows = []
    for g in eff.columns:
        a = eff.loc[aml, g].dropna().values
        b = eff.loc[comp_ids, g].dropna().values
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
        rows.append(dict(Column=g, AML_mean=hm, AML_dep=hdep, AML_n=len(a), AML_dep_frac=hfrac,
                         Comp_mean=cm, Comp_dep=cdep, Comp_n=len(b), Comp_dep_frac=cfrac,
                         Enrichment=hm - cm, FracDiff=hfrac - cfrac, GenomeWide_dep_frac=gwfrac,
                         PanEssential=bool(hfrac >= 0.70 and cfrac >= 0.70), p_welch=p))
    s = pd.DataFrame(rows)
    s['q_BH'] = bh(s['p_welch'].values)
    s['EnrichRank'] = s['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
    s['Gene'] = s['Column'].str.split(' (', regex=False).str[0]
    s = s.sort_values('Enrichment')
    strict = s[(s.Enrichment <= -0.10) & (s.AML_dep_frac >= 0.25) & (s.FracDiff >= 0.10) & (~s.PanEssential) & (s.GenomeWide_dep_frac <= 0.60)]
    relaxed = s[(s.Enrichment <= -0.05) & (s.AML_dep_frac >= 0.20) & (s.FracDiff >= 0.10) & (~s.PanEssential) & (s.GenomeWide_dep_frac <= 0.60)]
    print(f'--- {tag}: genes {len(s)}  q<0.05: {int((s.q_BH < 0.05).sum())}  selective q<0.05: {int(((s.q_BH < 0.05) & (s.Enrichment < 0)).sum())}')
    print(f'    strict pass: {len(strict)}  relaxed pass: {len(relaxed)}  pan-essential: {int(s.PanEssential.sum())}')
    print(f'    genes scored in <28 AML lines: {int((s.AML_n < 28).sum())}')
    print(s.head(20)[['Gene', 'AML_mean', 'AML_dep', 'AML_dep_frac', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'q_BH']].to_string(index=False))
    return s, strict, relaxed

scan, strict, relaxed = scan_against(comp, 'PAN-CANCER comparator')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
strict.to_csv(OUT + 'strict_bar.csv', index=False)
relaxed.to_csv(OUT + 'relaxed_bar.csv', index=False)

hscan, hstrict, hrelaxed = scan_against(haem, 'HAEMATOLOGY comparator')
hscan.to_csv(OUT + 'dependency_scan_haem.csv', index=False)
hstrict.to_csv(OUT + 'strict_bar_haem.csv', index=False)
hrelaxed.to_csv(OUT + 'relaxed_bar_haem.csv', index=False)

CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'ITGB1', 'ITGB5', 'RAC1', 'BRAF', 'MAPK1',
         'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'MYC', 'CCNE1', 'FOXM1', 'PAX8',
         'FLT3', 'IDH1', 'IDH2', 'NPM1', 'KMT2A', 'BCL2', 'CD33', 'CDK6', 'KLF5',
         'FOSL1', 'HNF1B', 'PARD6B', 'CXADR', 'DOT1L', 'MEN1', 'FLT3LG', 'MCL1',
         'BCL2L1', 'KIT', 'JAK2', 'TET2', 'DNMT3A', 'ASXL1', 'RUNX1', 'CEBPA', 'GATA2',
         'SPI1', 'MYB', 'HOXA9', 'MEIS1', 'PBX3', 'CD47', 'CD123', 'IL3RA', 'SIRPA']
for g in CHECK:
    r = scan[scan.Gene == g]
    h = hscan[hscan.Gene == g]
    if len(r):
        print(g, 'PAN enrich rank', int(r.iloc[0].EnrichRank), 'enrich', round(r.iloc[0].Enrichment, 3),
              'dep', f"{int(r.iloc[0].AML_dep)}/28", 'q', f"{r.iloc[0].q_BH:.2g}",
              '| HAEM enrich rank', int(h.iloc[0].EnrichRank), 'enrich', round(h.iloc[0].Enrichment, 3),
              'dep frac diff', round(h.iloc[0].FracDiff, 3), 'q', f"{h.iloc[0].q_BH:.2g}")
    else:
        print(g, 'NOT FOUND in scan')
