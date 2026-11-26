"""Cell-line RNA: AML vs all other lines (primary leg) and AML vs other
blood-cancer lines (sensitivity leg), DepMap 24Q4 expression matrix."""
import pandas as pd, numpy as np
from scipy import stats
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/aml/'
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('expr shape', expr.shape, 'missing', int(expr.isna().sum().sum()))
aml_ids = model.loc[model['OncotreePrimaryDisease'] == 'Acute Myeloid Leukemia', 'ModelID'].tolist()
aml = [x for x in aml_ids if x in expr.index]
comp = [x for x in expr.index if x not in set(aml_ids)]
BLOOD = ['Mature B-Cell Neoplasms', 'B-Lymphoblastic Leukemia/Lymphoma',
         'Myeloproliferative Neoplasms', 'Mature T and NK Neoplasms',
         'T-Lymphoblastic Leukemia/Lymphoma', 'Non-Hodgkin Lymphoma',
         'Hodgkin Lymphoma', 'Acute Leukemias of Ambiguous Lineage']
blood_ids = model.loc[model['OncotreePrimaryDisease'].isin(BLOOD), 'ModelID'].tolist()
haem = [x for x in blood_ids if x in expr.index and x not in set(aml_ids)]
nc = model.loc[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(expr.index)), 'ModelID'].tolist()
print('AML RNA n=', len(aml), 'pan-cancer RNA comparator n=', len(comp), 'of which Non-Cancerous:', len(nc))
print('haem RNA comparator n=', len(haem))
print('AML missing RNA:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in aml_ids if i not in expr.index])
print('AML lines in CRISPR but not RNA:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in aml_ids if i in pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', usecols=[0]).iloc[:, 0].tolist() and i not in expr.index])

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

rows = []
for g in expr.columns:
    a = expr.loc[aml, g].dropna().values
    b = expr.loc[comp, g].dropna().values
    h = expr.loc[haem, g].dropna().values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    hm = h.mean() if len(h) else np.nan
    d = am - bm
    dh = am - hm
    fc = (2**am - 1) / (2**bm - 1) if (2**bm - 1) != 0 else np.nan
    fch = (2**am - 1) / (2**hm - 1) if (2**hm - 1) != 0 else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    ph = stats.ttest_ind(a, h, equal_var=False)[1] if len(a) >= 2 and len(h) >= 2 else np.nan
    sym = g.split(' (')[0]
    rows.append(dict(Gene=sym, Column=g, AML_mean_log2TPM=am, AML_n=len(a),
                     Comp_mean_log2TPM=bm, Comp_n=len(b), DeltaLog2=d, TPM_FC=fc, p=p,
                     HaemComp_mean_log2TPM=hm, HaemComp_n=len(h), DeltaLog2_haem=dh, TPM_FC_haem=fch, p_haem=ph))
rna = pd.DataFrame(rows)
rna['q'] = bh(rna['p'].values)
rna['q_haem'] = bh(rna['p_haem'].values)
rna = rna.sort_values('DeltaLog2', ascending=False)
rna.to_csv(OUT + 'celline_rna.csv', index=False)
print('genes', len(rna))
bar = rna[(rna.TPM_FC >= 1.5) & (rna.DeltaLog2 >= 0.5) & (rna.p < 0.05)]
bar_h = rna[(rna.TPM_FC_haem >= 1.5) & (rna.DeltaLog2_haem >= 0.5) & (rna.p_haem < 0.05)]
print('cell-line RNA bar pass (pan-cancer):', len(bar), ' (haem):', len(bar_h))
print(rna.head(20)[['Gene', 'AML_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'q']].to_string(index=False))
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'MYB', 'CBFB', 'SPI1', 'RUNX1', 'CEBPA', 'CDK6', 'FLT3', 'IDH1', 'IDH2', 'NPM1', 'KMT2A', 'BCL2', 'CD33', 'IL3RA', 'KLF5', 'FOSL1', 'HNF1B', 'PARD6B', 'CXADR', 'JAK2', 'MCL1', 'BCL2L1', 'KIT', 'DOT1L', 'MEN1', 'MEIS1', 'HOXA9', 'GATA2', 'ZEB2', 'LMO2', 'CSF2RB']:
    r = rna[rna.Gene == g]
    print(g, r[['Gene', 'DeltaLog2', 'TPM_FC', 'p', 'DeltaLog2_haem', 'TPM_FC_haem', 'p_haem']].to_string(index=False) if len(r) else 'NOT FOUND')
