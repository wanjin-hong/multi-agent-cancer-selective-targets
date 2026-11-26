"""Cell-line RNA: RCC vs all other lines, DepMap 24Q4 expression matrix."""
import pandas as pd, numpy as np
from scipy import stats
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/renal-cell/'
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('expr shape', expr.shape, 'missing', int(expr.isna().sum().sum()))
rcc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Renal Cell Carcinoma', 'ModelID'].tolist()
rcc = [x for x in rcc_ids if x in expr.index]
comp = [x for x in expr.index if x not in set(rcc_ids)]
print('RCC RNA n=', len(rcc), 'comparator n=', len(comp))
print('RCC missing RNA:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in rcc_ids if i not in expr.index])

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

rows = []
for g in expr.columns:
    a = expr.loc[rcc, g].dropna().values
    b = expr.loc[comp, g].dropna().values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    d = am - bm
    fc = (2**am - 1) / (2**bm - 1) if (2**bm - 1) != 0 else np.nan
    if len(a) >= 2 and len(b) >= 2:
        t, p = stats.ttest_ind(a, b, equal_var=False)
    else:
        p = np.nan
    sym = g.split(' (')[0]
    rows.append(dict(Gene=sym, Column=g, RCC_mean_log2TPM=am, RCC_n=len(a), Comp_mean_log2TPM=bm, Comp_n=len(b),
                     DeltaLog2=d, TPM_FC=fc, p=p))
rna = pd.DataFrame(rows)
rna['q'] = bh(rna['p'].values)
rna = rna.sort_values('DeltaLog2', ascending=False)
rna.to_csv(OUT + 'celline_rna.csv', index=False)
print('genes', len(rna))
bar = rna[(rna.TPM_FC >= 1.5) & (rna.DeltaLog2 >= 0.5) & (rna.p < 0.05)]
print('cell-line RNA bar pass:', len(bar))
print(rna.head(20)[['Gene', 'RCC_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'q']].to_string())
CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'VHL', 'HIF1A', 'EPAS1', 'EGLN1', 'KLF5', 'FERMT2', 'EGFR', 'HNF1B',
         'FOSL1', 'SOX9', 'MT1E', 'ITGAV', 'TCF7L2', 'PARD3', 'PARD6B',
         'ARHGAP29', 'ITGB5', 'CDK6', 'TIPARP', 'KRAS', 'CA9', 'NDUFA4L2',
         'SLC2A1', 'VEGFA', 'PAX8', 'MET', 'KDR', 'FLT1', 'PDGFRA', 'FGFR1',
         'CFLAR', 'MDM2', 'CDH2', 'HAPSTR1', 'ZFP36L1', 'UNC50', 'SEPHS2']
for g in CHECK:
    r = rna[rna.Gene == g]
    print(g, r[['Gene', 'DeltaLog2', 'TPM_FC', 'p', 'q']].to_string(index=False) if len(r) else 'NOT FOUND')
