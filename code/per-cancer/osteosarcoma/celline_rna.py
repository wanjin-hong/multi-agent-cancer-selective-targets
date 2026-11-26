"""Cell-line RNA: osteosarcoma vs all other lines, DepMap 24Q4 expression matrix."""
import pandas as pd, numpy as np
from scipy import stats
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/osteosarcoma/'
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('expr shape', expr.shape, 'missing', int(expr.isna().sum().sum()))
os_ids = model.loc[model['OncotreePrimaryDisease'] == 'Osteosarcoma', 'ModelID'].tolist()
os_ = [x for x in os_ids if x in expr.index]
comp = [x for x in expr.index if x not in set(os_ids)]
print('OS RNA n=', len(os_), 'comparator n=', len(comp))
nc = model.loc[model['OncotreePrimaryDisease'] == 'Non-Cancerous', 'ModelID'].tolist()
print('non-cancerous lines in expression matrix:', len([x for x in nc if x in expr.index]),
      'of which in comparator:', len([x for x in nc if x in expr.index and x not in set(os_ids)]))
print('OS missing RNA:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in os_ids if i not in expr.index])
print('OS lines in RNA but not CRISPR:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in os_ids if i in expr.index and i not in pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', usecols=[0]).iloc[:, 0].tolist()])

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

rows = []
for g in expr.columns:
    a = expr.loc[os_, g].dropna().values
    b = expr.loc[comp, g].dropna().values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    d = am - bm
    fc = (2**am - 1) / (2**bm - 1) if (2**bm - 1) != 0 else np.nan
    if len(a) >= 2 and len(b) >= 2:
        t, p = stats.ttest_ind(a, b, equal_var=False)
    else:
        p = np.nan
    sym = g.split(' (')[0]
    rows.append(dict(Gene=sym, Column=g, OS_mean_log2TPM=am, OS_n=len(a), Comp_mean_log2TPM=bm, Comp_n=len(b),
                     DeltaLog2=d, TPM_FC=fc, p=p))
rna = pd.DataFrame(rows)
rna['q'] = bh(rna['p'].values)
rna = rna.sort_values('DeltaLog2', ascending=False)
rna.to_csv(OUT + 'celline_rna.csv', index=False)
print('genes', len(rna))
bar = rna[(rna.TPM_FC >= 1.5) & (rna.DeltaLog2 >= 0.5) & (rna.p < 0.05)]
print('cell-line RNA bar pass:', len(bar))
print(rna.head(20)[['Gene', 'OS_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'q']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'CCNE1', 'TP53', 'RB1', 'MDM2', 'CDK4', 'CDK6', 'FGFR1', 'PDGFRA', 'IGF1R', 'MET', 'SRC', 'RUNX2', 'SATB2', 'SP7', 'DLX5', 'KLF5', 'HNF1B', 'FOSL1', 'PAX8', 'CXADR', 'PARD6B', 'ITGB5', 'ITGB1']:
    r = rna[rna.Gene == g]
    print(g, r[['Gene', 'OS_n', 'DeltaLog2', 'TPM_FC', 'p', 'q']].to_string(index=False) if len(r) else 'NOT FOUND')
