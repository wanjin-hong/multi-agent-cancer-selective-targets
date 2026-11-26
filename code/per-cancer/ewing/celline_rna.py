"""Cell-line RNA: Ewing vs all other lines, DepMap 24Q4 expression matrix."""
import pandas as pd, numpy as np
from scipy import stats
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/ewing/'
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('expr shape', expr.shape, 'missing', int(expr.isna().sum().sum()))
ew_ids = model.loc[model['OncotreePrimaryDisease'] == 'Ewing Sarcoma', 'ModelID'].tolist()
ew = [x for x in ew_ids if x in expr.index]
comp = [x for x in expr.index if x not in set(ew_ids)]
print('Ewing RNA n=', len(ew), 'comparator n=', len(comp))
print('Ewing missing RNA:', [model.loc[model.ModelID == i, 'CellLineName'].iloc[0] for i in ew_ids if i not in expr.index])
nc = model[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(expr.index))]
print('non-cancerous lines in RNA comparator:', len(nc), nc['CellLineName'].tolist())
sub = model.set_index('ModelID')
print('Ewing RNA lines:', [sub.loc[i, 'CellLineName'] for i in ew])

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

rows = []
for g in expr.columns:
    a = expr.loc[ew, g].dropna().values
    b = expr.loc[comp, g].dropna().values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    d = am - bm
    fc = (2**am - 1) / (2**bm - 1) if (2**bm - 1) != 0 else np.nan
    if len(a) >= 2 and len(b) >= 2:
        t, p = stats.ttest_ind(a, b, equal_var=False)
    else:
        p = np.nan
    sym = g.split(' (')[0]
    rows.append(dict(Gene=sym, Column=g, Ewing_mean_log2TPM=am, Ewing_n=len(a), Comp_mean_log2TPM=bm, Comp_n=len(b),
                     DeltaLog2=d, TPM_FC=fc, p=p))
rna = pd.DataFrame(rows)
rna['q'] = bh(rna['p'].values)
rna = rna.sort_values('DeltaLog2', ascending=False)
rna.to_csv(OUT + 'celline_rna.csv', index=False)
print('genes', len(rna))
bar = rna[(rna.TPM_FC >= 1.5) & (rna.DeltaLog2 >= 0.5) & (rna.p < 0.05)]
print('cell-line RNA bar pass:', len(bar))
print(rna.head(25)[['Gene', 'Ewing_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'q']].to_string())
watch = ['TEAD1','TEAD2','TEAD3','TEAD4','YAP1','WWTR1','PTK2','PKN2','FLI1','EWSR1',
         'NKX2-2','NR0B1','EZH2','BCL11B','GLI1','STEAP1','CAV1','CD99','CCND1','CDK4',
         'EGFR','KLF5','FERMT2','HNF1B','FOSL1','SOX9','MT1E','ITGAV','TCF7L2','PARD3',
         'PARD6B','ARHGAP29','ITGB5','CDK6','TIPARP','KRAS','BRAF','MAPK1','TP63','MYBL2',
         'CXADR','SOX10','GATA3','PPARG','GRHL2','SKP2','ADAR','YEATS2','STAG2','STAG1',
         'TRIM8','E2F3','MYC','MYCN','IGF1R','KIT','PDGFRA','FGFR1','FGFR4','ERBB2','ERBB3',
         'PIK3CA','MDM2','MDM4','AURKB','PLK1','PARP1','BRCA1','BRCA2','RB1','TP53','CDKN2A',
         'ATP1A1','SLC7A11','GPX4','CHMP4B','UBE2N','GNG5','MARCHF5','TRIM8','SNAI2']
for g in watch:
    r = rna[rna.Gene == g]
    print(g, r[['Gene','DeltaLog2','TPM_FC','p','q']].to_string(index=False) if len(r) else 'NOT FOUND')
