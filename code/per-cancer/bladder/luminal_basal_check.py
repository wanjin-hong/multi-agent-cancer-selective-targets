"""Luminal/basal expression split check for the bladder CRISPR group.

Model.csv carries no luminal/basal annotation (all 35 lines carry subtype
'Bladder Urothelial Carcinoma'), so this is a descriptive expression check
only, not a second primary analysis and not a stratified scan: per-line mean
log2(TPM+1) of luminal markers (GATA3, FOXA1, PPARG, UPK1B, KRT20, UPK2,
UPK3A, FGFR3) minus basal markers (KRT5, KRT6A, KRT14, TP63, CD44, KRT16,
KRT17) across the 33 bladder lines with RNA data. The split is reported
descriptively; at 5 luminal-leaning vs 14 basal-leaning lines no stratified
scan is powered, and the marker set is an unvalidated proxy — the three
FGFR3-activating-mutant lines (J82, UM-UC-16, UM-UC-6) all score
basal-leaning, against the usual luminal-papillary association.
"""
import pandas as pd, numpy as np

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/bladder/'
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
blad_ids = model.loc[model['OncotreePrimaryDisease'] == 'Bladder Urothelial Carcinoma', 'ModelID'].tolist()
blad = [x for x in blad_ids if x in expr.index]
name = dict(zip(model['ModelID'], model['CellLineName']))

def col(g):
    c = [x for x in expr.columns if x.startswith(g + ' (')]
    return c[0] if c else None

lum = ['GATA3', 'FOXA1', 'PPARG', 'UPK1B', 'KRT20', 'UPK2', 'UPK3A', 'FGFR3']
bas = ['KRT5', 'KRT6A', 'KRT14', 'TP63', 'CD44', 'KRT16', 'KRT17']
rows = []
for mid in blad:
    l = np.mean([expr.loc[mid, col(g)] for g in lum if col(g)])
    b = np.mean([expr.loc[mid, col(g)] for g in bas if col(g)])
    rows.append(dict(Line=name[mid], ModelID=mid, Luminal=l, Basal=b, Diff=l - b))
df = pd.DataFrame(rows).sort_values('Diff', ascending=False)
df.to_csv(OUT + 'luminal_basal_check.csv', index=False)
print(df.to_string(index=False))
print('\nluminal-leaning (Diff>1):', int((df.Diff > 1).sum()),
      'basal-leaning (Diff<-1):', int((df.Diff < -1).sum()),
      'intermediate:', int((df.Diff.abs() <= 1).sum()))
mut = pd.read_csv(OUT + 'bladder_mutations.csv')
fg = set(mut[mut.HugoSymbol == 'FGFR3']['ModelID'])
print('FGFR3-mutant lines FGFR3 log2TPM:', {name[m]: round(float(expr.loc[m, col('FGFR3')]), 2) for m in fg if m in expr.index})
print('FGFR3 log2TPM all bladder lines median:', round(float(expr.loc[blad, col('FGFR3')].median()), 2))
