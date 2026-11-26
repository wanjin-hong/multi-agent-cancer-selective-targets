"""TP53 / RB1 mutation check in the 17 osteosarcoma CRISPR lines.

Reads OmicsSomaticMutations.csv (dedup: the file repeats each variant
2-3 times), filters to the 17 OS lines, reports TP53 and RB1 variant
status per line, and splits the group TP53-LoF vs TP53-other to ask
whether TP53 status splits WWTR1/PTK2/TEAD1 dependency.
"""
import pandas as pd, numpy as np
from scipy import stats
import glob

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/osteosarcoma/'
mut_path = glob.glob('/home/hatch/workspace/browser_downloads/sess-1145554613/*OmicsSomaticMutations.csv')[0]
print('mutation file:', mut_path)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
os_ids = model.loc[model['OncotreePrimaryDisease'] == 'Osteosarcoma', 'ModelID'].tolist()
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
os_ = [x for x in os_ids if x in eff.index]
names = {r.ModelID: r.CellLineName for r in model.itertuples()}

cols = ['HugoSymbol', 'ModelID', 'VariantType', 'VepImpact', 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF']
mut = pd.read_csv(mut_path, usecols=cols, low_memory=False)
print('rows read:', len(mut))
mut = mut[mut.ModelID.isin(os_)]
print('rows in OS lines:', len(mut))
before = len(mut)
mut = mut.drop_duplicates()
print('after dedup:', len(mut), '(removed', before - len(mut), 'duplicate rows)')
mut.to_csv(OUT + 'os_mutations_dedup.csv', index=False)

for gene in ['TP53', 'RB1', 'MDM2', 'CDK4', 'ATRX', 'DLG2']:
    g = mut[mut.HugoSymbol == gene]
    print(f'\n{gene}: {len(g)} variant rows in {g.ModelID.nunique()} lines')
    if len(g):
        print(g[['ModelID', 'VariantType', 'VepImpact', 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF']].to_string(index=False))

# TP53 LoF definition: LikelyLoF == True, or VepImpact HIGH with nonsense/frameshift/splice
tp53 = mut[mut.HugoSymbol == 'TP53'].copy()
lof_ids = set(tp53.loc[(tp53.LikelyLoF == True) | ((tp53.VepImpact == 'HIGH') & (tp53.VariantType.isin(['Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site']))) , 'ModelID'])
any_tp53 = set(tp53.ModelID)
print('\nTP53-mutant lines:', len(any_tp53), 'of', len(os_), [names[i] for i in sorted(any_tp53)])
print('TP53-LoF lines:', len(lof_ids), [names[i] for i in sorted(lof_ids)])
print('TP53-WT (no TP53 variant) lines:', len(os_) - len(any_tp53), [names[i] for i in os_ if i not in any_tp53])
rb1 = mut[mut.HugoSymbol == 'RB1']
print('RB1-mutant lines:', len(set(rb1.ModelID)), [names[i] for i in sorted(set(rb1.ModelID))])

# dependency split by TP53 status (descriptive at n=17)
def gene_vals(gene):
    col = [c for c in eff.columns if c.startswith(gene + ' (')][0]
    return eff.loc[os_, col]
for gene in ['WWTR1', 'PTK2', 'TEAD1', 'YAP1', 'FERMT2', 'ITGAV', 'FGFR1', 'PDGFRA', 'RUNX2', 'TP53', 'RB1', 'MDM2', 'CDK4']:
    v = gene_vals(gene)
    lof = v.loc[[i for i in os_ if i in lof_ids]].dropna()
    wt = v.loc[[i for i in os_ if i not in lof_ids]].dropna()
    p = stats.ttest_ind(lof, wt, equal_var=False)[1] if len(lof) >= 2 and len(wt) >= 2 else np.nan
    print(f'{gene}: TP53-LoF mean {lof.mean():.3f} (n={len(lof)}, dep {(lof < -0.5).sum()}) vs other mean {wt.mean():.3f} (n={len(wt)}, dep {(wt < -0.5).sum()}) Welch p={p:.3f}')
