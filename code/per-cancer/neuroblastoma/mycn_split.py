"""MYCN split: does MYCN dependency sit only in MYCN-high lines?

MYCN RNA (DepMap 24Q4 expression matrix) is cleanly bimodal across the
36 NB lines with RNA: 26 lines at log2(TPM+1) >= 8.14, a gap, then
SK-N-FI at 5.74 and 9 lines <= 2.72. Rule (stated before looking at
dependency): MYCN-high = MYCN log2(TPM+1) >= 8.0; MYCN-low = below.
This is an RNA-expression split, NOT a copy-number annotation — DepMap
24Q4 on disk carries no MYCN amplification call, so no line is labelled
"MYCN-amplified" here; the high group is consistent with the
literature-known amplified lines (SK-N-BE(2), KELLY, IMR-32, NGP,
CHP-134, NMB, NB-1643, CHP-212) but the label is expression-based.
6 of the 39 CRISPR lines have no RNA data (KP-N-YS, NB10, NB13, NB17,
NB6, NB7) and are unclassified.
"""
import pandas as pd, numpy as np
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/neuroblastoma/'
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
names = model.set_index('ModelID')['CellLineName'].to_dict()
nb_ids = model.loc[model['OncotreePrimaryDisease'] == 'Neuroblastoma', 'ModelID'].tolist()
nb = [x for x in nb_ids if x in eff.index]
mycn = expr['MYCN (4613)']
hi = [i for i in nb if i in mycn.index and mycn[i] >= 8.0]
lo = [i for i in nb if i in mycn.index and mycn[i] < 8.0]
un = [i for i in nb if i not in mycn.index]
print('MYCN-high CRISPR lines:', len(hi), [names[i] for i in hi])
print('MYCN-low CRISPR lines:', len(lo), [names[i] for i in lo])
print('unclassified (no RNA):', len(un), [names[i] for i in un])

def frac(gene, lines):
    a = eff.loc[lines, gene].dropna()
    return f"{(a < -0.5).sum()}/{len(a)}", a.mean()

for g in ['MYCN (4613)', 'MYC (4609)', 'MYCL (4610)', 'ALK (238)', 'PHOX2B (8929)',
          'ISL1 (3670)', 'GATA3 (2625)', 'HAND2 (9464)', 'PHOX2A (401)', 'CDK6 (1021)',
          'CCND1 (595)', 'AURKA (6790)', 'MDM2 (4193)', 'BCL2 (596)', 'TEAD1 (7003)',
          'YAP1 (10413)', 'WWTR1 (25937)', 'PTK2 (5747)', 'PKN2 (5586)', 'LIN28B (389421)',
          'SOX11 (6664)', 'HDAC2 (3066)', 'EZH2 (2146)', 'CHEK1 (1111)', 'TERT (7015)']:
    fh, mh = frac(g, hi)
    fl, ml = frac(g, lo)
    print(f"{g.split(' (')[0]:8s} MYCN-high dep {fh} mean {mh:.3f} | MYCN-low dep {fl} mean {ml:.3f}")

# which MYCN-dependent lines are high vs low
dep = eff.loc[nb, 'MYCN (4613)']
dep_lines = dep[dep < -0.5].index.tolist()
print('\nMYCN-dependent lines:', len(dep_lines))
print('  high:', sum(1 for i in dep_lines if i in hi), ' low:', sum(1 for i in dep_lines if i in lo),
      ' unclassified:', sum(1 for i in dep_lines if i in un))
print('  low/unclassified MYCN-dependent:', [names[i] for i in dep_lines if i not in hi])
print('MYCN Chronos in low lines:', {names[i]: round(float(dep[i]), 3) for i in lo})
print('MYCN Chronos in unclassified lines:', {names[i]: round(float(dep[i]), 3) for i in un})
# ALK-mutant split (activating mutations from nb_mutations_dedup.csv)
mut = pd.read_csv(OUT + 'nb_mutations_dedup.csv')
act = mut[(mut.HugoSymbol == 'ALK') & (mut.ProteinChange.str.contains('F1174|R1275|F1245', na=False))]
act_ids = model.loc[model.CellLineName.isin(act['CellLineName']), 'ModelID'].tolist()
act_crispr = [i for i in act_ids if i in eff.index]
print('\nALK-activating CRISPR lines:', len(act_crispr), [names[i] for i in act_crispr])
alk_dep = eff.loc[act_crispr, 'ALK (238)']
print('ALK dependent among ALK-mutant:', int((alk_dep < -0.5).sum()), '/', len(alk_dep),
      'mean', round(float(alk_dep.mean()), 3))
wt = [i for i in nb if i not in act_crispr]
alk_wt = eff.loc[wt, 'ALK (238)']
print('ALK dependent among other NB lines:', int((alk_wt < -0.5).sum()), '/', len(alk_wt),
      'mean', round(float(alk_wt.mean()), 3))
