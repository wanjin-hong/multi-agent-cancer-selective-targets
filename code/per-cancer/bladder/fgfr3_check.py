"""FGFR3 (and key bladder-driver) mutation check in the 30 bladder CRISPR lines.

DepMap 24Q4 OmicsSomaticMutations.csv: deduplicate (the file repeats each
variant 2-3 times), then report, per gene, which of the 30 bladder lines
carry a variant, with VariantType / VepImpact / ProteinChange / Hotspot /
AMClass / LikelyLoF. FGFR3 mutation is the luminal-papillary proxy, since
Model.csv carries no luminal/basal annotation. Also cross-tab FGFR3
dependency (Chronos < -0.5) against FGFR3 mutation status.
"""
import pandas as pd, numpy as np

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/bladder/'
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
blad_ids = model.loc[model['OncotreePrimaryDisease'] == 'Bladder Urothelial Carcinoma', 'ModelID'].tolist()
blad = [x for x in blad_ids if x in eff.index]
name = dict(zip(model['ModelID'], model['CellLineName']))

mut = pd.read_csv('/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv',
                  usecols=['HugoSymbol', 'ModelID', 'VariantType', 'VepImpact', 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF'],
                  low_memory=False)
print('mutation rows raw:', len(mut))
mut = mut.drop_duplicates()
print('after dedup:', len(mut))
mut_blad = mut[mut['ModelID'].isin(blad)].copy()
print('bladder-line mutation rows:', len(mut_blad))
mut_blad.to_csv(OUT + 'bladder_mutations.csv', index=False)

GENES = ['FGFR3', 'ERBB2', 'PIK3CA', 'TP53', 'RB1', 'CDKN2A', 'FGFR1', 'KRAS', 'HRAS', 'NRAS', 'TERT', 'KDM6A', 'ARID1A', 'STAG2', 'ELF3', 'CDKN1A', 'TSC1']
for g in GENES:
    sub = mut_blad[mut_blad['HugoSymbol'] == g]
    if len(sub) == 0:
        print(g, ': no variants in bladder lines')
        continue
    lines = sorted(set(sub['ModelID']))
    print(f"\n{g}: {len(lines)}/{len(blad)} lines with a variant")
    for _, r in sub.sort_values(['ModelID', 'ProteinChange']).iterrows():
        print(f"  {name.get(r['ModelID'], r['ModelID'])}: {r['VariantType']} {r['ProteinChange']} impact={r['VepImpact']} hotspot={r['Hotspot']} AMClass={r['AMClass']} LoF={r['LikelyLoF']}")

# FGFR3 dependency vs mutation
fgfr3_col = [c for c in eff.columns if c.startswith('FGFR3 (')][0]
fgfr3_mut_lines = set(mut_blad.loc[mut_blad['HugoSymbol'] == 'FGFR3', 'ModelID'])
print('\nFGFR3 Chronos by line:')
for mid in blad:
    print(f"  {name[mid]:12s} Chronos={eff.loc[mid, fgfr3_col]:+.3f} dep={eff.loc[mid, fgfr3_col] < -0.5} FGFR3-mut={mid in fgfr3_mut_lines}")
# activating FGFR3 mutations: hotspot missense (VariantType is 'SNV' in this
# file, not 'Missense_Mutation' — the earlier filter using VariantType
# returned nothing and was wrong)
act = mut_blad[(mut_blad['HugoSymbol'] == 'FGFR3') & (mut_blad['Hotspot'] == True) & (mut_blad['VepImpact'] == 'MODERATE')]
print('\nFGFR3 missense lines:', sorted(set(act['ModelID'])))
print('their FGFR3 Chronos:', {name[m]: round(float(eff.loc[m, fgfr3_col]), 3) for m in sorted(set(act['ModelID']))})
print('non-mutant lines FGFR3 Chronos mean:', round(float(eff.loc[[m for m in blad if m not in fgfr3_mut_lines], fgfr3_col].mean()), 3))
