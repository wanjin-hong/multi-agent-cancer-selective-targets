"""Extract somatic mutations for neuroblastoma lines from DepMap 24Q4
OmicsSomaticMutations.csv (642 MB). The file repeats each variant 2-3
times — deduplicate on (ModelID, HugoSymbol, ProteinChange, VariantType).
Focus genes: ALK, PHOX2B, ATRX, TP53, NRAS, MYCN, BRAF, PTPN11, NF1.
ALK activating mutations (F1174, R1275, F1245) define the ALK-mutant
split for the Q4 check.
"""
import pandas as pd, numpy as np
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/neuroblastoma/'
MUT = '/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv'
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
nb_ids = set(model.loc[model['OncotreePrimaryDisease'] == 'Neuroblastoma', 'ModelID'])
names = model.set_index('ModelID')['CellLineName'].to_dict()
FOCUS = ['ALK', 'PHOX2B', 'ATRX', 'TP53', 'NRAS', 'MYCN', 'BRAF', 'PTPN11', 'NF1', 'ARID1A', 'ARID1B']
cols = ['ModelID', 'HugoSymbol', 'VariantType', 'VepImpact', 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF']
rows = []
n_raw = 0
for chunk in pd.read_csv(MUT, usecols=cols, chunksize=500000, low_memory=False):
    n_raw += len(chunk)
    sub = chunk[chunk['ModelID'].isin(nb_ids) & chunk['HugoSymbol'].isin(FOCUS)]
    if len(sub):
        rows.append(sub)
print('raw rows scanned:', n_raw)
mut = pd.concat(rows, ignore_index=True)
print('NB focus-gene rows before dedup:', len(mut))
mut = mut.drop_duplicates(subset=['ModelID', 'HugoSymbol', 'ProteinChange', 'VariantType'])
print('after dedup:', len(mut))
mut['CellLineName'] = mut['ModelID'].map(names)
mut = mut.sort_values(['HugoSymbol', 'CellLineName'])
mut.to_csv(OUT + 'nb_mutations_dedup.csv', index=False)
print(mut[['CellLineName', 'HugoSymbol', 'VariantType', 'VepImpact', 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF']].to_string(index=False))
print()
print('ALK-mutant NB lines (any ALK variant):')
alk = mut[mut.HugoSymbol == 'ALK']
print(sorted(alk['CellLineName'].unique().tolist()))
print('ALK activating (F1174/R1275/F1245) lines:')
act = alk[alk['ProteinChange'].str.contains('F1174|R1275|F1245', na=False)]
print(sorted(act['CellLineName'].unique().tolist()))
print('counts per gene:')
print(mut.groupby('HugoSymbol')['CellLineName'].nunique().to_string())
