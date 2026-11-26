"""Ewing fusion / STAG2 checks.
1. Per-line FLI1 Chronos across the 27 Ewing lines (which 7 are non-dependent?).
2. Somatic mutations in the 27 lines from OmicsSomaticMutations.csv (dedup):
   STAG2 LoF vs STAG1 dependency (synthetic-lethality hypothesis), TP53,
   CDKN2A, and EWSR1/FLI1/ERG point mutations (fusions are not in this file).
Fusion annotation: Model.csv carries none, and the mutation file records
point mutations/indels, not fusions — no EWSR1-FLI1 vs EWSR1-ERG split is
possible from the data on disk. The literature-known split is noted in the
report (CADO-ES1 is the established EWSR1-ERG line in this panel).
"""
import pandas as pd, numpy as np, os

OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/ewing/'
BASE = '/home/hatch/workspace/research_notes/'
cols = ['FLI1 (2313)', 'EWSR1 (2130)', 'ERG (2078)', 'STAG1 (10274)', 'STAG2 (10735)']
hdr = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', nrows=0)
idx = [0] + [list(hdr.columns).index(c) for c in cols]
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0, usecols=idx)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
ew_ids = model.loc[model['OncotreePrimaryDisease'] == 'Ewing Sarcoma', 'ModelID'].tolist()
ew = [x for x in ew_ids if x in eff.index]
sub = model.set_index('ModelID')
rows = []
for i in ew:
    rows.append(dict(ModelID=i, Line=sub.loc[i, 'CellLineName'],
                     FLI1=eff.loc[i, 'FLI1 (2313)'], EWSR1=eff.loc[i, 'EWSR1 (2130)'],
                     ERG=eff.loc[i, 'ERG (2078)'], STAG1=eff.loc[i, 'STAG1 (10274)'],
                     STAG2=eff.loc[i, 'STAG2 (10735)']))
df = pd.DataFrame(rows)
df.to_csv(OUT + 'fusion_check.csv', index=False)
print(df.to_string(index=False))
print()
print('FLI1-dependent (<-0.5):', int((df.FLI1 < -0.5).sum()), '/ 27; non-dependent lines:',
      df.loc[df.FLI1 >= -0.5, 'Line'].tolist())
print('ERG-dependent:', int((df.ERG < -0.5).sum()))
print('STAG1-dependent:', int((df.STAG1 < -0.5).sum()), '/ 27; STAG2-dependent:',
      int((df.STAG2 < -0.5).sum()))

# ---- mutation file: extract rows for the 27 lines ----
mut_path = '/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv'
print('mutation file exists:', os.path.exists(mut_path), os.path.getsize(mut_path) if os.path.exists(mut_path) else '')
# grep the 27 ModelIDs (ACH- ids) — header first
ids = '|'.join(ew)
with open(OUT + 'ewing_mutations.csv', 'w') as out:
    with open(mut_path) as f:
        header = f.readline()
        out.write(header)
        n = 0
        for line in f:
            # Columns: Unnamed:0, SequencingID, ModelID, ... — ModelID is 3rd (index 2)
            parts = line.split(',', 3)
            if len(parts) >= 3 and parts[2] in set(ew):
                out.write(line); n += 1
print('mutation rows for 27 lines (with dupes):', n)
mut = pd.read_csv(OUT + 'ewing_mutations.csv', low_memory=False)
print('cols:', list(mut.columns))
print('rows:', len(mut))
mut_dedup = mut.drop_duplicates(subset=['HugoSymbol', 'ModelID', 'VariantType', 'ProteinChange'])
print('dedup rows:', len(mut_dedup))
print()
print('STAG2 rows:')
print(mut_dedup[mut_dedup.HugoSymbol == 'STAG2'][['HugoSymbol', 'ModelID', 'VariantType', 'VepImpact', 'ProteinChange', 'LikelyLoF']].to_string(index=False))
mut_dedup = mut_dedup.merge(df[['ModelID', 'Line', 'STAG1']], on='ModelID', how='left')
stag2_lines = mut_dedup[mut_dedup.HugoSymbol == 'STAG2']['ModelID'].unique().tolist()
print('lines with any STAG2 variant:', [sub.loc[i, 'CellLineName'] for i in stag2_lines])
lof = mut_dedup[(mut_dedup.HugoSymbol == 'STAG2') & (mut_dedup.LikelyLoF == True)]
print('lines with STAG2 LikelyLoF:', [sub.loc[i, 'CellLineName'] for i in lof['ModelID'].unique().tolist()])
print('STAG1 Chronos in STAG2-LoF lines:', df[df.ModelID.isin(lof['ModelID'].unique())][['Line', 'STAG1']].to_string(index=False))
print('STAG1 Chronos in other lines: mean', df[~df.ModelID.isin(lof['ModelID'].unique())]['STAG1'].mean(),
      'dependent:', int((df[~df.ModelID.isin(lof['ModelID'].unique())]['STAG1'] < -0.5).sum()), '/',
 int((~df.ModelID.isin(lof['ModelID'].unique())).sum()))
print()
for g in ['TP53', 'CDKN2A', 'EWSR1', 'FLI1', 'ERG', 'CDK4', 'CCND1', 'KRAS', 'BRAF', 'PIK3CA', 'PTEN', 'RB1', 'SMARCA4']:
    gg = mut_dedup[mut_dedup.HugoSymbol == g]
    print(g, 'variant rows (dedup):', len(gg),
          '| lines:', sorted(sub.loc[i, 'CellLineName'] for i in gg['ModelID'].unique().tolist()) if len(gg) else '')
mut_dedup.to_csv(OUT + 'ewing_mutations_dedup.csv', index=False)
