"""Mutation stratification for the AML repeat: which of the 28 CRISPR-panel AML
lines carry the known AML driver mutations, and does dependency track them?

Source: DepMap 24Q4 OmicsSomaticMutations.csv (642 MB; the file repeats each
variant 2-3 times — deduplicate on all columns). Columns used: HugoSymbol,
ModelID, VariantType, VepImpact, ProteinChange, Hotspot, AMClass, LikelyLoF.
KMT2A rearrangements are fusions and are NOT in this file — KMT2A-r status
below is the literature-annotated line list (flagged as such, not file-derived):
MV4-11 (KMT2A-AFF1), MOLM-13, MOLM-14, THP-1, NOMO-1, MONO-MAC-1
(all KMT2A-MLLT3). FLT3-ITD in this file: VariantType INS with a ProteinChange
containing 'dup' (tandem duplication in the juxtamembrane domain, exons 14-15).
"""
import pandas as pd, numpy as np
BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/aml/'
MUT = '/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv'
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
aml_ids = model.loc[model['OncotreePrimaryDisease'] == 'Acute Myeloid Leukemia', 'ModelID'].tolist()
aml = [x for x in aml_ids if x in eff.index]
name = dict(zip(model['ModelID'], model['CellLineName']))
print('AML CRISPR lines:', len(aml))

GENES = ['FLT3', 'NPM1', 'IDH1', 'IDH2', 'DNMT3A', 'TET2', 'RUNX1', 'TP53',
         'NRAS', 'KRAS', 'KIT', 'CEBPA', 'WT1', 'ASXL1', 'STAG2', 'SRSF2',
         'U2AF1', 'SF3B1', 'PTPN11', 'JAK2', 'CBL', 'BCOR', 'EZH2', 'GATA2']
mut = pd.read_csv(MUT, usecols=['HugoSymbol', 'ModelID', 'VariantType', 'VepImpact',
                                 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF'],
                  low_memory=False)
print('mutation rows (raw):', len(mut))
mut = mut.drop_duplicates()
print('mutation rows (dedup):', len(mut))
m = mut[(mut['ModelID'].isin(aml)) & (mut['HugoSymbol'].isin(GENES))].copy()
print('rows in AML panel x genes of interest:', len(m))
m.to_csv(OUT + 'aml_mutations.csv', index=False)

# FLT3 detail: all FLT3 variants in the panel
flt3 = m[m['HugoSymbol'] == 'FLT3']
print('\nFLT3 variants in the 28 AML lines:')
print(flt3[['ModelID', 'VariantType', 'ProteinChange', 'Hotspot', 'AMClass', 'VepImpact']].to_string(index=False))
flt3 = flt3.copy()
flt3['Line'] = flt3['ModelID'].map(name)
print(flt3[['Line', 'VariantType', 'ProteinChange', 'Hotspot']].to_string(index=False))

def dep(gene, lines):
    col = [c for c in eff.columns if c.startswith(gene + ' (')]
    if not col:
        return None
    s = eff.loc[lines, col[0]]
    return s

# FLT3-ITD lines: VariantType 'insertion' with an 'ins' ProteinChange in the
# juxtamembrane domain (exons 14-15); the three insertions are all Hotspot=True.
# (An earlier 'INS'+'dup' rule found none — corrected 21 Nov 2026.)
itd = sorted(set(flt3[(flt3['VariantType'] == 'insertion') & (flt3['ProteinChange'].str.contains('ins', na=False))]['ModelID'].tolist()))
flt3_mut = sorted(set(flt3['ModelID'].tolist()))
print('\nFLT3-ITD lines (INS+dup):', [name[i] for i in itd])
print('Any FLT3 mutation:', [name[i] for i in flt3_mut])
for gene in ['FLT3', 'IDH1', 'IDH2', 'NPM1', 'KMT2A', 'BCL2', 'CD33', 'DOT1L', 'MEN1', 'MCL1']:
    s = dep(gene, aml)
    if s is None:
        print(gene, 'not in CRISPR matrix'); continue
    print(f"\n{gene} Chronos in FLT3-ITD lines (n={len(itd)}):",
          [round(float(v), 2) for v in s.loc[itd].values] if itd else 'none',
          ' dependent:', int((s.loc[itd] < -0.5).sum()) if itd else 0)
    print(f"{gene} Chronos in other FLT3-mutant lines:",
          [round(float(v), 2) for v in s.loc[[i for i in flt3_mut if i not in itd]].values] if len(flt3_mut) > len(itd) else 'none')
    print(f"{gene} Chronos in FLT3-WT lines (n={len(aml) - len(flt3_mut)}):",
          'dependent:', int((s.loc[[i for i in aml if i not in flt3_mut]] < -0.5).sum()),
          'mean:', round(float(s.loc[[i for i in aml if i not in flt3_mut]].mean()), 3))

# per-gene mutation table for the report
print('\nMutation counts per gene in the 28-line panel (dedup, any variant):')
print(m.groupby('HugoSymbol')['ModelID'].nunique().sort_values(ascending=False).to_string())
print('\nNPM1 variants:')
print(m[m['HugoSymbol'] == 'NPM1'][['ModelID', 'VariantType', 'ProteinChange']].to_string(index=False))
print('\nIDH1/IDH2 variants:')
print(m[m['HugoSymbol'].isin(['IDH1', 'IDH2'])][['ModelID', 'VariantType', 'ProteinChange', 'Hotspot']].to_string(index=False))
print('\nDNMT3A / TET2 / RUNX1 / TP53 / NRAS / KRAS / KIT / CEBPA / WT1 lines:')
for g in ['DNMT3A', 'TET2', 'RUNX1', 'TP53', 'NRAS', 'KRAS', 'KIT', 'CEBPA', 'WT1', 'STAG2']:
    sub = m[m['HugoSymbol'] == g]
    print(g, ':', [name[i] for i in sorted(set(sub['ModelID']))], '|',
          sub['ProteinChange'].tolist() if len(sub) else '')

# KMT2A-r literature list: dependency of DOT1L / MEN1 / KMT2A in those lines
kmt2a_r = ['MV4-11', 'MOLM-13', 'MOLM-14', 'THP-1', 'NOMO-1', 'MONO-MAC-1']
inv = {v: k for k, v in name.items()}
k_ids = [inv[x] for x in kmt2a_r if x in inv]
k_ids = [i for i in k_ids if i in aml]
print('\nKMT2A-r lines in CRISPR panel (literature-annotated):', [name[i] for i in k_ids])
for gene in ['KMT2A', 'DOT1L', 'MEN1', 'MEIS1', 'HOXA9', 'FLT3', 'BCL2']:
    s = dep(gene, aml)
    if s is None:
        print(gene, 'not in CRISPR matrix'); continue
    a = s.loc[k_ids]; b = s.loc[[i for i in aml if i not in k_ids]]
    print(f"{gene}: KMT2A-r mean {a.mean():.3f} dep {int((a < -0.5).sum())}/{len(a)} | other mean {b.mean():.3f} dep {int((b < -0.5).sum())}/{len(b)}")
