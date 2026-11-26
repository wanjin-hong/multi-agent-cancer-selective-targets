import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/aml/'
scan = pd.read_csv(OUT + 'dependency_scan.csv')          # primary, pan-cancer comparator
hscan = pd.read_csv(OUT + 'dependency_scan_haem.csv')    # sensitivity, haematology comparator
rna = pd.read_csv(OUT + 'celline_rna.csv')               # Gene=Symbol, pan + haem RNA legs
rna = rna.rename(columns={'AML_n': 'AML_RNA_n', 'Comp_n': 'Comp_RNA_n', 'HaemComp_n': 'HaemComp_RNA_n'})
xc = pd.read_csv(OUT + 'tcga_rna.csv')                   # cross-cancer, THE TCGA leg (no normals exist for LAML)
m = scan.merge(rna[['Gene', 'AML_mean_log2TPM', 'AML_RNA_n', 'Comp_mean_log2TPM', 'Comp_RNA_n',
                    'DeltaLog2', 'TPM_FC', 'p', 'q',
                    'HaemComp_mean_log2TPM', 'HaemComp_RNA_n', 'DeltaLog2_haem', 'TPM_FC_haem', 'p_haem', 'q_haem']],
               on='Gene', how='left')
m = m.rename(columns={'p': 'CellRNA_p', 'q': 'CellRNA_q'})
xc = xc.rename(columns={'Gene': 'Symbol'})
m = m.merge(xc[['Symbol', 'TCGA_RNA_n_AML', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_AML', 'TCGA_RNA_mean_other',
                'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']],
            left_on='Gene', right_on='Symbol', how='left')
h = hscan[['Gene', 'Enrichment', 'FracDiff', 'AML_dep_frac', 'Comp_dep_frac', 'q_BH', 'EnrichRank', 'PanEssential']].copy()
h.columns = ['Gene', 'Haem_Enrichment', 'Haem_FracDiff', 'Haem_AML_dep_frac', 'Haem_Comp_dep_frac',
             'Haem_q', 'Haem_EnrichRank', 'Haem_PanEssential']
m = m.merge(h, on='Gene', how='left')
m = m.rename(columns={'p_welch': 'Dep_p', 'q_BH': 'Dep_q'})
# bars — TCGA leg is the cross-cancer comparison (design c); NO protein leg exists for AML
m['CellRNA_pass'] = (m['TPM_FC'] >= 1.5) & (m['DeltaLog2'] >= 0.5) & (m['CellRNA_p'] < 0.05)
m['TCGA_XC_pass'] = (m['TCGA_RNA_log2FC'] >= 0.5) & (m['TCGA_RNA_p'] < 0.05)
m['OverExpr_pass'] = m['CellRNA_pass'] | m['TCGA_XC_pass']
m['Strict_pass'] = (m['Enrichment'] <= -0.10) & (m['AML_dep_frac'] >= 0.25) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Relaxed_pass'] = (m['Enrichment'] <= -0.05) & (m['AML_dep_frac'] >= 0.20) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Haem_Strict_pass'] = (m['Haem_Enrichment'] <= -0.10) & (m['AML_dep_frac'] >= 0.25) & (m['Haem_FracDiff'] >= 0.10) & (~m['Haem_PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Haem_Relaxed_pass'] = (m['Haem_Enrichment'] <= -0.05) & (m['AML_dep_frac'] >= 0.20) & (m['Haem_FracDiff'] >= 0.10) & (~m['Haem_PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
print('strict:', int(m.Strict_pass.sum()), 'relaxed:', int(m.Relaxed_pass.sum()))
print('strict & overexpr:', int((m.Strict_pass & m.OverExpr_pass).sum()))
print('relaxed & overexpr:', int((m.Relaxed_pass & m.OverExpr_pass).sum()))
print('haem strict:', int(m.Haem_Strict_pass.sum()), 'haem relaxed:', int(m.Haem_Relaxed_pass.sum()))
print('haem relaxed & overexpr:', int((m.Haem_Relaxed_pass & m.OverExpr_pass).sum()))
print('strict & overexpr breakdown: cellRNA only', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & ~m.TCGA_XC_pass).sum()),
      '| XC only', int(((m.Strict_pass & m.OverExpr_pass) & m.TCGA_XC_pass & ~m.CellRNA_pass).sum()),
      '| both RNA', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & m.TCGA_XC_pass).sum()))
print('relaxed genes with no XC data:', int(m[m.Relaxed_pass]['TCGA_RNA_log2FC'].isna().sum()))
print('relaxed genes with no cell-line RNA data:', int(m[m.Relaxed_pass]['DeltaLog2'].isna().sum()))
# combined score across the relaxed set (pre-defined; TCGA component = cross-cancer log2FC; no protein bonus exists)
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel = m[m.Relaxed_pass].copy()
rel['Score'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['TCGA_RNA_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0))
rel['Tier'] = np.where(rel['Strict_pass'], 'A', 'B')
cand = rel[rel['OverExpr_pass']].copy().sort_values('Score', ascending=False)
print('candidates (pan-cancer primary):', len(cand))
print(cand[['Gene', 'Tier', 'AML_mean', 'AML_dep', 'AML_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
            'TPM_FC', 'DeltaLog2', 'TCGA_RNA_log2FC', 'Score']].head(30).to_string(index=False))
cand.to_csv(OUT + 'candidates.csv', index=False)
# haem-comparator candidate list (same score formula, haem dependency components)
hrel = m[m.Haem_Relaxed_pass].copy()
hrel['Haem_Score'] = z(-hrel['Haem_Enrichment']) + z(hrel['Haem_FracDiff']) + z(hrel['TCGA_RNA_log2FC'].fillna(0)) + z(hrel['DeltaLog2'].fillna(0))
hrel['Tier'] = np.where(hrel['Haem_Strict_pass'], 'A', 'B')
hcand = hrel[hrel['OverExpr_pass']].copy().sort_values('Haem_Score', ascending=False)
print('candidates (haem sensitivity):', len(hcand))
print(hcand[['Gene', 'Tier', 'AML_mean', 'AML_dep', 'AML_n', 'Haem_Enrichment', 'Haem_FracDiff',
             'TPM_FC', 'DeltaLog2', 'TCGA_RNA_log2FC', 'Haem_Score']].head(30).to_string(index=False))
hcand.to_csv(OUT + 'candidates_haem.csv', index=False)
m.to_csv(OUT + 'overexpression.csv', index=False)
print('ranks 21-25:', cand.iloc[20:25]['Gene'].tolist() if len(cand) > 20 else 'n/a')
print('haem ranks 21-25:', hcand.iloc[20:25]['Gene'].tolist() if len(hcand) > 20 else 'n/a')
print('small-denominator in top20 (AML_n<28):')
print(cand.head(20)[cand.head(20).AML_n < 28][['Gene', 'AML_n', 'AML_dep']].to_string(index=False))
print('overlap top20 pan vs haem:', len(set(cand.head(20)['Gene']) & set(hcand.head(20)['Gene'])))
print('pan top20 not in haem top20:', sorted(set(cand.head(20)['Gene']) - set(hcand.head(20)['Gene'])))
print('haem top20 not in pan top20:', sorted(set(hcand.head(20)['Gene']) - set(cand.head(20)['Gene'])))
print('pan candidates not haem candidates:', sorted(set(cand['Gene']) - set(hcand['Gene']))[:50])
print('haem candidates not pan candidates:', sorted(set(hcand['Gene']) - set(cand['Gene']))[:50])
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'MYB', 'CBFB',
          'SPI1', 'RUNX1', 'CEBPA', 'CDK6', 'FLT3', 'IDH1', 'IDH2', 'NPM1', 'KMT2A', 'BCL2', 'CD33',
          'IL3RA', 'KLF5', 'FOSL1', 'HNF1B', 'PARD6B', 'CXADR', 'JAK2', 'MCL1', 'BCL2L1', 'KIT',
          'DOT1L', 'MEN1', 'MEIS1', 'HOXA9', 'GATA2', 'ZEB2', 'LMO2', 'CSF2RB', 'GFI1', 'DNM2', 'NAMPT']:
    r = m[m.Gene == g]
    if len(r) == 0:
        print(g, 'NOT IN CRISPR SCAN (absent from CRISPRGeneEffect.csv)')
        continue
    cr = cand[cand.Gene == g]
    hr = hcand[hcand.Gene == g]
    print(g, 'enrich rank', int(r.iloc[0].EnrichRank),
          'strict', bool(r.iloc[0].Strict_pass), 'relaxed', bool(r.iloc[0].Relaxed_pass),
          'overexpr', bool(r.iloc[0].OverExpr_pass),
          'cellRNA FC', round(r.iloc[0].TPM_FC, 2) if pd.notna(r.iloc[0].TPM_FC) else None,
          'XC', round(r.iloc[0].TCGA_RNA_log2FC, 2) if pd.notna(r.iloc[0].TCGA_RNA_log2FC) else None,
          'haem rank', int(r.iloc[0].Haem_EnrichRank),
          'haem relaxed', bool(r.iloc[0].Haem_Relaxed_pass),
          'candidate rank', (int(cand.index.get_loc(cr.index[0])) + 1 if len(cr) else None),
          'score', (round(float(cr.iloc[0].Score), 2) if len(cr) else None),
          'haem cand rank', (int(hcand.index.get_loc(hr.index[0])) + 1 if len(hr) else None))
