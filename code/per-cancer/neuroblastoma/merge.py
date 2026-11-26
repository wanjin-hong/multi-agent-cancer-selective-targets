import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/'
scan = pd.read_csv(OUT + 'dependency_scan_raw.csv')
scan['Symbol'] = scan['Gene'].str.split(' (', regex=False).str[0]
scan['Entrez'] = scan['Gene'].str.extract(r'\((\d+)\)$')[0].astype('Int64')
rna = pd.read_csv(OUT + 'celline_rna.csv')  # Gene=Symbol, TPM_FC, DeltaLog2, p, q
rna = rna.rename(columns={'p': 'CellRNA_p', 'q': 'CellRNA_q'})
xc = pd.read_csv(OUT + 'tcga_cross_cancer.csv')  # cross-cancer patient RNA = the TCGA leg (design c)
m = scan.merge(rna[['Gene', 'NB_mean_log2TPM', 'NB_n', 'Comp_mean_log2TPM', 'Comp_n', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q']],
               left_on='Symbol', right_on='Gene', how='left', suffixes=('', '_rna'))
xc = xc.rename(columns={'Gene': 'Symbol'})
m = m.merge(xc[['Symbol', 'TCGA_XC_n_NBL', 'TCGA_XC_n_other', 'TCGA_XC_mean_NBL', 'TCGA_XC_mean_other', 'TCGA_XC_log2FC', 'TCGA_XC_p', 'TCGA_XC_q']],
            on='Symbol', how='left')
m = m.rename(columns={'p': 'Dep_p', 'q': 'Dep_q'})
# bars — no protein leg exists for neuroblastoma (no NBL RPPA/CPTAC profile),
# so the +1.0 protein bonus is dropped (glioma precedent) and the
# over-expression bar is cell-line RNA OR patient cross-cancer RNA.
m['CellRNA_pass'] = (m['TPM_FC'] >= 1.5) & (m['DeltaLog2'] >= 0.5) & (m['CellRNA_p'] < 0.05)
m['TCGA_XC_pass'] = (m['TCGA_XC_log2FC'] >= 0.5) & (m['TCGA_XC_p'] < 0.05)
m['OverExpr_pass'] = m['CellRNA_pass'] | m['TCGA_XC_pass']
m['Strict_pass'] = (m['Enrichment'] <= -0.10) & (m['NB_dep_frac'] >= 0.25) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Relaxed_pass'] = (m['Enrichment'] <= -0.05) & (m['NB_dep_frac'] >= 0.20) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
print('strict:', int(m.Strict_pass.sum()), 'relaxed:', int(m.Relaxed_pass.sum()))
print('strict & overexpr:', int((m.Strict_pass & m.OverExpr_pass).sum()))
print('relaxed & overexpr:', int((m.Relaxed_pass & m.OverExpr_pass).sum()))
print('strict & overexpr breakdown: cellRNA only', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & ~m.TCGA_XC_pass).sum()),
      '| XC only', int(((m.Strict_pass & m.OverExpr_pass) & m.TCGA_XC_pass & ~m.CellRNA_pass).sum()),
      '| both RNA', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & m.TCGA_XC_pass).sum()))
print('relaxed genes with no XC data:', int(m[m.Relaxed_pass]['TCGA_XC_log2FC'].isna().sum()))
print('relaxed genes with no cell-line RNA data:', int(m[m.Relaxed_pass]['DeltaLog2'].isna().sum()))
# combined score across the relaxed set (pre-defined; TCGA component = cross-cancer log2FC)
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel = m[m.Relaxed_pass].copy()
rel['Score'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['TCGA_XC_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0))
rel['Tier'] = np.where(rel['Strict_pass'], 'A', 'B')
cand = rel[rel['OverExpr_pass']].copy().sort_values('Score', ascending=False)
print('candidates:', len(cand))
print(cand[['Symbol', 'Tier', 'NB_mean', 'NB_dep', 'NB_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'TPM_FC', 'DeltaLog2', 'TCGA_XC_log2FC', 'Score']].head(30).to_string())
cand.to_csv(OUT + 'candidates.csv', index=False)
m.to_csv(OUT + 'overexpression.csv', index=False)
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist() if len(cand) > 20 else 'n/a (list shorter than 20)')
print('small-denominator in top20 (NB_n<39):')
print(cand.head(20)[cand.head(20).NB_n < 39][['Symbol', 'NB_n', 'NB_dep']].to_string())
CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'MYCN', 'MYC', 'MYCL', 'ALK', 'PHOX2B', 'ATRX', 'KLF5', 'FERMT2', 'EGFR', 'HNF1B',
         'FOSL1', 'SOX9', 'MT1E', 'ITGAV', 'TCF7L2', 'PARD3', 'PARD6B',
         'ARHGAP29', 'ITGB5', 'CDK6', 'TIPARP', 'KRAS', 'HAND2', 'GATA3',
         'ISL1', 'ASCL1', 'BCL2', 'MDM2', 'AURKA', 'CHEK1', 'EZH2', 'LIN28B',
         'TERT', 'SOX11', 'DLK1', 'DBH', 'TH', 'CHGA', 'SYP', 'RET',
         'NTRK1', 'NTRK2', 'NGFR']
for g in CHECK:
    r = m[m.Symbol == g]
    cr = cand[cand.Symbol == g]
    print(g, 'enrich rank', int((m.Enrichment < r.iloc[0].Enrichment).sum()) + 1,
          'strict', bool(r.iloc[0].Strict_pass), 'relaxed', bool(r.iloc[0].Relaxed_pass),
          'overexpr', bool(r.iloc[0].OverExpr_pass),
          'cellRNA FC', round(r.iloc[0].TPM_FC, 2) if pd.notna(r.iloc[0].TPM_FC) else None,
          'XC', round(r.iloc[0].TCGA_XC_log2FC, 2) if pd.notna(r.iloc[0].TCGA_XC_log2FC) else None,
          'candidate rank', (int(cand.index.get_loc(cr.index[0])) + 1 if len(cr) else None),
          'score', (round(float(cr.iloc[0].Score), 2) if len(cr) else None))
