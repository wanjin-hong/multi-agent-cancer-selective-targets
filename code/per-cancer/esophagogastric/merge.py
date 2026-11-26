import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/'
scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan['Symbol'] = scan['Gene'].str.split(' (', regex=False).str[0]
scan['Entrez'] = scan['Gene'].str.extract(r'\((\d+)\)$')[0].astype('Int64')
rna = pd.read_csv(OUT + 'celline_rna.csv')  # Gene=Symbol, TPM_FC, DeltaLog2, p, q
rna = rna.rename(columns={'p': 'CellRNA_p', 'q': 'CellRNA_q'})
tvn = pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv')
xc = pd.read_csv(OUT + 'tcga_rna.csv')  # cross-cancer, context only
rp = pd.read_csv(OUT + 'tcga_rppa.csv')
m = scan.merge(rna[['Gene', 'ESGA_mean_log2TPM', 'ESGA_n', 'Comp_mean_log2TPM', 'Comp_n', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
                    'Gastric_mean_log2TPM', 'Gastric_n', 'EsoGEJ_mean_log2TPM', 'EsoGEJ_n']],
               left_on='Symbol', right_on='Gene', how='left', suffixes=('', '_rna'))
tvn = tvn.rename(columns={'Gene': 'Symbol'})
m = m.merge(tvn[['Symbol', 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
                   'GDC_TvN_STAD_log2FC', 'GDC_TvN_STAD_p', 'GDC_TvN_ESCA_log2FC', 'GDC_TvN_ESCA_p']],
            on='Symbol', how='left')
xc = xc.rename(columns={'Gene': 'Symbol'})
m = m.merge(xc[['Symbol', 'TCGA_RNA_n_ESGA', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_ESGA', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']],
            on='Symbol', how='left')
rp = rp.rename(columns={'Gene': 'Symbol'})
m = m.merge(rp[['Symbol', 'RPPA_n_ESGA', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p']],
            on='Symbol', how='left')
m = m.rename(columns={'p': 'Dep_p', 'q': 'Dep_q'})
# bars — cross-cancer is context only, NOT part of the bar (HCC/PDAC/HNSC design)
m['CellRNA_pass'] = (m['TPM_FC'] >= 1.5) & (m['DeltaLog2'] >= 0.5) & (m['CellRNA_p'] < 0.05)
m['TCGA_TvN_pass'] = (m['GDC_TvN_log2FC'] >= 0.5) & (m['GDC_TvN_p'] < 0.05)
m['TCGA_XC_pass'] = (m['TCGA_RNA_log2FC'] >= 0.5) & (m['TCGA_RNA_p'] < 0.05)
m['RPPA_pass'] = (m['RPPA_delta'] >= 0.5) & (m['RPPA_p'] < 0.05)
m['OverExpr_pass'] = m['CellRNA_pass'] | m['TCGA_TvN_pass'] | m['RPPA_pass']
m['Strict_pass'] = (m['Enrichment'] <= -0.10) & (m['ESGA_dep_frac'] >= 0.25) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Relaxed_pass'] = (m['Enrichment'] <= -0.05) & (m['ESGA_dep_frac'] >= 0.20) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
print('strict:', int(m.Strict_pass.sum()), 'relaxed:', int(m.Relaxed_pass.sum()))
print('strict & overexpr:', int((m.Strict_pass & m.OverExpr_pass).sum()))
print('relaxed & overexpr:', int((m.Relaxed_pass & m.OverExpr_pass).sum()))
print('strict & overexpr breakdown: cellRNA only', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & ~m.TCGA_TvN_pass & ~m.RPPA_pass).sum()),
      '| TvN only', int(((m.Strict_pass & m.OverExpr_pass) & m.TCGA_TvN_pass & ~m.CellRNA_pass & ~m.RPPA_pass).sum()),
      '| both RNA', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & m.TCGA_TvN_pass).sum()),
      '| RPPA', int(((m.Strict_pass & m.OverExpr_pass) & m.RPPA_pass).sum()))
print('relaxed genes with no TvN data:', int(m[m.Relaxed_pass]['GDC_TvN_log2FC'].isna().sum()))
print('relaxed genes with no cell-line RNA data:', int(m[m.Relaxed_pass]['DeltaLog2'].isna().sum()))
print('relaxed genes with no XC data:', int(m[m.Relaxed_pass]['TCGA_RNA_log2FC'].isna().sum()))
# combined score across the relaxed set (pre-defined; TCGA component = pooled tumour-vs-normal log2FC)
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel = m[m.Relaxed_pass].copy()
rel['Score'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['GDC_TvN_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0)) + np.where(rel['RPPA_pass'].fillna(False), 1.0, 0.0)
rel['Tier'] = np.where(rel['Strict_pass'], 'A', 'B')
cand = rel[rel['OverExpr_pass']].copy().sort_values('Score', ascending=False)
print('candidates:', len(cand))
print(cand[['Symbol', 'Tier', 'ESGA_mean', 'ESGA_dep', 'ESGA_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'TPM_FC', 'DeltaLog2', 'GDC_TvN_log2FC', 'TCGA_RNA_log2FC', 'RPPA_delta', 'Score']].head(30).to_string())
cand.to_csv(OUT + 'candidates.csv', index=False)
m.to_csv(OUT + 'overexpression.csv', index=False)
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist() if len(cand) > 20 else 'n/a (list shorter than 20)')
print('small-denominator in top20 (ESGA_n<45):')
print(cand.head(20)[cand.head(20).ESGA_n < 45][['Symbol', 'ESGA_n', 'ESGA_dep']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'KLF5', 'SOX9', 'HNF1B', 'FOSL1', 'CDX2', 'ERBB2', 'CDH17', 'GATA6', 'CTNNB1', 'TCF7L2', 'GFPT1']:
    r = m[m.Symbol == g]
    cr = cand[cand.Symbol == g]
    print(g, 'enrich rank', int((m.Enrichment < r.iloc[0].Enrichment).sum()) + 1,
          'strict', bool(r.iloc[0].Strict_pass), 'relaxed', bool(r.iloc[0].Relaxed_pass),
          'overexpr', bool(r.iloc[0].OverExpr_pass),
          'cellRNA FC', round(r.iloc[0].TPM_FC, 2) if pd.notna(r.iloc[0].TPM_FC) else None,
          'TvN', round(r.iloc[0].GDC_TvN_log2FC, 2) if pd.notna(r.iloc[0].GDC_TvN_log2FC) else None,
          'XC', round(r.iloc[0].TCGA_RNA_log2FC, 2) if pd.notna(r.iloc[0].TCGA_RNA_log2FC) else None,
          'RPPA', round(r.iloc[0].RPPA_delta, 2) if pd.notna(r.iloc[0].RPPA_delta) else None,
          'candidate rank', (int(cand.index.get_loc(cr.index[0])) + 1 if len(cr) else None),
          'score', (round(float(cr.iloc[0].Score), 2) if len(cr) else None))
