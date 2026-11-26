import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/osteosarcoma/'
scan = pd.read_csv(OUT + 'dependency_scan.csv')  # Gene='SYM (Entrez)', Symbol plain
rna = pd.read_csv(OUT + 'celline_rna.csv')  # Gene=Symbol, TPM_FC, DeltaLog2, p, q
rna = rna.rename(columns={'p': 'CellRNA_p', 'q': 'CellRNA_q'})
xc = pd.read_csv(OUT + 'tcga_rna.csv')  # cross-cancer TARGET-OS vs pooled TCGA, THE TCGA leg
m = scan.merge(rna[['Gene', 'OS_mean_log2TPM', 'OS_n', 'Comp_mean_log2TPM', 'Comp_n', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q']],
               left_on='Symbol', right_on='Gene', how='left', suffixes=('', '_rna'))
xc = xc.rename(columns={'Gene': 'Symbol'})
m = m.merge(xc[['Symbol', 'TCGA_RNA_n_OS', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_OS', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']],
            on='Symbol', how='left')
m = m.rename(columns={'p': 'Dep_p', 'q': 'Dep_q'})
# bars — TCGA leg is the cross-cancer TARGET-vs-TCGA comparison (design c);
# NO protein leg exists (no CPTAC OS study, no OS RPPA profile), so no
# protein bar and no +1.0 bonus for any gene.
m['CellRNA_pass'] = (m['TPM_FC'] >= 1.5) & (m['DeltaLog2'] >= 0.5) & (m['CellRNA_p'] < 0.05)
m['TCGA_XC_pass'] = (m['TCGA_RNA_log2FC'] >= 0.5) & (m['TCGA_RNA_p'] < 0.05)
m['OverExpr_pass'] = m['CellRNA_pass'] | m['TCGA_XC_pass']
m['Strict_pass'] = (m['Enrichment'] <= -0.10) & (m['OS_dep_frac'] >= 0.25) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Relaxed_pass'] = (m['Enrichment'] <= -0.05) & (m['OS_dep_frac'] >= 0.20) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
print('strict:', int(m.Strict_pass.sum()), 'relaxed:', int(m.Relaxed_pass.sum()))
print('strict & overexpr:', int((m.Strict_pass & m.OverExpr_pass).sum()))
print('relaxed & overexpr:', int((m.Relaxed_pass & m.OverExpr_pass).sum()))
print('strict & overexpr breakdown: cellRNA only', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & ~m.TCGA_XC_pass).sum()),
      '| XC only', int(((m.Strict_pass & m.OverExpr_pass) & m.TCGA_XC_pass & ~m.CellRNA_pass).sum()),
      '| both RNA', int(((m.Strict_pass & m.OverExpr_pass) & m.CellRNA_pass & m.TCGA_XC_pass).sum()))
print('relaxed genes with no XC data:', int(m[m.Relaxed_pass]['TCGA_RNA_log2FC'].isna().sum()))
print('relaxed genes with no cell-line RNA data:', int(m[m.Relaxed_pass]['DeltaLog2'].isna().sum()))
# combined score across the relaxed set (pre-defined; TCGA component =
# cross-cancer log2FC; NO protein bonus — no protein leg exists)
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel = m[m.Relaxed_pass].copy()
rel['Score'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['TCGA_RNA_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0))
rel['Tier'] = np.where(rel['Strict_pass'], 'A', 'B')
cand = rel[rel['OverExpr_pass']].copy().sort_values('Score', ascending=False)
print('candidates:', len(cand))
print(cand[['Symbol', 'Tier', 'OS_mean', 'OS_dep', 'OS_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'TPM_FC', 'DeltaLog2', 'TCGA_RNA_log2FC', 'Score']].head(30).to_string())
cand.to_csv(OUT + 'candidates.csv', index=False)
m.to_csv(OUT + 'overexpression.csv', index=False)
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist() if len(cand) > 20 else 'n/a (list shorter than 20)')
print('small-denominator in top20 (OS_n<17):')
print(cand.head(20)[cand.head(20).OS_n < 17][['Symbol', 'OS_n', 'OS_dep']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'CCNE1', 'TP53', 'RB1', 'MDM2', 'CDK4', 'CDK6', 'FGFR1', 'PDGFRA', 'IGF1R', 'MET', 'SRC', 'RUNX2', 'SATB2', 'SP7', 'DLX5', 'KLF5', 'HNF1B', 'FOSL1', 'PAX8', 'CXADR', 'PARD6B', 'ITGB5', 'ITGB1']:
    r = m[m.Symbol == g]
    cr = cand[cand.Symbol == g]
    print(g, 'enrich rank', int(r.iloc[0].EnrichRank),
          'strict', bool(r.iloc[0].Strict_pass), 'relaxed', bool(r.iloc[0].Relaxed_pass),
          'overexpr', bool(r.iloc[0].OverExpr_pass),
          'cellRNA FC', round(r.iloc[0].TPM_FC, 2) if pd.notna(r.iloc[0].TPM_FC) else None,
          'XC', round(r.iloc[0].TCGA_RNA_log2FC, 2) if pd.notna(r.iloc[0].TCGA_RNA_log2FC) else None,
          'candidate rank', (int(cand.index.get_loc(cr.index[0])) + 1 if len(cr) else None),
          'score', (round(float(cr.iloc[0].Score), 2) if len(cr) else None))
