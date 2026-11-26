import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/small-cell-lung/'
scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan['Symbol'] = scan['Gene'].str.split(' (', regex=False).str[0]
scan['Entrez'] = scan['Gene'].str.extract(r'\((\d+)\)$')[0].astype('Int64')
rna = pd.read_csv(OUT + 'celline_rna.csv')  # Gene=Symbol, TPM_FC, DeltaLog2, p, q
rna = rna.rename(columns={'p': 'CellRNA_p', 'q': 'CellRNA_q'})
xc = pd.read_csv(OUT + 'tcga_rna.csv')  # fetched cross-cancer values — RETAINED for transparency but NOT used as a bar (scale check failed: ucologne median raw 12.6 vs TCGA RSEM median 904.7; see report Methods)
m = scan.merge(rna[['Gene', 'SCLC_mean_log2TPM', 'SCLC_n', 'Comp_mean_log2TPM', 'Comp_n', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q']],
               left_on='Symbol', right_on='Gene', how='left', suffixes=('', '_rna'))
xc = xc.rename(columns={'Gene': 'Symbol'})
m = m.merge(xc[['Symbol', 'TCGA_RNA_n_SCLC', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_SCLC', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']],
            on='Symbol', how='left')
m = m.rename(columns={'p': 'Dep_p', 'q': 'Dep_q'})
# bars — the patient-RNA leg is INVALID for SCLC (cross-study scale trap, see
# report): ucologne values are on a different scale from TCGA RSEM, so no
# TCGA bar is applied and the TCGA score component is 0 for every gene
# (playbook §5: missing TCGA values are z-scored as 0). No protein leg
# exists for SCLC (no RPPA/CPTAC profile), so no protein bonus fires.
m['CellRNA_pass'] = (m['TPM_FC'] >= 1.5) & (m['DeltaLog2'] >= 0.5) & (m['CellRNA_p'] < 0.05)
m['TCGA_XC_pass'] = False
m['RPPA_pass'] = False
m['OverExpr_pass'] = m['CellRNA_pass']
m['Strict_pass'] = (m['Enrichment'] <= -0.10) & (m['SCLC_dep_frac'] >= 0.25) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
m['Relaxed_pass'] = (m['Enrichment'] <= -0.05) & (m['SCLC_dep_frac'] >= 0.20) & (m['FracDiff'] >= 0.10) & (~m['PanEssential']) & (m['GenomeWide_dep_frac'] <= 0.60)
print('strict:', int(m.Strict_pass.sum()), 'relaxed:', int(m.Relaxed_pass.sum()))
print('strict & overexpr:', int((m.Strict_pass & m.OverExpr_pass).sum()))
print('relaxed & overexpr:', int((m.Relaxed_pass & m.OverExpr_pass).sum()))
print('relaxed genes with no cell-line RNA data:', int(m[m.Relaxed_pass]['DeltaLog2'].isna().sum()))
print('cell-line RNA bar pass (genome-wide):', int(m.CellRNA_pass.sum()))
# combined score across the relaxed set: z(-enrichment) + z(frac diff) +
# z(cell-line DeltaLog2); TCGA component = 0 for every gene (leg unavailable)
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel = m[m.Relaxed_pass].copy()
rel['Score'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['DeltaLog2'].fillna(0))
rel['Tier'] = np.where(rel['Strict_pass'], 'A', 'B')
cand = rel[rel['OverExpr_pass']].copy().sort_values('Score', ascending=False)
print('candidates:', len(cand))
print(cand[['Symbol', 'Tier', 'SCLC_mean', 'SCLC_dep', 'SCLC_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'TPM_FC', 'DeltaLog2', 'TCGA_RNA_log2FC', 'Score']].head(30).to_string())
cand.to_csv(OUT + 'candidates.csv', index=False)
m.to_csv(OUT + 'overexpression.csv', index=False)
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist() if len(cand) > 20 else 'n/a (list shorter than 20)')
print('small-denominator in top20 (SCLC_n<26):')
print(cand.head(20)[cand.head(20).SCLC_n < 26][['Symbol', 'SCLC_n', 'SCLC_dep']].to_string())
print('low-n (2-line) genes among candidates:')
print(cand[cand.SCLC_n <= 2][['Symbol', 'SCLC_n', 'SCLC_dep', 'Score']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'ITGB5', 'ITGB1', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'MYCL', 'MYCN', 'ASCL1', 'NEUROD1', 'POU2F3', 'INSM1', 'DLL3', 'BCL2', 'EZH2', 'KDM1A', 'RB1', 'TP53', 'SOX10', 'KLF5', 'FOSL1', 'NEUROG2', 'CHGA', 'SYP', 'FOXM1', 'CCNE1', 'CCNE2', 'E2F3', 'SKP2', 'SLC2A1', 'NAMPT', 'SMARCA4', 'EP300', 'DOT1L', 'TIMM50', 'RNASEK']:
    r = m[m.Symbol == g]
    cr = cand[cand.Symbol == g]
    print(g, 'enrich rank', int((m.Enrichment < r.iloc[0].Enrichment).sum()) + 1,
          'strict', bool(r.iloc[0].Strict_pass), 'relaxed', bool(r.iloc[0].Relaxed_pass),
          'overexpr', bool(r.iloc[0].OverExpr_pass),
          'cellRNA FC', round(r.iloc[0].TPM_FC, 2) if pd.notna(r.iloc[0].TPM_FC) else None,
          'XC(log2FC, NOT a bar)', round(r.iloc[0].TCGA_RNA_log2FC, 2) if pd.notna(r.iloc[0].TCGA_RNA_log2FC) else None,
          'candidate rank', (int(cand.index.get_loc(cr.index[0])) + 1 if len(cr) else None),
          'score', (round(float(cr.iloc[0].Score), 2) if len(cr) else None))
