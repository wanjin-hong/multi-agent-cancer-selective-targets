import pandas as pd, numpy as np, re
scan = pd.read_csv('dependency_scan.csv')
scan['Symbol'] = scan['Gene'].str.split(' (', regex=False).str[0]
scan['Entrez'] = scan['Gene'].str.extract(r'\((\d+)\)$')[0].astype('Int64')
rna = pd.read_csv('celline_rna.csv')  # Gene=Symbol, TPM_FC, DeltaLog2, p, q
rna = rna.rename(columns={'p':'CellRNA_p','q':'CellRNA_q'})
tvn = pd.read_csv('tcga_tumour_vs_normal.csv')
xc = pd.read_csv('tcga_rna.csv')  # cross-cancer
m = scan.merge(rna[['Gene','HCC_mean_log2TPM','HCC_n','Comp_mean_log2TPM','Comp_n','DeltaLog2','TPM_FC','CellRNA_p','CellRNA_q']],
               left_on='Symbol', right_on='Gene', how='left', suffixes=('','_rna'))
tvn = tvn.rename(columns={'Gene':'Symbol'})
m = m.merge(tvn[['Symbol','GDC_TvN_n_tumour','GDC_TvN_n_normal','GDC_TvN_mean_tumour','GDC_TvN_mean_normal','GDC_TvN_log2FC','GDC_TvN_p','GDC_TvN_q']],
            on='Symbol', how='left')
xc = xc.rename(columns={'Gene':'Symbol'})
m = m.merge(xc[['Symbol','TCGA_RNA_n_LIHC','TCGA_RNA_n_other','TCGA_RNA_mean_LIHC','TCGA_RNA_mean_other','TCGA_RNA_log2FC','TCGA_RNA_p','TCGA_RNA_q','TCGA_refNormal_z_mean','RPPA_n_LIHC','RPPA_n_other','RPPA_delta','RPPA_p']],
            on='Symbol', how='left')
m = m.rename(columns={'p':'Dep_p','q':'Dep_q'})
# bars
m['CellRNA_pass'] = (m['TPM_FC']>=1.5)&(m['DeltaLog2']>=0.5)&(m['CellRNA_p']<0.05)
m['TCGA_TvN_pass'] = (m['GDC_TvN_log2FC']>=0.5)&(m['GDC_TvN_p']<0.05)
m['TCGA_XC_pass'] = (m['TCGA_RNA_log2FC']>=0.5)&(m['TCGA_RNA_p']<0.05)
m['RPPA_pass'] = (m['RPPA_delta']>=0.5)&(m['RPPA_p']<0.05)
m['OverExpr_pass'] = m['CellRNA_pass']|m['TCGA_TvN_pass']|m['RPPA_pass']
m['Strict_pass'] = (m['Enrichment']<=-0.10)&(m['HCC_dep_frac']>=0.25)&(m['FracDiff']>=0.10)&(~m['PanEssential'])&(m['GenomeWide_dep_frac']<=0.60)
m['Relaxed_pass'] = (m['Enrichment']<=-0.05)&(m['HCC_dep_frac']>=0.20)&(m['FracDiff']>=0.10)&(~m['PanEssential'])&(m['GenomeWide_dep_frac']<=0.60)
print('strict:', m.Strict_pass.sum(), 'relaxed:', m.Relaxed_pass.sum())
print('strict & overexpr:', (m.Strict_pass & m.OverExpr_pass).sum())
print('relaxed & overexpr:', (m.Relaxed_pass & m.OverExpr_pass).sum())
print('strict & overexpr & no TCGA data:', ((m.Strict_pass & m.OverExpr_pass) & m['GDC_TvN_log2FC'].isna()).sum())
# combined score across relaxed genes passing overexpr bar
sub = m[m.Relaxed_pass & m.OverExpr_pass].copy()
def z(s):
    s = s.astype(float)
    return (s - s.mean())/s.std(ddof=1)
sub['Score'] = z(-sub['Enrichment']) + z(sub['FracDiff']) + z(sub['GDC_TvN_log2FC'].fillna(0)) + z(sub['DeltaLog2'].fillna(0)) + np.where(sub['RPPA_pass'].fillna(False), 1.0, 0.0)
# NOTE: TNBC template z-scored TCGA/cell-line with missing=0 across the relaxed set; here score is across the relaxed&overexpr set? 
# TNBC: score computed for each of the 152 relaxed genes, then top 20 taken from the 25 passing both. Follow exactly:
rel = m[m.Relaxed_pass].copy()
rel['Score'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['GDC_TvN_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0)) + np.where(rel['RPPA_pass'].fillna(False), 1.0, 0.0)
rel['Tier'] = np.where(rel['Strict_pass'], 'A', 'B')
cand = rel[rel['OverExpr_pass']].copy().sort_values('Score', ascending=False)
print('candidates:', len(cand))
print(cand[['Symbol','Tier','HCC_mean','HCC_dep','HCC_n','Comp_dep_frac','Enrichment','FracDiff','TPM_FC','DeltaLog2','GDC_TvN_log2FC','TCGA_RNA_log2FC','Score']].head(30).to_string())
cand.to_csv('candidates.csv', index=False)
m.to_csv('overexpression.csv', index=False)
# top-20 table columns
top20 = cand.head(20)
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
print('small-denominator in top20 (HCC_n<21):', top20[top20.HCC_n<21][['Symbol','HCC_n','HCC_dep']].to_string())
# headline genes full rows
for g in ['TEAD1','YAP1','WWTR1','PTK2','PKN2','GRB2','ITGAV','FERMT2','RAC1']:
    r = m[m.Symbol==g]
    print(g, 'enrich rank', int((m.Enrichment < r.iloc[0].Enrichment).sum())+1,
          'strict', bool(r.iloc[0].Strict_pass), 'relaxed', bool(r.iloc[0].Relaxed_pass),
          'overexpr', bool(r.iloc[0].OverExpr_pass),
          'cellRNA FC', round(r.iloc[0].TPM_FC,2) if pd.notna(r.iloc[0].TPM_FC) else None,
          'TvN', round(r.iloc[0].GDC_TvN_log2FC,2) if pd.notna(r.iloc[0].GDC_TvN_log2FC) else None,
          'XC', round(r.iloc[0].TCGA_RNA_log2FC,2) if pd.notna(r.iloc[0].TCGA_RNA_log2FC) else None,
          'z', round(r.iloc[0].TCGA_refNormal_z_mean,2) if pd.notna(r.iloc[0].TCGA_refNormal_z_mean) else None)
