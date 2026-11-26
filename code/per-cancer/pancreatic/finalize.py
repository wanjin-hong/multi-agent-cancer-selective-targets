import pandas as pd, numpy as np
m = pd.read_csv('overexpression.csv')
cand = pd.read_csv('candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol','Gene','PDAC_mean','PDAC_dep','PDAC_n','PDAC_dep_frac','Comp_mean','Comp_dep','Comp_n','Comp_dep_frac',
          'Enrichment','FracDiff','GenomeWide_dep_frac','Dep_p','Dep_q','PanEssential']].copy()
scan = scan.rename(columns={'Gene':'Column','Dep_p':'p_welch','Dep_q':'q_BH','Symbol':'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv('dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'KLF5': 'Transcription factor; no direct inhibitor; preclinical KLF5 degraders reported only',
 'FOSL1': 'AP-1 transcription factor (FRA1); no direct inhibitor',
 'SPRR1B': 'Small proline-rich epithelial differentiation protein; no drug',
 'DOCK5': 'Rac GEF (DOCK family); no drug; preclinical tool inhibitors only',
 'NCKAP1': 'WAVE complex subunit; no drug',
 'SLC2A1': 'GLUT1 glucose transporter; no approved inhibitor; GLUT1 inhibitors (e.g. BAY-876) preclinical only',
 'SOX9': 'Transcription factor; no direct inhibitor',
 'PARD6B': 'Cell polarity protein; no drug',
 'JUNB': 'AP-1 transcription factor; no direct inhibitor',
 'FZD5': 'Frizzled-5 Wnt receptor; no FZD5-selective drug; broad Wnt/porcupine inhibitors (ETC-159, RXC004) investigational',
 'BCAR1': 'p130Cas adaptor/scaffold; no drug',
 'CFLAR': 'c-FLIP, apoptosis regulator; no direct inhibitor; preclinical c-FLIP inhibitors only',
 'EGFR': 'Approved EGFR TKIs (erlotinib, gefitinib, osimertinib) and antibodies (cetuximab, panitumumab); erlotinib + gemcitabine approved in PDAC with marginal benefit',
 'HNF1B': 'Transcription factor; no direct inhibitor; pancreatic lineage marker',
 'MT1E': 'Metallothionein; undruggable as a target; culture-sensitive',
 'DLG5': 'Scaffold/polarity protein; no drug',
 'TCF7L2': 'Wnt transcription factor (TCF4); no direct inhibitor',
 'TIPARP': 'PARP7 mono-ADP-ribosyltransferase; RBN-2397 (PARP7 inhibitor) in Phase 1',
 'PXN': 'Paxillin, focal-adhesion adaptor; no drug',
 'CXADR': 'Coxsackievirus/adenovirus receptor, adhesion; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank','Symbol','Tier','PDAC_mean','PDAC_dep','PDAC_n','Comp_mean','Comp_dep','Comp_n','Comp_dep_frac','Enrichment','FracDiff',
        'GenomeWide_dep_frac','TPM_FC','DeltaLog2','CellRNA_p','CellRNA_q','GDC_TvN_log2FC','GDC_TvN_p','GDC_TvN_q',
        'TCGA_RNA_log2FC','TCGA_RNA_p','TCGA_RNA_q','RPPA_delta','RPPA_p','CPTAC_n_tumour','CPTAC_median','Score','Druggability']
top20[cols].to_csv('top20.csv', index=False)
print(top20[['Rank','Symbol','Tier','Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol','Gene','Entrez','PDAC_mean','PDAC_dep','PDAC_n','PDAC_dep_frac','Comp_mean','Comp_dep','Comp_n','Comp_dep_frac',
 'Enrichment','FracDiff','GenomeWide_dep_frac','PanEssential','Dep_p','Dep_q',
 'PDAC_mean_log2TPM','Comp_mean_log2TPM','DeltaLog2','TPM_FC','CellRNA_p','CellRNA_q',
 'GDC_TvN_n_tumour','GDC_TvN_n_normal','GDC_TvN_mean_tumour','GDC_TvN_mean_normal','GDC_TvN_log2FC','GDC_TvN_p','GDC_TvN_q',
 'TCGA_RNA_n_PAAD','TCGA_RNA_n_other','TCGA_RNA_mean_PAAD','TCGA_RNA_mean_other','TCGA_RNA_log2FC','TCGA_RNA_p','TCGA_RNA_q',
 'RPPA_n_PAAD','RPPA_n_other','RPPA_delta','RPPA_p',
 'CPTAC_n_tumour','CPTAC_median','CPTAC_mean',
 'CellRNA_pass','TCGA_TvN_pass','TCGA_XC_pass','RPPA_pass','OverExpr_pass','Strict_pass','Relaxed_pass']
m[keep].to_csv('overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
