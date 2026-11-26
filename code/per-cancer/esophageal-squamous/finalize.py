import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophageal-squamous/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'ESCC_mean', 'ESCC_dep', 'ESCC_n', 'ESCC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'TP63': 'Transcription factor (p63, squamous lineage master regulator); no direct inhibitor',
 'EGFR': 'Approved EGFR TKIs (erlotinib, gefitinib, osimertinib) and antibodies (cetuximab) in lung / HNSC / colorectal cancer — none ESCC-labelled; EGFR TKIs have been tested in ESCC trials without an approved indication',
 'SOX2': 'Transcription factor (squamous stemness); no direct inhibitor',
 'ADAR': 'ADAR1 RNA-editing enzyme; no clinical inhibitor; preclinical ADAR1 inhibitors in development',
 'KLF5': 'Transcription factor; no direct inhibitor; preclinical KLF5 degraders reported only',
 'TIMM17A': 'Mitochondrial inner-membrane translocase subunit; no drug',
 'CTTN': 'Cortactin, actin-binding scaffold (11q13 amplicon with CCND1); no direct inhibitor',
 'MYBL2': 'Transcription factor (B-Myb), cell-cycle driver; no direct inhibitor',
 'STX4': 'Syntaxin-4, SNARE membrane-fusion protein; no drug',
 'CCDC51': 'Mitochondrial inner-membrane protein; no drug',
 'RNF168': 'E3 ubiquitin ligase (DNA-damage response); no drug',
 'SLC2A1': 'GLUT1 glucose transporter; no clinical inhibitor; preclinical GLUT1 inhibitors only (e.g. BAY-876)',
 'FADD': 'Fas-associated death-domain adaptor (11q13 amplicon); no drug',
 'SREBF1': 'SREBP1 transcription factor, lipogenesis; no clinical inhibitor; fatostatin is preclinical only',
 'SLC39A10': 'ZIP10 zinc transporter; no drug',
 'PTDSS1': 'Phosphatidylserine synthase 1; no drug',
 'UBA6': 'Ubiquitin-activating enzyme E1 (UBE1L2); no drug',
 'MRPL18': 'Mitochondrial ribosomal protein L18; no drug',
 'ERBB3': 'HER3; no approved direct inhibitor; HER3-directed ADC patritumab deruxtecan in Phase 2/3 in other cancers',
 'ERBB2': 'HER2; approved antibodies and ADCs (trastuzumab, trastuzumab deruxtecan) in HER2+ cancers; HER2 amplification is uncommon in ESCC',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'ESCC_mean', 'ESCC_dep', 'ESCC_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q', 'RPPA_delta', 'RPPA_p', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'ESCC_mean', 'ESCC_dep', 'ESCC_n', 'ESCC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'ESCC_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
 'TCGA_RNA_n_ESCC', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_ESCC', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_ESCC', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p',
 'CellRNA_pass', 'TCGA_TvN_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
