import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/head-and-neck/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'HNSC_mean', 'HNSC_dep', 'HNSC_n', 'HNSC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'TP63': 'Transcription factor (p63, squamous lineage master regulator); no direct inhibitor',
 'EGFR': 'Approved EGFR TKIs (erlotinib, gefitinib, osimertinib) and antibodies; cetuximab is FDA-approved in HNSC (with radiation / with platinum chemotherapy)',
 'SNAI2': 'SLUG, EMT transcription factor; no direct inhibitor',
 'ITGB1': 'Integrin beta-1; no beta-1-selective drug; pan-integrin agents have failed clinically',
 'KLF5': 'Transcription factor; no direct inhibitor; preclinical KLF5 degraders reported only',
 'ITGA3': 'Integrin alpha-3; no selective drug',
 'ADAR': 'ADAR1 RNA-editing enzyme; no clinical inhibitor; preclinical ADAR1 inhibitors in development',
 'EFR3A': 'Phosphoinositide / membrane trafficking regulator; no drug',
 'EXT1': 'Heparan sulfate biosynthesis enzyme; no drug',
 'CDK6': 'CDK4/6 inhibitors approved (palbociclib, ribociclib, abemaciclib) in HR+/HER2- breast cancer; not HNSC-labelled',
 'FOSL1': 'AP-1 transcription factor (FRA1); no direct inhibitor',
 'TRAF2': 'TNF-receptor adaptor; no drug',
 'CKAP2L': 'Mitotic spindle protein; no drug',
 'KIF22': 'Kinesin motor protein; no drug',
 'ENO1': 'Glycolytic enzyme (alpha-enolase); no clinical inhibitor; preclinical ENO1 inhibitors only',
 'CCDC51': 'Mitochondrial inner-membrane protein; no drug',
 'ERBB2': 'HER2; approved antibodies and ADCs (trastuzumab, trastuzumab deruxtecan) in HER2+ cancers; HER2 amplification is rare in HNSC (~3-5%)',
 'TIPARP': 'PARP7 mono-ADP-ribosyltransferase; RBN-2397 (PARP7 inhibitor) in Phase 1',
 'RRAGC': 'RagC GTPase, mTORC1 amino-acid sensing; no drug',
 'CAP1': 'Adenylyl cyclase-associated protein, actin dynamics; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'HNSC_mean', 'HNSC_dep', 'HNSC_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q', 'RPPA_delta', 'RPPA_p', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'HNSC_mean', 'HNSC_dep', 'HNSC_n', 'HNSC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'HNSC_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
 'TCGA_RNA_n_HNSC', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_HNSC', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_HNSC', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p',
 'CellRNA_pass', 'TCGA_TvN_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
