import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/renal-cell/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'RCC_mean', 'RCC_dep', 'RCC_n', 'RCC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'HNF1B': 'Transcription factor (hepatocyte nuclear factor 1-beta, kidney lineage); no direct inhibitor',
 'PAX8': 'Transcription factor (paired-box 8, renal lineage master regulator); no direct inhibitor',
 'CFLAR': 'c-FLIP, anti-apoptotic regulator of caspase-8; no clinical inhibitor; preclinical FLIP inhibitors only',
 'CDH2': 'N-cadherin; no selective drug; ADH-1 (exherin) was discontinued after early trials',
 'MDM2': 'MDM2-p53 interaction inhibitors in Phase 1/2 (idasanutlin, siremadlin, milademetan); no MDM2 inhibitor approved; activity requires wild-type TP53',
 'ITGAV': 'Integrin alpha-v; cilengitide failed Phase 3 in glioblastoma; no approved alpha-v drug',
 'WWTR1': 'TAZ; no direct inhibitor — drugged indirectly via pan-TEAD inhibitors in clinic (VT3989, ISM6331, ODM-212)',
 'FERMT2': 'Kindlin-2, integrin adhesion scaffold; no drug',
 'PARD6B': 'Cell-polarity protein (Par6B); no drug',
 'UBC': 'Ubiquitin C; no drug — and scored in only 3 lines (low-n flag), so the rank is insecure',
 'PARD3': 'Cell-polarity protein (Par3); no drug',
 'ZFP36L1': 'RNA-binding protein (tristetraprolin-like 1); no drug',
 'MZT2A': 'Mitotic spindle-organising protein (MOZART2); no drug',
 'BRAT1': 'BRCA1-associated ATM activator; no drug',
 'KIF2C': 'Kinesin-13 (MCAK), mitotic microtubule depolymerase; no clinical inhibitor',
 'EPAS1': 'HIF-2alpha; belzutifan (Welireg) FDA-approved — Aug 2021 for VHL-disease-associated RCC, Dec 2023 for advanced RCC after PD-1/PD-L1 inhibitor and VEGF-TKI; the only gene in the list with a drug approved in this cancer',
 'TUBA3D': 'Tubulin alpha-3D; no selective drug (microtubule agents are non-selective)',
 'FAF2': 'UBXD8, ER-associated degradation adaptor; no drug',
 'UBA5': 'E1 enzyme of the UFMylation pathway; no clinical UBA5 inhibitor (TAK-243 inhibits UBA1, not UBA5)',
 'UFM1': 'Ubiquitin-fold modifier 1; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'RCC_mean', 'RCC_dep', 'RCC_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q', 'RPPA_delta', 'RPPA_p', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'RCC_mean', 'RCC_dep', 'RCC_n', 'RCC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'RCC_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
 'TCGA_RNA_n_KIRC', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_KIRC', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_KIRC', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p',
 'CellRNA_pass', 'TCGA_TvN_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
