import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/small-cell-lung/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'SCLC_mean', 'SCLC_dep', 'SCLC_n', 'SCLC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'ASCL1': 'Lineage master transcription factor (SCLC-A); no direct inhibitor',
 'INSM1': 'Neuroendocrine transcription factor (zinc finger); no drug',
 'SKP2': 'SCF ubiquitin-ligase subunit (targets p27 for degradation); preclinical inhibitors only (e.g. SZL-P1-41), no clinical drug',
 'E2F3': 'E2F transcription factor; no direct inhibitor',
 'PRKRA': 'PACT, dsRNA-binding PKR activator; no drug',
 'FOXA2': 'Forkhead transcription factor (endoderm/neuroendocrine lineage); no direct inhibitor',
 'CCNE2': 'Cyclin E2; no direct inhibitor — CCNE-amplified tumours are sensitive to CDK2 inhibitors (PF-07104091, BLU-222, Phase 1) and WEE1/PKMYT1 inhibition, indirect routes only (same route as ovarian CCNE1)',
 'CXADR': 'Coxsackie/adenovirus receptor, adhesion; no drug — scored in only 2 of 26 lines (DMS 273, NCI-H1882), low-n flag',
 'MIB1': 'Mindbomb E3 ubiquitin ligase (Notch ligand ubiquitination); no drug',
 'ANKRD36': 'Uncharacterised ankyrin-repeat protein; no drug — scored in only 2 of 26 lines (DMS 273, NCI-H1882), low-n flag',
 'TSPYL5': 'Testis-specific protein Y-encoded-like 5 (nucleosome assembly); no drug',
 'EPC2': 'Enhancer of Polycomb homolog 2 (Polycomb repressive machinery); no drug',
 'CIT': 'Citron Rho-interacting kinase (cytokinesis); preclinical inhibitors only, no clinical drug',
 'TFDP1': 'E2F dimerisation partner DP1; no direct inhibitor',
 'YEATS2': 'Histone acetylation reader (ATAC complex); no drug',
 'ASPM': 'Mitotic spindle protein; no drug',
 'PNISR': 'Splicing regulator (PNN-interacting serine/arginine-rich protein); no drug',
 'SQLE': 'Squalene epoxidase (cholesterol biosynthesis); no oncology drug — terbinafine is an approved antifungal SQLE inhibitor, not a cancer therapy',
 'KIF14': 'Kinesin (mitotic); no drug',
 'WSB2': 'SOCS-box E3 ligase adaptor; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'SCLC_mean', 'SCLC_dep', 'SCLC_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'TCGA_RNA_log2FC', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order (keep fetched TCGA columns for transparency; TCGA_XC_pass is False for all — leg invalid, see report)
keep = ['Symbol', 'Gene', 'Entrez', 'SCLC_mean', 'SCLC_dep', 'SCLC_n', 'SCLC_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'SCLC_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'TCGA_RNA_n_SCLC', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_SCLC', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'CellRNA_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
