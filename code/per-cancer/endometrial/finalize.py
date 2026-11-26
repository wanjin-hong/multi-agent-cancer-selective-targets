import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/endometrial/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'ENDO_mean', 'ENDO_dep', 'ENDO_n', 'ENDO_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'PAX8': 'Paired-box transcription factor (Mullerian lineage master regulator); no direct inhibitor',
 'HNF1B': 'Transcription factor (hepatocyte nuclear factor 1-beta, Mullerian/renal lineage); no direct inhibitor',
 'CCNE1': 'Cyclin E1; no direct inhibitor — CCNE1-amplified tumours are sensitive to CDK2 inhibitors (PF-07104091, BLU-222, Phase 1) and WEE1/PKMYT1 inhibition, indirect routes only',
 'HS6ST1': 'Heparan sulfate 6-O-sulfotransferase; no drug',
 'MECOM': 'Transcription factor (EVI1, MDS1-EVI1); no direct inhibitor',
 'ARHGAP29': 'RhoGAP; no drug',
 'MT1X': 'Metallothionein 1X; undruggable, culture-sensitive (same caveat as MT1E)',
 'PDZK1IP1': 'MAP17, small membrane protein; no drug',
 'PARD6B': 'Cell-polarity protein (Par6B); no drug',
 'TACC3': 'Transforming acidic coiled-coil 3, centrosome/mitotic; no direct inhibitor (FGFR inhibitors act on the FGFR3-TACC3 fusion partner in fusion cancers, not applicable here)',
 'TMSB10': 'Thymosin beta-10, actin-binding protein; no drug',
 'CXADR': 'Coxsackievirus/adenovirus receptor, adhesion; no drug — scored in only 7 of 28 lines, low-n flag',
 'BRCA1': 'DNA-repair protein; no direct inhibitor — PARP inhibitors exploit BRCA loss, not BRCA dependency, so this entry does not imply a PARP route',
 'XPR1': 'Phosphate exporter; no drug',
 'RPL22L1': 'Ribosomal protein L22-like 1; no drug',
 'NFKBIE': 'IkappaB epsilon (NF-kB inhibitor); no drug',
 'NME2': 'Nucleoside diphosphate kinase 2; no drug',
 'EPOP': 'PRC2-associated protein (Elongin BC and Polycomb repressive complex 2-associated protein); no drug',
 'BRCA2': 'DNA-repair protein; no direct inhibitor — PARP inhibitors exploit BRCA loss, not BRCA dependency, so this entry does not imply a PARP route',
 'UBE2S': 'Ubiquitin-conjugating enzyme E2 S; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'ENDO_mean', 'ENDO_dep', 'ENDO_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q', 'RPPA_delta', 'RPPA_p', 'CPTAC_median', 'CPTAC_n_tumour', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'ENDO_mean', 'ENDO_dep', 'ENDO_n', 'ENDO_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'ENDO_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
 'TCGA_RNA_n_UCEC', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_UCEC', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_UCEC', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p', 'CPTAC_n_tumour', 'CPTAC_median', 'CPTAC_mean',
 'CellRNA_pass', 'TCGA_TvN_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
