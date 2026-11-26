import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/bladder/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'Bladder_mean', 'Bladder_dep', 'Bladder_n', 'Bladder_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'PPARG': 'PPAR-gamma, luminal lineage transcription factor; approved drugs are AGONISTS (thiazolidinediones: pioglitazone, rosiglitazone) — the wrong direction for a dependency; no clinical PPARG antagonist',
 'GRHL2': 'Grainyhead-like 2 transcription factor (epithelial identity); no drug',
 'TP63': 'Transcription factor (p63, basal/squamous lineage master regulator); no direct inhibitor',
 'KLF5': 'Transcription factor; no direct inhibitor; preclinical KLF5 degraders reported only',
 'SOX4': 'Transcription factor (SOX family); no direct inhibitor',
 'SKP2': 'S-phase kinase-associated protein 2, SCF E3 ubiquitin ligase; no clinical inhibitor (preclinical SKP2 inhibitors only)',
 'FAM50A': 'XAP5-like protein, spliceosome-associated; no drug',
 'TRIP13': 'AAA ATPase, mitotic spindle checkpoint; no clinical inhibitor (preclinical TRIP13 inhibitors reported)',
 'RXRA': 'Retinoid X receptor alpha; bexarotene is an RXR AGONIST (approved in CTCL) — the wrong direction for a dependency; no clinical RXRA antagonist',
 'E2F3': 'Transcription factor; no direct inhibitor',
 'KATNB1': 'Katanin regulatory subunit, microtubule severing; no drug',
 'SLC35B2': 'Nucleotide-sugar transporter (PAPS transporter); no drug',
 'PTK2': 'FAK; defactinib + avutometinib approved in combination (May 2025, KRAS-mutant low-grade serous ovarian cancer); FAK inhibitors in trials; not bladder-labelled',
 'WDR83': 'WD-repeat scaffold protein; no drug',
 'EGFR': 'Approved EGFR TKIs and antibodies in lung/HNSC/colon; no established bladder indication; EGFR inhibitors are not bladder-labelled',
 'FOSL1': 'AP-1 transcription factor (FRA1); no direct inhibitor',
 'TEAD3': 'No TEAD3-selective drug; pan-TEAD inhibitors in clinic (VT3989, ISM6331, ODM-212)',
 'TRAF2': 'TNF-receptor adaptor; no drug',
 'TUBB4B': 'Beta-tubulin; microtubule-targeting cytotoxics exist but are not gene-selective',
 'ENO1': 'Glycolytic enzyme (alpha-enolase); no clinical inhibitor; preclinical ENO1 inhibitors only',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'Bladder_mean', 'Bladder_dep', 'Bladder_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q', 'RPPA_delta', 'RPPA_p', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'Bladder_mean', 'Bladder_dep', 'Bladder_n', 'Bladder_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'Bladder_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
 'TCGA_RNA_n_BLCA', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_BLCA', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_BLCA', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p',
 'CellRNA_pass', 'TCGA_TvN_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
