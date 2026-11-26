import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/ovarian/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'OV_mean', 'OV_dep', 'OV_n', 'OV_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'PAX8': 'Paired-box transcription factor (Mullerian lineage master regulator); no direct inhibitor',
 'ARHGAP29': 'Rho GTPase-activating protein (RhoGAP); no drug',
 'FERMT2': 'Kindlin-2, integrin adhesion scaffold; no drug',
 'CCNE1': 'Cyclin E1; no direct inhibitor — CCNE1-amplified tumours are sensitive to CDK2 inhibitors (PF-07104091, BLU-222, Phase 1) and WEE1/PKMYT1 inhibition, indirect routes only',
 'PARD6B': 'Cell-polarity protein (Par6); no drug',
 'XPR1': 'Phosphate exporter; no drug',
 'PARD3': 'Cell-polarity protein (Par3); no drug',
 'KLF5': 'Transcription factor; no direct inhibitor; preclinical KLF5 degraders reported only',
 'WWTR1': 'TAZ, Hippo co-activator; drugged only indirectly, via pan-TEAD inhibitors in clinic (VT3989, ISM6331, ODM-212)',
 'C10orf95': 'Uncharacterised open reading frame; no drug — scored in only 28 of 58 lines (low-n flag)',
 'SOCS3': 'Suppressor of cytokine signalling 3 (JAK-STAT negative regulator); no drug',
 'ZFP36L1': 'RNA-binding zinc-finger protein; no drug',
 'HS6ST1': 'Heparan sulfate 6-O-sulfotransferase; no drug',
 'FAAP20': 'Fanconi anaemia-associated protein 20; no drug',
 'H2AX': 'Histone H2AX (DNA-damage marker); no drug',
 'BTF3L4': 'Basic transcription factor 3-like 4; no drug',
 'RELA': 'NF-kappaB p65 subunit; no direct inhibitor (IKK inhibitors act upstream)',
 'RBM10': 'RNA-binding motif protein 10; no drug',
 'ITGAV': 'Integrin alpha-v; cilengitide (pan-alpha-v) already failed Phase 3 in GBM',
 'FAAP24': 'Fanconi anaemia-associated protein 24; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'OV_mean', 'OV_dep', 'OV_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
        'RPPA_delta', 'RPPA_p', 'CPTAC_n_tumour', 'CPTAC_median', 'CPTAC_mean', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-28:', cand.iloc[20:28]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'OV_mean', 'OV_dep', 'OV_n', 'OV_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'OV_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'TCGA_RNA_n_OV', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_OV', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_OV', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p',
 'CPTAC_n_tumour', 'CPTAC_median', 'CPTAC_mean',
 'CellRNA_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
