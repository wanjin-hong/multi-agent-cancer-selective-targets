import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol', 'Gene', 'ESGA_mean', 'ESGA_dep', 'ESGA_n', 'ESGA_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential',
          'Gastric_dep', 'Gastric_n', 'EsoGEJ_dep', 'EsoGEJ_n']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'KLF5': 'Transcription factor (Kruppel-like factor 5, intestinal/gastric lineage); no direct inhibitor; preclinical KLF5 degraders reported only',
 'TCF7L2': 'Wnt transcription factor (TCF4); no direct inhibitor',
 'RAB10': 'Small GTPase, vesicle trafficking; no drug',
 'CTNNB1': 'Beta-catenin; no approved direct inhibitor; tegavivint (beta-catenin/TBL1 inhibitor) in Phase 1/2 and PRI-724 (CBP/beta-catenin) in Phase 1/2 — clinical-stage Wnt/beta-catenin route, not approved',
 'MYBL2': 'Transcription factor (B-Myb); no direct inhibitor',
 'GART': 'Purine biosynthesis trifunctional enzyme; no GART-selective drug — antifolates (pemetrexed, methotrexate) inhibit folate-dependent purine synthesis at the pathway level',
 'YAP1': 'Hippo co-activator; no direct inhibitor — drugged indirectly via pan-TEAD inhibitors in clinic (VT3989, ISM6331, ODM-212)',
 'PAICS': 'Purine biosynthesis enzyme; no clinical inhibitor',
 'LRPPRC': 'Mitochondrial RNA-binding protein; no drug',
 'GFPT1': 'Hexosamine pathway rate-limiting enzyme (glutamine-fructose-6-phosphate transaminase 1); no clinical inhibitor',
 'HNF1B': 'Transcription factor (hepatocyte nuclear factor 1-beta, pancreatic/gastric lineage); no direct inhibitor',
 'CAD': 'Pyrimidine biosynthesis trifunctional enzyme; no direct inhibitor — DHODH inhibitors (brequinar, BAY2402234) act downstream in the same pathway',
 'MRPS35': 'Mitochondrial ribosomal protein S35; no drug',
 'PFAS': 'Purine biosynthesis enzyme (phosphoribosylformylglycinamidine synthase); no drug',
 'STXBP3': 'Vesicle trafficking / syntaxin-binding protein; no drug',
 'ATP2B1': 'Plasma membrane calcium ATPase (PMCA1); no drug',
 'ADSS2': 'Purine biosynthesis enzyme (adenylosuccinate synthetase 2); no drug',
 'UMPS': 'Pyrimidine biosynthesis enzyme; converts 5-FU to active FUMP — a 5-FU sensitivity determinant, not a drug target',
 'AHCY': 'S-adenosylhomocysteine hydrolase; no clinical inhibitor',
 'GFM1': 'Mitochondrial translation elongation factor G1; no drug',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'ESGA_mean', 'ESGA_dep', 'ESGA_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q', 'RPPA_delta', 'RPPA_p', 'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order (rename the RNA-merge n columns to avoid collision with dependency n)
m = m.rename(columns={'ESGA_n_rna': 'ESGA_RNA_n', 'Comp_n_rna': 'Comp_RNA_n',
                      'Gastric_n_rna': 'Gastric_RNA_n', 'EsoGEJ_n_rna': 'EsoGEJ_RNA_n'})
keep = ['Symbol', 'Gene', 'Entrez', 'ESGA_mean', 'ESGA_dep', 'ESGA_n', 'ESGA_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'Gastric_dep', 'Gastric_n', 'EsoGEJ_dep', 'EsoGEJ_n',
 'ESGA_mean_log2TPM', 'ESGA_RNA_n', 'Comp_mean_log2TPM', 'Comp_RNA_n', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'Gastric_mean_log2TPM', 'Gastric_RNA_n', 'EsoGEJ_mean_log2TPM', 'EsoGEJ_RNA_n',
 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_mean_tumour', 'GDC_TvN_mean_normal', 'GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q',
 'GDC_TvN_STAD_log2FC', 'GDC_TvN_STAD_p', 'GDC_TvN_ESCA_log2FC', 'GDC_TvN_ESCA_p',
 'TCGA_RNA_n_ESGA', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_ESGA', 'TCGA_RNA_mean_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
 'RPPA_n_ESGA', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p',
 'CellRNA_pass', 'TCGA_TvN_pass', 'TCGA_XC_pass', 'RPPA_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
cols_present = [c for c in keep if c in m.columns]
m[cols_present].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
print('overexpression columns kept:', len(cols_present), 'of', len(keep))
