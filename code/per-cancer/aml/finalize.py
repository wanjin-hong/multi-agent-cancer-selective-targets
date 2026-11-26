import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/aml/'
m = pd.read_csv(OUT + 'overexpression.csv')
m['Entrez'] = m['Column'].str.extract(r'\((\d+)\)$')[0].astype('Int64')
cand = pd.read_csv(OUT + 'candidates.csv')
hcand = pd.read_csv(OUT + 'candidates_haem.csv')
# --- dependency_scan.csv in programme column style (plain Symbol first) ---
scan = m[['Gene', 'Column', 'AML_mean', 'AML_dep', 'AML_n', 'AML_dep_frac', 'Comp_mean', 'Comp_dep',
          'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q',
          'PanEssential', 'EnrichRank']].copy()
scan = scan.rename(columns={'Dep_p': 'p_welch', 'Dep_q': 'q_BH'})
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
hscan = pd.read_csv(OUT + 'dependency_scan_haem.csv')
hscan = hscan[['Gene', 'Column', 'AML_mean', 'AML_dep', 'AML_n', 'AML_dep_frac', 'Comp_mean', 'Comp_dep',
               'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'p_welch', 'q_BH',
               'PanEssential', 'EnrichRank']].copy()
hscan = hscan.sort_values('Enrichment')
hscan.to_csv(OUT + 'dependency_scan_haem.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'MYB': 'Myb transcription factor; no direct inhibitor (preclinical MYB inhibitors only)',
 'SPI1': 'PU.1 transcription factor; no direct inhibitor',
 'ZEB2': 'Transcription factor; no direct inhibitor',
 'LMO2': 'LIM-domain transcription factor/scaffold; no drug',
 'CBFB': 'CBF-beta subunit; no direct inhibitor — dependency is in wild-type-CBFB lines (no inv(16) line in the CRISPR panel; ME-1 lacks CRISPR data)',
 'GFI1': 'Transcriptional repressor; no direct inhibitor (LSD1 inhibition derepresses GFI1 — indirect)',
 'CDK6': 'CDK4/6 inhibitors approved in HR+/HER2- breast cancer (palbociclib, ribociclib, abemaciclib); CDK6-selective inhibitors in development; palbociclib tested in AML Phase 1/2',
 'RUNX1': 'Transcription factor; no direct inhibitor',
 'CEBPA': 'Transcription factor; no direct inhibitor',
 'CSF2RB': 'Common beta chain of IL-3/IL-5/GM-CSF receptors; no beta-chain-selective drug',
 'FLI1': 'ETS transcription factor; no direct inhibitor (TK216 targets EWSR1-FLI1 in Ewing, not wild-type FLI1)',
 'CCND3': 'Cyclin D3; no direct inhibitor — partner CDK6 is druggable (CDK4/6 inhibitors act on the CDK6-CCND3 complex)',
 'HHEX': 'Transcription factor; no direct inhibitor',
 'ZFP36L2': 'RNA-binding zinc-finger protein; no drug',
 'LDB1': 'LIM-domain-binding scaffold; no drug',
 'WDR49': 'WD-repeat protein; no drug — scored in only 5 of 28 lines (low-n flag)',
 'TYMS': 'Thymidylate synthase; inhibited by approved antifolates (5-fluorouracil, pemetrexed) — drug at the gene',
 'MEF2C': 'Transcription factor; no direct inhibitor',
 'CITED2': 'Transcriptional co-activator; no drug',
 'TAFAZZIN': 'Tafazzin, cardiolipin remodelling; no drug; no cell-line RNA data (candidate via TCGA leg only)',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Symbol'] = top20['Gene']
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'AML_mean', 'AML_dep', 'AML_n', 'Comp_mean', 'Comp_dep', 'Comp_n',
        'Comp_dep_frac', 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac',
        'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
        'Haem_Enrichment', 'Haem_FracDiff', 'Haem_EnrichRank',
        'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Gene'].tolist())
# haem top 20 (sensitivity list)
htop = hcand.head(20).copy()
htop['Rank'] = range(1, 21)
htop['Symbol'] = htop['Gene']
hcols = ['Rank', 'Symbol', 'Tier', 'AML_mean', 'AML_dep', 'AML_n', 'Haem_Enrichment', 'Haem_FracDiff',
         'Haem_EnrichRank', 'TPM_FC', 'DeltaLog2', 'TCGA_RNA_log2FC', 'Haem_Score']
htop[hcols].to_csv(OUT + 'top20_haem.csv', index=False)
print('haem top 20:', htop['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Gene', 'Column', 'Entrez', 'AML_mean', 'AML_dep', 'AML_n', 'AML_dep_frac', 'Comp_mean', 'Comp_dep',
        'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
        'EnrichRank', 'AML_mean_log2TPM', 'AML_RNA_n', 'Comp_mean_log2TPM', 'Comp_RNA_n', 'DeltaLog2', 'TPM_FC',
        'CellRNA_p', 'CellRNA_q', 'HaemComp_mean_log2TPM', 'HaemComp_RNA_n', 'DeltaLog2_haem', 'TPM_FC_haem',
        'p_haem', 'q_haem', 'TCGA_RNA_n_AML', 'TCGA_RNA_n_other', 'TCGA_RNA_mean_AML', 'TCGA_RNA_mean_other',
        'TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q',
        'Haem_Enrichment', 'Haem_FracDiff', 'Haem_AML_dep_frac', 'Haem_Comp_dep_frac', 'Haem_q', 'Haem_EnrichRank',
        'Haem_PanEssential', 'CellRNA_pass', 'TCGA_XC_pass', 'OverExpr_pass',
        'Strict_pass', 'Relaxed_pass', 'Haem_Strict_pass', 'Haem_Relaxed_pass']
keep = [c for c in keep if c in m.columns]
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; overexpression rows:', len(m), 'top20 rows:', len(top20))
