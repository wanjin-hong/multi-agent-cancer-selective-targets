import pandas as pd, numpy as np
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/'
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
# --- dependency_scan.csv in TNBC column style: plain Gene symbol + EnrichRank ---
scan = m[['Symbol', 'Gene', 'NB_mean', 'NB_dep', 'NB_n', 'NB_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
          'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'Dep_p', 'Dep_q', 'PanEssential']].copy()
scan = scan.rename(columns={'Gene': 'Column', 'Dep_p': 'p_welch', 'Dep_q': 'q_BH', 'Symbol': 'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv(OUT + 'dependency_scan.csv', index=False)
# --- top20.csv with druggability (filled in after the merge; edit DRUG below) ---
drug = {
 'ISL1': 'Transcription factor (ISL LIM homeobox 1, sympathoadrenal lineage); no direct inhibitor',
 'PHOX2A': 'Transcription factor (paired-like homeobox 2A / ARIX, autonomic-neuron lineage); no direct inhibitor',
 'MYCN': 'Transcription factor (N-Myc, the disease driver); no direct MYCN inhibitor — indirect routes only: Aurora A inhibition destabilises MYCN protein (alisertib, Phase 2 in neuroblastoma), BET inhibition (preclinical/Phase 1), ODC1 inhibition with DFMO (difluoromethylornithine, in neuroblastoma trials)',
 'HAND2': 'Transcription factor (heart and neural crest derivatives-expressed 2, sympathoadrenal lineage); no direct inhibitor',
 'GATA3': 'Transcription factor (GATA-binding factor 3, sympathoadrenal lineage); no direct inhibitor',
 'LDB1': 'Transcriptional cofactor (LIM-domain binding 1, scaffolds ISL1/HAND2 complexes); no drug',
 'SOX11': 'Transcription factor (SRY-box 11, neuronal lineage); no direct inhibitor',
 'IRS2': 'Adaptor protein (insulin receptor substrate 2, IGF1R/PI3K signalling); no selective drug — upstream IGF1R inhibitors exist but none approved',
 'ASCL1': 'Transcription factor (achaete-scute homolog 1, neuroendocrine lineage); no direct inhibitor',
 'TUBE1': 'Tubulin epsilon 1 (centriole/microtubule protein); no selective drug (microtubule agents are non-selective)',
 'CDIN1': 'CDAN1-interacting nuclease 1 (chromatin/cell-cycle protein); no drug',
 'CDK6': 'CDK4/6 inhibitors approved in HR+/HER2- breast cancer (palbociclib, ribociclib, abemaciclib) — not neuroblastoma-labelled; third cancer for CDK6 in this programme (TNBC, HNSC, neuroblastoma)',
 'RTTN': 'Rotatin (centrosome/cilia protein, mitotic spindle positioning); no drug',
 'SMARCA4': 'BRG1, SWI/SNF chromatin-remodelling ATPase; no direct inhibitor (BRG1/BRM degraders preclinical; SMARCA4-deficient contexts sensitised to EZH2 inhibition — indirect only)',
 'PPM1D': 'Phosphatase WIP1 (p53-regulating PP2C phosphatase); no clinical inhibitor (GSK2830371 and analogues preclinical only)',
 'RCOR2': 'REST corepressor 2 (neuronal-gene repression complex); no drug',
 'CCNL1': 'Cyclin L1 (CDK11 partner); no selective drug (CDK11 inhibitors preclinical only)',
 'TRA2B': 'Splicing factor (transformer-2 beta, RNA processing); no drug',
 'MSI2': 'Musashi-2, RNA-binding protein (stem-cell regulator); no clinical inhibitor (preclinical MSI2 inhibitors only)',
 'MDM4': 'MDMX, p53 negative regulator; no approved MDM4 inhibitor — dual MDM2/MDM4 inhibitors preclinical/early Phase 1 only',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, len(top20) + 1)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank', 'Symbol', 'Tier', 'NB_mean', 'NB_dep', 'NB_n', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac', 'Enrichment', 'FracDiff',
        'GenomeWide_dep_frac', 'TPM_FC', 'DeltaLog2', 'CellRNA_p', 'CellRNA_q', 'TCGA_XC_log2FC', 'TCGA_XC_p', 'TCGA_XC_q',
        'Score', 'Druggability']
top20[cols].to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist() if len(cand) > 20 else 'n/a')
# tidy overexpression.csv column order
keep = ['Symbol', 'Gene', 'Entrez', 'NB_mean', 'NB_dep', 'NB_n', 'NB_dep_frac', 'Comp_mean', 'Comp_dep', 'Comp_n', 'Comp_dep_frac',
 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac', 'PanEssential', 'Dep_p', 'Dep_q',
 'NB_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC', 'CellRNA_p', 'CellRNA_q',
 'TCGA_XC_n_NBL', 'TCGA_XC_n_other', 'TCGA_XC_mean_NBL', 'TCGA_XC_mean_other', 'TCGA_XC_log2FC', 'TCGA_XC_p', 'TCGA_XC_q',
 'CellRNA_pass', 'TCGA_XC_pass', 'OverExpr_pass', 'Strict_pass', 'Relaxed_pass']
m[keep].to_csv(OUT + 'overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
print('top20 missing druggability notes:', top20[top20.Druggability.isna()]['Symbol'].tolist())
