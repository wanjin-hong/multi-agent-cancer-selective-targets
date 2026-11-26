import pandas as pd, numpy as np
m = pd.read_csv('overexpression.csv')
cand = pd.read_csv('candidates.csv')
# --- dependency_scan.csv in TNBC column style ---
scan = m[['Symbol','Gene','HCC_mean','HCC_dep','HCC_n','HCC_dep_frac','Comp_mean','Comp_dep','Comp_n','Comp_dep_frac',
          'Enrichment','FracDiff','GenomeWide_dep_frac','Dep_p','Dep_q','PanEssential']].copy()
scan = scan.rename(columns={'Gene':'Column','Dep_p':'p_welch','Dep_q':'q_BH','Symbol':'Gene'})
scan['EnrichRank'] = scan['Enrichment'].rank(method='min', na_option='bottom').astype('Int64')
scan = scan.sort_values('Enrichment')
scan.to_csv('dependency_scan.csv', index=False)
# --- top20.csv with druggability ---
drug = {
 'ITGAV': 'No approved integrin-alpha-v drug; cilengitide (alpha-v beta3/beta5) failed Phase 3 in glioblastoma; abituzumab discontinued',
 'HNF1B': 'Transcription factor; no direct inhibitor; HCC/liver lineage marker',
 'XYLT2': 'Xylosyltransferase 2; no drug',
 'PTK2': 'FAK inhibitors investigational only: defactinib, ifebemtinib (IN10018), avutometinib (VS-6766); none approved',
 'FGFR1': 'Approved pan-FGFR inhibitors: erdafitinib (FGFR1-4), pemigatinib (FGFR1-3), futibatinib (FGFR1-4) - labelled for urothelial/cholangiocarcinoma, not HCC; HCC-relevant axis is FGF19-FGFR4, not FGFR1',
 'PPRC1': 'Transcriptional coactivator; no drug',
 'FERMT2': 'Kindlin-2; no drug',
 'ITGB5': 'Integrin beta5 (pairs with alpha-v); covered by cilengitide historically; no beta5-selective drug',
 'SLC35B2': 'Nucleotide-sugar transporter; no drug',
 'ELMO2': 'No drug',
 'ACTB': 'Beta-actin; undruggable as a target (essential cytoskeleton)',
 'NCKAP1': 'WAVE complex subunit; no drug',
 'PARD3': 'Cell polarity protein; no drug',
 'PDCL3': 'Phosducin-like 3 chaperone; no drug',
 'STT3A': 'Oligosaccharyltransferase catalytic subunit; OST inhibitor NGI-1 preclinical only',
 'ARHGAP29': 'RhoGAP; no drug',
 'FAF2': 'UBXD8, ER-associated degradation; no drug',
 'ACTR3': 'ARP3, Arp2/3 complex; CK-666 preclinical tool inhibitor only',
 'LMNA': 'Lamin A/C; no drug',
 'UBA5': 'E1 enzyme for UFMylation; preclinical UBA5 inhibitors only',
}
top20 = cand.head(20).copy()
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(drug)
cols = ['Rank','Symbol','Tier','HCC_mean','HCC_dep','HCC_n','Comp_mean','Comp_dep','Comp_n','Comp_dep_frac','Enrichment','FracDiff',
        'GenomeWide_dep_frac','TPM_FC','DeltaLog2','CellRNA_p','CellRNA_q','GDC_TvN_log2FC','GDC_TvN_p','GDC_TvN_q',
        'TCGA_RNA_log2FC','TCGA_RNA_p','TCGA_RNA_q','TCGA_refNormal_z_mean','RPPA_delta','RPPA_p','Score','Druggability']
top20[cols].to_csv('top20.csv', index=False)
print(top20[['Rank','Symbol','Tier','Score']].to_string(index=False))
print('ranks 21-25:', cand.iloc[20:25]['Symbol'].tolist())
# tidy overexpression.csv column order
keep = ['Symbol','Gene','Entrez','HCC_mean','HCC_dep','HCC_n','HCC_dep_frac','Comp_mean','Comp_dep','Comp_n','Comp_dep_frac',
 'Enrichment','FracDiff','GenomeWide_dep_frac','PanEssential','Dep_p','Dep_q',
 'HCC_mean_log2TPM','Comp_mean_log2TPM','DeltaLog2','TPM_FC','CellRNA_p','CellRNA_q',
 'GDC_TvN_n_tumour','GDC_TvN_n_normal','GDC_TvN_mean_tumour','GDC_TvN_mean_normal','GDC_TvN_log2FC','GDC_TvN_p','GDC_TvN_q',
 'TCGA_RNA_n_LIHC','TCGA_RNA_n_other','TCGA_RNA_mean_LIHC','TCGA_RNA_mean_other','TCGA_RNA_log2FC','TCGA_RNA_p','TCGA_RNA_q',
 'TCGA_refNormal_z_mean','RPPA_n_LIHC','RPPA_n_other','RPPA_delta','RPPA_p',
 'CellRNA_pass','TCGA_TvN_pass','TCGA_XC_pass','RPPA_pass','OverExpr_pass','Strict_pass','Relaxed_pass']
m[keep].to_csv('overexpression.csv', index=False)
print('files written; rows:', len(m), len(top20))
