"""Ewing merge + finalize: intersect relaxed-bar genes with the over-expression
bar and rank by combined score.

Over-expression legs available for Ewing:
  (a) cell-line RNA (DepMap 24Q4) — the ONLY leg. TPM FC >= 1.5 AND
      DeltaLog2 >= 0.5 AND p < 0.05.
  (b) patient RNA — UNAVAILABLE. No TCGA Ewing study exists; no GDC Ewing
      project exists; the three Ewing studies in cBioPortal
      (es_dfarber_broad_2014, es_iocurie_2014, es_dsrct_msk_2023) carry
      mutations / structural variants / CNA only — no mRNA profile. The
      tumour-vs-normal and cross-cancer patient-RNA legs therefore do not
      exist for this cancer, and no patient-RNA column is imputed.
  (c) protein — UNAVAILABLE. No CPTAC sarcoma/Ewing study and no Ewing RPPA
      profile exist in cBioPortal. No +1.0 protein bonus fires for any gene.

Combined score (adapted per playbook §5, missing legs z-scored as 0):
  score = z(-enrichment) + z(frac diff) + z(cell-line DeltaLog2) + z(TCGA log2FC = 0)
  each z computed across the relaxed-bar genes. No protein bonus.
Scores are comparable within this cancer only.
"""
import pandas as pd, numpy as np

OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/ewing/'
scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan['GeneSym'] = scan['Gene'].str.split(' \\(').str[0]
rna = pd.read_csv(OUT + 'celline_rna.csv')

# ---- overexpression.csv: all genes, all evidence columns ----
# NB: both frames carry 'p'/'q' — suffixes put the RNA ones in p_rna/q_rna;
# the scan's dependency p/q stay as p/q. Do NOT rename 'p' itself.
m = scan.merge(rna[['Gene', 'Ewing_mean_log2TPM', 'Ewing_n', 'Comp_mean_log2TPM',
                    'Comp_n', 'DeltaLog2', 'TPM_FC', 'p', 'q']],
               left_on='GeneSym', right_on='Gene', how='left', suffixes=('', '_rna'))
m = m.rename(columns={'p': 'Dependency_p', 'q': 'Dependency_q',
                      'p_rna': 'CellLineRNA_p', 'q_rna': 'CellLineRNA_q'})
m['CellLineRNA_pass'] = (m['TPM_FC'] >= 1.5) & (m['DeltaLog2'] >= 0.5) & (m['CellLineRNA_p'] < 0.05)
m['TCGA_log2FC'] = np.nan
m['TCGA_p'] = np.nan
m['TCGA_pass'] = False
m['TCGA_note'] = 'no Ewing patient-RNA profile exists (see report Methods)'
m['Protein_delta'] = np.nan
m['Protein_p'] = np.nan
m['Protein_pass'] = False
m['Protein_note'] = 'no Ewing CPTAC/RPPA profile exists (see report Methods)'
m['OverExpr_pass'] = m['CellLineRNA_pass']  # only leg available
m.to_csv(OUT + 'overexpression.csv', index=False)
print('overexpression rows:', len(m), 'cell-line RNA bar pass:', int(m.CellLineRNA_pass.sum()))

# ---- bars + intersection ----
strict = scan[(scan.Enrichment <= -0.10) & (scan.Ewing_dep_frac >= 0.25) &
              (scan.FracDiff >= 0.10) & (~scan.PanEssential) &
              (scan.GenomeWide_dep_frac <= 0.60)].copy()
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.Ewing_dep_frac >= 0.20) &
               (scan.FracDiff >= 0.10) & (~scan.PanEssential) &
               (scan.GenomeWide_dep_frac <= 0.60)].copy()
print('strict:', len(strict), 'relaxed:', len(relaxed))

def z(s):
    return (s - s.mean()) / s.std(ddof=1)

cand = relaxed.merge(m[['GeneSym', 'DeltaLog2', 'TPM_FC', 'CellLineRNA_p',
                        'CellLineRNA_pass', 'OverExpr_pass']],
                     on='GeneSym', how='left')
cand = cand[cand['OverExpr_pass'] == True].copy()
print('relaxed + OE:', len(cand))
cand_strict = cand[cand['Gene'].isin(set(strict['Gene']))]
print('strict + OE:', len(cand_strict))

# combined score across the RELAXED-bar genes (playbook §5: z across relaxed-bar genes)
rel = relaxed.merge(m[['GeneSym', 'DeltaLog2']], on='GeneSym', how='left')
rel['z_neg_enrich'] = z(-rel['Enrichment'])
rel['z_fracdiff'] = z(rel['FracDiff'])
rel['z_celline'] = z(rel['DeltaLog2'].fillna(0))
rel['z_tcga'] = 0.0  # no patient-RNA leg exists for Ewing
rel['Score'] = rel['z_neg_enrich'] + rel['z_fracdiff'] + rel['z_celline'] + rel['z_tcga']
scores = rel[['GeneSym', 'Score']].set_index('GeneSym')['Score']
cand['Score'] = cand['GeneSym'].map(scores)
cand['Tier'] = np.where(cand['Gene'].isin(set(strict['Gene'])), 'A', 'B')
cand = cand.sort_values('Score', ascending=False).reset_index(drop=True)
cand['CandidateRank'] = cand.index + 1
cand.to_csv(OUT + 'candidates.csv', index=False)
top20 = cand.head(20).copy()
top20.to_csv(OUT + 'top20.csv', index=False)
print(cand.head(30)[['CandidateRank', 'GeneSym', 'Tier', 'Ewing_mean', 'Ewing_dep',
                     'Enrichment', 'FracDiff', 'TPM_FC', 'DeltaLog2', 'Score']].to_string(index=False))
print()
print('ranks 21-25:')
print(cand.iloc[20:25][['CandidateRank', 'GeneSym', 'Tier', 'Ewing_dep', 'Enrichment', 'TPM_FC', 'Score']].to_string(index=False))
print()
print('low-n candidates (Ewing_n < 27):')
print(cand[cand.Ewing_n < 27][['CandidateRank', 'GeneSym', 'Tier', 'Ewing_n', 'Ewing_dep', 'Enrichment', 'Score']].to_string(index=False))
print()
print('low-n relaxed-bar genes NOT candidates:')
rel_nc = relaxed[~relaxed.Gene.isin(set(cand.Gene))]
print(rel_nc[rel_nc.Ewing_n < 27][['GeneSym', 'Ewing_n', 'Ewing_dep', 'Enrichment']].to_string(index=False))
