"""Process cBioPortal patient-RNA data for osteosarcoma.

- Cross-cancer TPM (THE TCGA leg of the over-expression bar, design (c)):
  TARGET osteosarcoma primary tumours (os_target_gdc, 88 with TPM data)
  vs pooled primary tumours of 6 other cancers (BRCA, LUAD, LUSC,
  COADREAD, STAD, LIHC) from the GDC-harmonised TCGA studies
  (*_tcga_gdc_mrna_seq_tpm). All values are GDC STAR TPM, log2(TPM+1)
  here. Pre-registered design: there is no TCGA osteosarcoma study and
  no adjacent normals, so tumour-vs-normal is impossible — same design
  as melanoma/colon/glioma/ovarian, with TARGET as the tumour source.
  Scale check: median log2(TPM+1) of the two pooled sets must be within
  0.5 of each other, else the leg is rejected (SCLC precedent).
- Protein: NO protein leg. No CPTAC osteosarcoma study exists in
  cBioPortal and os_target_gdc carries no RPPA profile. The +1.0 bonus
  fires for nobody; recorded, not imputed.
Primary-tumour samples only, by TCGA barcode sample-type 01.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/osteosarcoma/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/osteosarcoma/'
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
inv = {v: k for k, v in entrez.items()}

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

def load_profile(fname):
    return pd.DataFrame(json.load(open(os.path.join(DL, fname))))

def is_primary_tumour(sid):
    # TCGA barcodes: 'TCGA-3C-AAAU-01' (sample type = last field, '01');
    # TARGET barcodes: 'TARGET-40-0A4HLD-01A' (sample type = last field,
    # '01A'). The old sid[13:15] rule only fits TCGA.
    return isinstance(sid, str) and sid.split('-')[-1].startswith('01')

os_ = load_profile('tcga_os_tpm.json')
print('OS TPM rows:', len(os_), 'samples:', os_['sampleId'].nunique(),
      'primary tumours:', os_[os_['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
others = pd.concat([load_profile(f'tcga_{c}_tpm.json') for c in ['brca', 'luad', 'lusc', 'coad', 'stad', 'hcc']])
print('other-cancer TPM rows:', len(others), 'samples:', others['sampleId'].nunique(),
      'primary tumours:', others[others['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
os_ = os_[os_['sampleId'].apply(is_primary_tumour)].copy()
others = others[others['sampleId'].apply(is_primary_tumour)].copy()
os_['value'] = np.log2(os_['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
# scale check: pooled medians of log2(TPM+1)
print('SCALE CHECK median log2(TPM+1): OS %.3f vs other %.3f, diff %.3f'
      % (os_['value'].median(), others['value'].median(), abs(os_['value'].median() - others['value'].median())))
print('SCALE CHECK mean log2(TPM+1): OS %.3f vs other %.3f'
      % (os_['value'].mean(), others['value'].mean()))
rows = []
for eid, g in inv.items():
    a = os_.loc[os_['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=eid, TCGA_RNA_n_OS=len(a), TCGA_RNA_n_other=len(b),
                     TCGA_RNA_mean_OS=am, TCGA_RNA_mean_other=bm,
                     TCGA_RNA_log2FC=am - bm, TCGA_RNA_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_RNA_q'] = bh(xc['TCGA_RNA_p'].values)
xc = xc.sort_values('TCGA_RNA_log2FC', ascending=False)
xc.to_csv(OUT + 'tcga_rna.csv', index=False)
bar = xc[(xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)]
print('cross-cancer TCGA bar pass:', len(bar))
print(xc.head(15)[['Gene', 'TCGA_RNA_n_OS', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_q']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'CCNE1', 'TP53', 'RB1', 'MDM2', 'CDK4', 'CDK6', 'FGFR1', 'PDGFRA', 'IGF1R', 'MET', 'SRC', 'RUNX2', 'SATB2', 'SP7', 'DLX5', 'KLF5', 'HNF1B', 'FOSL1', 'PAX8', 'CXADR', 'PARD6B', 'ITGB5', 'ITGB1']:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']].to_string(index=False) if len(r) else 'NOT FOUND')
