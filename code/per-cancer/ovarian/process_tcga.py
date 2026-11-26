"""Process cBioPortal TCGA data for ovarian cancer.

- Cross-cancer mRNA (THE TCGA leg of the over-expression bar, design (c)):
  OV primary tumours vs pooled primary tumours of 6 other cancers
  (BRCA, LUAD, LUSC, COADREAD, STAD, LIHC), RSEM log2-transformed here.
  This is the pre-registered TCGA design: TCGA-OV has ZERO adjacent normals
  in GDC (426 primary tumours, 8 recurrent, 0 normal), so tumour-vs-normal
  is impossible — same design as melanoma/colon/glioma.
- RPPA protein (the protein bar): OV vs pooled other cancers.
  LIHC has no RPPA profile (404) — pooled others for RPPA are the 5 remaining
  studies (BRCA, LUAD, LUSC, COADREAD, STAD).
- CPTAC protein (ov_tcga_protein_quantification): sample composition checked
  by TCGA barcode; if no normal samples exist, values are reported
  descriptively (tumour medians) only — no bar, no bonus.
Primary-tumour samples only, by TCGA barcode sample-type 01.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/ovarian/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/ovarian/'
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
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'

def sample_type(sid):
    return sid[13:15] if isinstance(sid, str) and len(sid) >= 15 else '??'

# ---------- cross-cancer mRNA (the TCGA leg) ----------
ov = load_profile('tcga_ov_mrna.json')
print('OV mRNA rows:', len(ov), 'samples:', ov['sampleId'].nunique(),
      'primary tumours:', ov[ov['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
others = pd.concat([load_profile(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
print('other-cancer mRNA rows:', len(others), 'samples:', others['sampleId'].nunique(),
      'primary tumours:', others[others['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
ov = ov[ov['sampleId'].apply(is_primary_tumour)].copy()
others = others[others['sampleId'].apply(is_primary_tumour)].copy()
ov['value'] = np.log2(ov['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
rows = []
for eid, g in inv.items():
    a = ov.loc[ov['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=eid, TCGA_RNA_n_OV=len(a), TCGA_RNA_n_other=len(b),
                     TCGA_RNA_mean_OV=am, TCGA_RNA_mean_other=bm,
                     TCGA_RNA_log2FC=am - bm, TCGA_RNA_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_RNA_q'] = bh(xc['TCGA_RNA_p'].values)
xc = xc.sort_values('TCGA_RNA_log2FC', ascending=False)
xc.to_csv(OUT + 'tcga_rna.csv', index=False)
bar = xc[(xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)]
print('cross-cancer TCGA bar pass:', len(bar))
print(xc.head(15)[['Gene', 'TCGA_RNA_n_OV', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_q']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'CCNE1', 'FOXM1', 'PAX8', 'WT1']:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']].to_string(index=False) if len(r) else 'NOT FOUND')

# ---------- RPPA (protein bar) ----------
ov_r = load_profile('tcga_ov_rppa.json')
print('\nOV RPPA rows:', len(ov_r), 'samples:', ov_r['sampleId'].nunique(),
      'primary tumours:', ov_r[ov_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
other_r = pd.concat([load_profile(f'tcga_{c}_rppa.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
print('other RPPA rows:', len(other_r), 'samples:', other_r['sampleId'].nunique(),
      'primary tumours:', other_r[other_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
ov_r = ov_r[ov_r['sampleId'].apply(is_primary_tumour)].copy()
other_r = other_r[other_r['sampleId'].apply(is_primary_tumour)].copy()
rrows = []
for eid, g in inv.items():
    a = ov_r.loc[ov_r['entrezGeneId'] == eid, 'value'].values
    b = other_r.loc[other_r['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rrows.append(dict(Gene=g, Entrez=eid, RPPA_n_OV=len(a), RPPA_n_other=len(b),
                      RPPA_mean_OV=am, RPPA_mean_other=bm, RPPA_delta=am - bm, RPPA_p=p))
rp = pd.DataFrame(rrows)
rp = rp.sort_values('RPPA_delta', ascending=False)
rp.to_csv(OUT + 'tcga_rppa.csv', index=False)
print(rp[rp.RPPA_n_OV > 0][['Gene', 'RPPA_n_OV', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p']].to_string())
bar_r = rp[(rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)]
print('RPPA bar pass:', len(bar_r), bar_r['Gene'].tolist() if len(bar_r) else '')

# ---------- CPTAC (descriptive; check sample composition) ----------
cp = load_profile('tcga_ov_cptac.json')
print('\nCPTAC rows:', len(cp), 'samples:', cp['sampleId'].nunique())
types = cp['sampleId'].apply(sample_type).value_counts(dropna=False).to_dict()
print('CPTAC sample types (barcode 01=primary tumour, 11=normal):', types)
crows = []
for eid, g in inv.items():
    a = cp.loc[cp['entrezGeneId'] == eid, 'value'].astype(float).values
    crows.append(dict(Gene=g, Entrez=eid, CPTAC_n_tumour=len(a),
                      CPTAC_median=np.median(a) if len(a) else np.nan,
                      CPTAC_mean=np.mean(a) if len(a) else np.nan))
cpt = pd.DataFrame(crows)
cpt.to_csv(OUT + 'tcga_cptac.csv', index=False)
print(cpt.sort_values('CPTAC_median', ascending=False).head(15).to_string(index=False))
