"""Process cBioPortal TCGA data for endometrial (UCEC).

- Cross-cancer mRNA (context only, NOT part of the over-expression bar):
  UCEC primary tumours vs pooled primary tumours of 6 other cancers
  (BRCA, LUAD, LUSC, COADREAD, STAD, LIHC), RSEM log2-transformed here.
  All UCEC histologies are used — the patient cohort is the study as a whole,
  as in HNSC/BLCA.
- RPPA protein (the protein bar): UCEC primary tumours vs pooled
  other cancers. LIHC has no RPPA profile (404) — pooled others for RPPA are
  the 5 remaining studies (BRCA, LUAD, LUSC, COADREAD, STAD).
- CPTAC protein (ucec_cptac_2020_protein_quantification): descriptive only;
  sample composition checked by TCGA barcode — if no normals, no bar.
Primary-tumour samples only, by TCGA barcode sample-type 01.
The over-expression bar's TCGA leg is the GDC tumour-vs-normal axis
(tcga_tumour_vs_normal.csv), not this cross-cancer file.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/endometrial/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/endometrial/'
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
fetch_genes = json.load(open(os.path.join(DL, 'fetch_genes.json')))
inv = {entrez[g]: g for g in fetch_genes}

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
    return sid[13:15] if isinstance(sid, str) and len(sid) >= 15 else None

# ---------- cross-cancer mRNA ----------
ucec = load_profile('tcga_ucec_mrna.json')
print('UCEC mRNA rows:', len(ucec), 'samples:', ucec['sampleId'].nunique(),
      'primary tumours:', ucec[ucec['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
others = pd.concat([load_profile(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
print('other-cancer mRNA rows:', len(others), 'samples:', others['sampleId'].nunique(),
      'primary tumours:', others[others['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
ucec = ucec[ucec['sampleId'].apply(is_primary_tumour)].copy()
others = others[others['sampleId'].apply(is_primary_tumour)].copy()
ucec['value'] = np.log2(ucec['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
rows = []
for eid, g in inv.items():
    a = ucec.loc[ucec['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=eid, TCGA_RNA_n_UCEC=len(a), TCGA_RNA_n_other=len(b),
                     TCGA_RNA_mean_UCEC=am, TCGA_RNA_mean_other=bm,
                     TCGA_RNA_log2FC=am - bm, TCGA_RNA_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_RNA_q'] = bh(xc['TCGA_RNA_p'].values)
xc = xc.sort_values('TCGA_RNA_log2FC', ascending=False)
xc.to_csv(OUT + 'tcga_rna.csv', index=False)
bar = xc[(xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)]
print('cross-cancer bar pass:', len(bar))
print(xc.head(15)[['Gene', 'TCGA_RNA_n_UCEC', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_q']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'PAX8', 'HNF1B', 'CCNE1', 'PARD6B', 'ARHGAP29', 'PARD3', 'KLF5', 'FOSL1', 'SOX9', 'ERBB2', 'CDK6', 'ITGAV', 'FERMT2', 'ITGB1', 'ITGB5', 'ITGA3', 'MYBL2', 'WT1', 'SOX17']:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']].to_string(index=False) if len(r) else 'NOT FOUND')

# ---------- RPPA ----------
ucec_r = load_profile('tcga_ucec_rppa.json')
print('\nUCEC RPPA rows:', len(ucec_r), 'samples:', ucec_r['sampleId'].nunique(),
      'primary tumours:', ucec_r[ucec_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
other_r = pd.concat([load_profile(f'tcga_{c}_rppa.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
print('other RPPA rows:', len(other_r), 'samples:', other_r['sampleId'].nunique(),
      'primary tumours:', other_r[other_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
ucec_r = ucec_r[ucec_r['sampleId'].apply(is_primary_tumour)].copy()
other_r = other_r[other_r['sampleId'].apply(is_primary_tumour)].copy()
rrows = []
for eid, g in inv.items():
    a = ucec_r.loc[ucec_r['entrezGeneId'] == eid, 'value'].values
    b = other_r.loc[other_r['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rrows.append(dict(Gene=g, Entrez=eid, RPPA_n_UCEC=len(a), RPPA_n_other=len(b),
                      RPPA_mean_UCEC=am, RPPA_mean_other=bm, RPPA_delta=am - bm, RPPA_p=p))
rp = pd.DataFrame(rrows)
rp = rp.sort_values('RPPA_delta', ascending=False)
rp.to_csv(OUT + 'tcga_rppa.csv', index=False)
print(rp[rp.RPPA_n_UCEC > 0][['Gene', 'RPPA_n_UCEC', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p']].to_string())
bar_r = rp[(rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)]
print('RPPA bar pass:', len(bar_r), bar_r['Gene'].tolist() if len(bar_r) else '')

# ---------- CPTAC (descriptive; check sample composition) ----------
cp = load_profile('tcga_ucec_cptac.json')
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
