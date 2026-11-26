"""Process cBioPortal TCGA data for ESGA.

- Cross-cancer mRNA (context only, NOT part of the over-expression bar):
  STAD + ESCA-adenocarcinoma primary tumours vs pooled primary tumours of 6
  other cancers (BRCA, LUAD, LUSC, COADREAD, LIHC, HNSC), RSEM
  log2-transformed here. ESCA tumours are filtered to adenocarcinoma by TCGA
  barcode prefix, using the adenocarcinoma case list from the GDC cases API
  (gdc_esca_cases.json) — ESCA is mixed squamous + adenocarcinoma.
- RPPA protein (the protein bar): pooled STAD + ESCA-adenocarcinoma primary
  tumours vs pooled other-cancer primary tumours (BRCA, LUAD, LUSC, COADREAD,
  HNSC — LIHC has no RPPA profile).
Primary-tumour samples only, by TCGA barcode sample-type 01.
The over-expression bar's TCGA leg is the GDC tumour-vs-normal axis
(tcga_tumour_vs_normal.csv), not this cross-cancer file.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/'
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
inv = {v: k for k, v in entrez.items()}

# ESCA adenocarcinoma sample prefixes (TCGA barcode first 12 chars)
cases = json.load(open(os.path.join(DL, 'gdc_esca_cases.json')))
adeno_prefixes = set()
for c in cases:
    diags = [d.get('primary_diagnosis', '') for d in c.get('diagnoses', [])]
    if any('adenocarcinoma' in d.lower() for d in diags):
        adeno_prefixes.add(c['submitter_id'][:12])
print('ESCA adenocarcinoma case prefixes:', len(adeno_prefixes))

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

def is_esca_adeno(sid):
    return isinstance(sid, str) and len(sid) >= 12 and sid[:12] in adeno_prefixes

# ---------- cross-cancer mRNA ----------
stad = load_profile('tcga_stad_mrna.json')
esca = load_profile('tcga_esca_mrna.json')
print('STAD mRNA rows:', len(stad), 'samples:', stad['sampleId'].nunique(),
      'primary tumours:', stad[stad['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
print('ESCA mRNA rows:', len(esca), 'samples:', esca['sampleId'].nunique(),
      'primary tumours:', esca[esca['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique(),
      'adeno primary tumours:', esca[esca['sampleId'].apply(lambda s: is_primary_tumour(s) and is_esca_adeno(s))]['sampleId'].nunique())
print('ESCA mRNA genes missing from fetch:', sorted(set(entrez.values()) - set(esca['entrezGeneId'].unique())))
others = pd.concat([load_profile(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'lihc', 'hnsc']])
print('other-cancer mRNA rows:', len(others), 'samples:', others['sampleId'].nunique(),
      'primary tumours:', others[others['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
esga_m = pd.concat([stad[stad['sampleId'].apply(is_primary_tumour)],
                    esca[esca['sampleId'].apply(lambda s: is_primary_tumour(s) and is_esca_adeno(s))]]).copy()
others = others[others['sampleId'].apply(is_primary_tumour)].copy()
esga_m['value'] = np.log2(esga_m['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
rows = []
for eid, g in inv.items():
    a = esga_m.loc[esga_m['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=eid, TCGA_RNA_n_ESGA=len(a), TCGA_RNA_n_other=len(b),
                     TCGA_RNA_mean_ESGA=am, TCGA_RNA_mean_other=bm,
                     TCGA_RNA_log2FC=am - bm, TCGA_RNA_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_RNA_q'] = bh(xc['TCGA_RNA_p'].values)
xc = xc.sort_values('TCGA_RNA_log2FC', ascending=False)
xc.to_csv(OUT + 'tcga_rna.csv', index=False)
bar = xc[(xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)]
print('cross-cancer bar pass:', len(bar))
print(xc.head(15)[['Gene', 'TCGA_RNA_n_ESGA', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_q']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'KLF5', 'SOX9', 'HNF1B', 'FOSL1', 'CDX2', 'ERBB2', 'CDH17', 'GATA6', 'CTNNB1', 'TCF7L2', 'GFPT1']:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']].to_string(index=False) if len(r) else 'NOT FOUND')

# ---------- RPPA ----------
stad_r = load_profile('tcga_stad_rppa.json')
esca_r = load_profile('tcga_esca_rppa.json')
print('\nSTAD RPPA rows:', len(stad_r), 'samples:', stad_r['sampleId'].nunique(),
      'primary tumours:', stad_r[stad_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
print('ESCA RPPA rows:', len(esca_r), 'samples:', esca_r['sampleId'].nunique(),
      'primary tumours:', esca_r[esca_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique(),
      'adeno primary tumours:', esca_r[esca_r['sampleId'].apply(lambda s: is_primary_tumour(s) and is_esca_adeno(s))]['sampleId'].nunique())
print('RPPA genes measured (STAD):', sorted(stad_r['entrezGeneId'].map(lambda e: inv.get(e, e)).unique()))
print('RPPA genes measured (ESCA):', sorted(esca_r['entrezGeneId'].map(lambda e: inv.get(e, e)).unique()))
other_r = pd.concat([load_profile(f'tcga_{c}_rppa.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'hnsc']])
print('other RPPA rows:', len(other_r), 'samples:', other_r['sampleId'].nunique(),
      'primary tumours:', other_r[other_r['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
esga_r = pd.concat([stad_r[stad_r['sampleId'].apply(is_primary_tumour)],
                    esca_r[esca_r['sampleId'].apply(lambda s: is_primary_tumour(s) and is_esca_adeno(s))]]).copy()
other_r = other_r[other_r['sampleId'].apply(is_primary_tumour)].copy()
rrows = []
for eid, g in inv.items():
    a = esga_r.loc[esga_r['entrezGeneId'] == eid, 'value'].values
    b = other_r.loc[other_r['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rrows.append(dict(Gene=g, Entrez=eid, RPPA_n_ESGA=len(a), RPPA_n_other=len(b),
                      RPPA_mean_ESGA=am, RPPA_mean_other=bm, RPPA_delta=am - bm, RPPA_p=p))
rp = pd.DataFrame(rrows)
rp = rp.sort_values('RPPA_delta', ascending=False)
rp.to_csv(OUT + 'tcga_rppa.csv', index=False)
print(rp[rp.RPPA_n_ESGA > 0][['Gene', 'RPPA_n_ESGA', 'RPPA_n_other', 'RPPA_delta', 'RPPA_p']].to_string())
bar_r = rp[(rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)]
print('RPPA bar pass:', len(bar_r), bar_r['Gene'].tolist() if len(bar_r) else '')
