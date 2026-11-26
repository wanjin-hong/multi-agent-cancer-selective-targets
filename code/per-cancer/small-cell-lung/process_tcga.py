"""Process cBioPortal patient RNA for SCLC.

Cross-cancer mRNA (THE patient leg of the over-expression bar, design (c)):
SCLC primary tumours from sclc_ucologne_2015 (George et al., Nature 2015;
81 tumours, RNA-seq RSEM — log2-transformed here) vs pooled primary tumours
of 6 other TCGA cancers (BRCA, LUAD, LUSC, COADREAD, STAD, LIHC), RSEM
log2-transformed here. Pre-registered design: there is no TCGA SCLC study
and ucologne has tumours only (no normals), so tumour-vs-normal is
impossible — same design as melanoma/colon/glioma/ovarian.
Primary-tumour samples only, by TCGA barcode sample-type 01 (ucologne
sample IDs are not TCGA barcodes; all 81 are primary tumours per the study).
Protein: no SCLC RPPA/CPTAC profile exists — no protein leg, no bonus.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/small-cell-lung/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/small-cell-lung/'
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

sclc = load_profile('tcga_sclc_mrna.json')
print('SCLC mRNA rows:', len(sclc), 'samples:', sclc['sampleId'].nunique(),
      'genes:', sclc['entrezGeneId'].nunique())
others = pd.concat([load_profile(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
print('other-cancer mRNA rows:', len(others), 'samples:', others['sampleId'].nunique(),
      'primary tumours:', others[others['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
others = others[others['sampleId'].apply(is_primary_tumour)].copy()
sclc['value'] = np.log2(sclc['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
fetch_genes = json.load(open(os.path.join(DL, 'fetch_genes.json')))
print('SCLC genes missing from fetch:', [g for g in fetch_genes if entrez.get(g) not in set(sclc['entrezGeneId'])])
rows = []
for g in fetch_genes:
    eid = entrez.get(g)
    a = sclc.loc[sclc['entrezGeneId'] == eid, 'value'].values if eid is not None else np.array([])
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values if eid is not None else np.array([])
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=eid, TCGA_RNA_n_SCLC=len(a), TCGA_RNA_n_other=len(b),
                     TCGA_RNA_mean_SCLC=am, TCGA_RNA_mean_other=bm,
                     TCGA_RNA_log2FC=am - bm, TCGA_RNA_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_RNA_q'] = bh(xc['TCGA_RNA_p'].values)
xc = xc.sort_values('TCGA_RNA_log2FC', ascending=False)
xc.to_csv(OUT + 'tcga_rna.csv', index=False)
bar = xc[(xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)]
print('cross-cancer TCGA bar pass (of fetched):', len(bar))
print(xc.head(20)[['Gene', 'TCGA_RNA_n_SCLC', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_q']].to_string())
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'MYCL', 'MYCN', 'ASCL1', 'NEUROD1', 'POU2F3', 'INSM1', 'DLL3', 'BCL2', 'EZH2', 'KDM1A', 'RB1', 'TP53', 'SOX10', 'KLF5', 'FOSL1', 'NEUROG2', 'CHGA', 'SYP', 'FOXM1', 'CCNE1', 'CCNE2', 'E2F3', 'SKP2', 'SLC2A1', 'NAMPT', 'SMARCA4', 'EP300', 'DOT1L', 'TIMM50', 'RNASEK']:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']].to_string(index=False) if len(r) else 'NOT FOUND')
