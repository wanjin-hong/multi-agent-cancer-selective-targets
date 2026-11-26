"""TARGET-NBL vs pooled other-cancer tumours, GDC STAR counts.

Design (c), pre-registered in gdc_fetch.py: all 153 TARGET-NBL primary
tumours vs 150 pooled primary tumours (25 each BRCA/LUAD/LUSC/COAD/
STAD/LIHC, seed 20261121), TPM, log2(TPM+1). Genes restricted to the
19,193 protein-coding genes of the DepMap expression matrix so the
table aligns with the scan. Scale check first: both cohorts are GDC
STAR-count TPM, so the medians should agree; the comparison is
rejected if the cohort medians of per-gene mean log2(TPM+1) differ
by >1.0 log2 unit (the SCLC ucologne rejection rule).
"""
import json, glob, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/'
subset = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
nbl_ids = {h['file_id'] for h in subset['nbl']}
pool_ids = {h['file_id'] for h in subset['pooled']}
expr_cols = pd.read_csv('/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', nrows=0).columns[1:]
genes = set(c.split(' (')[0] for c in expr_cols)
print('protein-coding genes in expression matrix:', len(genes))

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

data = {}
n_nbl = n_pool = 0
missing = set()
for d in glob.glob(os.path.join(DL, 'gdc_subset/*/')):
    fid = os.path.basename(d.rstrip('/'))
    if fid in nbl_ids:
        cohort = 'nbl'
    elif fid in pool_ids:
        cohort = 'pool'
    else:
        continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv:
        continue
    if cohort == 'nbl':
        n_nbl += 1
    else:
        n_pool += 1
    df = pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name', 'tpm_unstranded'])
    df = df[df['gene_name'].isin(genes)]
    missing |= (genes - set(df['gene_name']))
    vals = dict(zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)))
    for g, v in vals.items():
        data.setdefault(g, {'nbl': [], 'pool': []})
        data[g][cohort].append(v)
print('tumour files read: NBL', n_nbl, 'pooled', n_pool)
print('genes never seen in any file:', sorted(missing))
rows = []
for g in sorted(genes):
    a = np.array(data.get(g, {}).get('nbl', []), dtype=float)
    b = np.array(data.get(g, {}).get('pool', []), dtype=float)
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    fc = am - bm
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, TCGA_XC_n_NBL=len(a), TCGA_XC_n_other=len(b),
                     TCGA_XC_mean_NBL=am, TCGA_XC_mean_other=bm,
                     TCGA_XC_log2FC=fc, TCGA_XC_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_XC_q'] = bh(xc['TCGA_XC_p'].values)
xc = xc.sort_values('TCGA_XC_log2FC', ascending=False)
# scale check: cohort medians of per-gene mean log2(TPM+1)
med_nbl = float(np.nanmedian(xc['TCGA_XC_mean_NBL']))
med_pool = float(np.nanmedian(xc['TCGA_XC_mean_other']))
print(f'SCALE CHECK: median per-gene mean log2(TPM+1): NBL {med_nbl:.3f} vs pooled {med_pool:.3f}, diff {med_nbl - med_pool:.3f}')
print('scale check:', 'PASS (same GDC pipeline)' if abs(med_nbl - med_pool) <= 1.0 else 'FAIL — reject comparison')
xc.to_csv(OUT + 'tcga_cross_cancer.csv', index=False)
bar = xc[(xc.TCGA_XC_log2FC >= 0.5) & (xc.TCGA_XC_p < 0.05)]
print('cross-cancer bar pass:', len(bar))
print(xc.head(25)[['Gene', 'TCGA_XC_n_NBL', 'TCGA_XC_n_other', 'TCGA_XC_log2FC', 'TCGA_XC_q']].to_string())
CHECK = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
         'MYCN', 'MYC', 'MYCL', 'ALK', 'PHOX2B', 'ATRX', 'KLF5', 'FERMT2', 'EGFR', 'HNF1B',
         'FOSL1', 'SOX9', 'MT1E', 'ITGAV', 'TCF7L2', 'PARD3', 'PARD6B',
         'ARHGAP29', 'ITGB5', 'CDK6', 'TIPARP', 'KRAS', 'HAND2', 'GATA3',
         'ISL1', 'ASCL1', 'BCL2', 'MDM2', 'AURKA', 'CHEK1', 'EZH2', 'LIN28B',
         'TERT', 'SOX11', 'DLK1', 'DBH', 'TH', 'NCAM1', 'CHGA', 'SYP', 'RET',
         'NTRK1', 'NTRK2', 'NGFR', 'HAPSTR1']
for g in CHECK:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_XC_log2FC', 'TCGA_XC_p', 'TCGA_XC_q']].to_string(index=False) if len(r) else 'NOT FOUND')
