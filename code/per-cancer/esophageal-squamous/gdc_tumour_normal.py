"""TCGA ESCC tumour vs adjacent normal from GDC STAR counts.

Design (b), pre-registered before seeing results: ALL 95 ESCA squamous-cell
primary tumours vs ALL 13 solid-tissue normals in GDC. No subsampling — the
study is small, so every squamous tumour is used. 13 normals is fewer than
HNSC's 44 but far more than PDAC's 4: the axis can discriminate, weakly, and
the bar is applied exactly as pre-defined.
"""
import json, glob, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophageal-squamous/downloads'
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
genes = set(entrez)
files = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

data = {}
n_tum = n_norm = 0
missing_genes = set()
for d in glob.glob(os.path.join(DL, 'gdc_subset/*/')):
    fid = os.path.basename(d.rstrip('/'))
    if fid not in fid2type: continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv: continue
    st = fid2type[fid]
    if st == 'Primary Tumor': n_tum += 1
    elif st == 'Solid Tissue Normal': n_norm += 1
    else: continue
    df = pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name', 'tpm_unstranded'])
    df = df[df['gene_name'].isin(genes)]
    got = set(df['gene_name'])
    missing_genes |= (genes - got)
    vals = dict(zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)))
    for g, v in vals.items():
        data.setdefault(g, {'tumour': [], 'normal': []})
        data[g]['tumour' if st == 'Primary Tumor' else 'normal'].append(v)
print('tumour files', n_tum, 'normal files', n_norm)
print('genes never seen in any file:', sorted(missing_genes))
rows = []
for g in sorted(entrez):
    a = np.array(data.get(g, {}).get('tumour', []), dtype=float)
    b = np.array(data.get(g, {}).get('normal', []), dtype=float)
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    fc = am - bm
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=entrez[g], GDC_TvN_n_tumour=len(a), GDC_TvN_n_normal=len(b),
                     GDC_TvN_mean_tumour=am, GDC_TvN_mean_normal=bm, GDC_TvN_log2FC=fc, GDC_TvN_p=p))
tvn = pd.DataFrame(rows)
tvn['GDC_TvN_q'] = bh(tvn['GDC_TvN_p'].values)
tvn = tvn.sort_values('GDC_TvN_log2FC', ascending=False)
tvn.to_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophageal-squamous/tcga_tumour_vs_normal.csv', index=False)
print(tvn.head(25)[['Gene', 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_log2FC', 'GDC_TvN_q']].to_string())
bar = tvn[(tvn.GDC_TvN_log2FC >= 0.5) & (tvn.GDC_TvN_p < 0.05)]
print('GDC tumour-vs-normal bar pass:', len(bar), sorted(bar.Gene.tolist()))
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA']:
    r = tvn[tvn.Gene == g]
    print(g, r[['GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q']].to_string(index=False) if len(r) else 'NOT FOUND')
