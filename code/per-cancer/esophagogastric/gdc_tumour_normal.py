"""TCGA ESGA tumour vs adjacent normal from GDC STAR counts.

Pooled design (pre-registered): 150 random STAD primary tumours (seed 20261121)
+ ALL 89 ESCA adenocarcinoma primary tumours vs ALL 49 adjacent normals
(36 STAD + 13 ESCA). ESCA squamous tumours excluded at the file-selection stage
via the GDC cases API (primary_diagnosis).
Also writes per-study TvN (STAD-only, ESCA-adeno-only) as a sensitivity check,
decided before seeing results.
"""
import json, glob, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/'
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
genes = set(entrez)
files = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
# project per file: rebuild from per-project dumps
proj_of = {}
for proj, fn in [('TCGA-STAD', 'gdc_tcga-stad_files.json'), ('TCGA-ESCA', 'gdc_tcga-esca_files.json')]:
    for h in json.load(open(os.path.join(DL, fn))):
        proj_of[h['file_id']] = proj
fid2type = {h['file_id']: (h['cases'][0]['samples'][0]['sample_type'], proj_of[h['file_id']]) for h in files}

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

data = {}   # gene -> {'tumour': [...], 'normal': [...], 'tum_stad': [...], 'norm_stad': [...], 'tum_esca': [...], 'norm_esca': [...]}
n_tum = n_norm = n_tum_stad = n_tum_esca = n_norm_stad = n_norm_esca = 0
missing_genes = set()
for d in glob.glob(os.path.join(DL, '*/')):  # GDC /data tarball extracts flat UUID dirs into DL
    fid = os.path.basename(d.rstrip('/'))
    if fid not in fid2type: continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv: continue
    st, proj = fid2type[fid]
    if st == 'Primary Tumor':
        n_tum += 1
        if proj == 'TCGA-STAD': n_tum_stad += 1
        else: n_tum_esca += 1
    elif st == 'Solid Tissue Normal':
        n_norm += 1
        if proj == 'TCGA-STAD': n_norm_stad += 1
        else: n_norm_esca += 1
    else:
        continue
    df = pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name', 'tpm_unstranded'])
    df = df[df['gene_name'].isin(genes)]
    got = set(df['gene_name'])
    missing_genes |= (genes - got)
    vals = dict(zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)))
    for g, v in vals.items():
        e = data.setdefault(g, {'tumour': [], 'normal': [], 'tum_stad': [], 'norm_stad': [], 'tum_esca': [], 'norm_esca': []})
        if st == 'Primary Tumor':
            e['tumour'].append(v)
            e['tum_stad' if proj == 'TCGA-STAD' else 'tum_esca'].append(v)
        else:
            e['normal'].append(v)
            e['norm_stad' if proj == 'TCGA-STAD' else 'norm_esca'].append(v)
print('tumour files', n_tum, '(STAD', n_tum_stad, 'ESCA-adeno', n_tum_esca, ')',
      'normal files', n_norm, '(STAD', n_norm_stad, 'ESCA', n_norm_esca, ')')
print('genes never seen in any file:', sorted(missing_genes))

def welch(a, b):
    a = np.array(a, dtype=float); b = np.array(b, dtype=float)
    if len(a) < 2 or len(b) < 2: return np.nan, np.nan, np.nan
    am, bm = a.mean(), b.mean()
    return am, bm, stats.ttest_ind(a, b, equal_var=False)[1]

rows = []
for g in sorted(entrez):
    e = data.get(g, {})
    t, n = e.get('tumour', []), e.get('normal', [])
    am, bm, p = welch(t, n)
    ts, ns, ps = welch(e.get('tum_stad', []), e.get('norm_stad', []))
    te, ne, pe = welch(e.get('tum_esca', []), e.get('norm_esca', []))
    rows.append(dict(Gene=g, Entrez=entrez[g],
                     GDC_TvN_n_tumour=len(t), GDC_TvN_n_normal=len(n),
                     GDC_TvN_mean_tumour=am, GDC_TvN_mean_normal=bm,
                     GDC_TvN_log2FC=(am - bm) if am is not None and not np.isnan(am) and not np.isnan(bm) else np.nan,
                     GDC_TvN_p=p,
                     GDC_TvN_STAD_log2FC=(ts - ns) if ts is not None and not np.isnan(ts) and not np.isnan(ns) else np.nan,
                     GDC_TvN_STAD_p=ps,
                     GDC_TvN_ESCA_log2FC=(te - ne) if te is not None and not np.isnan(te) and not np.isnan(ne) else np.nan,
                     GDC_TvN_ESCA_p=pe))
tvn = pd.DataFrame(rows)
tvn['GDC_TvN_q'] = bh(tvn['GDC_TvN_p'].values)
tvn = tvn.sort_values('GDC_TvN_log2FC', ascending=False)
tvn.to_csv(OUT + 'tcga_tumour_vs_normal.csv', index=False)
print(tvn.head(25)[['Gene', 'GDC_TvN_n_tumour', 'GDC_TvN_n_normal', 'GDC_TvN_log2FC', 'GDC_TvN_q']].to_string())
bar = tvn[(tvn.GDC_TvN_log2FC >= 0.5) & (tvn.GDC_TvN_p < 0.05)]
print('GDC pooled tumour-vs-normal bar pass:', len(bar), sorted(bar.Gene.tolist()))
print('STAD-only bar pass:', int(((tvn.GDC_TvN_STAD_log2FC >= 0.5) & (tvn.GDC_TvN_STAD_p < 0.05)).sum()))
print('ESCA-only bar pass:', int(((tvn.GDC_TvN_ESCA_log2FC >= 0.5) & (tvn.GDC_TvN_ESCA_p < 0.05)).sum()))
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2', 'EGFR', 'PIK3CA', 'KLF5', 'SOX9', 'HNF1B', 'FOSL1', 'CDX2', 'ERBB2', 'CDH17', 'GATA6', 'CTNNB1', 'TCF7L2', 'GFPT1']:
    r = tvn[tvn.Gene == g]
    print(g, r[['GDC_TvN_log2FC', 'GDC_TvN_p', 'GDC_TvN_q', 'GDC_TvN_STAD_log2FC', 'GDC_TvN_ESCA_log2FC']].to_string(index=False) if len(r) else 'NOT FOUND')
