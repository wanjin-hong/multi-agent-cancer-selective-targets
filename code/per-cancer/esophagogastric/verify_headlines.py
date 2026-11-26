"""Independent verification of ESGA headline numbers, recomputed from the raw
files by a separate code path (no shared functions with scan.py / merge.py).

Checks: group counts; Chronos means / dependent counts / enrichments /
enrichment ranks for KLF5, YAP1, TEAD1-4, WWTR1/TAZ, PTK2, PKN2, ITGAV,
FERMT2, CTNNB1, TCF7L2, PIK3CA, KRAS, EGFR, MET, ERBB2; bar counts
(strict 143, relaxed 241, candidates 161, strict+OE 91); cell-line RNA
fold changes (verified formula); GDC tumour-vs-normal log2FC from the raw
STAR-count TSVs; RPPA deltas from the raw cBioPortal JSON; top-20 order and
combined scores recomputed from the merged components; q<0.05 counts;
site-split dependent counts.
"""
import json, glob, os, numpy as np, pandas as pd
from scipy import stats

OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/'
BASE = '/home/hatch/workspace/research_notes/'
DL = OUT + 'downloads/'
fails, passes = [], []

def check(name, got, want, tol=0.0):
    ok = (got == want) if tol == 0 else (abs(got - want) <= tol)
    (passes if ok else fails).append((name, got, want))
    if not ok:
        print(f'  FAIL {name}: got {got} want {want}')

# ---- raw files, fresh reads ----
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
gids = model.loc[model['OncotreePrimaryDisease'] == 'Esophagogastric Adenocarcinoma', 'ModelID'].tolist()
grp = [x for x in gids if x in eff.index]
cmp_ = [x for x in eff.index if x not in set(gids)]
check('CRISPR group n', len(grp), 45)
check('CRISPR comparator n', len(cmp_), 1133)
check('total CRISPR lines', len(eff), 1178)
grna = [x for x in gids if x in expr.index]
cmprna = [x for x in expr.index if x not in set(gids)]
check('RNA group n', len(grna), 52)
check('RNA comparator n', len(cmprna), 1621)

def col(gene):
    hits = [c for c in eff.columns if c.startswith(gene + ' (')]
    assert len(hits) == 1, (gene, hits)
    return hits[0]

def enrich(gene):
    a = eff.loc[grp, col(gene)].dropna()
    b = eff.loc[cmp_, col(gene)].dropna()
    return a.mean(), (a < -0.5).sum(), len(a), b.mean(), (b < -0.5).sum(), a.mean() - b.mean()

# enrichment ranks recomputed from scratch
means_grp = eff.loc[grp].mean()
means_cmp = eff.loc[cmp_].mean()
enr = (means_grp - means_cmp).sort_values()
rank = {c.split(' (')[0]: int((enr < enr[c]).sum()) + 1 for c in enr.index}

scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan_idx = scan.set_index('Gene')
for gene, (m, d, n, cm, cd, e) in {
    'KLF5': enrich('KLF5'), 'YAP1': enrich('YAP1'), 'TEAD1': enrich('TEAD1'),
    'TEAD2': enrich('TEAD2'), 'TEAD3': enrich('TEAD3'), 'TEAD4': enrich('TEAD4'),
    'WWTR1': enrich('WWTR1'), 'PTK2': enrich('PTK2'), 'PKN2': enrich('PKN2'),
    'ITGAV': enrich('ITGAV'), 'FERMT2': enrich('FERMT2'), 'CTNNB1': enrich('CTNNB1'),
    'TCF7L2': enrich('TCF7L2'), 'PIK3CA': enrich('PIK3CA'), 'KRAS': enrich('KRAS'),
    'EGFR': enrich('EGFR'), 'MET': enrich('MET'), 'ERBB2': enrich('ERBB2'),
    'CCND1': enrich('CCND1'), 'GRB2': enrich('GRB2'), 'RAB10': enrich('RAB10'),
}.items():
    row = scan_idx.loc[gene]
    check(f'{gene} ESGA mean', round(row['ESGA_mean'], 4), round(m, 4))
    check(f'{gene} ESGA dep', int(row['ESGA_dep']), int(d))
    check(f'{gene} enrichment', round(row['Enrichment'], 4), round(e, 4))
    check(f'{gene} enrich rank', int(row['EnrichRank']), rank[gene])
    check(f'{gene} comp mean', round(row['Comp_mean'], 4), round(cm, 4))

# bar counts from the merged table
ov = pd.read_csv(OUT + 'overexpression.csv')
check('strict bar', int(ov.Strict_pass.sum()), 143)
check('relaxed bar', int(ov.Relaxed_pass.sum()), 241)
check('strict+OE', int(((ov.Strict_pass) & (ov.OverExpr_pass)).sum()), 91)
check('relaxed+OE (candidates)', int(((ov.Relaxed_pass) & (ov.OverExpr_pass)).sum()), 161)
check('q<0.05', int((ov.Dep_q < 0.05).sum()), 212)
check('q<0.05 selective', int(((ov.Dep_q < 0.05) & (ov.Enrichment < 0)).sum()), 130)
check('pan-essential', int(ov.PanEssential.sum()), 1309)

# cell-line RNA fold changes (verified formula)
def rnacheck(gene):
    a = expr.loc[grna, col(gene)].dropna()
    b = expr.loc[cmprna, col(gene)].dropna()
    fc = (2 ** a.mean() - 1) / (2 ** b.mean() - 1)
    row = ov[ov.Symbol == gene].iloc[0]
    check(f'{gene} RNA FC', round(row['TPM_FC'], 3), round(fc, 3))
    check(f'{gene} RNA dlog2', round(row['DeltaLog2'], 3), round(a.mean() - b.mean(), 3))
for g in ['KLF5', 'YAP1', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'WWTR1', 'PTK2', 'PKN2', 'CTNNB1', 'TCF7L2', 'HNF1B', 'SOX9', 'CDX2', 'EGFR', 'ERBB2']:
    rnacheck(g)

# GDC TvN from raw TSVs (independent of gdc_tumour_normal.py)
entrez = json.load(open(DL + 'entrez_map.json'))
files = json.load(open(DL + 'gdc_subset_files.json'))
proj_of = {}
for proj, fn in [('TCGA-STAD', 'gdc_tcga-stad_files.json'), ('TCGA-ESCA', 'gdc_tcga-esca_files.json')]:
    for h in json.load(open(DL + fn)):
        proj_of[h['file_id']] = proj
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
vals = {}
for d in glob.glob(DL + '*/'):
    fid = os.path.basename(d.rstrip('/'))
    if fid not in fid2type:
        continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv:
        continue
    st = fid2type[fid]
    if st not in ('Primary Tumor', 'Solid Tissue Normal'):
        continue
    df = pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name', 'tpm_unstranded'])
    for gn, tv in zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)):
        if gn in entrez:
            vals.setdefault(gn, {'tumour': [], 'normal': []})
            vals[gn]['tumour' if st == 'Primary Tumor' else 'normal'].append(tv)
check('GDC tumour files', sum(1 for h in files if fid2type[h['file_id']] == 'Primary Tumor'), 239)
check('GDC normal files', sum(1 for h in files if fid2type[h['file_id']] == 'Solid Tissue Normal'), 49)
for g in ['KLF5', 'YAP1', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'WWTR1', 'PTK2', 'PKN2', 'CTNNB1', 'TCF7L2', 'HNF1B', 'SOX9', 'CDX2', 'EGFR', 'ERBB2', 'FERMT2', 'ITGAV', 'KRAS', 'MET']:
    t = np.array(vals[g]['tumour']); n = np.array(vals[g]['normal'])
    row = ov[ov.Symbol == g].iloc[0]
    check(f'{g} TvN log2FC', round(row['GDC_TvN_log2FC'], 3), round(t.mean() - n.mean(), 3))
    check(f'{g} TvN n_tum', int(row['GDC_TvN_n_tumour']), len(t))
    check(f'{g} TvN n_norm', int(row['GDC_TvN_n_normal']), len(n))

# RPPA deltas from raw cBioPortal JSON
def rppa(gene):
    eid = entrez[gene]
    a = []
    for fn in ['tcga_stad_rppa.json', 'tcga_esca_rppa.json']:
        d = pd.DataFrame(json.load(open(DL + fn)))
        d = d[d['entrezGeneId'] == eid]
        # primary tumours only; ESCA filtered to adenocarcinoma by barcode prefix
        d = d[d['sampleId'].str[13:15] == '01']
        if 'esca' in fn:
            cases = json.load(open(DL + 'gdc_esca_cases.json'))
            pre = {c['submitter_id'][:12] for c in cases
                   if any('adenocarcinoma' in x.get('primary_diagnosis', '').lower() for x in c.get('diagnoses', []))}
            d = d[d['sampleId'].str[:12].isin(pre)]
        a.extend(d['value'].tolist())
    b = []
    for fn in ['tcga_brca_rppa.json', 'tcga_luad_rppa.json', 'tcga_lusc_rppa.json', 'tcga_coadread_rppa.json', 'tcga_hnsc_rppa.json']:
        d = pd.DataFrame(json.load(open(DL + fn)))
        d = d[(d['entrezGeneId'] == eid) & (d['sampleId'].str[13:15] == '01')]
        b.extend(d['value'].tolist())
    return np.mean(a) - np.mean(b), len(a), len(b)
for g in ['RICTOR', 'MET', 'YAP1', 'WWTR1', 'EGFR', 'ERBB2', 'CTNNB1', 'MAPK1']:
    delta, na, nb = rppa(g)
    row = ov[ov.Symbol == g].iloc[0]
    check(f'{g} RPPA delta', round(row['RPPA_delta'], 3), round(delta, 3))
    check(f'{g} RPPA n_ESGA', int(row['RPPA_n_ESGA']), na)

# top-20 order + scores recomputed from merged components
cand = pd.read_csv(OUT + 'candidates.csv')
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 order', top20['Symbol'].tolist(), cand.head(20)['Symbol'].tolist())
rel = ov[ov.Relaxed_pass].copy()
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel['Score2'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['GDC_TvN_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0)) + np.where(rel['RPPA_pass'].fillna(False), 1.0, 0.0)
for _, r in top20.iterrows():
    s2 = rel[rel.Symbol == r['Symbol']].iloc[0]['Score2']
    check(f"score {r['Symbol']}", round(r['Score'], 3), round(s2, 3))
check('top20 tiers', top20['Tier'].tolist(), cand.head(20)['Tier'].tolist())
check('ranks 21-25', cand.iloc[20:25]['Symbol'].tolist(), ['MRPS18B', 'SUZ12', 'SLC7A5', 'DHODH', 'ETNK1'])

# site-split dependent counts for headline genes
sub = model.set_index('ModelID').loc[grp, ['CellLineName', 'OncotreeSubtype']]
gastric = sub[sub['OncotreeSubtype'].str.contains('Stomach', na=False)].index.tolist()
eso = sub[sub['OncotreeSubtype'].str.contains('Esophageal|Gastroesophageal', na=False)].index.tolist()
check('gastric CRISPR n', len(gastric), 35)
check('esoGEJ CRISPR n', len(eso), 10)
for g in ['KLF5', 'YAP1', 'CTNNB1', 'TEAD1', 'SOX9', 'ERBB2']:
    row = ov[ov.Symbol == g].iloc[0]
    check(f'{g} gastric dep', int(row['Gastric_dep']), int((eff.loc[gastric, col(g)] < -0.5).sum()))
    check(f'{g} esoGEJ dep', int(row['EsoGEJ_dep']), int((eff.loc[eso, col(g)] < -0.5).sum()))

print(f'\n{len(passes)} checks passed, {len(fails)} failed')
if fails:
    for name, got, want in fails:
        print('FAILED:', name, 'got', got, 'want', want)
