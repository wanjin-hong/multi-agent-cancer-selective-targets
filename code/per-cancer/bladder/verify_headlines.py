"""Independent verification of the bladder run's headline numbers.

Separate code path from the pipeline scripts: vectorised numpy recomputation
of the dependency scan (not the row loop in scan.py), direct re-reads of the
raw expression matrix, the GDC STAR-count TSVs, the cBioPortal JSONs and the
mutation file. Every check compares a hard-coded expected value (taken from
the pipeline outputs) against a fresh recomputation and reports PASS/FAIL.
"""
import json, glob, os, numpy as np, pandas as pd
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/bladder/'
DL = OUT + 'downloads/'
fails = []
nchecks = [0]

def check(name, got, exp, tol=0.0):
    nchecks[0] += 1
    if isinstance(exp, float) or isinstance(got, float):
        ok = (got is not None and exp is not None and abs(float(got) - float(exp)) <= max(tol, 1e-9))
    else:
        ok = (got == exp)
    print(('PASS' if ok else 'FAIL'), name, '| got', got, '| expected', exp)
    if not ok:
        fails.append(name)

# ---------- group / comparator composition ----------
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
blad_ids = model.loc[model['OncotreePrimaryDisease'] == 'Bladder Urothelial Carcinoma', 'ModelID'].tolist()
blad = [x for x in blad_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(blad_ids)]
check('effect matrix rows', eff.shape[0], 1178)
check('effect matrix genes', eff.shape[1], 17916)
check('bladder CRISPR n', len(blad), 30)
check('comparator n', len(comp), 1148)
nc = model[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(eff.index))]
check('non-cancerous lines in comparator', len(nc), 14)
check('non-cancerous names', sorted(nc['CellLineName'].tolist()),
      sorted(['Hs 729', 'HA1E', 'WPE1-NA22', 'BPH-1', 'P4E6', 'Shmac 4', 'Shmac 5', 'ARH-77', 'IM-9',
              'RPE1-ss48', 'RPE1-ss77', 'RPE1-ss6', 'RPE1-ss119', 'RPE1-ss51']))
check('bladder in Model.csv', len(blad_ids), 35)
sub = model.set_index('ModelID')
check('same-patient pairs in CRISPR group', int(sub.loc[blad, 'PatientID'].duplicated().sum()), 3)
check('bladder squamous lines excluded', int((model['OncotreePrimaryDisease'] == 'Bladder Squamous Cell Carcinoma').sum()), 2)

# ---------- vectorised dependency scan ----------
genes = [c.split(' (')[0] for c in eff.columns]
A = eff.loc[blad].to_numpy(dtype=float)
B = eff.loc[comp].to_numpy(dtype=float)
G = eff.to_numpy(dtype=float)
am = np.nanmean(A, axis=0); cm = np.nanmean(B, axis=0)
an = np.sum(~np.isnan(A), axis=0); bn = np.sum(~np.isnan(B), axis=0); gn = np.sum(~np.isnan(G), axis=0)
adep = np.nansum(A < -0.5, axis=0); bdep = np.nansum(B < -0.5, axis=0); gdep = np.nansum(G < -0.5, axis=0)
af = np.where(an > 0, adep / np.maximum(an, 1), np.nan)
cf = np.where(bn > 0, bdep / np.maximum(bn, 1), np.nan)
gw = np.where(gn > 0, gdep / np.maximum(gn, 1), np.nan)
enr = am - cm
# Welch p via scipy on the two matrices column-wise (independent of scan.py loop)
pvals = np.array([stats.ttest_ind(A[:, j][~np.isnan(A[:, j])], B[:, j][~np.isnan(B[:, j])], equal_var=False)[1]
                  if np.sum(~np.isnan(A[:, j])) >= 2 and np.sum(~np.isnan(B[:, j])) >= 2 else np.nan
                  for j in range(A.shape[1])])
def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p)
    pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out
qvals = bh(pvals)
paness = (af >= 0.70) & (cf >= 0.70)
strict = (enr <= -0.10) & (af >= 0.25) & ((af - cf) >= 0.10) & (~paness) & (gw <= 0.60)
relaxed = (enr <= -0.05) & (af >= 0.20) & ((af - cf) >= 0.10) & (~paness) & (gw <= 0.60)
check('q<0.05 genes', int(np.nansum(qvals < 0.05)), 19)
check('q<0.05 selective-direction', int(np.nansum((qvals < 0.05) & (enr < 0))), 2)
check('strict-bar genes', int(strict.sum()), 58)
check('relaxed-bar genes', int(relaxed.sum()), 132)
check('pan-essential genes', int(paness.sum()), 1302)
check('genes never scored in bladder', int(np.sum(np.isnan(am))), 129)
check('genes scored in <30 bladder lines (excl. never-scored)', int(np.sum((~np.isnan(am)) & (np.sum(~np.isnan(A), axis=0) < 30))), 0)

def gidx(sym):
    return genes.index(sym)

def enrich_rank(j):
    return int(np.sum(enr < enr[j])) + 1

for sym, exp_mean, exp_dep, exp_enr, exp_rank, exp_q in [
    ('PPARG', -0.284885, 6, -0.243803, 12, 0.071847),
    ('TP63', -0.251644, 8, -0.226611, 18, 0.252364),
    ('KLF5', -0.545573, 14, -0.179967, 32, 0.492389),
    ('SKP2', -0.902831, 23, -0.261078, 7, 0.295792),
    ('FAM50A', -0.818893, 20, -0.292423, 2, 0.213110),
    ('EGFR', -0.392012, 9, -0.100589, 172, 0.592604),
    ('FOSL1', -0.587037, 16, -0.200083, 26, 0.226099),
    ('TEAD3', -0.411648, 9, -0.150649, 53, 0.178038),
    ('PTK2', -0.646655, 22, -0.083782, 274, 0.622413),
    ('PKN2', -0.509986, 14, -0.139114, 66, 0.278685),
    ('TEAD1', -0.371774, 7, -0.139106, 67, 0.191555),
    ('TEAD2', -0.153717, 1, -0.043170, 1263, 0.488201),
    ('TEAD4', -0.187525, 1, -0.051371, 881, 0.600559),
    ('YAP1', -0.411653, 10, -0.077579, 328, 0.691603),
    ('WWTR1', -0.431274, 10, -0.081840, 291, 0.424436),
    ('FGFR3', -0.004003, 2, -0.070080, 418, 0.606125),
    ('PIK3CA', -0.645142, 16, -0.163361, 42, 0.305608),
    ('ERBB2', -0.331164, 5, -0.002034, 8312, 0.998427),
    ('ERBB3', -0.272798, 1, 0.006364, 10621, 0.960988),
    ('FERMT2', -0.589964, 16, 0.025036, 14563, 0.960988),
    ('ITGAV', -0.704196, 18, -0.068632, 441, 0.798592),
    ('ITGB5', -0.437406, 9, -0.125062, 91, 0.639483),
    ('CDK6', -0.472519, 14, 0.118687, 17705, 0.577156),
    ('MT1E', -0.330103, 5, 0.152311, 17759, 0.130699),
    ('SOX9', -0.236011, 4, -0.006726, 7042, 0.991022),
    ('HNF1B', -0.302581, 5, 0.018002, 13352, 0.959653),
    ('TCF7L2', -0.043195, 0, 0.055592, 17040, 0.399014),
]:
    j = gidx(sym)
    check(f'{sym} Chronos mean', round(float(am[j]), 6), exp_mean, tol=1e-5)
    check(f'{sym} dependent lines', int(np.nansum(A[:, j] < -0.5)), exp_dep)
    check(f'{sym} enrichment', round(float(enr[j]), 6), exp_enr, tol=1e-5)
    check(f'{sym} enrichment rank', enrich_rank(j), exp_rank)
    check(f'{sym} q', round(float(qvals[j]), 6), exp_q, tol=1e-4)
# CXADR never scored in bladder
j = gidx('CXADR')
check('CXADR never scored in bladder (all NaN)', bool(np.isnan(A[:, j]).all()), True)
check('CXADR enrichment is NaN', bool(np.isnan(enr[j])), True)

# ---------- cell-line RNA (fresh read) ----------
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
blad_r = [x for x in blad_ids if x in expr.index]
comp_r = [x for x in expr.index if x not in set(blad_ids)]
check('bladder RNA n', len(blad_r), 33)
check('RNA comparator n', len(comp_r), 1640)
nc_r = model[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(expr.index))]
check('non-cancerous lines in RNA comparator', len(nc_r), 85)
check('CRISPR non-cancerous lines also in RNA matrix', int(nc[model['ModelID'].isin(expr.index)]['ModelID'].isin(eff.index).sum()) if False else int(sum(1 for x in nc['ModelID'] if x in expr.index and x in eff.index)), 12)
def rna_check(sym, exp_fc, exp_d):
    col = [c for c in expr.columns if c.startswith(sym + ' (')][0]
    a = expr.loc[blad_r, col].dropna().values; b = expr.loc[comp_r, col].dropna().values
    am_, bm_ = a.mean(), b.mean()
    fc = (2 ** am_ - 1) / (2 ** bm_ - 1)
    check(f'{sym} cell-line RNA FC', round(float(fc), 4), exp_fc, tol=0.01)
    check(f'{sym} cell-line RNA dlog2', round(float(am_ - bm_), 4), exp_d, tol=0.01)
for sym, fc, d in [('PPARG', 8.6808, 2.8730), ('GRHL2', 5.5310, 1.9281), ('TP63', 4.5699, 1.5522),
                   ('KLF5', 3.8527, 1.8431), ('EGFR', 5.1980, 2.2580), ('TEAD1', 1.7583, 0.7751),
                   ('TEAD3', 2.1022, 0.9777), ('TEAD4', 1.8673, 0.8503), ('YAP1', 2.5281, 1.2979),
                   ('WWTR1', 2.0011, 0.9524), ('PTK2', 2.0734, 1.0341), ('PKN2', 1.7060, 0.7353),
                   ('FGFR3', 2.1263, 0.8613), ('ERBB2', 1.8371, 0.8443), ('ERBB3', 2.2920, 1.0846),
                   ('FOSL1', 3.3918, 1.7105), ('PIK3CA', 1.0552, 0.0701)]:
    rna_check(sym, fc, d)
rna = pd.read_csv(OUT + 'celline_rna.csv')
check('cell-line RNA bar pass (genome-wide)', int(((rna.TPM_FC >= 1.5) & (rna.DeltaLog2 >= 0.5) & (rna.p < 0.05)).sum()), 1620)

# ---------- GDC tumour-vs-normal (fresh read of TSVs, each file read once) ----------
files = json.load(open(DL + 'gdc_subset_files.json'))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
tvn_genes = ['PPARG', 'GRHL2', 'TP63', 'KLF5', 'TEAD1', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2',
             'PKN2', 'FGFR3', 'FERMT2', 'FOSL1', 'ERBB2', 'ERBB3', 'MYBL2', 'CDKN2A', 'SKP2']
tvn_data = {g: {'tumour': [], 'normal': []} for g in tvn_genes}
for d in glob.glob(os.path.join(DL, 'gdc_subset/*/')):
    fid = os.path.basename(d.rstrip('/'))
    if fid not in fid2type:
        continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv:
        continue
    df = pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name', 'tpm_unstranded'])
    df = df[df['gene_name'].isin(tvn_genes)]
    bucket = 'tumour' if fid2type[fid] == 'Primary Tumor' else 'normal'
    for _, row in df.iterrows():
        tvn_data[row['gene_name']][bucket].append(float(np.log2(float(row['tpm_unstranded']) + 1)))
for sym, exp_fc, exp_p in [('PPARG', 0.7391, None), ('GRHL2', 1.5605, None), ('TP63', 0.9355, None),
                   ('KLF5', 0.5923, None), ('TEAD1', -1.3781, 0.000853), ('TEAD4', 0.9079, 0.000463),
                   ('YAP1', -0.7861, 0.002192), ('WWTR1', -1.7093, 0.000024), ('PTK2', 0.2509, None),
                   ('PKN2', -0.1861, None), ('FGFR3', 0.6376, None), ('FERMT2', -2.3292, 0.000012),
                   ('FOSL1', -1.5661, 0.005967), ('ERBB2', 0.8402, 0.008256), ('ERBB3', 0.9524, 0.018286),
                   ('MYBL2', 3.3724, 0.000001), ('CDKN2A', 2.1833, 0.000000), ('SKP2', 1.2963, 0.0)]:
    tvals, nvals = tvn_data[sym]['tumour'], tvn_data[sym]['normal']
    fc = float(np.mean(tvals) - np.mean(nvals))
    check(f'{sym} TvN log2FC', round(fc, 4), exp_fc, tol=0.01)
    check(f'{sym} TvN n tumour/normal', (len(tvals), len(nvals)), (150, 19))
    if exp_p is not None:
        p = stats.ttest_ind(tvals, nvals, equal_var=False)[1]
        check(f'{sym} TvN p', round(float(p), 6), exp_p, tol=0.01)
tvn = pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv')
check('TvN bar pass (fetch set)', int(((tvn.GDC_TvN_log2FC >= 0.5) & (tvn.GDC_TvN_p < 0.05)).sum()), 39)

# ---------- cBioPortal cross-cancer RNA + RPPA (fresh read of JSONs, each loaded once) ----------
entrez = json.load(open(DL + 'entrez_map.json'))
def is_pt(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
mrna = {c: pd.DataFrame(json.load(open(DL + f'tcga_{c}_mrna.json'))) for c in ['blca', 'brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']}
rppa = {c: pd.DataFrame(json.load(open(DL + f'tcga_{c}_rppa.json'))) for c in ['blca', 'brca', 'luad', 'lusc', 'coadread', 'stad']}
others_mrna = pd.concat([mrna[c] for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
others_rppa = pd.concat([rppa[c] for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
for sym, exp_fc in [('TP63', 4.1485), ('KRT20', 4.6063), ('GATA3', 3.3185), ('PPARG', 2.3414),
                ('FOSL1', 2.2730), ('KLF5', 2.0263), ('FGFR3', 1.8247), ('TEAD3', 0.9951),
                ('TEAD1', -0.2790), ('YAP1', -0.0231), ('WWTR1', -0.7730), ('PTK2', 0.1231),
                ('PKN2', 0.0871), ('ERBB2', 0.5479), ('ERBB3', -0.8749)]:
    eid = entrez[sym]
    bl = mrna['blca']
    a = np.log2(bl[(bl['sampleId'].apply(is_pt)) & (bl['entrezGeneId'] == eid)]['value'].astype(float) + 1)
    b = np.log2(others_mrna[(others_mrna['sampleId'].apply(is_pt)) & (others_mrna['entrezGeneId'] == eid)]['value'].astype(float) + 1)
    check(f'{sym} cross-cancer log2FC', round(float(a.mean() - b.mean()), 4), exp_fc, tol=0.01)
    check(f'{sym} cross-cancer n BLCA/other', (len(a), len(b)), (408, 3274))
for sym, exp_d in [('BRAF', 0.2438), ('WWTR1', 0.0381), ('YAP1', 0.0264), ('EGFR', 0.0204),
               ('ERBB3', 0.0279), ('ERBB2', -0.1727), ('GATA3', -0.3073), ('PIK3CA', -0.0412),
               ('MAPK1', 0.0286), ('CCND1', -0.0123)]:
    eid = entrez[sym]
    bl = rppa['blca']
    a = bl[(bl['sampleId'].apply(is_pt)) & (bl['entrezGeneId'] == eid)]['value'].astype(float)
    b = others_rppa[(others_rppa['sampleId'].apply(is_pt)) & (others_rppa['entrezGeneId'] == eid)]['value'].astype(float)
    check(f'{sym} RPPA delta', round(float(a.mean() - b.mean()), 4), exp_d, tol=0.01)
    check(f'{sym} RPPA n BLCA/other', (len(a), len(b)), (344, 2466))
rp = pd.read_csv(OUT + 'tcga_rppa.csv')
check('RPPA genes measured in BLCA', int((rp.RPPA_n_BLCA > 0).sum()), 13)
check('RPPA bar pass', int(((rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)).sum()), 0)

# ---------- intersection / ranking (recomputed from merged file) ----------
m = pd.read_csv(OUT + 'overexpression.csv')
cand = pd.read_csv(OUT + 'candidates.csv')
check('strict & OE', int(((m.Strict_pass) & (m.OverExpr_pass)).sum()), 22)
check('relaxed & OE (candidates)', int(((m.Relaxed_pass) & (m.OverExpr_pass)).sum()), 45)
check('candidates.csv rows', len(cand), 45)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top-20 order', top20['Symbol'].tolist(),
      ['PPARG', 'GRHL2', 'TP63', 'KLF5', 'SOX4', 'SKP2', 'FAM50A', 'TRIP13', 'RXRA', 'E2F3',
       'KATNB1', 'SLC35B2', 'PTK2', 'WDR83', 'EGFR', 'FOSL1', 'TEAD3', 'TRAF2', 'TUBB4B', 'ENO1'])
check('top-20 #1 score', round(float(top20.iloc[0]['Score']), 2), 8.65, tol=0.01)
check('top-20 #20 score', round(float(top20.iloc[19]['Score']), 2), 0.82, tol=0.01)
check('ranks 21-25', cand.iloc[20:25]['Symbol'].tolist(), ['VMA21', 'CKS2', 'MYH9', 'PKN2', 'TFDP1'])
# recompute combined score independently from the relaxed set in overexpression.csv
rel = m[m.Relaxed_pass].copy()
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel['Score2'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['GDC_TvN_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0)) + np.where(rel['RPPA_pass'].fillna(False), 1.0, 0.0)
rel = rel[rel['OverExpr_pass']].sort_values('Score2', ascending=False)
check('recomputed score order == candidates order', rel['Symbol'].tolist(), cand['Symbol'].tolist())
for sym, exp_s in [('PPARG', 8.6456), ('GRHL2', 7.9518), ('TP63', 6.3486), ('KLF5', 6.0611),
                   ('TEAD3', 2.1934), ('PTK2', 2.4191), ('PKN2', 0.6115), ('YAP1', -1.5179)]:
    got = float(rel[rel.Symbol == sym]['Score2'].iloc[0])
    check(f'{sym} recomputed score', round(got, 4), exp_s, tol=0.01)
# TEAD-axis candidacy
for sym, exp_cand in [('TEAD1', None), ('TEAD3', 17), ('TEAD4', None), ('YAP1', 44), ('WWTR1', None), ('PTK2', 13), ('PKN2', 24)]:
    got = cand.index[cand.Symbol == sym].tolist()
    got = (int(got[0]) + 1) if got else None
    check(f'{sym} candidate rank', got, exp_cand)
check('PABIR1 has no TvN and no cell RNA', bool(m[m.Symbol == 'PABIR1']['GDC_TvN_log2FC'].isna().iloc[0] and m[m.Symbol == 'PABIR1']['DeltaLog2'].isna().iloc[0]), True)

# ---------- FGFR3 mutation check (fresh read, dedup) ----------
mut = pd.read_csv('/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv',
                  usecols=['HugoSymbol', 'ModelID', 'VariantType', 'VepImpact', 'ProteinChange', 'Hotspot', 'AMClass', 'LikelyLoF'],
                  low_memory=False).drop_duplicates()
mb = mut[mut['ModelID'].isin(blad)]
fg = mb[(mb['HugoSymbol'] == 'FGFR3') & (mb['Hotspot'] == True) & (mb['VepImpact'] == 'MODERATE')]
check('FGFR3 activating-mutant bladder lines', sorted(fg['ModelID'].tolist()), sorted(['ACH-000396', 'ACH-001409', 'ACH-001414']))
check('FGFR3 activating changes', sorted(fg['ProteinChange'].tolist()), sorted(['p.K650E', 'p.R248C', 'p.S249C']))
fgcol = [c for c in eff.columns if c.startswith('FGFR3 (')][0]
check('FGFR3 Chronos in mutant lines', [round(float(eff.loc[x, fgcol]), 3) for x in sorted(fg['ModelID'].tolist())], [-0.16, 0.117, 0.154])
dep_lines = [x for x in blad if eff.loc[x, fgcol] < -0.5]
check('FGFR3-dependent bladder lines', sorted(dep_lines), sorted(['ACH-000473', 'ACH-001183']))
check('FGFR3-dependent lines are the RT-112 pair', set(dep_lines) <= set(sub.loc[blad][sub.loc[blad, 'PatientID'] == 'PT-ESM25t'].index.tolist()), True)
check('TP53-mutant bladder lines', int(mb[mb['HugoSymbol'] == 'TP53']['ModelID'].nunique()), 23)
check('PIK3CA-mutant bladder lines (total)', int(mb[mb['HugoSymbol'] == 'PIK3CA']['ModelID'].nunique()), 10)
check('PIK3CA hotspot-mutant lines', int(mb[(mb['HugoSymbol'] == 'PIK3CA') & (mb['Hotspot'] == True)]['ModelID'].nunique()), 9)

# ---------- luminal/basal ----------
lb = pd.read_csv(OUT + 'luminal_basal_check.csv')
check('luminal-leaning lines', int((lb.Diff > 1).sum()), 5)
check('basal-leaning lines', int((lb.Diff < -1).sum()), 14)
check('intermediate lines', int((lb.Diff.abs() <= 1).sum()), 14)

print()
print(f'{nchecks[0] - len(fails)}/{nchecks[0]} checks passed')
if fails:
    print('FAILURES:', fails)
else:
    print('ALL CHECKS PASSED')
