"""Independent verification of the HNSC headline numbers, separate code path
from the pipeline: recompute from the RAW files (CRISPRGeneEffect.csv,
expression matrix, GDC STAR-count TSVs, cBioPortal RPPA JSON) and compare
against dependency_scan.csv / overexpression.csv / top20.csv.
Prints MISMATCH for any disagreement; exits 0 only if all checks pass.
"""
import json, glob, os, numpy as np, pandas as pd

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/head-and-neck/'
DL = OUT + 'downloads/'
fails = []

def check(name, got, want, tol=1e-6):
    ok = (got == want) if isinstance(want, (int, str, bool)) else (abs(got - want) <= tol * max(1, abs(want)))
    print(('OK   ' if ok else 'MISMATCH'), name, 'got', got, 'want', want)
    if not ok:
        fails.append(name)

# --- raw group counts ---
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
hnsc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Head and Neck Squamous Cell Carcinoma', 'ModelID'].tolist()
hn_c = [x for x in hnsc_ids if x in eff.index]
hn_r = [x for x in hnsc_ids if x in expr.index]
comp_c = [x for x in eff.index if x not in set(hnsc_ids)]
check('HNSC CRISPR n', len(hn_c), 72)
check('HNSC RNA n', len(hn_r), 71)
check('comparator CRISPR n', len(comp_c), 1106)

scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan = scan.set_index('Gene')

def col(gene):
    return [c for c in eff.columns if c.startswith(gene + ' (')][0]

# --- dependency headline numbers, recomputed from the raw matrix ---
for gene, want_mean, want_dep, want_enr in [
        ('TP63', -0.499495, 40, -0.499202),
        ('EGFR', -0.704061, 49, -0.436773),
        ('TEAD1', -0.327667, 20, -0.097411),
        ('YAP1', -0.476607, 29, -0.149708),
        ('WWTR1', -0.399213, 21, -0.0508),
        ('PTK2', -0.421704, 25, 0.152631),
        ('PKN2', -0.488637, 35, -0.121658)]:
    c = col(gene)
    a = eff.loc[hn_c, c].dropna()
    b = eff.loc[comp_c, c].dropna()
    check(f'{gene} HNSC mean', round(a.mean(), 6), want_mean, 1e-4)
    check(f'{gene} HNSC dep', int((a < -0.5).sum()), want_dep)
    check(f'{gene} enrichment', round(a.mean() - b.mean(), 6), want_enr, 1e-4)

# --- enrichment ranks, recomputed (dependency_scan.csv is in TNBC column style:
# Gene = symbol, Column = 'SYMBOL (Entrez)') ---
enr = scan['Enrichment']
for gene, want_rank in [('TP63', 1), ('EGFR', 3), ('TEAD1', 400), ('YAP1', 160), ('WWTR1', 1279), ('PTK2', 17845), ('PKN2', 241), ('ITGAV', 17914), ('FERMT2', 17913)]:
    check(f'{gene} enrichment rank', int((enr < enr[gene]).sum()) + 1, want_rank)

# --- bar counts, recomputed from the scan file ---
check('strict pass', int(((scan.Enrichment <= -0.10) & (scan.HNSC_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 125)
check('relaxed pass', int(((scan.Enrichment <= -0.05) & (scan.HNSC_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 203)
check('q<0.05 count', int((scan.q_BH < 0.05).sum()), 4037)

# --- cell-line RNA, recomputed from the raw expression matrix ---
def expr_col(gene):
    return [c for c in expr.columns if c.startswith(gene + ' (')][0]

for gene, want_fc, want_d in [('TP63', 36.846666, 4.219137), ('EGFR', 6.59384, 2.588687), ('YAP1', 1.965094, 0.941775), ('TEAD1', 0.7612, None)]:
    c = expr_col(gene)
    a = expr.loc[hn_r, c].dropna()
    b = expr.loc[[x for x in expr.index if x not in set(hnsc_ids)], c].dropna()
    am, bm = a.mean(), b.mean()
    fc = (2 ** am - 1) / (2 ** bm - 1)
    check(f'{gene} cell-line RNA FC', round(fc, 4), want_fc, 1e-3)
    if want_d is not None:
        check(f'{gene} cell-line RNA DeltaLog2', round(am - bm, 4), want_d, 1e-3)

# --- GDC tumour-vs-normal, recomputed from the raw TSVs ---
files = json.load(open(DL + 'gdc_subset_files.json'))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
want_genes = ['TP63', 'EGFR', 'TEAD1', 'TEAD2', 'YAP1', 'WWTR1', 'PTK2', 'ITGAV']
data = {g: {'tumour': [], 'normal': []} for g in want_genes}
n_tum = n_norm = 0
for d in glob.glob(DL + 'gdc_subset/*/'):
    fid = os.path.basename(d.rstrip('/'))
    if fid not in fid2type:
        continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv:
        continue
    st = fid2type[fid]
    if st == 'Primary Tumor':
        n_tum += 1
    elif st == 'Solid Tissue Normal':
        n_norm += 1
    else:
        continue
    df = pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name', 'tpm_unstranded'])
    df = df[df['gene_name'].isin(want_genes)]
    for g, v in zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)):
        data[g]['tumour' if st == 'Primary Tumor' else 'normal'].append(v)
check('GDC tumour files', n_tum, 150)
check('GDC normal files', n_norm, 44)
tvn = pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv').set_index('Gene')
for g in want_genes:
    fc = np.mean(data[g]['tumour']) - np.mean(data[g]['normal'])
    check(f'{g} GDC TvN log2FC', round(fc, 4), round(tvn.loc[g, 'GDC_TvN_log2FC'], 4), 1e-3)

# --- RPPA, recomputed from the raw JSON ---
emap = json.load(open(DL + 'entrez_map.json'))
inv = {v: k for k, v in emap.items()}
hnsc_r = pd.DataFrame(json.load(open(DL + 'tcga_hnsc_rppa.json')))
other_r = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_rppa.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
def is_primary(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
hnsc_r = hnsc_r[hnsc_r['sampleId'].apply(is_primary)]
other_r = other_r[other_r['sampleId'].apply(is_primary)]
rp = pd.read_csv(OUT + 'tcga_rppa.csv').set_index('Gene')
for g in ['TP63', 'ERBB2', 'WWTR1', 'EGFR', 'CCND1']:
    eid = emap[g]
    a = hnsc_r.loc[hnsc_r['entrezGeneId'] == eid, 'value'].astype(float).values
    b = other_r.loc[other_r['entrezGeneId'] == eid, 'value'].astype(float).values
    check(f'{g} RPPA delta', round(a.mean() - b.mean(), 4), round(rp.loc[g, 'RPPA_delta'], 4), 1e-3)

# --- intersection counts ---
oe = pd.read_csv(OUT + 'overexpression.csv')
check('strict & overexpr', int(((oe.Strict_pass) & (oe.OverExpr_pass)).sum()), 70)
check('relaxed & overexpr', int(((oe.Relaxed_pass) & (oe.OverExpr_pass)).sum()), 109)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 #1', top20.iloc[0]['Symbol'], 'TP63')
check('top20 #1 score', round(top20.iloc[0]['Score'], 2), 19.53, 1e-2)
check('top20 #2', top20.iloc[1]['Symbol'], 'EGFR')
check('top20 #20', top20.iloc[19]['Symbol'], 'CAP1')

print()
if fails:
    print('FAILURES:', fails)
    raise SystemExit(1)
print('ALL CHECKS PASSED')
