"""Independent verification of the ESCC headline numbers, separate code path
from the pipeline: recompute from the RAW files (CRISPRGeneEffect.csv,
expression matrix, GDC STAR-count TSVs, cBioPortal RPPA JSON) and compare
against dependency_scan.csv / overexpression.csv / top20.csv.
Prints MISMATCH for any disagreement; exits 0 only if all checks pass.
"""
import json, glob, os, numpy as np, pandas as pd

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/esophageal-squamous/'
DL = OUT + 'downloads/'
ESGA_DL = BASE + 'multi-cancer-selective-targets-20261120/esophagogastric/downloads/'
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
escc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Esophageal Squamous Cell Carcinoma', 'ModelID'].tolist()
es_c = [x for x in escc_ids if x in eff.index]
es_r = [x for x in escc_ids if x in expr.index]
comp_c = [x for x in eff.index if x not in set(escc_ids)]
comp_r = [x for x in expr.index if x not in set(escc_ids)]
check('ESCC CRISPR n', len(es_c), 24)
check('ESCC RNA n', len(es_r), 25)
check('comparator CRISPR n', len(comp_c), 1154)
check('comparator RNA n', len(comp_r), 1648)
check('non-cancerous lines in comparator', int(model[model.ModelID.isin(comp_c)]['OncotreePrimaryDisease'].eq('Non-Cancerous').sum()), 14)
check('ESCC lines in Model.csv', len(escc_ids), 29)
check('same-patient pairs in ESCC group', int(model[model.ModelID.isin(es_c)]['PatientID'].duplicated().sum()), 0)

scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan = scan.set_index('Gene')

def col(gene):
    return [c for c in eff.columns if c.startswith(gene + ' (')][0]

# --- dependency headline numbers, recomputed from the raw matrix ---
for gene, want_mean, want_dep, want_enr in [
        ('TP63', -0.328209, 6, -0.303589),
        ('EGFR', -0.516295, 13, -0.226934),
        ('SOX2', -0.327263, 7, -0.252182),
        ('KLF5', -0.542771, 13, -0.176171),
        ('TEAD1', -0.146728, 1, 0.091343),
        ('YAP1', -0.231586, 3, 0.106636),
        ('WWTR1', -0.343342, 5, 0.008346),
        ('PTK2', -0.435128, 9, 0.132579),
        ('PKN2', -0.347068, 3, 0.027915),
        ('ITGAV', -0.321796, 2, 0.322077),
        ('FERMT2', -0.284545, 1, 0.336676),
        ('ITGB1', -0.393599, 9, -0.188209),
        ('CCND1', -1.378841, 19, -0.405181)]:
    c = col(gene)
    a = eff.loc[es_c, c].dropna()
    b = eff.loc[comp_c, c].dropna()
    check(f'{gene} ESCC mean', round(a.mean(), 6), want_mean, 1e-4)
    check(f'{gene} ESCC dep', int((a < -0.5).sum()), want_dep)
    check(f'{gene} enrichment', round(a.mean() - b.mean(), 6), want_enr, 1e-4)

# --- enrichment ranks, recomputed ---
enr = scan['Enrichment']
for gene, want_rank in [('TP63', 23), ('EGFR', 76), ('SOX2', 56), ('KLF5', 162), ('TEAD1', 17551),
                        ('YAP1', 17680), ('WWTR1', 11115), ('PTK2', 17784), ('PKN2', 14369),
                        ('ITGAV', 17912), ('FERMT2', 17913), ('FOSL1', 12536), ('ITGB1', 136),
                        ('CCND1', 6), ('GRB2', 37), ('ITGB5', 17893), ('CDK6', 11511)]:
    check(f'{gene} enrichment rank', int((enr < enr[gene]).sum()) + 1, want_rank)

# --- bar counts, recomputed from the scan file ---
check('strict pass', int(((scan.Enrichment <= -0.10) & (scan.ESCC_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 153)
check('relaxed pass', int(((scan.Enrichment <= -0.05) & (scan.ESCC_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 230)
check('q<0.05 count', int((scan.q_BH < 0.05).sum()), 570)
check('q<0.05 selective-direction', int(((scan.q_BH < 0.05) & (scan.Enrichment < 0)).sum()), 249)
check('pan-essential count', int(scan.PanEssential.sum()), 1309)

# --- cell-line RNA, recomputed from the raw expression matrix ---
def expr_col(gene):
    return [c for c in expr.columns if c.startswith(gene + ' (')][0]

for gene, want_fc, want_d in [('TP63', 38.800902, 4.394909), ('EGFR', 8.722961, 2.992724),
                              ('SOX2', 7.905075, 2.328654), ('YAP1', 2.281415, 1.152923),
                              ('TEAD1', 0.969948, -0.041286), ('KLF5', 6.512168, 2.584798),
                              ('WWTR1', 1.356475, 0.41489), ('PTK2', 1.693555, 0.74604)]:
    c = expr_col(gene)
    a = expr.loc[es_r, c].dropna()
    b = expr.loc[comp_r, c].dropna()
    am, bm = a.mean(), b.mean()
    fc = (2 ** am - 1) / (2 ** bm - 1)
    check(f'{gene} cell-line RNA FC', round(fc, 4), want_fc, 1e-3)
    check(f'{gene} cell-line RNA DeltaLog2', round(am - bm, 4), want_d, 1e-3)

# --- GDC tumour-vs-normal, recomputed from the raw TSVs ---
files = json.load(open(DL + 'gdc_subset_files.json'))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
want_genes = ['TP63', 'EGFR', 'SOX2', 'KLF5', 'TEAD1', 'TEAD2', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'ITGAV', 'FOSL1']
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
check('GDC tumour files', n_tum, 95)
check('GDC normal files', n_norm, 13)
tvn = pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv').set_index('Gene')
for g in want_genes:
    fc = np.mean(data[g]['tumour']) - np.mean(data[g]['normal'])
    check(f'{g} GDC TvN log2FC', round(fc, 4), round(tvn.loc[g, 'GDC_TvN_log2FC'], 4), 1e-3)

# --- RPPA, recomputed from the raw JSON (squamous-filtered ESCA side) ---
emap = json.load(open(DL + 'entrez_map.json'))
cases = json.load(open(os.path.join(ESGA_DL, 'gdc_esca_cases.json')))
sq_sub = set()
for c in cases:
    diags = [d.get('primary_diagnosis', '') for d in c.get('diagnoses', [])]
    if any('squamous' in d.lower() for d in diags):
        sq_sub.add(c['submitter_id'])
def is_primary(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
def is_sq(sid):
    return isinstance(sid, str) and sid[:12] in sq_sub
esca_r = pd.DataFrame(json.load(open(DL + 'tcga_esca_rppa.json')))
other_r = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_rppa.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
esca_r = esca_r[esca_r['sampleId'].apply(lambda s: is_primary(s) and is_sq(s))]
other_r = other_r[other_r['sampleId'].apply(is_primary)]
check('RPPA squamous ESCC samples', int(esca_r['sampleId'].nunique()), 78)
rp = pd.read_csv(OUT + 'tcga_rppa.csv').set_index('Gene')
for g in ['EGFR', 'WWTR1', 'ERBB2', 'CCND1', 'YAP1']:
    eid = emap[g]
    a = esca_r.loc[esca_r['entrezGeneId'] == eid, 'value'].astype(float).values
    b = other_r.loc[other_r['entrezGeneId'] == eid, 'value'].astype(float).values
    check(f'{g} RPPA delta', round(a.mean() - b.mean(), 4), round(rp.loc[g, 'RPPA_delta'], 4), 1e-3)

# --- cross-cancer RNA, recomputed from the raw JSON (squamous-filtered ESCA side) ---
esca_m = pd.DataFrame(json.load(open(DL + 'tcga_esca_mrna.json')))
other_m = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_mrna.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
esca_m = esca_m[esca_m['sampleId'].apply(lambda s: is_primary(s) and is_sq(s))].copy()
other_m = other_m[other_m['sampleId'].apply(is_primary)].copy()
esca_m['value'] = np.log2(esca_m['value'].astype(float) + 1)
other_m['value'] = np.log2(other_m['value'].astype(float) + 1)
xc = pd.read_csv(OUT + 'tcga_rna.csv').set_index('Gene')
for g in ['TP63', 'EGFR', 'SOX2', 'KLF5']:
    eid = emap[g]
    a = esca_m.loc[esca_m['entrezGeneId'] == eid, 'value'].values
    b = other_m.loc[other_m['entrezGeneId'] == eid, 'value'].values
    check(f'{g} cross-cancer RNA log2FC', round(a.mean() - b.mean(), 4), round(xc.loc[g, 'TCGA_RNA_log2FC'], 4), 1e-3)

# --- intersection counts and top 20 ---
oe = pd.read_csv(OUT + 'overexpression.csv')
check('strict & overexpr', int(((oe.Strict_pass) & (oe.OverExpr_pass)).sum()), 89)
check('relaxed & overexpr', int(((oe.Relaxed_pass) & (oe.OverExpr_pass)).sum()), 131)
check('RPPA bar pass (any gene)', int(oe.RPPA_pass.sum()), 0)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 #1', top20.iloc[0]['Symbol'], 'TP63')
check('top20 #1 score', round(top20.iloc[0]['Score'], 2), 17.31, 1e-2)
check('top20 #2', top20.iloc[1]['Symbol'], 'EGFR')
check('top20 #2 score', round(top20.iloc[1]['Score'], 2), 9.85, 1e-2)
check('top20 #5', top20.iloc[4]['Symbol'], 'KLF5')
check('top20 #20', top20.iloc[19]['Symbol'], 'ERBB2')
check('top20 #20 score', round(top20.iloc[19]['Score'], 2), 3.14, 1e-2)
check('top20 Tier A count', int((top20['Tier'] == 'A').sum()), 20)

print()
if fails:
    print('FAILURES:', fails)
    raise SystemExit(1)
print('ALL CHECKS PASSED')
