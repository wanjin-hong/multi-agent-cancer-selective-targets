"""Independent verification of the endometrial headline numbers, separate code path
from the pipeline: recompute from the RAW files (CRISPRGeneEffect.csv,
expression matrix, GDC STAR-count tarball subset, cBioPortal mRNA/RPPA/CPTAC
JSON) and compare against dependency_scan.csv / overexpression.csv / top20.csv.
Prints MISMATCH for any disagreement; exits 0 only if all checks pass.
"""
import json, os, glob, numpy as np, pandas as pd

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/endometrial/'
DL = OUT + 'downloads/'
fails = []

def check(name, got, want, tol=1e-6):
    ok = (got == want) if isinstance(want, (int, str, bool, list)) else (abs(got - want) <= tol * max(1, abs(want)))
    print(('OK   ' if ok else 'MISMATCH'), name, 'got', got, 'want', want)
    if not ok:
        fails.append(name)

# --- raw group counts ---
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
en_ids = model.loc[model['OncotreePrimaryDisease'] == 'Endometrial Carcinoma', 'ModelID'].tolist()
en_c = [x for x in en_ids if x in eff.index]
en_r = [x for x in en_ids if x in expr.index]
comp_c = [x for x in eff.index if x not in set(en_ids)]
comp_r = [x for x in expr.index if x not in set(en_ids)]
check('ENDO CRISPR n', len(en_c), 28)
check('ENDO RNA n', len(en_r), 36)
check('comparator CRISPR n', len(comp_c), 1150)
check('comparator RNA n', len(comp_r), 1637)
nc_c = model.loc[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(eff.index)), 'ModelID'].tolist()
nc_r = model.loc[(model['OncotreePrimaryDisease'] == 'Non-Cancerous') & (model['ModelID'].isin(expr.index)), 'ModelID'].tolist()
check('non-cancerous in effect matrix', len(nc_c), 14)
check('non-cancerous in expression matrix', len(nc_r), 85)
check('non-cancerous in CRISPR comparator', len(set(nc_c) & set(comp_c)), 14)
check('non-cancerous in RNA comparator', len(set(nc_r) & set(comp_r)), 85)
# same-patient pairs inside the CRISPR panel
inc = model[model.ModelID.isin(en_c)]
dups = inc[inc.duplicated('PatientID', keep=False)].sort_values('PatientID')
check('same-patient pairs in CRISPR panel', len(dups), 4)
check('same-patient pair names', sorted(dups['CellLineName'].tolist()), ['HEC-1', 'HEC-1-B', 'HHUA', 'SNG-M'])

scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan = scan.set_index('Gene')

def col(gene):
    return [c for c in eff.columns if c.startswith(gene + ' (')][0]

# --- dependency headline numbers, recomputed from the raw matrix ---
for gene, want_mean, want_dep, want_enr in [
        ('PAX8', -0.740523, 18, -0.548546),
        ('HNF1B', -0.652167, 13, -0.340127),
        ('CCNE1', -0.513020, 10, -0.323059),
        ('WRN', -0.543166, 11, -0.445605),
        ('WWTR1', -0.330387, 5, 0.021645),
        ('YAP1', -0.309167, 5, 0.027537),
        ('TEAD1', -0.197319, 2, 0.039838),
        ('PTK2', -0.604973, 16, -0.040941),
        ('PKN2', -0.427426, 12, -0.054303),
        ('FERMT2', -0.779492, 19, -0.169151),
        ('PARD6B', -0.398910, 9, -0.234061),
        ('ARHGAP29', -0.457274, 12, -0.183501)]:
    c = col(gene)
    a = eff.loc[en_c, c].dropna()
    b = eff.loc[comp_c, c].dropna()
    check(f'{gene} ENDO mean', round(a.mean(), 6), want_mean, 1e-4)
    check(f'{gene} ENDO dep', int((a < -0.5).sum()), want_dep)
    check(f'{gene} enrichment', round(a.mean() - b.mean(), 6), want_enr, 1e-3)

# --- enrichment ranks, recomputed ---
enr = scan['Enrichment']
for gene, want_rank in [('PAX8', 1), ('WRN', 2), ('HNF1B', 4), ('CCNE1', 7), ('PARD6B', 20),
                        ('ARHGAP29', 52), ('FERMT2', 67), ('KLF5', 16128), ('FOSL1', 17572),
                        ('TEAD1', 15576), ('YAP1', 14160), ('WWTR1', 13260), ('PTK2', 2552),
                        ('PKN2', 1518), ('EGFR', 5516), ('CDK6', 17915), ('CXADR', 284)]:
    check(f'{gene} enrichment rank', int((enr < enr[gene]).sum()) + 1, want_rank)

# --- bar counts, recomputed from the scan file ---
check('strict pass', int(((scan.Enrichment <= -0.10) & (scan.ENDO_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 54)
check('relaxed pass', int(((scan.Enrichment <= -0.05) & (scan.ENDO_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 164)
check('q<0.05 count', int((scan.q_BH < 0.05).sum()), 209)
check('CXADR scored lines', int(scan.loc['CXADR', 'ENDO_n']), 7)

# --- cell-line RNA, recomputed from the raw expression matrix ---
def expr_col(gene):
    return [c for c in expr.columns if c.startswith(gene + ' (')][0]

for gene, want_fc, want_d in [('PAX8', 21.140828, 3.894049), ('HNF1B', 9.623616, 2.200730),
                              ('CCNE1', 2.003436, 0.939410), ('WWTR1', 2.3554, 1.1799),
                              ('YAP1', 1.69, 0.73), ('TEAD1', 1.63, 0.67),
                              ('FERMT2', 2.0337, 0.9725), ('SOX17', 22.133927, 2.15579)]:
    c = expr_col(gene)
    a = expr.loc[en_r, c].dropna()
    b = expr.loc[comp_r, c].dropna()
    am, bm = a.mean(), b.mean()
    fc = (2 ** am - 1) / (2 ** bm - 1)
    check(f'{gene} cell-line RNA FC', round(fc, 4), want_fc, 1e-3)
    check(f'{gene} cell-line RNA DeltaLog2', round(am - bm, 4), want_d, 1e-3)

# --- GDC tumour-vs-normal, recomputed from the raw STAR-count files ---
emap = json.load(open(DL + 'entrez_map.json'))
files = json.load(open(DL + 'gdc_subset_files.json'))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
genes = ['PAX8', 'CCNE1', 'HNF1B', 'YAP1', 'WWTR1', 'TEAD1', 'FERMT2', 'ITGAV', 'PKN2', 'TEAD2', 'TEAD4']
data = {}
n_tum = n_norm = 0
for d in glob.glob(os.path.join(DL, 'gdc_subset/*/')):
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
    df = df[df['gene_name'].isin(genes)]
    for g, v in zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)):
        data.setdefault(g, {'tumour': [], 'normal': []})
        data[g]['tumour' if st == 'Primary Tumor' else 'normal'].append(v)
check('GDC tumour files', n_tum, 150)
check('GDC normal files', n_norm, 35)
tvn = pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv').set_index('Gene')
for g in genes:
    a = np.array(data[g]['tumour']); b = np.array(data[g]['normal'])
    check(f'{g} GDC TvN log2FC', round(a.mean() - b.mean(), 4), round(tvn.loc[g, 'GDC_TvN_log2FC'], 4), 1e-3)

# --- cross-cancer TCGA RNA, recomputed from the raw JSON ---
def is_primary(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
uc_m = pd.DataFrame(json.load(open(DL + 'tcga_ucec_mrna.json')))
other_m = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_mrna.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
uc_m = uc_m[uc_m['sampleId'].apply(is_primary)].copy()
other_m = other_m[other_m['sampleId'].apply(is_primary)].copy()
uc_m['value'] = np.log2(uc_m['value'].astype(float) + 1)
other_m['value'] = np.log2(other_m['value'].astype(float) + 1)
check('UCEC mRNA primary tumours', int(uc_m['sampleId'].nunique()), 176)
xc = pd.read_csv(OUT + 'tcga_rna.csv').set_index('Gene')
for g in ['PAX8', 'SOX17', 'WWTR1', 'TEAD1', 'YAP1', 'CCNE1', 'HNF1B', 'WT1']:
    eid = emap[g]
    a = uc_m.loc[uc_m['entrezGeneId'] == eid, 'value'].values
    b = other_m.loc[other_m['entrezGeneId'] == eid, 'value'].values
    check(f'{g} TCGA XC log2FC', round(a.mean() - b.mean(), 4), round(xc.loc[g, 'TCGA_RNA_log2FC'], 4), 1e-3)

# --- RPPA, recomputed from the raw JSON ---
uc_rp = pd.DataFrame(json.load(open(DL + 'tcga_ucec_rppa.json')))
other_rp = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_rppa.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
uc_rp = uc_rp[uc_rp['sampleId'].apply(is_primary)]
other_rp = other_rp[other_rp['sampleId'].apply(is_primary)]
rp = pd.read_csv(OUT + 'tcga_rppa.csv').set_index('Gene')
for g in ['CCNE1', 'WWTR1', 'YAP1', 'EGFR', 'BRCA2', 'ERBB2', 'PXN']:
    eid = emap[g]
    a = uc_rp.loc[uc_rp['entrezGeneId'] == eid, 'value'].astype(float).values
    b = other_rp.loc[other_rp['entrezGeneId'] == eid, 'value'].astype(float).values
    check(f'{g} RPPA delta', round(a.mean() - b.mean(), 4), round(rp.loc[g, 'RPPA_delta'], 4), 1e-3)
check('RPPA bar pass count', int(((rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)).sum()), 0)
check('RPPA measured genes', int((rp.RPPA_n_UCEC > 0).sum()), 24)

# --- CPTAC medians, recomputed from the raw JSON ---
cp = pd.DataFrame(json.load(open(DL + 'tcga_ucec_cptac.json')))
cpt = pd.read_csv(OUT + 'tcga_cptac.csv').set_index('Gene')
check('CPTAC samples', int(cp['sampleId'].nunique()), 95)
for g in ['PAX8', 'WWTR1', 'TEAD1', 'CCNE1', 'PDZK1IP1', 'MECOM']:
    eid = emap[g]
    a = cp.loc[cp['entrezGeneId'] == eid, 'value'].astype(float).values
    check(f'{g} CPTAC median', round(float(np.median(a)), 4), round(cpt.loc[g, 'CPTAC_median'], 4), 1e-3)

# --- intersection counts ---
oe = pd.read_csv(OUT + 'overexpression.csv')
check('strict & overexpr', int(((oe.Strict_pass) & (oe.OverExpr_pass)).sum()), 22)
check('relaxed & overexpr', int(((oe.Relaxed_pass) & (oe.OverExpr_pass)).sum()), 41)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 #1', top20.iloc[0]['Symbol'], 'PAX8')
check('top20 #1 score', round(top20.iloc[0]['Score'], 2), 23.16, 1e-2)
check('top20 #2', top20.iloc[1]['Symbol'], 'HNF1B')
check('top20 #3', top20.iloc[2]['Symbol'], 'CCNE1')
check('top20 #20', top20.iloc[19]['Symbol'], 'UBE2S')
check('top20 #20 score', round(top20.iloc[19]['Score'], 2), 1.90, 1e-2)

print()
if fails:
    print('FAILURES:', fails)
    raise SystemExit(1)
print('ALL CHECKS PASSED')
