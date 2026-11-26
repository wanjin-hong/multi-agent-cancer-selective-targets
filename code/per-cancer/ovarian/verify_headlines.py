"""Independent verification of the ovarian headline numbers, separate code path
from the pipeline: recompute from the RAW files (CRISPRGeneEffect.csv,
expression matrix, cBioPortal mRNA/RPPA JSON) and compare against
dependency_scan.csv / overexpression.csv / top20.csv.
Prints MISMATCH for any disagreement; exits 0 only if all checks pass.
"""
import json, os, numpy as np, pandas as pd

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/ovarian/'
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
ov_ids = model.loc[model['OncotreePrimaryDisease'] == 'Ovarian Epithelial Tumor', 'ModelID'].tolist()
ov_c = [x for x in ov_ids if x in eff.index]
ov_r = [x for x in ov_ids if x in expr.index]
comp_c = [x for x in eff.index if x not in set(ov_ids)]
check('OV CRISPR n', len(ov_c), 58)
check('OV RNA n', len(ov_r), 64)
check('comparator CRISPR n', len(comp_c), 1120)

scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan = scan.set_index('Gene')

def col(gene):
    return [c for c in eff.columns if c.startswith(gene + ' (')][0]

# --- dependency headline numbers, recomputed from the raw matrix ---
for gene, want_mean, want_dep, want_enr in [
        ('PAX8', -0.763616, 33, -0.587528),
        ('FERMT2', -0.917741, 46, -0.319090),
        ('CCNE1', -0.456771, 18, -0.272550),
        ('WWTR1', -0.483047, 24, -0.138341),
        ('YAP1', -0.407611, 18, -0.075267),
        ('TEAD1', -0.275069, 11, -0.040872),
        ('PTK2', -0.651728, 37, -0.091213),
        ('PKN2', -0.440525, 22, -0.069535),
        ('KLF5', -0.501726, 26, -0.138349)]:
    c = col(gene)
    a = eff.loc[ov_c, c].dropna()
    b = eff.loc[comp_c, c].dropna()
    check(f'{gene} OV mean', round(a.mean(), 6), want_mean, 1e-4)
    check(f'{gene} OV dep', int((a < -0.5).sum()), want_dep)
    check(f'{gene} enrichment', round(a.mean() - b.mean(), 6), want_enr, 1e-4)

# --- enrichment ranks, recomputed ---
enr = scan['Enrichment']
for gene, want_rank in [('PAX8', 1), ('FERMT2', 5), ('CCNE1', 8), ('WWTR1', 61), ('YAP1', 347), ('TEAD1', 1382), ('PTK2', 212), ('PKN2', 422), ('ITGAV', 133), ('WT1', 47), ('KLF5', 60)]:
    check(f'{gene} enrichment rank', int((enr < enr[gene]).sum()) + 1, want_rank)

# --- bar counts, recomputed from the scan file ---
check('strict pass', int(((scan.Enrichment <= -0.10) & (scan.OV_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 38)
check('relaxed pass', int(((scan.Enrichment <= -0.05) & (scan.OV_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)).sum()), 94)
check('q<0.05 count', int((scan.q_BH < 0.05).sum()), 1033)

# --- cell-line RNA, recomputed from the raw expression matrix ---
def expr_col(gene):
    return [c for c in expr.columns if c.startswith(gene + ' (')][0]

for gene, want_fc, want_d in [('PAX8', 20.196553, 3.802787), ('WWTR1', 2.353267, 1.177908), ('YAP1', 2.15377, 1.070998), ('TEAD1', 1.475372, 0.531901), ('FERMT2', 2.014742, 0.9589), ('CCNE1', 1.795895, 0.78865)]:
    c = expr_col(gene)
    a = expr.loc[ov_r, c].dropna()
    b = expr.loc[[x for x in expr.index if x not in set(ov_ids)], c].dropna()
    am, bm = a.mean(), b.mean()
    fc = (2 ** am - 1) / (2 ** bm - 1)
    check(f'{gene} cell-line RNA FC', round(fc, 4), want_fc, 1e-3)
    if want_d is not None:
        check(f'{gene} cell-line RNA DeltaLog2', round(am - bm, 4), want_d, 1e-3)

# --- cross-cancer TCGA RNA, recomputed from the raw JSON ---
emap = json.load(open(DL + 'entrez_map.json'))
def is_primary(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
ov_m = pd.DataFrame(json.load(open(DL + 'tcga_ov_mrna.json')))
other_m = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_mrna.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
ov_m = ov_m[ov_m['sampleId'].apply(is_primary)].copy()
other_m = other_m[other_m['sampleId'].apply(is_primary)].copy()
ov_m['value'] = np.log2(ov_m['value'].astype(float) + 1)
other_m['value'] = np.log2(other_m['value'].astype(float) + 1)
check('OV mRNA primary tumours', int(ov_m['sampleId'].nunique()), 303)
xc = pd.read_csv(OUT + 'tcga_rna.csv').set_index('Gene')
for g in ['PAX8', 'WT1', 'WWTR1', 'TEAD1', 'YAP1', 'CCNE1', 'FERMT2']:
    eid = emap[g]
    a = ov_m.loc[ov_m['entrezGeneId'] == eid, 'value'].values
    b = other_m.loc[other_m['entrezGeneId'] == eid, 'value'].values
    check(f'{g} TCGA XC log2FC', round(a.mean() - b.mean(), 4), round(xc.loc[g, 'TCGA_RNA_log2FC'], 4), 1e-3)

# --- RPPA, recomputed from the raw JSON ---
ov_rp = pd.DataFrame(json.load(open(DL + 'tcga_ov_rppa.json')))
other_rp = pd.concat([pd.DataFrame(json.load(open(DL + f'tcga_{c}_rppa.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
ov_rp = ov_rp[ov_rp['sampleId'].apply(is_primary)]
other_rp = other_rp[other_rp['sampleId'].apply(is_primary)]
rp = pd.read_csv(OUT + 'tcga_rppa.csv').set_index('Gene')
for g in ['CCNE1', 'WWTR1', 'YAP1', 'EGFR', 'PARP1']:
    eid = emap[g]
    a = ov_rp.loc[ov_rp['entrezGeneId'] == eid, 'value'].astype(float).values
    b = other_rp.loc[other_rp['entrezGeneId'] == eid, 'value'].astype(float).values
    check(f'{g} RPPA delta', round(a.mean() - b.mean(), 4), round(rp.loc[g, 'RPPA_delta'], 4), 1e-3)
check('RPPA bar pass count', int(((rp.RPPA_delta >= 0.5) & (rp.RPPA_p < 0.05)).sum()), 0)

# --- intersection counts ---
oe = pd.read_csv(OUT + 'overexpression.csv')
check('strict & overexpr', int(((oe.Strict_pass) & (oe.OverExpr_pass)).sum()), 15)
check('relaxed & overexpr', int(((oe.Relaxed_pass) & (oe.OverExpr_pass)).sum()), 28)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 #1', top20.iloc[0]['Symbol'], 'PAX8')
check('top20 #1 score', round(top20.iloc[0]['Score'], 2), 25.05, 1e-2)
check('top20 #2', top20.iloc[1]['Symbol'], 'ARHGAP29')
check('top20 #9', top20.iloc[8]['Symbol'], 'WWTR1')
check('top20 #20', top20.iloc[19]['Symbol'], 'FAAP24')

print()
if fails:
    print('FAILURES:', fails)
    raise SystemExit(1)
print('ALL CHECKS PASSED')
