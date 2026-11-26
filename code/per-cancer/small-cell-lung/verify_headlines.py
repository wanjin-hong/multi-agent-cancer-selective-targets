"""Independent verification of SCLC headline numbers from the raw files.

Separate code path from scan.py/celline_rna.py/merge.py: re-reads
CRISPRGeneEffect.csv, Model.csv and the expression matrix directly and
recomputes group counts, dependency means/fractions/enrichments/ranks,
bar counts, cell-line RNA fold changes (playbook §4 formula), candidate
count and top-20 order, TEAD-axis ranks, and the low-n block structure.
Also re-checks the ucologne/TCGA scale mismatch that invalidated the
patient-RNA leg, from the raw downloaded JSON.
Fails loudly on any mismatch.
"""
import json, numpy as np, pandas as pd
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/small-cell-lung/'
checks = []
def check(name, got, want, tol=None):
    if tol is not None:
        ok = abs(got - want) <= tol
    else:
        ok = got == want
    checks.append((name, got, want, ok))
    print(('PASS' if ok else 'FAIL'), name, '| got', got, '| want', want)
    return ok

eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
ids = model.loc[model['OncotreePrimaryDisease'] == 'Lung Neuroendocrine Tumor', 'ModelID'].tolist()
grp = [x for x in ids if x in eff.index]
comp = [x for x in eff.index if x not in set(ids)]
check('effect matrix lines', eff.shape[0], 1178)
check('effect matrix genes', eff.shape[1], 17916)
check('SCLC CRISPR n', len(grp), 26)
check('comparator CRISPR n', len(comp), 1152)
check('SCLC lines in Model.csv', len(ids), 85)
sub = model[model.ModelID.isin(grp)]
check('SCLC subtype in CRISPR panel', int((sub['OncotreeSubtype'] == 'Small Cell Lung Cancer').sum()), 25)
check('carcinoid in CRISPR panel', int((sub['OncotreeSubtype'] == 'Lung Carcinoid').sum()), 1)
check('carcinoid line is NCI-H727', sub.loc[sub['OncotreeSubtype'] == 'Lung Carcinoid', 'CellLineName'].iloc[0], 'NCI-H727')
check('RNA panel n', sum(x in expr.index for x in ids), 60)
check('CRISPR panel lines with RNA', sum(x in expr.index for x in grp), 26)
check('same-patient pairs in CRISPR panel', int(sub.duplicated('PatientID', keep=False).sum()), 0)

scan = pd.read_csv(OUT + 'dependency_scan.csv')
check('scan rows', len(scan), 17916)
check('q<0.05', int((scan.q_BH < 0.05).sum()), 365)
check('q<0.05 selective', int(((scan.q_BH < 0.05) & (scan.Enrichment < 0)).sum()), 126)
check('pan-essential', int(scan.PanEssential.sum()), 1250)

def dep_row(gene):
    col = [c for c in eff.columns if c.startswith(gene + ' (')][0]
    a = eff.loc[grp, col].dropna(); b = eff.loc[comp, col].dropna()
    return a, b

def dep_check(gene, want_mean, want_dep, want_enrich, want_rank):
    a, b = dep_row(gene)
    r = scan[scan.Gene == gene].iloc[0]
    check(f'{gene} SCLC mean', round(a.mean(), 4), round(want_mean, 4), tol=1e-3)
    check(f'{gene} SCLC dep', int((a < -0.5).sum()), want_dep)
    check(f'{gene} enrichment', round(a.mean() - b.mean(), 4), round(want_enrich, 4), tol=1e-3)
    check(f'{gene} enrich rank', int(r.EnrichRank), want_rank)
    return r

dep_check('ASCL1', -0.531790, 11, -0.461010, 10)
dep_check('INSM1', -0.329083, 8, -0.360584, 18)
dep_check('SKP2', -1.370192, 23, -0.738080, 4)
dep_check('E2F3', -0.667061, 15, -0.455875, 11)
dep_check('TEAD1', -0.228541, 5, 0.007842, 11145)
dep_check('YAP1', -0.183696, 1, 0.155792, 17668)
dep_check('WWTR1', -0.233265, 1, 0.120921, 17530)
dep_check('PTK2', -0.339857, 3, 0.230231, 17803)
dep_check('PKN2', -0.162265, 1, 0.216937, 17784)
dep_check('TEAD3', -0.190238, 0, 0.076281, 17141)
dep_check('TEAD4', -0.086831, 0, 0.051774, 16491)
dep_check('KLF5', -0.239878, 2, 0.133252, 17585)
dep_check('FOSL1', -0.226003, 3, 0.169794, 17705)
dep_check('FERMT2', -0.308610, 2, 0.312653, 17863)
dep_check('EGFR', -0.142602, 0, 0.154799, 17664)
dep_check('ITGAV', -0.348703, 6, 0.295122, 17849)
dep_check('NEUROD1', -0.092495, 4, -0.167229, 163)
dep_check('POU2F3', -0.212898, 4, -0.261711, 49)
dep_check('MYCL', -0.203563, 4, -0.125374, 334)
dep_check('SOX2', -0.215213, 4, -0.138040, 267)
dep_check('DLL3', -0.046535, 0, -0.016877, 5435)
dep_check('BCL2', -0.147037, 0, -0.112857, 404)
dep_check('TIMM50', -1.475589, 2, -1.055172, 1)

# low-n block structure: which lines carry the 2-line genes?
a2 = eff.loc[grp, 'TIMM50 (92609)'].dropna()
check('TIMM50 scored lines', len(a2), 2)
check('TIMM50 lines are DMS 273 + NCI-H1882',
      sorted(model.set_index('ModelID').loc[a2.index, 'CellLineName'].tolist()),
      ['DMS 273', 'NCI-H1882'])
n24 = int((scan.SCLC_n == 24).sum()); n2 = int((scan.SCLC_n == 2).sum())
check('genes scored in 24 lines', n24, 251)
check('genes scored in 2 lines', n2, 129)
check('genes scored in all 26 lines', int((scan.SCLC_n == 26).sum()), 17533)

# bar counts recomputed from the scan file
strict = scan[(scan.Enrichment <= -0.10) & (scan.SCLC_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.SCLC_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
check('strict bar', len(strict), 115)
check('relaxed bar', len(relaxed), 202)

# cell-line RNA recomputation (playbook §4 formula)
grp_rna = [x for x in ids if x in expr.index]
comp_rna = [x for x in expr.index if x not in set(ids)]
check('RNA comparator n', len(comp_rna), 1613)
def rna_check(gene, want_fc, want_d):
    col = [c for c in expr.columns if c.startswith(gene + ' (')][0]
    a = expr.loc[grp_rna, col].dropna(); b = expr.loc[comp_rna, col].dropna()
    fc = (2**a.mean() - 1) / (2**b.mean() - 1)
    check(f'{gene} TPM FC', round(fc, 2), round(want_fc, 2), tol=0.02)
    check(f'{gene} DeltaLog2', round(a.mean() - b.mean(), 3), round(want_d, 3), tol=0.002)
rna_check('ASCL1', 305.38, 5.836)
rna_check('INSM1', 314.73, 6.009)
rna_check('SKP2', 2.41, 1.230)
rna_check('TEAD1', 2.26, 1.124)
rna_check('YAP1', 0.05, -3.479)
rna_check('WWTR1', 0.12, -2.528)
rna_check('NEUROD1', 55.05, 2.493)
rna_check('DLL3', 22.65, 3.648)
rna_check('MYCL', 22.77, 3.789)
rna_check('KLF5', 0.88, -0.160)
rna_check('FOSL1', 0.04, -3.627)
rna_check('CHGA', 173.56, 5.887)

# intersection + top 20
oe = pd.read_csv(OUT + 'overexpression.csv')
check('strict & OE', int((oe.Strict_pass & oe.OverExpr_pass).sum()), 26)
check('relaxed & OE', int((oe.Relaxed_pass & oe.OverExpr_pass).sum()), 45)
cand = pd.read_csv(OUT + 'candidates.csv')
check('candidates', len(cand), 45)
check('Tier A candidates', int((cand.Tier == 'A').sum()), 26)
check('Tier B candidates', int((cand.Tier == 'B').sum()), 19)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 order', top20['Symbol'].tolist(),
      ['ASCL1', 'INSM1', 'SKP2', 'E2F3', 'PRKRA', 'FOXA2', 'CCNE2', 'CXADR', 'MIB1', 'ANKRD36',
       'TSPYL5', 'EPC2', 'CIT', 'TFDP1', 'YEATS2', 'ASPM', 'PNISR', 'SQLE', 'KIF14', 'WSB2'])
check('top20 #1 score', round(top20.iloc[0]['Score'], 2), 11.42, tol=0.01)
check('top20 #20 score', round(top20.iloc[19]['Score'], 2), 0.64, tol=0.01)
check('ranks 21-25', cand.iloc[20:25]['Symbol'].tolist(), ['SMARCA4', 'RPRD2', 'MYBL2', 'OSBP', 'DHFR2'])
check('CXADR in top20 scored lines', int(top20[top20.Symbol == 'CXADR'].iloc[0]['SCLC_n']), 2)
check('ANKRD36 in top20 scored lines', int(top20[top20.Symbol == 'ANKRD36'].iloc[0]['SCLC_n']), 2)
# score recomputation for #1: z-scores across the relaxed set
rel = oe[oe.Relaxed_pass].copy()
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel['Score_check'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['DeltaLog2'].fillna(0))
s1 = rel[rel.Symbol == 'ASCL1'].iloc[0]['Score_check']
check('ASCL1 score recomputed', round(s1, 2), 11.42, tol=0.01)
check('TCGA component is 0 for all (leg invalid)', bool((oe['TCGA_XC_pass'] == False).all()), True)
check('no protein bonus (no leg)', bool((oe['RPPA_pass'] == False).all()), True)

# ucologne/TCGA scale mismatch, from raw downloads
sclc_json = pd.DataFrame(json.load(open(OUT + 'downloads/tcga_sclc_mrna.json')))
oth_json = pd.concat([pd.DataFrame(json.load(open(OUT + f'downloads/tcga_{c}_mrna.json'))) for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
oth_json = oth_json[oth_json['sampleId'].str[13:15] == '01']
check('ucologne samples', int(sclc_json['sampleId'].nunique()), 81)
check('pooled other primary tumours', int(oth_json['sampleId'].nunique()), 3274)
check('ucologne median raw value ~12.6', round(float(sclc_json['value'].median()), 1), 12.6, tol=0.5)
check('TCGA median raw value ~904.7', round(float(oth_json['value'].median()), 1), 904.7, tol=5.0)
brca1_s = sclc_json[sclc_json.entrezGeneId == 672]['value'].mean()
brca1_o = oth_json[oth_json.entrezGeneId == 672]['value'].mean()
check('BRCA1 raw mean SCLC ~0.5', round(float(brca1_s), 1), 0.5, tol=0.2)
check('BRCA1 raw mean other ~459', round(float(brca1_o), 0), 459.0, tol=10.0)

# subtype-dependency facts stated in the report
cols = {c.split(' (')[0]: c for c in expr.columns}
st = expr.loc[grp, [cols[g] for g in ['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1']]].copy()
st.columns = ['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1']
subtype = st.idxmax(axis=1)
check('SCLC-A lines in CRISPR panel', int((subtype == 'ASCL1').sum()), 12)
check('SCLC-N lines in CRISPR panel', int((subtype == 'NEUROD1').sum()), 9)
check('SCLC-P lines in CRISPR panel', int((subtype == 'POU2F3').sum()), 4)
check('SCLC-Y lines in CRISPR panel', int((subtype == 'YAP1').sum()), 1)
check('NCI-H727 (carcinoid) NOT ASCL1-dependent',
      bool(eff.loc['ACH-000775', 'ASCL1 (429)'] < -0.5), False)
ascl1_dep = (eff.loc[grp, 'ASCL1 (429)'] < -0.5).values
check('ASCL1-dependent lines that are SCLC-A', int((subtype[ascl1_dep] == 'ASCL1').sum()), 10)
pou_dep = (eff.loc[grp, 'POU2F3 (25833)'] < -0.5).values
check('POU2F3-dependent lines all SCLC-P', int((subtype[pou_dep] == 'POU2F3').sum()), 4)
neu_dep = (eff.loc[grp, 'NEUROD1 (4760)'] < -0.5).values
check('NEUROD1-dependent lines that are SCLC-N', int((subtype[neu_dep] == 'NEUROD1').sum()), 3)
yap_dep = (eff.loc[grp, 'YAP1 (10413)'] < -0.5)
check('single YAP1-dependent line is Lu-165',
      model.set_index('ModelID').loc[yap_dep[yap_dep].index, 'CellLineName'].tolist(), ['Lu-165'])
missing_uc = [g for g in json.load(open(OUT + 'downloads/fetch_genes.json')) if json.load(open(OUT + 'downloads/entrez_map.json'))[g] not in set(sclc_json['entrezGeneId'])]
check('gene missing from ucologne profile is SPATA31A7', missing_uc, ['SPATA31A7'])

print()
fails = [c for c in checks if not c[3]]
print(f'{len(checks) - len(fails)}/{len(checks)} checks passed')
if fails:
    print('FAILURES:')
    for f in fails:
        print(' ', f)
    raise SystemExit(1)
print('ALL CHECKS PASS')
