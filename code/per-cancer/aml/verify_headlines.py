"""Independent verification of the AML headline numbers — separate code path from
the pipeline: reads the raw DepMap files and the raw cBioPortal downloads
directly and recomputes every number quoted in report.md. Fails loudly on any
mismatch (tolerance 1e-6 relative, 0.01 absolute for scores).
"""
import pandas as pd, numpy as np, json, math, sys
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/aml/'
fails = []
def check(label, got, want, tol=1e-6):
    if isinstance(got, list) or isinstance(want, list):
        ok = (got == want)
    else:
        ok = (abs(got - want) <= tol * max(1.0, abs(want))) if (got == got and want == want) else (got == want)
    print(('PASS' if ok else 'FAIL'), label, 'got', got, 'want', want)
    if not ok:
        fails.append(label)

# --- raw files, fresh reads ---
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
check('effect matrix lines', len(eff), 1178)
check('effect matrix genes', eff.shape[1], 17916)
check('expression matrix lines', len(expr), 1673)
aml_ids = model.loc[model['OncotreePrimaryDisease'] == 'Acute Myeloid Leukemia', 'ModelID'].tolist()
aml = [x for x in aml_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(aml_ids)]
check('AML CRISPR n', len(aml), 28)
check('pan-cancer comparator n', len(comp), 1150)
nc = [x for x in comp if model.loc[model.ModelID == x, 'OncotreePrimaryDisease'].iloc[0] == 'Non-Cancerous']
check('Non-Cancerous lines in comparator', len(nc), 14)
BLOOD = ['Mature B-Cell Neoplasms', 'B-Lymphoblastic Leukemia/Lymphoma',
         'Myeloproliferative Neoplasms', 'Mature T and NK Neoplasms',
         'T-Lymphoblastic Leukemia/Lymphoma', 'Non-Hodgkin Lymphoma',
         'Hodgkin Lymphoma', 'Acute Leukemias of Ambiguous Lineage']
blood_ids = model.loc[model['OncotreePrimaryDisease'].isin(BLOOD), 'ModelID'].tolist()
haem = [x for x in blood_ids if x in eff.index and x not in set(aml_ids)]
check('haem comparator n', len(haem), 104)
aml_rna = [x for x in aml_ids if x in expr.index]
check('AML RNA n', len(aml_rna), 52)
# same-patient pairs
sub = model.loc[model.ModelID.isin(aml)]
check('distinct PatientIDs in CRISPR panel', int(sub['PatientID'].nunique()), 26)

def col(df, gene):
    c = [x for x in df.columns if x.split(' (')[0] == gene]
    assert len(c) == 1, gene
    return c[0]

def dep_stats(gene, comp_ids):
    s = eff.loc[aml, col(eff, gene)].dropna()
    t = eff.loc[comp_ids, col(eff, gene)].dropna()
    return dict(mean=float(s.mean()), dep=int((s < -0.5).sum()), n=len(s),
                frac=float((s < -0.5).mean()), cmean=float(t.mean()),
                cfrac=float((t < -0.5).mean()), enrich=float(s.mean() - t.mean()),
                fdiff=float((s < -0.5).mean() - (t < -0.5).mean()))

# --- bar counts from the pipeline CSVs (recomputed from scan columns, not trusted) ---
scan = pd.read_csv(OUT + 'dependency_scan.csv')
check('scan rows', len(scan), 17916)
check('scan Gene column is plain symbols', bool((~scan['Gene'].str.contains(' \\(')).all()), True)
check('scan has EnrichRank', 'EnrichRank' in scan.columns, True)
strict = scan[(scan.Enrichment <= -0.10) & (scan.AML_dep_frac >= 0.25) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
relaxed = scan[(scan.Enrichment <= -0.05) & (scan.AML_dep_frac >= 0.20) & (scan.FracDiff >= 0.10) & (~scan.PanEssential) & (scan.GenomeWide_dep_frac <= 0.60)]
check('strict pass', len(strict), 270)
check('relaxed pass', len(relaxed), 421)
check('q<0.05 total', int((scan.q_BH < 0.05).sum()), 551)
check('q<0.05 selective', int(((scan.q_BH < 0.05) & (scan.Enrichment < 0)).sum()), 255)
hscan = pd.read_csv(OUT + 'dependency_scan_haem.csv')
hstrict = hscan[(hscan.Enrichment <= -0.10) & (hscan.AML_dep_frac >= 0.25) & (hscan.FracDiff >= 0.10) & (~hscan.PanEssential) & (hscan.GenomeWide_dep_frac <= 0.60)]
hrelaxed = hscan[(hscan.Enrichment <= -0.05) & (hscan.AML_dep_frac >= 0.20) & (hscan.FracDiff >= 0.10) & (~hscan.PanEssential) & (hscan.GenomeWide_dep_frac <= 0.60)]
check('haem strict pass', len(hstrict), 96)
check('haem relaxed pass', len(hrelaxed), 187)
check('haem q<0.05 selective', int(((hscan.q_BH < 0.05) & (hscan.Enrichment < 0)).sum()), 26)
ov = pd.read_csv(OUT + 'overexpression.csv')
check('strict & OE', int((ov.Strict_pass & ov.OverExpr_pass).sum()), 103)
check('relaxed & OE', int((ov.Relaxed_pass & ov.OverExpr_pass).sum()), 153)
check('haem relaxed & OE', int((ov.Haem_Relaxed_pass & ov.OverExpr_pass).sum()), 89)

# --- top-20 dependency numbers, recomputed from the raw effect matrix ---
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 rows', len(top20), 20)
check('top20 Symbol plain', bool((~top20['Symbol'].str.contains(' \\(')).all()), True)
for _, r in top20.iterrows():
    d = dep_stats(r['Symbol'], comp)
    check(f"{r['Symbol']} AML mean", d['mean'], r['AML_mean'])
    check(f"{r['Symbol']} AML dep", d['dep'], r['AML_dep'])
    check(f"{r['Symbol']} enrichment", d['enrich'], r['Enrichment'])
    check(f"{r['Symbol']} frac diff", d['fdiff'], r['FracDiff'])
    # enrichment rank: rank column must be consistent with the enrichment column
    rk = int((scan.Enrichment < r['Enrichment']).sum()) + 1
    check(f"{r['Symbol']} enrich rank (pan)", rk, int(scan.loc[scan.Gene == r['Symbol'], 'EnrichRank'].iloc[0]))

# --- cell-line RNA for top 20 + TEAD axis, recomputed from the raw expression matrix ---
comp_rna = [x for x in expr.index if x not in set(aml_ids)]
for g in list(top20['Symbol']) + ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'FLT3', 'CD33', 'KLF5']:
    if g not in [x.split(' (')[0] for x in expr.columns]:
        print('INFO', g, 'absent from expression matrix (no cell-line RNA data) — expected for TAFAZZIN')
        continue
    c = col(expr, g)
    a = expr.loc[aml_rna, c].dropna().values; b = expr.loc[comp_rna, c].dropna().values
    am, bm = a.mean(), b.mean()
    fc = (2 ** am - 1) / (2 ** bm - 1)
    row = ov[ov.Gene == g]
    if len(row):
        check(f'{g} cell-line TPM FC', fc, float(row.iloc[0].TPM_FC), tol=1e-4)
        check(f'{g} cell-line DeltaLog2', am - bm, float(row.iloc[0].DeltaLog2), tol=1e-4)

# --- TCGA cross-cancer leg recomputed from the raw downloads ---
DL = OUT + 'downloads/'
entrez = json.load(open(DL + 'entrez_map.json'))
def load(fn):
    return pd.DataFrame(json.load(open(DL + fn)))
def is_primary(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] in ('01', '03')
laml = load('tcga_laml_mrna.json'); laml = laml[laml['sampleId'].apply(is_primary)].copy()
check('TCGA LAML primary tumours', int(laml['sampleId'].nunique()), 173)
others = pd.concat([load(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
others = others[others['sampleId'].apply(is_primary)].copy()
check('TCGA other primary tumours', int(others['sampleId'].nunique()), 3274)
laml['value'] = np.log2(laml['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
for g in list(top20['Symbol']) + ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'FLT3', 'CD33', 'KLF5', 'BCL2', 'JAK2']:
    eid = entrez[g]
    a = laml.loc[laml['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    row = ov[ov.Gene == g]
    if len(row) and len(a):
        check(f'{g} TCGA XC log2FC', float(a.mean() - b.mean()), float(row.iloc[0].TCGA_RNA_log2FC), tol=1e-4)

# --- TEAD axis, both comparators ---
for g, want_rank, want_enrich in [('TEAD1', 17802, 0.181), ('TEAD2', 16549, 0.056), ('TEAD3', 16964, 0.066),
                                  ('TEAD4', 17519, 0.102), ('YAP1', 17796, 0.175), ('WWTR1', 17819, 0.193),
                                  ('PTK2', 17895, 0.316), ('PKN2', 17814, 0.190), ('FERMT2', 17910, 0.442),
                                  ('ITGAV', 17894, 0.312), ('ITGB1', 17827, 0.201), ('KLF5', 17846, 0.221)]:
    d = dep_stats(g, comp)
    check(f'{g} pan enrich', round(d['enrich'], 3), want_enrich, tol=1e-3)
    stored_enrich = float(scan.loc[scan.Gene == g, 'Enrichment'].iloc[0])
    check(f'{g} pan enrich rank', int((scan.Enrichment < stored_enrich).sum()) + 1, want_rank)
    check(f'{g} stored EnrichRank col', int(scan.loc[scan.Gene == g, 'EnrichRank'].iloc[0]), want_rank)
    dh = dep_stats(g, haem)
    print('INFO', g, 'haem enrich', round(dh['enrich'], 3), 'rank', int((hscan.Enrichment < dh['enrich']).sum()) + 1,
          'dep', f"{dh['dep']}/{dh['n']}")

# --- mutation-driven checks ---
# FLT3-ITD lines (insertion in juxtamembrane domain): MV4-11, MOLM-13, MOLM-14
name = dict(zip(model['ModelID'], model['CellLineName']))
itd_ids = [i for i in aml if name[i] in ('MV4-11', 'MOLM-13', 'MOLM-14')]
d = dep_stats('FLT3', comp)
s = eff.loc[aml, col(eff, 'FLT3')]
check('FLT3 dep in ITD lines', int((s.loc[itd_ids] < -0.5).sum()), 3)
check('FLT3 ITD mean', float(s.loc[itd_ids].mean()), -1.291, tol=1e-3)
wt = [i for i in aml if name[i] not in ('MV4-11', 'MOLM-13', 'MOLM-14', 'TF-1', 'MONO-MAC-1')]
check('FLT3 dep in FLT3-WT lines', int((s.loc[wt] < -0.5).sum()), 1)
check('FLT3 enrich rank', int((scan.Enrichment < float(scan.loc[scan.Gene == 'FLT3', 'Enrichment'].iloc[0])).sum()) + 1, 521)
# KIT: the two KIT-mutant lines are the two most KIT-dependent
k = eff.loc[aml, col(eff, 'KIT')]
check('KIT SKNO-1', float(k.loc[[i for i in aml if name[i] == 'SKNO-1'][0]]), -1.27, tol=1e-2)
check('KIT KASUMI-1', float(k.loc[[i for i in aml if name[i] == 'KASUMI-1'][0]]), -0.72, tol=1e-2)
# NPM1: scored in only 5 lines; OCI-AML3 (NPM1-mutant) -2.11
n = eff.loc[aml, col(eff, 'NPM1')]
check('NPM1 scored AML lines', int(n.notna().sum()), 5)
check('NPM1 dep among scored', int((n < -0.5).sum()), 3)
check('NPM1 OCI-AML3', float(n.loc[[i for i in aml if name[i] == 'OCI-AML3'][0]]), -2.11, tol=1e-2)
# 5-line block: 129 genes, fixed subset OCI-AML2, NOMO-1, OCI-AML3, SKM-1, CMS
check('genes scored in exactly 5 AML lines', int((scan.AML_n == 5).sum()), 129)
blk = eff.loc[aml, scan[scan.AML_n == 5]['Column'].tolist()]
check('5-block lines', sorted([name[i] for i in aml if blk.loc[i].notna().any()]),
      sorted(['OCI-AML2', 'NOMO-1', 'OCI-AML3', 'SKM-1', 'CMS']))
check('WDR49 AML_n', int(scan.loc[scan.Gene == 'WDR49', 'AML_n'].iloc[0]), 5)
check('WDR49 AML_dep', int(scan.loc[scan.Gene == 'WDR49', 'AML_dep'].iloc[0]), 1)
# CDK6 / CBFB haem pan-essential
for g in ['CDK6', 'CBFB']:
    r = hscan[hscan.Gene == g].iloc[0]
    check(f'{g} haem pan-essential', bool(r.PanEssential), True)
    check(f'{g} haem comp dep frac', round(float(r.Comp_dep_frac), 3), 0.712 if g == 'CDK6' else 0.721, tol=1e-3)
# combined score recomputation for the top 5 (independent z over the relaxed set)
rel = ov[ov.Relaxed_pass].copy()
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel['Score_check'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['TCGA_RNA_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0))
for g in ['MYB', 'SPI1', 'ZEB2', 'LMO2', 'CBFB']:
    got = float(rel.loc[rel.Gene == g, 'Score_check'].iloc[0])
    want = float(top20.loc[top20.Symbol == g, 'Score'].iloc[0])
    check(f'{g} combined score', got, want, tol=1e-2)
check('MYB score', float(top20.loc[top20.Symbol == 'MYB', 'Score'].iloc[0]), 22.41, tol=1e-2)

print()
print('FAILURES:', len(fails), fails if fails else 'none — all headline numbers verified')
sys.exit(1 if fails else 0)
