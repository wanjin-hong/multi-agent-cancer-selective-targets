"""Independent verification of the osteosarcoma headlines, recomputed from
the raw files by a separate code path (no reuse of scan.py / merge.py
intermediates except the final CSVs being checked).

Checks: group counts (CRISPR 17, RNA 25, comparator sizes, non-cancerous
lines in each matrix), dependency means/fractions/enrichments/ranks for
the top-20 and TEAD-axis genes, bar counts (strict 100, relaxed 198,
strict+OE 21, relaxed+OE 33, candidates 33), cell-line RNA fold changes
(verified formula), TCGA cross-cancer log2FCs, combined scores for the
top 20, and the TP53/RB1 mutation counts. Prints PASS/FAIL per check.
"""
import csv, json, math, re, statistics
import pandas as pd

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/osteosarcoma/'
fails, passes = [], []

def check(name, got, want, tol=1e-6):
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    (passes if ok else fails).append(name)
    if not ok:
        print(f'FAIL {name}: got {got!r} want {want!r}')

# ---- group membership, read straight from Model.csv ----
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
os_ids = model.loc[model['OncotreePrimaryDisease'] == 'Osteosarcoma', 'ModelID'].tolist()
check('Model.csv OS lines', len(os_ids), 28)
with open(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv') as f:
    hdr = next(csv.reader(f))
    eff_ids = [row[0] for row in csv.reader(f)]
check('effect matrix lines', len(eff_ids), 1178)
check('effect matrix genes', len(hdr) - 1, 17916)
os_cr = [i for i in os_ids if i in set(eff_ids)]
check('OS CRISPR n', len(os_cr), 17)
check('comparator CRISPR n', len(eff_ids) - len(os_cr), 1161)
nc = model.loc[model['OncotreePrimaryDisease'] == 'Non-Cancerous', 'ModelID'].tolist()
check('non-cancerous in effect matrix', len([i for i in nc if i in set(eff_ids)]), 14)
with open(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv') as f:
    ehdr = next(csv.reader(f))
    ex_ids = [row[0] for row in csv.reader(f)]
check('expression matrix lines', len(ex_ids), 1673)
check('expression matrix genes', len(ehdr) - 1, 19193)
os_rna = [i for i in os_ids if i in set(ex_ids)]
check('OS RNA n', len(os_rna), 25)
check('comparator RNA n', len(ex_ids) - len(os_rna), 1648)
check('non-cancerous in expression matrix', len([i for i in nc if i in set(ex_ids)]), 85)
check('OS lines in both panels', len([i for i in os_cr if i in set(ex_ids)]), 17)

# ---- dependency values, recomputed with the csv module ----
eff = {}
with open(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv') as f:
    r = csv.reader(f)
    h = next(r)
    col = {name: j for j, name in enumerate(h)}
    for row in r:
        eff[row[0]] = row

def col_for(sym):
    return next(c for c in hdr[1:] if c.startswith(sym + ' ('))

def dep_stats(sym):
    c = col_for(sym)
    a = [float(eff[i][col[c]]) for i in os_cr if eff[i][col[c]] != '']
    b = [float(eff[i][col[c]]) for i in eff_ids if i not in set(os_cr) and eff[i][col[c]] != '']
    gw = [float(eff[i][col[c]]) for i in eff_ids if eff[i][col[c]] != '']
    return dict(mean=sum(a) / len(a), n=len(a), dep=sum(1 for x in a if x < -0.5),
                cmean=sum(b) / len(b), cfrac=sum(1 for x in b if x < -0.5) / len(b),
                enr=sum(a) / len(a) - sum(b) / len(b),
                fracdiff=sum(1 for x in a if x < -0.5) / len(a) - sum(1 for x in b if x < -0.5) / len(b),
                gwfrac=sum(1 for x in gw if x < -0.5) / len(gw))

scan = pd.read_csv(OUT + 'dependency_scan.csv')
def row(sym):
    return scan[scan.Symbol == sym].iloc[0]

for sym in ['RUNX2', 'FGFR1', 'PDGFRA', 'FERMT2', 'GPX8', 'EXT1', 'CBFB', 'IRS1',
            'WWTR1', 'SURF4', 'SEC23A', 'SLC35B2', 'ITGAV', 'E2F3', 'SCNM1', 'UBE2S',
            'EXT2', 'FBXO42', 'MT2A', 'H2AC11', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4',
            'YAP1', 'PTK2', 'PKN2', 'TP53', 'RB1', 'CDK2', 'KLF5', 'HNF1B', 'FOSL1']:
    d = dep_stats(sym)
    s = row(sym)
    check(f'{sym} OS mean', round(d['mean'], 6), round(s['OS_mean'], 6))
    check(f'{sym} OS n', d['n'], int(s['OS_n']))
    check(f'{sym} OS dep', d['dep'], int(s['OS_dep']))
    check(f'{sym} enrichment', round(d['enr'], 6), round(s['Enrichment'], 6))
    check(f'{sym} frac diff', round(d['fracdiff'], 6), round(s['FracDiff'], 6))
    check(f'{sym} gw frac', round(d['gwfrac'], 6), round(s['GenomeWide_dep_frac'], 6))
    check(f'{sym} enrich rank', int((scan.Enrichment < s['Enrichment']).sum()) + 1, int(s['EnrichRank']))

# bar counts from the scan CSV
def bar(df, strict):
    if strict:
        return df[(df.Enrichment <= -0.10) & (df.OS_dep_frac >= 0.25) & (df.FracDiff >= 0.10)
                  & (~df.PanEssential) & (df.GenomeWide_dep_frac <= 0.60)]
    return df[(df.Enrichment <= -0.05) & (df.OS_dep_frac >= 0.20) & (df.FracDiff >= 0.10)
              & (~df.PanEssential) & (df.GenomeWide_dep_frac <= 0.60)]
check('strict bar count', len(bar(scan, True)), 100)
check('relaxed bar count', len(bar(scan, False)), 198)
check('q<0.05 total', int((scan.q < 0.05).sum()), 101)
check('q<0.05 selective', int(((scan.q < 0.05) & (scan.Enrichment < 0)).sum()), 14)
check('pan-essential count', int(scan.PanEssential.sum()), 1282)
check('genes scored in 3 OS lines', int((scan.OS_n == 3).sum()), 129)

# ---- cell-line RNA, recomputed with the csv module ----
expr = {}
with open(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv') as f:
    r = csv.reader(f)
    eh = next(r)
    ecol = {name: j for j, name in enumerate(eh)}
    for row in r:
        expr[row[0]] = row
rna = pd.read_csv(OUT + 'celline_rna.csv')
def rna_stats(sym):
    c = next(cc for cc in eh[1:] if cc.startswith(sym + ' ('))
    a = [float(expr[i][ecol[c]]) for i in os_rna if expr[i][ecol[c]] != '']
    b = [float(expr[i][ecol[c]]) for i in ex_ids if i not in set(os_rna) and expr[i][ecol[c]] != '']
    am, bm = sum(a) / len(a), sum(b) / len(b)
    return dict(d=am - bm, fc=(2 ** am - 1) / (2 ** bm - 1))
for sym in ['RUNX2', 'FGFR1', 'PDGFRA', 'FERMT2', 'WWTR1', 'TEAD1', 'TEAD2', 'TEAD3',
            'YAP1', 'PTK2', 'PKN2', 'ITGAV', 'TP53', 'RB1', 'SP7', 'DLX5', 'SATB2']:
    d = rna_stats(sym)
    s = rna[rna.Gene == sym].iloc[0]
    check(f'{sym} cell-line dlog2', round(d['d'], 6), round(s['DeltaLog2'], 6))
    check(f'{sym} cell-line TPM FC', round(d['fc'], 6), round(s['TPM_FC'], 6))
check('cell-line RNA bar pass', int(((rna.TPM_FC >= 1.5) & (rna.DeltaLog2 >= 0.5) & (rna.p < 0.05)).sum()), 1700)

# ---- TCGA cross-cancer leg, recomputed from the downloaded JSON ----
osj = json.load(open(OUT + 'downloads/tcga_os_tpm.json'))
oth = []
for c in ['brca', 'luad', 'lusc', 'coad', 'stad', 'hcc']:
    oth += json.load(open(OUT + f'downloads/tcga_{c}_tpm.json'))
def primary(sid):
    return sid.split('-')[-1].startswith('01')
os_s = {x['sampleId'] for x in osj if primary(x['sampleId'])}
ot_s = {x['sampleId'] for x in oth if primary(x['sampleId'])}
check('TARGET-OS primary tumours', len(os_s), 88)
check('pooled other-cancer primary tumours', len(ot_s), 3347)
entrez = json.load(open(OUT + 'downloads/entrez_map.json'))
xc = pd.read_csv(OUT + 'tcga_rna.csv')
for sym in ['RUNX2', 'FGFR1', 'PDGFRA', 'FERMT2', 'WWTR1', 'TEAD1', 'TEAD3', 'TEAD4',
            'YAP1', 'PTK2', 'PKN2', 'SP7', 'DLX5', 'SATB2', 'CCNE1', 'GPX8']:
    eid = entrez[sym]
    a = [math.log2(float(x['value']) + 1) for x in osj if x['entrezGeneId'] == eid and primary(x['sampleId'])]
    b = [math.log2(float(x['value']) + 1) for x in oth if x['entrezGeneId'] == eid and primary(x['sampleId'])]
    fc = sum(a) / len(a) - sum(b) / len(b)
    s = xc[xc.Gene == sym].iloc[0]
    check(f'{sym} TCGA XC log2FC', round(fc, 6), round(s['TCGA_RNA_log2FC'], 6))
    check(f'{sym} TCGA XC n OS', len(a), int(s['TCGA_RNA_n_OS']))
check('TCGA XC bar pass', int(((xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)).sum()), 37)

# ---- intersection + scores ----
oe = pd.read_csv(OUT + 'overexpression.csv')
check('strict+OE', int(((oe.Strict_pass) & (oe.OverExpr_pass)).sum()), 21)
check('relaxed+OE', int(((oe.Relaxed_pass) & (oe.OverExpr_pass)).sum()), 33)
cand = pd.read_csv(OUT + 'candidates.csv')
check('candidates', len(cand), 33)
top20 = pd.read_csv(OUT + 'top20.csv')
check('top20 order', top20['Symbol'].tolist(),
      ['RUNX2', 'FGFR1', 'PDGFRA', 'FERMT2', 'GPX8', 'EXT1', 'CBFB', 'IRS1', 'WWTR1',
       'SURF4', 'SEC23A', 'SLC35B2', 'ITGAV', 'E2F3', 'SCNM1', 'UBE2S', 'EXT2',
       'FBXO42', 'MT2A', 'H2AC11'])
# recompute combined scores from the relaxed set in overexpression.csv
rel = oe[oe.Relaxed_pass].copy()
def z(s):
    s = s.astype(float)
    return (s - s.mean()) / s.std(ddof=1)
rel['Score2'] = z(-rel['Enrichment']) + z(rel['FracDiff']) + z(rel['TCGA_RNA_log2FC'].fillna(0)) + z(rel['DeltaLog2'].fillna(0))
for sym in ['RUNX2', 'FGFR1', 'PDGFRA', 'FERMT2', 'WWTR1', 'H2AC11']:
    got = float(rel[rel.Symbol == sym]['Score2'].iloc[0])
    want = float(cand[cand.Symbol == sym]['Score'].iloc[0])
    check(f'{sym} combined score', round(got, 4), round(want, 4))

# ---- mutations ----
mut = pd.read_csv(OUT + 'os_mutations_dedup.csv')
check('OS mutation rows after dedup', len(mut), 2883)
check('TP53-mutant lines', mut[mut.HugoSymbol == 'TP53'].ModelID.nunique(), 3)
check('RB1-mutant lines', mut[mut.HugoSymbol == 'RB1'].ModelID.nunique(), 5)
check('TP53 LoF lines', int(((mut.HugoSymbol == 'TP53') & (mut.LikelyLoF == True)).sum()), 3)

print(f'\n{len(passes)} checks passed, {len(fails)} failed')
if fails:
    print('FAILED:', fails)
