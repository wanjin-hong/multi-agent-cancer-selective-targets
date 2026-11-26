"""Independent verification of every headline number in the Ewing report.

Separate code path from the pipeline: reads the raw CSVs with the csv module
(no pandas merge), recomputes means/fractions/enrichments/RNA fold changes
with the playbook formula, re-derives the bar counts and the top-20 order
from the recomputed z-scores, and re-checks the STAG2/STAG1, low-n-block,
comparator-composition and no-patient-leg claims. Exits non-zero on any
mismatch. Run: python3 verify_headlines.py
"""
import csv, json, math, statistics, sys
from collections import Counter

BASE = '/home/hatch/workspace/research_notes/'
EW = BASE + 'multi-cancer-selective-targets-20261120/ewing/'
fails = []
def check(name, got, want, tol=0.0):
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print(('PASS' if ok else 'FAIL'), name, 'got', got, 'want', want)
    if not ok:
        fails.append(name)

# ---------- Model.csv: group ----------
with open(BASE + 'vrk1-codependency-20261108/Model.csv') as f:
    model = list(csv.DictReader(f))
ew_model = [r for r in model if r['OncotreePrimaryDisease'] == 'Ewing Sarcoma']
check('Ewing lines in Model.csv', len(ew_model), 51)

# ---------- CRISPRGeneEffect.csv (csv module, header gene labels) ----------
with open(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv') as f:
    rdr = csv.reader(f)
    hdr = next(rdr)
    eff = {row[0]: row[1:] for row in rdr}
gene_of = {h.split(' (')[0]: h for h in hdr[1:]}
check('CRISPR genes', len(hdr) - 1, 17916)
check('CRISPR lines', len(eff), 1178)
check('CD99 absent from CRISPR scan', 'CD99' in gene_of, False)
ew_cr = [r['ModelID'] for r in ew_model if r['ModelID'] in eff]
check('Ewing CRISPR n', len(ew_cr), 27)
comp = [m for m in eff if m not in set(r['ModelID'] for r in ew_model)]
check('comparator CRISPR n', len(comp), 1151)
nc = [r for r in model if r['OncotreePrimaryDisease'] == 'Non-Cancerous' and r['ModelID'] in eff]
check('non-cancerous lines in CRISPR comparator', len(nc), 14)
dups = Counter(r['PatientID'] for r in ew_model if r['ModelID'] in eff)
check('same-patient pairs in CRISPR group', sum(v - 1 for v in dups.values() if v > 1), 1)

def colvals(gene, ids):
    j = hdr.index(gene_of[gene]) - 1
    return [float(eff[i][j]) for i in ids if eff[i][j] != '']

def depstats(gene):
    a = colvals(gene, ew_cr); b = colvals(gene, comp)
    gw = [float(v) for row in eff.values() for v in [row[hdr.index(gene_of[gene]) - 1]] if v != '']
    return dict(n=len(a), mean=statistics.fmean(a), dep=sum(v < -0.5 for v in a),
                cmean=statistics.fmean(b), cdep=sum(v < -0.5 for v in b) / len(b),
                enrich=statistics.fmean(a) - statistics.fmean(b),
                fracdiff=sum(v < -0.5 for v in a) / len(a) - sum(v < -0.5 for v in b) / len(b),
                gwfrac=sum(v < -0.5 for v in gw) / len(gw))

for g, want_mean, want_dep, want_enrich in [
        ('FLI1', -0.934, 20, -0.931), ('EWSR1', -1.792, 27, -0.950),
        ('ZEB2', -0.673, 21, -0.353), ('STAG1', -0.770, 15, -0.635),
        ('CDK4', -1.530, 24, -0.783), ('CCND1', -1.467, 24, -0.496),
        ('TRIM8', -0.742, 19, -0.654), ('IGF1R', -0.730, 19, -0.425),
        ('TEAD1', -0.062, 0, 0.179), ('WWTR1', -0.163, 0, 0.193),
        ('PTK2', -0.279, 3, 0.293), ('PKN2', -0.202, 0, 0.177),
        ('YAP1', -0.216, 2, 0.123), ('KLF5', -0.188, 2, 0.186),
        ('EGFR', -0.100, 0, 0.199), ('CDK6', -0.048, 0, 0.553),
        ('BCL11B', -0.246, 6, -0.268), ('MT1E', -0.676, 16, -0.202),
        ('MYBL2', -0.525, 15, -0.056), ('ERG', -0.046, 0, 0.007)]:
    d = depstats(g)
    check(f'{g} Ewing mean', round(d['mean'], 3), want_mean, 0.001)
    check(f'{g} Ewing dep', d['dep'], want_dep)
    check(f'{g} enrichment', round(d['enrich'], 3), want_enrich, 0.001)

# enrichment ranks from the scan file (recomputed ordering check for key genes)
scan = list(csv.DictReader(open(EW + 'dependency_scan.csv')))
order = sorted(scan, key=lambda r: float(r['Enrichment']))
rank = {r['Gene'].split(' (')[0]: i + 1 for i, r in enumerate(order)}
for g, want_rank in [('EWSR1', 1), ('FLI1', 2), ('CDK4', 3), ('TRIM8', 7),
                     ('STAG1', 8), ('CCND1', 12), ('IGF1R', 22), ('BCL11B', 84),
                     ('TEAD2', 2657), ('TEAD1', 17759), ('WWTR1', 17794),
                     ('PTK2', 17882), ('CDK6', 17916), ('KLF5', 17779),
                     ('EGFR', 17807), ('FOSL1', 17816), ('FERMT2', 17893),
                     ('ITGAV', 17868), ('YAP1', 17547), ('PKN2', 17755)]:
    check(f'{g} enrichment rank', rank[g], want_rank)

# bar counts recomputed from the scan file
def passes(r, strict):
    e, f, d, gw = float(r['Enrichment']), float(r['Ewing_dep_frac']), float(r['FracDiff']), float(r['GenomeWide_dep_frac'])
    pe = r['PanEssential'] == 'True'
    if strict:
        return e <= -0.10 and f >= 0.25 and d >= 0.10 and not pe and gw <= 0.60
    return e <= -0.05 and f >= 0.20 and d >= 0.10 and not pe and gw <= 0.60
check('strict-bar genes', sum(passes(r, True) for r in scan), 159)
check('relaxed-bar genes', sum(passes(r, False) for r in scan), 272)
check('q<0.05 total', sum(float(r['q']) < 0.05 for r in scan), 749)
check('q<0.05 selective-direction', sum(float(r['q']) < 0.05 and float(r['Enrichment']) < 0 for r in scan), 302)
check('pan-essential genes', sum(r['PanEssential'] == 'True' for r in scan), 1293)
ncount = Counter(int(r['Ewing_n']) for r in scan)
check('genes scored in 12 Ewing lines', ncount[12], 129)
check('genes scored in 19 Ewing lines', ncount[19], 580)
check('genes scored in 18 Ewing lines', ncount[18], 2)
check('genes scored in all 27 Ewing lines', ncount[27], 17205)

# ---------- expression matrix ----------
with open(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv') as f:
    rdr = csv.reader(f)
    ehdr = next(rdr)
    expr = {row[0]: row[1:] for row in rdr}
egene = {h.split(' (')[0]: h for h in ehdr[1:]}
ew_rna = [r['ModelID'] for r in ew_model if r['ModelID'] in expr]
check('Ewing RNA n', len(ew_rna), 24)
check('Ewing lines in BOTH CRISPR and RNA panels',
      len(set(ew_rna) & set(ew_cr)), 17)
check('CRISPR lines missing RNA', len(set(ew_cr) - set(ew_rna)), 10)
check('RNA lines missing CRISPR', len(set(ew_rna) - set(ew_cr)), 7)
comp_rna = [m for m in expr if m not in set(r['ModelID'] for r in ew_model)]
check('comparator RNA n', len(comp_rna), 1649)
nc_rna = [r for r in model if r['OncotreePrimaryDisease'] == 'Non-Cancerous' and r['ModelID'] in expr]
check('non-cancerous lines in RNA comparator', len(nc_rna), 85)

def rnafc(gene):
    j = ehdr.index(egene[gene]) - 1
    a = [float(expr[i][j]) for i in ew_rna if expr[i][j] != '']
    b = [float(expr[i][j]) for i in comp_rna if expr[i][j] != '']
    am, bm = statistics.fmean(a), statistics.fmean(b)
    return am - bm, (2 ** am - 1) / (2 ** bm - 1)
for g, want_d, want_fc in [('FLI1', 2.090, 5.71), ('NKX2-2', 4.368, 97.04),
                           ('NR0B1', 4.914, 97.80), ('BCL11B', 4.237, 64.69),
                           ('CDH11', 4.582, 34.31), ('ZEB2', 3.037, 9.77),
                           ('TEAD2', 3.048, 9.31), ('TEAD1', 0.309, 1.25),
                           ('YAP1', -0.145, 0.90), ('WWTR1', -0.235, 0.84),
                           ('PTK2', 1.125, 2.21), ('CD99', 3.176, 9.11),
                           ('CCND1', 3.338, 10.30), ('IGF1R', -0.507, 0.66),
                           ('TRIM8', -0.934, 0.50), ('KLF5', -2.246, 0.13),
                           ('EGFR', -2.788, 0.06), ('MT1E', -3.175, 0.05)]:
    d, fc = rnafc(g)
    check(f'{g} RNA DeltaLog2', round(d, 3), want_d, 0.001)
    check(f'{g} RNA TPM FC', round(fc, 2), want_fc, 0.01)

# ---------- intersection + top-20 order, recomputed ----------
oe = {r['GeneSym']: r for r in csv.DictReader(open(EW + 'overexpression.csv'))}
relaxed = [r for r in scan if passes(r, False)]
cand = [r for r in relaxed if oe[r['Gene'].split(' (')[0]]['OverExpr_pass'] == 'True']
check('relaxed + OE candidates', len(cand), 40)
check('strict + OE candidates', sum(passes(r, True) for r in cand), 28)
def z(vals):
    m, s = statistics.fmean(vals), statistics.pstdev(vals) if len(vals) > 1 else 1
    # pipeline used sample SD (ddof=1); recompute with sample SD
    s = statistics.stdev(vals)
    return [(v - m) / s for v in vals]
ze = z([-float(r['Enrichment']) for r in relaxed])
zf = z([float(r['FracDiff']) for r in relaxed])
zc = z([float(oe[r['Gene'].split(' (')[0]]['DeltaLog2'] or 0) for r in relaxed])
score = {relaxed[i]['Gene'].split(' (')[0]: ze[i] + zf[i] + zc[i] for i in range(len(relaxed))}
ranked = sorted(cand, key=lambda r: -score[r['Gene'].split(' (')[0]])
top20 = list(csv.DictReader(open(EW + 'top20.csv')))
check('top-20 order matches recomputed scores',
      [r['Gene'].split(' (')[0] for r in ranked[:20]],
      [r['GeneSym'] for r in top20])
check('#1 gene', top20[0]['GeneSym'], 'FLI1')
check('#1 score', round(float(top20[0]['Score']), 2), 14.40, 0.01)
check('#2 gene', top20[1]['GeneSym'], 'ZEB2')
check('top20 Tier A count', sum(r['Tier'] == 'A' for r in top20), 17)
check('top20 Tier B count', sum(r['Tier'] == 'B' for r in top20), 3)
lr = [r for r in top20 if r['GeneSym'] == 'LRRC37A3'][0]
check('LRRC37A3 scored lines', int([r for r in scan if r['Gene'].split(' (')[0] == 'LRRC37A3'][0]['Ewing_n']), 19)
check('MYBL2 candidate rank', int([r for r in csv.DictReader(open(EW + 'candidates.csv')) if r['GeneSym'] == 'MYBL2'][0]['CandidateRank']), 37)

# ---------- STAG2 / STAG1 ----------
mut = list(csv.DictReader(open(EW + 'ewing_mutations_dedup.csv')))
stag2_lof = sorted({r['ModelID'] for r in mut if r['HugoSymbol'] == 'STAG2' and r['LikelyLoF'] == 'True'})
check('STAG2 LikelyLoF lines', len(stag2_lof), 6)
fli1 = {r['ModelID']: float(r['FLI1']) for r in csv.DictReader(open(EW + 'fusion_check.csv'))}
stag1 = {r['ModelID']: float(r['STAG1']) for r in csv.DictReader(open(EW + 'fusion_check.csv'))}
check('STAG1-dependent in STAG2-LoF lines', sum(stag1[i] < -0.5 for i in stag2_lof), 6)
check('STAG1-dependent in other lines',
      sum(stag1[i] < -0.5 for i in stag1 if i not in set(stag2_lof)), 9)
check('FLI1-dependent lines', sum(v < -0.5 for v in fli1.values()), 20)
check('CADO-ES1 FLI1 Chronos', round(fli1['ACH-000210'], 3), -0.294, 0.001)
check('ERG-dependent lines', sum(1 for r in csv.DictReader(open(EW + 'fusion_check.csv')) if float(r['ERG']) < -0.5), 0)

# ---------- no patient-RNA / protein leg ----------
studies = json.load(open(EW + 'downloads/cbioportal_studies.json'))
ew_studies = [s['studyId'] for s in studies if 'ewing' in s['name'].lower()]
check('Ewing studies in cBioPortal', sorted(ew_studies),
      ['es_dfarber_broad_2014', 'es_dsrct_msk_2023', 'es_iocurie_2014'])
for s in ew_studies:
    profs = json.load(open(EW + f'downloads/{s}_profiles.json'))
    types = [p['molecularAlterationType'] for p in profs]
    check(f'{s}: no mRNA profile', 'MRNA_EXPRESSION' in types, False)
    check(f'{s}: no protein profile', 'PROTEIN_LEVEL' in types, False)
check('no CPTAC Ewing study', any('cptac' in s['studyId'] and 'ewing' in s['name'].lower() for s in studies), False)
check('no TCGA Ewing study', any(s['studyId'].startswith('sarc_tcga') and 'ewing' in s['name'].lower() for s in studies), False)

print()
if fails:
    print('FAILURES:', fails); sys.exit(1)
print('ALL CHECKS PASSED')
