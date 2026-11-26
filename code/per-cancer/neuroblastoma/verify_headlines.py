"""Independent verification of the neuroblastoma headline numbers.

Separate code path from the pipeline: recomputes every headline number
directly from the raw DepMap/GDC/mutation files (fresh reads, plain
numpy, no pipeline intermediates) and compares against the values in
report.md / the delivered CSVs. Any mismatch is printed as FAIL.
"""
import csv, json, glob, os, math
import numpy as np

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/neuroblastoma/'
DL = OUT + 'downloads/'
fails = []
n_checks = 0

def check(name, got, want, tol=1e-6):
    global n_checks
    n_checks += 1
    if isinstance(want, float) or isinstance(got, float):
        ok = got is not None and want is not None and abs(float(got) - float(want)) <= tol
    else:
        ok = got == want
    if not ok:
        fails.append((name, got, want))
        print(f'FAIL {name}: got {got!r}, want {want!r}')
    return ok

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

# ---- raw loads (independent of pipeline files) ----
model = {}
with open(BASE + 'vrk1-codependency-20261108/Model.csv') as f:
    for r in csv.DictReader(f):
        model[r['ModelID']] = r
nb_ids = [m for m, r in model.items() if r['OncotreePrimaryDisease'] == 'Neuroblastoma']
print('Model.csv NB rows:', len(nb_ids))
check('Model.csv NB rows', len(nb_ids), 51)

eff = {}
eff_genes = None
with open(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv') as f:
    rd = csv.reader(f)
    eff_genes = next(rd)[1:]
    for row in rd:
        eff[row[0]] = row[1:]
check('effect matrix lines', len(eff), 1178)
check('effect matrix genes', len(eff_genes), 17916)
nb_crispr = [i for i in nb_ids if i in eff]
comp_crispr = [i for i in eff if i not in set(nb_ids)]
check('NB CRISPR n', len(nb_crispr), 39)
check('comparator CRISPR n', len(comp_crispr), 1139)
nc_comp = sum(1 for i in comp_crispr if model[i]['OncotreePrimaryDisease'] == 'Non-Cancerous')
check('Non-Cancerous in CRISPR comparator', nc_comp, 14)
# same-patient pair
pids = [model[i]['PatientID'] for i in nb_crispr]
check('same-patient duplicate PatientIDs in CRISPR group', len(pids) - len(set(pids)), 1)

def gene_col(sym):
    for j, g in enumerate(eff_genes):
        if g.startswith(sym + ' ('):
            return j
    return None

def dep_stats(sym):
    j = gene_col(sym)
    a = np.array([float(eff[i][j]) for i in nb_crispr if eff[i][j] != ''], dtype=float)
    b = np.array([float(eff[i][j]) for i in comp_crispr if eff[i][j] != ''], dtype=float)
    return dict(mean=float(a.mean()), n=len(a), dep=int((a < -0.5).sum()),
                comp_mean=float(b.mean()), comp_dep_frac=float((b < -0.5).mean()),
                enrich=float(a.mean() - b.mean()),
                fracdiff=float((a < -0.5).mean() - (b < -0.5).mean()))

# enrichment ranks (recompute whole scan ranks with plain numpy + scipy t-test)
from scipy import stats as st
enrich = {}
qs = {}
ps = []
syms = []
for gi, g in enumerate(eff_genes):
    sym = g.split(' (')[0]
    syms.append(sym)
    a = np.array([float(eff[i][gi]) for i in nb_crispr if eff[i][gi] != ''], dtype=float)
    b = np.array([float(eff[i][gi]) for i in comp_crispr if eff[i][gi] != ''], dtype=float)
    enrich[sym] = float(a.mean() - b.mean())
    ps.append(st.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan)
qvals = bh(np.array(ps, dtype=float))
qmap = dict(zip(syms, qvals))
ranked = sorted(enrich.items(), key=lambda kv: kv[1])
rankmap = {s: i + 1 for i, (s, _) in enumerate(ranked)}

for sym, want in [('ISL1', (1, -0.731895, 25, 1.4e-06)),
                  ('MYCN', (4, -0.635375, 22, 2.5e-04)),
                  ('PHOX2A', (14, -0.415045, 20, 9.6e-05)),
                  ('HAND2', (9, -0.477168, 12, 1.7e-08)),
                  ('GATA3', (6, -0.560327, 20, 2.0e-05)),
                  ('TEAD1', (17381, 0.097322, 3, 0.088)),
                  ('YAP1', (17628, 0.127052, 3, 0.076)),
                  ('WWTR1', (17766, 0.163093, 1, 6.3e-07)),
                  ('PTK2', (6829, -0.009447, 20, 0.93)),
                  ('PKN2', (17784, 0.168163, 6, 0.0025)),
                  ('KLF5', (17844, 0.205362, 0, 2.8e-10)),
                  ('ALK', (105, -0.205751, 6, 0.032)),
                  ('PHOX2B', (28, -0.331309, 4, 2.2e-05)),
                  ('CDK6', (35, -0.309432, 30, 0.017)),
                  ('MT1E', (146, -0.184549, 25, 0.046)),
                  ('MYC', (17916, 1.304408, 11, 2.6e-09))]:
    s = dep_stats(sym)
    check(f'{sym} enrich rank', rankmap[sym], want[0])
    check(f'{sym} enrichment', round(s['enrich'], 6), want[1], tol=1e-5)
    check(f'{sym} NB dep count', s['dep'], want[2])
    check(f'{sym} q', qmap[sym], want[3], tol=max(want[3] * 0.05, 1e-9))
check('q<0.05 total', int((np.array(qvals) < 0.05).sum()), 1846)
check('q<0.05 selective direction', int(sum(1 for s, q in zip(syms, qvals) if q < 0.05 and enrich[s] < 0)), 839)

# bar counts (recomputed from dep_stats-like rules using the scan CSV is NOT
# independent — recompute bars from raw effect values here)
sym_to_col = {g.split(' (')[0]: k for k, g in enumerate(eff_genes)}

def passes(sym, strict):
    j = sym_to_col[sym]
    a = np.array([float(eff[i][j]) for i in nb_crispr if eff[i][j] != ''], dtype=float)
    b = np.array([float(eff[i][j]) for i in comp_crispr if eff[i][j] != ''], dtype=float)
    gw = np.array([float(eff[i][j]) for i in eff if eff[i][j] != ''], dtype=float)
    enr = a.mean() - b.mean()
    hf = (a < -0.5).mean(); cf = (b < -0.5).mean()
    pan = hf >= 0.70 and cf >= 0.70
    gwf = (gw < -0.5).mean()
    if strict:
        return enr <= -0.10 and hf >= 0.25 and (hf - cf) >= 0.10 and not pan and gwf <= 0.60
    return enr <= -0.05 and hf >= 0.20 and (hf - cf) >= 0.10 and not pan and gwf <= 0.60
n_strict = n_relaxed = 0
for s in syms:
    if passes(s, True):
        n_strict += 1
    if passes(s, False):
        n_relaxed += 1
check('strict bar count', n_strict, 149)
check('relaxed bar count', n_relaxed, 290)
# low-n block
col_index = {g: k for k, g in enumerate(eff_genes)}
n18 = sum(1 for g in eff_genes if sum(1 for i in nb_crispr if eff[i][col_index[g]] != '') == 18)
check('genes scored in only 18 NB lines', n18, 129)

# ---- cell-line RNA (independent read) ----
expr = {}
expr_genes = None
with open(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv') as f:
    rd = csv.reader(f)
    expr_genes = next(rd)[1:]
    for row in rd:
        expr[row[0]] = row[1:]
nb_rna = [i for i in nb_ids if i in expr]
comp_rna = [i for i in expr if i not in set(nb_ids)]
check('NB RNA n', len(nb_rna), 36)
check('comparator RNA n', len(comp_rna), 1637)
nc_rna = sum(1 for i in comp_rna if model[i]['OncotreePrimaryDisease'] == 'Non-Cancerous')
check('Non-Cancerous in RNA comparator', nc_rna, 85)

def rna_stats(sym):
    j = next(k for k, g in enumerate(expr_genes) if g.startswith(sym + ' ('))
    a = np.array([float(expr[i][j]) for i in nb_rna], dtype=float)
    b = np.array([float(expr[i][j]) for i in comp_rna], dtype=float)
    am, bm = a.mean(), b.mean()
    fc = (2**am - 1) / (2**bm - 1)
    return dict(dlog2=float(am - bm), fc=float(fc),
                p=float(st.ttest_ind(a, b, equal_var=False)[1]))
for sym, want_fc, want_d in [('MYCN', 287.52, 6.832), ('PHOX2B', 3811.87, 6.072),
                             ('PHOX2A', 2494.98, 6.917), ('HAND2', 569.66, 7.525),
                             ('YAP1', 0.129, -2.588), ('WWTR1', 0.105, -2.635),
                             ('TEAD1', 2.388, 1.202), ('ALK', 46.52, 3.372),
                             ('KLF5', 0.089, -2.568), ('MT1E', 0.10, -2.673)]:
    r = rna_stats(sym)
    check(f'{sym} cell-line RNA FC', round(r['fc'], 2), want_fc, tol=max(want_fc * 0.01, 0.01))
    check(f'{sym} cell-line RNA dlog2', round(r['dlog2'], 3), want_d, tol=0.01)
# MYCN bimodality
mj = next(k for k, g in enumerate(expr_genes) if g.startswith('MYCN ('))
mycn_vals = sorted(((float(expr[i][mj]), model[i]['CellLineName']) for i in nb_rna), reverse=True)
check('MYCN-high lines (>=8.0)', sum(1 for v, _ in mycn_vals if v >= 8.0), 26)
check('MYCN RNA top line', mycn_vals[0][1], 'SK-N-BE(2)')
check('MYCN RNA top value', round(mycn_vals[0][0], 2), 11.40)
check('MYCN RNA gap: 26th value', round(mycn_vals[25][0], 2), 8.14)
check('MYCN RNA 27th value (SK-N-FI)', round(mycn_vals[26][0], 2), 5.74)

# ---- GDC cross-cancer (independent read of the downloaded TSVs) ----
subset = json.load(open(DL + 'gdc_subset_files.json'))
nbl_ids = {h['file_id'] for h in subset['nbl']}
pool_ids = {h['file_id'] for h in subset['pooled']}
check('GDC NBL files', len(nbl_ids), 153)
check('GDC pooled files', len(pool_ids), 150)
data = {}
n_nbl = n_pool = 0
for d in glob.glob(DL + 'gdc_subset/*/'):
    fid = os.path.basename(d.rstrip('/'))
    cohort = 'nbl' if fid in nbl_ids else ('pool' if fid in pool_ids else None)
    if cohort is None:
        continue
    tsv = glob.glob(d + '*.tsv')
    if not tsv:
        continue
    if cohort == 'nbl':
        n_nbl += 1
    else:
        n_pool += 1
    with open(tsv[0]) as f:
        lines = [ln for ln in f if not ln.startswith('#')]
    rd = csv.DictReader(lines, delimiter='\t')
    for row in rd:
        if row['gene_name'] in ('MYCN', 'ISL1', 'YAP1', 'WWTR1', 'TEAD1', 'PHOX2B', 'KLF5', 'MT1E', 'ALK', 'CDK6'):
            data.setdefault(row['gene_name'], {'nbl': [], 'pool': []})
            data[row['gene_name']][cohort].append(math.log2(float(row['tpm_unstranded']) + 1))
check('GDC NBL files read', n_nbl, 153)
check('GDC pooled files read', n_pool, 150)
for sym, want in [('MYCN', 4.99), ('ISL1', 6.87), ('YAP1', -2.35), ('WWTR1', -1.19),
                  ('TEAD1', 1.28), ('PHOX2B', 8.74), ('KLF5', -4.17), ('MT1E', -2.15),
                  ('ALK', 4.67), ('CDK6', 2.30)]:
    a = np.array(data[sym]['nbl']); b = np.array(data[sym]['pool'])
    check(f'{sym} XC log2FC', round(float(a.mean() - b.mean()), 2), want, tol=0.01)

# ---- mutations (independent read of dedup file vs raw file counts) ----
import pandas as pd
mut = pd.read_csv(OUT + 'nb_mutations_dedup.csv')
check('dedup mutation rows', len(mut), 54)
alk = mut[(mut.HugoSymbol == 'ALK') & (mut.ProteinChange.str.contains('F1174|R1275|F1245', na=False))]
check('ALK activating-mutant lines', alk['CellLineName'].nunique(), 9)
act_crispr = [i for i in nb_crispr if model[i]['CellLineName'] in set(alk['CellLineName'])]
check('ALK activating-mutant CRISPR lines', len(act_crispr), 8)
alk_dep_mut = sum(1 for i in act_crispr if float(eff[i][gene_col('ALK')]) < -0.5)
check('ALK-dependent among ALK-mutant', alk_dep_mut, 4)
others = [i for i in nb_crispr if i not in act_crispr]
check('ALK-dependent among other NB lines', sum(1 for i in others if float(eff[i][gene_col('ALK')]) < -0.5), 2)

# ---- MYCN split (independent) ----
hi = [i for i in nb_crispr if i in expr and float(expr[i][mj]) >= 8.0]
lo = [i for i in nb_crispr if i in expr and float(expr[i][mj]) < 8.0]
un = [i for i in nb_crispr if i not in expr]
check('MYCN-high CRISPR lines', len(hi), 23)
check('MYCN-low CRISPR lines', len(lo), 9)
check('MYCN-unclassified CRISPR lines', len(un), 7)
j_mycn = gene_col('MYCN'); j_myc = gene_col('MYC')
check('MYCN dep in high', sum(1 for i in hi if float(eff[i][j_mycn]) < -0.5), 14)
check('MYCN dep in low', sum(1 for i in lo if float(eff[i][j_mycn]) < -0.5), 2)
check('MYC dep in low', sum(1 for i in lo if float(eff[i][j_myc]) < -0.5), 7)
check('MYC dep in high', sum(1 for i in hi if float(eff[i][j_myc]) < -0.5), 3)

# ---- delivered CSVs: top-20 order/scores and candidate counts ----
top = pd.read_csv(OUT + 'top20.csv')
check('top20 order', top['Symbol'].tolist(),
      ['ISL1', 'PHOX2A', 'MYCN', 'HAND2', 'GATA3', 'LDB1', 'SOX11', 'IRS2', 'ASCL1', 'TUBE1',
       'CDIN1', 'CDK6', 'RTTN', 'SMARCA4', 'PPM1D', 'RCOR2', 'CCNL1', 'TRA2B', 'MSI2', 'MDM4'])
check('top20 #1 score', round(float(top.iloc[0]['Score']), 2), 18.98, tol=0.01)
check('top20 #3 score (MYCN)', round(float(top.iloc[2]['Score']), 2), 18.41, tol=0.01)
check('top20 tiers', top['Tier'].tolist(), ['A'] * 8 + ['B'] + ['A'] * 6 + ['B'] + ['A'] * 4)
cand = pd.read_csv(OUT + 'candidates.csv')
check('candidate count', len(cand), 126)
check('strict+OE count', int(((cand['Tier'] == 'A')).sum()), 71)
oe = pd.read_csv(OUT + 'overexpression.csv')
check('overexpression rows', len(oe), 17916)
check('strict+OE in overexpression', int(((oe['Strict_pass']) & (oe['OverExpr_pass'])).sum()), 71)
check('relaxed+OE in overexpression', int(((oe['Relaxed_pass']) & (oe['OverExpr_pass'])).sum()), 126)
check('top20 all NB_n=39', int((top['NB_n'] == 39).all()), 1)
check('dependency_scan rows', len(pd.read_csv(OUT + 'dependency_scan.csv')), 17916)
check('dependency_scan Gene plain symbols', pd.read_csv(OUT + 'dependency_scan.csv')['Gene'].str.contains(' \\(').sum(), 0)

print(f'\n{n_checks} checks, {len(fails)} failures')
if fails:
    print('FAILURES:', fails)
else:
    print('ALL CHECKS PASS')
