"""Independent verification of the PDAC headline numbers.

Recomputes from the raw files with a separate code path from the pipeline:
- dependency means / dep counts from CRISPRGeneEffect.csv
- cell-line RNA delta / TPM FC from the expression matrix (playbook formula)
- TCGA tumour-vs-normal log2FC from the downloaded GDC STAR-count files
- cross-cancer log2FC and RPPA delta from the cBioPortal JSON downloads
- KRAS-mutant-only (42/45 lines) sensitivity check on the top enrichments
- which lines the 24-line genes (SPRR1B, CXADR, GGT1) were scored in
"""
import pandas as pd, numpy as np, json, glob, os, csv
from scipy import stats

BASE = '/home/hatch/workspace/research_notes/'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/pancreatic/'
DL = OUT + 'downloads/'
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
pdac_ids = model.loc[model['OncotreePrimaryDisease'] == 'Pancreatic Adenocarcinoma', 'ModelID'].tolist()
pdac = [x for x in pdac_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(pdac_ids)]
name = dict(zip(model.ModelID, model.CellLineName))

scan = pd.read_csv(OUT + 'dependency_scan.csv')
top20 = pd.read_csv(OUT + 'candidates.csv').head(20)
check = list(top20['Symbol']) + ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2',
                                  'KRAS', 'GRB2', 'ITGAV', 'FERMT2', 'RAC1', 'BRAF', 'MAPK1', 'CCND1']
print('== dependency re-check (raw CRISPRGeneEffect.csv) ==')
bad = 0
for g in check:
    col = scan.loc[scan['Gene'].str.startswith(g + ' ('), 'Gene'].iloc[0]
    a = eff.loc[pdac, col].dropna()
    b = eff.loc[comp, col].dropna()
    row = scan[scan.Gene == col].iloc[0]
    ok = (abs(a.mean() - row.PDAC_mean) < 1e-9 and (a < -0.5).sum() == row.PDAC_dep
          and abs((a < -0.5).mean() - row.PDAC_dep_frac) < 1e-9
          and abs(a.mean() - b.mean() - row.Enrichment) < 1e-9)
    if not ok:
        bad += 1
        print('MISMATCH', g, a.mean(), row.PDAC_mean, (a < -0.5).sum(), row.PDAC_dep)
print('dependency mismatches:', bad, 'of', len(check))

print('== cell-line RNA re-check (raw expression matrix) ==')
bad = 0
rna = pd.read_csv(OUT + 'celline_rna.csv')
pdac_r = [x for x in pdac_ids if x in expr.index]
comp_r = [x for x in expr.index if x not in set(pdac_ids)]
for g in check:
    col = [c for c in expr.columns if c.startswith(g + ' (')][0]
    a = expr.loc[pdac_r, col].dropna(); b = expr.loc[comp_r, col].dropna()
    fc = (2 ** a.mean() - 1) / (2 ** b.mean() - 1)
    row = rna[rna.Gene == g].iloc[0]
    if abs(a.mean() - b.mean() - row.DeltaLog2) > 1e-9 or abs(fc - row.TPM_FC) > 1e-9:
        bad += 1; print('MISMATCH', g, a.mean() - b.mean(), row.DeltaLog2, fc, row.TPM_FC)
print('RNA mismatches:', bad, 'of', len(check))
print('PDAC RNA lines in matrix:', len(pdac_r), ' (one of 52 has all-NaN rows?)')
print('lines with RNA but missing for KLF5:', [name[i] for i in pdac_r if pd.isna(expr.loc[i, [c for c in expr.columns if c.startswith("KLF5 (")][0]])])

print('== TCGA tumour-vs-normal re-check (raw GDC files, KLF5/FOSL1/SPRR1B/WWTR1/CCND1) ==')
files = json.load(open(DL + 'gdc_paad_files.json'))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
for g in ['KLF5', 'FOSL1', 'SPRR1B', 'WWTR1', 'CCND1', 'TEAD1', 'KRAS']:
    t, n = [], []
    for d in glob.glob(DL + 'gdc_subset/*/'):
        fid = os.path.basename(d.rstrip('/'))
        if fid not in fid2type: continue
        tsv = glob.glob(d + '*.tsv')
        if not tsv: continue
        with open(tsv[0]) as f:
            lines = [ln for ln in f if not ln.startswith('#')]
        rd = csv.DictReader(lines, delimiter='\t')
        for row in rd:
            if row['gene_name'] == g:
                v = np.log2(float(row['tpm_unstranded']) + 1)
                (t if fid2type[fid] == 'Primary Tumor' else n).append(v)
                break
    tvn = pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv')
    row = tvn[tvn.Gene == g].iloc[0]
    print(g, 'tum n=%d norm n=%d log2FC=%.4f (file %.4f) p=%.4g (file %.4g)' %
          (len(t), len(n), np.mean(t) - np.mean(n), row.GDC_TvN_log2FC,
           stats.ttest_ind(t, n, equal_var=False)[1], row.GDC_TvN_p))

print('== cross-cancer + RPPA re-check (cBioPortal JSONs) ==')
entrez = json.load(open(DL + 'entrez_map.json'))
xc = pd.read_csv(OUT + 'tcga_rna.csv')
rp = pd.read_csv(OUT + 'tcga_rppa.csv')
def is_pt(sid): return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
for g in ['KLF5', 'WWTR1', 'HNF1B', 'ITGAV', 'KRAS']:
    eid = entrez[g]
    paad = json.load(open(DL + 'tcga_paad_mrna.json'))
    a = [np.log2(x['value'] + 1) for x in paad if x['entrezGeneId'] == eid and is_pt(x['sampleId'])]
    b = []
    for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']:
        d = json.load(open(DL + f'tcga_{c}_mrna.json'))
        b += [np.log2(x['value'] + 1) for x in d if x['entrezGeneId'] == eid and is_pt(x['sampleId'])]
    row = xc[xc.Gene == g].iloc[0]
    print(g, 'XC n=%d/%d log2FC=%.4f (file %.4f)' % (len(a), len(b), np.mean(a) - np.mean(b), row.TCGA_RNA_log2FC))
for g in ['WWTR1', 'YAP1', 'PXN']:
    eid = entrez[g]
    paad = json.load(open(DL + 'tcga_paad_rppa.json'))
    a = [x['value'] for x in paad if x['entrezGeneId'] == eid and is_pt(x['sampleId'])]
    b = []
    for c, prof in [('brca', 'brca'), ('luad', 'luad'), ('lusc', 'lusc'), ('coadread', 'coadread'),
                    ('stad', 'stad'), ('lihc', 'lihc')]:
        d = json.load(open(DL + f'tcga_{c}_rppa.json'))
        b += [x['value'] for x in d if x['entrezGeneId'] == eid and is_pt(x['sampleId'])]
    row = rp[rp.Gene == g].iloc[0]
    print(g, 'RPPA n=%d/%d delta=%.4f (file %.4f)' % (len(a), len(b), np.mean(a) - np.mean(b), row.RPPA_delta))

print('== KRAS-mutant-only sensitivity (42 lines) ==')
# KRAS hotspot-mutant lines from the earlier mutation-file check
wt_names = {'BxPC-3', 'PACADD-137', 'PACADD-165'}
kras_mut = [i for i in pdac if name[i] not in wt_names]
print('KRAS-mutant CRISPR n =', len(kras_mut))
for g in list(top20['Symbol']) + ['WWTR1', 'TEAD1', 'PKN2', 'PTK2', 'KRAS']:
    col = scan.loc[scan['Gene'].str.startswith(g + ' ('), 'Gene'].iloc[0]
    a = eff.loc[kras_mut, col].dropna(); b = eff.loc[comp, col].dropna()
    row = scan[scan.Gene == col].iloc[0]
    print('%-8s all45 enrich %+.3f (%d/%d dep) | KRAS-mut42 enrich %+.3f (%d/%d dep)' %
          (g, row.Enrichment, row.PDAC_dep, row.PDAC_n, a.mean() - b.mean(), (a < -0.5).sum(), len(a)))

print('== 24-line genes: which lines were scored? ==')
for g in ['SPRR1B', 'CXADR', 'GGT1']:
    col = scan.loc[scan['Gene'].str.startswith(g + ' ('), 'Gene'].iloc[0]
    scored = [i for i in pdac if pd.notna(eff.loc[i, col])]
    print(g, 'scored in', len(scored), 'lines:', sorted(name[i] for i in scored))
print('genes scored in <45 PDAC lines:', int((scan.PDAC_n < 45).sum()),
      '| scored in <30:', int((scan.PDAC_n < 30).sum()),
      '| scored in exactly 24:', int((scan.PDAC_n == 24).sum()))
print('relaxed genes with no TvN data:', scan.merge(pd.read_csv(OUT + 'tcga_tumour_vs_normal.csv'), left_on=scan['Gene'].str.split(' (', regex=False).str[0], right_on='Gene', how='left') is not None)
ov = pd.read_csv(OUT + 'overexpression.csv')
rel = ov[ov.Relaxed_pass]
print('relaxed, no TvN:', rel[rel.GDC_TvN_log2FC.isna()]['Symbol'].tolist())
print('relaxed, no cell-line RNA:', rel[rel.DeltaLog2.isna()]['Symbol'].tolist())
print('relaxed, no CPTAC:', rel[rel.CPTAC_n_tumour.isna() | (rel.CPTAC_n_tumour == 0)]['Symbol'].tolist())
