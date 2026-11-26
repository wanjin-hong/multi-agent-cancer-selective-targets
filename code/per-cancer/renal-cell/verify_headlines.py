"""Independent verification of the RCC headline numbers, separate code path
from the pipeline: recompute from the RAW files (CRISPRGeneEffect.csv,
expression matrix, GDC STAR-count TSVs, cBioPortal mRNA/RPPA JSON, mutation
file) and compare against dependency_scan.csv / overexpression.csv / top20.csv.
Prints MISMATCH for any disagreement; exits 0 only if all checks pass.
"""
import json, glob, os, sys, numpy as np, pandas as pd

BASE = '/home/hatch/workspace/research_notes/'
OUT = BASE + 'multi-cancer-selective-targets-20261120/renal-cell/'
DL = OUT + 'downloads/'
fails = []

def check(name, got, want, tol=1e-6):
    if isinstance(want, (list, tuple)):
        ok = (list(got) == list(want))
    elif isinstance(want, (int, str, bool)):
        ok = (got == want)
    else:
        ok = (abs(got - want) <= tol * max(1, abs(want)))
    print(('OK   ' if ok else 'MISMATCH'), name, 'got', got, 'want', want)
    if not ok:
        fails.append(name)

# --- raw group counts ---
eff = pd.read_csv(BASE + 'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
expr = pd.read_csv(BASE + 'glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv', index_col=0)
model = pd.read_csv(BASE + 'vrk1-codependency-20261108/Model.csv', low_memory=False)
rcc_ids = model.loc[model['OncotreePrimaryDisease'] == 'Renal Cell Carcinoma', 'ModelID'].tolist()
rc_c = [x for x in rcc_ids if x in eff.index]
rc_r = [x for x in rcc_ids if x in expr.index]
comp_c = [x for x in eff.index if x not in set(rcc_ids)]
check('RCC CRISPR n', len(rc_c), 30)
check('RCC RNA n', len(rc_r), 37)
check('RCC lines in Model', len(rcc_ids), 57)
check('comparator CRISPR n', len(comp_c), 1148)

scan = pd.read_csv(OUT + 'dependency_scan.csv')
scan = scan.set_index('Gene')
over = pd.read_csv(OUT + 'overexpression.csv')
over = over.set_index('Symbol')
top20 = pd.read_csv(OUT + 'top20.csv')
cand = pd.read_csv(OUT + 'candidates.csv')

def col(gene):
    return [c for c in eff.columns if c.startswith(gene + ' (')][0]

# --- dependency headline numbers, recomputed from the raw matrix ---
for gene in ['HNF1B', 'PAX8', 'CFLAR', 'CDH2', 'MDM2', 'ITGAV', 'WWTR1', 'FERMT2',
             'PARD6B', 'UBC', 'PARD3', 'EPAS1', 'TEAD1', 'YAP1', 'PTK2', 'PKN2',
             'VHL', 'HIF1A', 'EGLN1', 'KLF5', 'CDK6', 'FGFR1', 'GRB2']:
    c = col(gene)
    a = eff.loc[rc_c, c].dropna()
    b = eff.loc[comp_c, c].dropna()
    row = scan.loc[gene]
    check(f'{gene} RCC mean', round(a.mean(), 6), round(row['RCC_mean'], 6), 1e-4)
    check(f'{gene} RCC dep', int((a < -0.5).sum()), int(row['RCC_dep']))
    check(f'{gene} RCC n', len(a), int(row['RCC_n']))
    check(f'{gene} enrichment', round(a.mean() - b.mean(), 6), round(row['Enrichment'], 6), 1e-4)
    check(f'{gene} enrich rank', int((scan['Enrichment'] < row['Enrichment']).sum()) + 1, int(row['EnrichRank']))

# --- bar counts ---
check('strict bar count', int((over['Strict_pass']).sum()), 122)
check('relaxed bar count', int((over['Relaxed_pass']).sum()), 228)
check('strict & OE count', int(((over['Strict_pass']) & (over['OverExpr_pass'])).sum()), 34)
check('relaxed & OE count (candidates)', len(cand), 54)
check('q<0.05 count', int((scan['q_BH'] < 0.05).sum()), 480)
check('q<0.05 selective count', int(((scan['q_BH'] < 0.05) & (scan['Enrichment'] < 0)).sum()), 133)
check('pan-essential count', int(scan['PanEssential'].sum()), 1280)
check('top20 order', top20['Symbol'].tolist(),
      ['HNF1B', 'PAX8', 'CFLAR', 'CDH2', 'MDM2', 'ITGAV', 'WWTR1', 'FERMT2', 'PARD6B', 'UBC',
       'PARD3', 'ZFP36L1', 'MZT2A', 'BRAT1', 'KIF2C', 'EPAS1', 'TUBA3D', 'FAF2', 'UBA5', 'UFM1'])
check('ranks 21-25', cand.iloc[20:25]['Symbol'].tolist(), ['CDC25B', 'FGFR1', 'TRAPPC2L', 'ATF4', 'CDK6'])

# --- cell-line RNA, recomputed from the raw matrix ---
def rna_col(gene):
    return [c for c in expr.columns if c.split(' (')[0] == gene][0]
for gene, want_fc in [('HNF1B', 48.239334), ('PAX8', 34.611746), ('WWTR1', 7.461765),
                      ('TEAD1', 1.548576), ('PTK2', 1.464363), ('EPAS1', 1.819193),
                      ('YAP1', 1.516722), ('FERMT2', 1.701139), ('EGFR', 5.562995),
                      ('CA9', 5.410638), ('VHL', 0.519724), ('HIF1A', 0.509992)]:
    c = rna_col(gene)
    a = expr.loc[rc_r, c].dropna().values
    b = expr.loc[[x for x in expr.index if x not in set(rcc_ids)], c].dropna().values
    fc = (2 ** a.mean() - 1) / (2 ** b.mean() - 1)
    row = over.loc[gene]
    check(f'{gene} cell-line TPM FC', round(fc, 4), round(row['TPM_FC'], 4), 1e-3)
    check(f'{gene} cell-line DeltaLog2', round(a.mean() - b.mean(), 4), round(row['DeltaLog2'], 4), 1e-3)

# --- GDC tumour-vs-normal, recomputed from the raw TSVs ---
files = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
tvn_raw = {}
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
    for g, v in zip(df['gene_name'], np.log2(df['tpm_unstranded'].astype(float) + 1)):
        tvn_raw.setdefault(g, {'tumour': [], 'normal': []})
        tvn_raw[g]['tumour' if st == 'Primary Tumor' else 'normal'].append(v)
check('GDC tumour files', n_tum, 150)
check('GDC normal files', n_norm, 72)
for gene in ['CA9', 'NDUFA4L2', 'EPAS1', 'WWTR1', 'HNF1B', 'PAX8', 'TEAD1', 'YAP1',
             'PTK2', 'VHL', 'HIF1A', 'EGFR', 'MET', 'FERMT2', 'PARD6B', 'KLF5', 'TEAD4']:
    a = np.array(tvn_raw[gene]['tumour'])
    b = np.array(tvn_raw[gene]['normal'])
    row = over.loc[gene]
    check(f'{gene} TvN log2FC', round(a.mean() - b.mean(), 4), round(row['GDC_TvN_log2FC'], 4), 1e-3)
    check(f'{gene} TvN n_tumour', len(a), int(row['GDC_TvN_n_tumour']))
    check(f'{gene} TvN n_normal', len(b), int(row['GDC_TvN_n_normal']))

# --- cross-cancer RNA + RPPA, recomputed from the raw cBioPortal JSONs ---
def load_profile(fname):
    return pd.DataFrame(json.load(open(os.path.join(DL, fname))))
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
inv = {v: k for k, v in entrez.items()}
def is_primary(sid):
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] == '01'
kirc = load_profile('tcga_kirc_mrna.json')
others = pd.concat([load_profile(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
kirc = kirc[kirc['sampleId'].apply(is_primary)].copy()
others = others[others['sampleId'].apply(is_primary)].copy()
kirc['value'] = np.log2(kirc['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
for gene in ['CA9', 'PAX8', 'HNF1B', 'WWTR1', 'EPAS1', 'VEGFA', 'MET', 'NDUFA4L2']:
    eid = entrez[gene]
    a = kirc.loc[kirc['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    row = over.loc[gene]
    check(f'{gene} cross-cancer log2FC', round(a.mean() - b.mean(), 4), round(row['TCGA_RNA_log2FC'], 4), 1e-3)
kirc_r = load_profile('tcga_kirc_rppa.json')
other_r = pd.concat([load_profile(f'tcga_{c}_rppa.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad']])
kirc_r = kirc_r[kirc_r['sampleId'].apply(is_primary)].copy()
other_r = other_r[other_r['sampleId'].apply(is_primary)].copy()
for gene in ['WWTR1', 'CDH2', 'CCNE1', 'YAP1', 'EGFR', 'MTOR']:
    eid = entrez[gene]
    a = kirc_r.loc[kirc_r['entrezGeneId'] == eid, 'value'].values
    b = other_r.loc[other_r['entrezGeneId'] == eid, 'value'].values
    row = over.loc[gene]
    if len(a) and len(b):
        check(f'{gene} RPPA delta', round(a.mean() - b.mean(), 4), round(row['RPPA_delta'], 4), 1e-3)
    else:
        check(f'{gene} RPPA n_KIRC', len(a), int(row['RPPA_n_KIRC']))

# --- VHL mutation count from the raw mutation file (deduped) ---
cols = "idx,SequencingID,ModelID,ModelConditionID,IsDefaultEntryForModel,IsDefaultEntryForMC,Chrom,Pos,Ref,Alt,AF,DP,RefCount,AltCount,GT,PS,VariantType,VariantInfo,DNAChange,ProteinChange,HugoSymbol,Exon,Intron,EnsemblGeneID,EnsemblFeatureID,HgncName,HgncFamily,UniprotID,DbsnpRsID,GcContent,VepImpact,VepBiotype,VepHgncID,VepExistingVariation,VepManeSelect,VepENSP,VepSwissprot,Sift,Polyphen,GnomadeAF,GnomadgAF,VepClinSig,VepSomatic,VepPliGeneValue,VepLofTool,OncogeneHighImpact,TumorSuppressorHighImpact,CivicID,CivicDescription,CivicScore,LikelyLoF,AMClass,AMPathogenicity,Rescue,RescueReason,Hotspot,EntrezGeneID".split(',')
mut = pd.read_csv('/tmp/rcc_mutations.csv', low_memory=False, header=None, names=cols)
mut = mut.drop_duplicates(subset=['ModelID', 'HugoSymbol', 'ProteinChange', 'VariantType'])
vhl = mut[(mut.HugoSymbol == 'VHL') & ((mut.LikelyLoF == True) | (mut.VepImpact.isin(['HIGH', 'MODERATE'])) | (mut.Hotspot == True))]
check('VHL damaging-variant lines', len(set(vhl.ModelID)), 23)
hif1a = mut[(mut.HugoSymbol == 'HIF1A') & ((mut.LikelyLoF == True) | (mut.VepImpact.isin(['HIGH', 'MODERATE'])) | (mut.Hotspot == True))]
check('HIF1A LoF lines', len(set(hif1a.ModelID)), 5)
pbrm1 = mut[(mut.HugoSymbol == 'PBRM1') & ((mut.LikelyLoF == True) | (mut.VepImpact.isin(['HIGH', 'MODERATE'])) | (mut.Hotspot == True))]
check('PBRM1 damaging-variant lines', len(set(pbrm1.ModelID)), 14)

print()
print('FAILS:', fails if fails else 'none — all checks pass')
sys.exit(1 if fails else 0)
