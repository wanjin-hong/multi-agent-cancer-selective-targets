import pandas as pd, numpy as np
BASE='/home/hatch/workspace/research_notes/'
eff = pd.read_csv(BASE+'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE+'vrk1-codependency-20261108/Model.csv', low_memory=False)
hcc_ids = model.loc[model['OncotreePrimaryDisease']=='Hepatocellular Carcinoma','ModelID'].tolist()
hcc = [x for x in hcc_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(hcc_ids)]
print('shape', eff.shape, 'HCC', len(hcc), 'comp', len(comp), 'missing', int(eff.isna().sum().sum()))
def enrich(col, h, c):
    a = eff.loc[h, col]; b = eff.loc[c, col]
    return a.mean()-b.mean(), (a<-0.5).mean(), (b<-0.5).mean(), len(a.dropna()), (eff[col]<-0.5).mean()
for g in ['GRB2 (2885)','ITGAV (3685)','PTK2 (5747)','FERMT2 (10979)','YAP1 (10413)','WWTR1 (25937)','PKN2 (5586)','TEAD1 (7003)','RAC1 (5879)']:
    e, hf, cf, n, gw = enrich(g, hcc, comp)
    print(f"{g}: enrich {e:.4f} HCCdep {hf:.3f} compdep {cf:.3f} HCCn {n} gw {gw:.3f}")
# SK-HEP-1 sensitivity: exclude ACH-000361
hcc2 = [x for x in hcc if x != 'ACH-000361']
print('SK-HEP-1 excluded, n=', len(hcc2))
for g in ['GRB2 (2885)','ITGAV (3685)','PTK2 (5747)','FERMT2 (10979)','HNF1B (6928)','FGFR1 (2260)']:
    e, hf, cf, n, gw = enrich(g, hcc2, comp)
    print(f"  {g}: enrich {e:.4f} HCCdep {hf:.3f}")
# counts from dependency_scan.csv
s = pd.read_csv('dependency_scan.csv')
print('scan rows', len(s), 'strict count check:', ((s.Enrichment<=-0.10)&(s.HCC_dep_frac>=0.25)&(s.FracDiff>=0.10)&(~s.PanEssential)&(s.GenomeWide_dep_frac<=0.60)).sum())
print('relaxed count check:', ((s.Enrichment<=-0.05)&(s.HCC_dep_frac>=0.20)&(s.FracDiff>=0.10)&(~s.PanEssential)&(s.GenomeWide_dep_frac<=0.60)).sum())
print('pan-essential:', s.PanEssential.sum(), 'q<0.05:', (s.q_BH<0.05).sum())
o = pd.read_csv('overexpression.csv')
print('overexpr rows', len(o), 'strict&overexpr:', ((o.Strict_pass)&(o.OverExpr_pass)).sum(), 'relaxed&overexpr:', ((o.Relaxed_pass)&(o.OverExpr_pass)).sum())
t20 = pd.read_csv('top20.csv')
print('top20 scores match candidates.csv:', np.allclose(t20.Score.values, pd.read_csv('candidates.csv').head(20).Score.values))
# PTK2 TvN recompute from GDC files directly (independent read of 2 genes)
import json, glob, os
files = json.load(open('downloads/gdc_subset_files.json'))
fid2type = {h['file_id']: h['cases'][0]['samples'][0]['sample_type'] for h in files}
for g in ['PTK2','ITGAV','FERMT2']:
    vals={'Primary Tumor':[], 'Solid Tissue Normal':[]}
    for d in glob.glob('downloads/gdc_subset/*/'):
        fid=os.path.basename(d.rstrip('/'))
        if fid not in fid2type: continue
        tsv=glob.glob(d+'*.tsv')
        if not tsv: continue
        df=pd.read_csv(tsv[0], sep='\t', comment='#', usecols=['gene_name','tpm_unstranded'])
        row=df[df.gene_name==g]
        if len(row): vals[fid2type[fid]].append(np.log2(float(row.iloc[0].tpm_unstranded)+1))
    print(f"VERIFY {g}: tumour mean {np.mean(vals['Primary Tumor']):.3f} (n={len(vals['Primary Tumor'])}), normal mean {np.mean(vals['Solid Tissue Normal']):.3f} (n={len(vals['Solid Tissue Normal'])}), log2FC {np.mean(vals['Primary Tumor'])-np.mean(vals['Solid Tissue Normal']):.3f}")
