import pandas as pd, numpy as np
from scipy import stats

BASE='/home/hatch/workspace/research_notes/'
eff = pd.read_csv(BASE+'vrk1-codependency-20261108/CRISPRGeneEffect.csv', index_col=0)
model = pd.read_csv(BASE+'vrk1-codependency-20261108/Model.csv', low_memory=False)
print('effect shape', eff.shape, 'missing', int(eff.isna().sum().sum()))
# PDAC = Pancreatic Adenocarcinoma only (excludes the single Pancreatic Neuroendocrine Tumor line QGP-1)
pdac_ids = model.loc[model['OncotreePrimaryDisease']=='Pancreatic Adenocarcinoma','ModelID'].tolist()
pdac = [x for x in pdac_ids if x in eff.index]
comp = [x for x in eff.index if x not in set(pdac_ids)]
print('PDAC CRISPR n=', len(pdac), 'comparator n=', len(comp))
print('PDAC lines:', [model.loc[model.ModelID==i,'CellLineName'].iloc[0] for i in pdac])
print('PDAC missing CRISPR:', [model.loc[model.ModelID==i,'CellLineName'].iloc[0] for i in pdac_ids if i not in eff.index])

def bh(p):
    p=np.asarray(p,dtype=float); n=len(p); out=np.full(n,np.nan)
    mask=~np.isnan(p)
    pp=p[mask]; order=np.argsort(pp); ranked=pp[order]
    q=ranked*n/np.arange(1,len(pp)+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    out[mask]=np.minimum(q,1.0)[np.argsort(order)]
    return out

rows=[]
for g in eff.columns:
    a = eff.loc[pdac, g].dropna().values
    b = eff.loc[comp, g].dropna().values
    gw = eff[g].dropna().values
    hm = a.mean() if len(a) else np.nan
    cm = b.mean() if len(b) else np.nan
    hdep = (a < -0.5).sum(); cdep=(b < -0.5).sum()
    hfrac = hdep/len(a) if len(a) else np.nan
    cfrac = cdep/len(b) if len(b) else np.nan
    gwfrac = (gw < -0.5).mean() if len(gw) else np.nan
    if len(a)>=2 and len(b)>=2:
        t,p = stats.ttest_ind(a,b,equal_var=False)
    else:
        p=np.nan
    rows.append(dict(Gene=g, PDAC_mean=hm, PDAC_n=len(a), PDAC_dep=hdep, PDAC_dep_frac=hfrac,
        Comp_mean=cm, Comp_n=len(b), Comp_dep=cdep, Comp_dep_frac=cfrac,
        Enrichment=hm-cm, FracDiff=hfrac-cfrac, GenomeWide_dep_frac=gwfrac,
        PanEssential=bool(hfrac>=0.70 and cfrac>=0.70), p=p))
scan=pd.DataFrame(rows)
scan['q']=bh(scan['p'].values)
scan=scan.sort_values('Enrichment')
scan.to_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/pancreatic/dependency_scan.csv', index=False)
print('genes', len(scan))
print('q<0.05:', int((scan.q<0.05).sum()), ' q<0.1:', int((scan.q<0.1).sum()))
print('q<0.05 selective-direction (Enrichment<0):', int(((scan.q<0.05)&(scan.Enrichment<0)).sum()))
print(scan[scan.q<0.05][['Gene','Enrichment','PDAC_dep_frac','Comp_dep_frac','q']].to_string())
strict = scan[(scan.Enrichment<=-0.10)&(scan.PDAC_dep_frac>=0.25)&(scan.FracDiff>=0.10)&(~scan.PanEssential)&(scan.GenomeWide_dep_frac<=0.60)]
relaxed = scan[(scan.Enrichment<=-0.05)&(scan.PDAC_dep_frac>=0.20)&(scan.FracDiff>=0.10)&(~scan.PanEssential)&(scan.GenomeWide_dep_frac<=0.60)]
print('strict pass:', len(strict), 'relaxed pass:', len(relaxed))
print('pan-essential:', int(scan.PanEssential.sum()))
print('never scored in any PDAC line:', int(scan.PDAC_n.eq(0).sum()))
print(scan.head(25)[['Gene','PDAC_mean','PDAC_dep','PDAC_dep_frac','Comp_mean','Comp_dep_frac','Enrichment','FracDiff','q']].to_string())
strict.to_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/pancreatic/strict_bar.csv', index=False)
relaxed.to_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/pancreatic/relaxed_bar.csv', index=False)
# headline genes
for g in ['TEAD1','TEAD2','TEAD3','TEAD4','YAP1','WWTR1','PTK2','PKN2','KRAS','GRB2','ITGAV','FERMT2','RAC1','BRAF','MAPK1','CCND1']:
    r=scan[scan.Gene.str.startswith(g+' (')]
    print(g, r[['Gene','PDAC_mean','PDAC_dep','PDAC_dep_frac','Comp_dep_frac','Enrichment','GenomeWide_dep_frac','q']].to_string(index=False) if len(r) else 'NOT FOUND')
    if len(r):
        print('  enrichment rank:', int((scan.Enrichment < r.iloc[0].Enrichment).sum())+1)
