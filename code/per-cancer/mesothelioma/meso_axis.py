import pandas as pd, numpy as np
from scipy.stats import ttest_ind

eff = pd.read_csv('/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv')
mod = pd.read_csv('/home/hatch/workspace/research_notes/vrk1-codependency-20261108/Model.csv')
meso_ids = set(mod.loc[mod['OncotreePrimaryDisease']=='Pleural Mesothelioma','ModelID'])
print('Model.csv mesothelioma lines:', len(meso_ids))
eff_ids = set(eff.iloc[:,0])
print('Effect matrix lines:', len(eff_ids))
group_ids = sorted(meso_ids & eff_ids)
print('Mesothelioma lines WITH CRISPR data:', len(group_ids))
print('Mesothelioma lines WITHOUT CRISPR data:', sorted(meso_ids - eff_ids))
comp_ids = sorted(eff_ids - meso_ids)
print('Comparator (non-mesothelioma) lines:', len(comp_ids))

g = eff[eff.iloc[:,0].isin(group_ids)].set_index(eff.columns[0])
c = eff[eff.iloc[:,0].isin(comp_ids)].set_index(eff.columns[0])

genes = ['TEAD1','TEAD2','TEAD3','TEAD4','YAP1','WWTR1','PTK2','PKN2']
# genome-wide scan for ranks and BH
rows=[]
for col in eff.columns[1:]:
    gv = g[col].dropna(); cv = c[col].dropna()
    enr = gv.mean() - cv.mean()
    rows.append((col.split(' (')[0], gv.mean(), cv.mean(), enr, (gv<-0.5).mean(), (cv<-0.5).mean(), len(gv), len(cv)))
scan = pd.DataFrame(rows, columns=['Gene','GroupMean','CompMean','Enrichment','GroupDepFrac','CompDepFrac','GroupN','CompN'])
scan['FracDiff'] = scan['GroupDepFrac'] - scan['CompDepFrac']
# Welch p for all genes
ps=[]
for col in eff.columns[1:]:
    gv = g[col].dropna(); cv = c[col].dropna()
    if len(gv)<2 or len(cv)<2 or gv.std()==0 and cv.std()==0:
        ps.append(np.nan)
    else:
        ps.append(ttest_ind(gv, cv, equal_var=False).pvalue)
scan['p_welch'] = ps
# BH
p = scan['p_welch'].values
n = len(p); q = np.full(n, np.nan)
valid = ~np.isnan(p)
pv = p[valid]; order = np.argsort(pv); ranked = pv[order]
qvals = ranked * n / (np.arange(len(ranked))+1)
qvals = np.minimum.accumulate(qvals[::-1])[::-1]
q[valid] = qvals[np.argsort(order)]
scan['q_BH'] = np.clip(q, 0, 1)
scan['EnrichRank'] = scan['Enrichment'].rank(method='min').astype(int)
scan = scan.sort_values('Enrichment')
scan.to_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/mesothelioma/tead_axis_full_scan.csv', index=False)
print('scan rows:', len(scan))

sub = scan[scan['Gene'].isin(genes)].set_index('Gene').loc[genes]
print(sub[['GroupMean','CompMean','Enrichment','GroupDepFrac','CompDepFrac','FracDiff','GroupN','p_welch','q_BH','EnrichRank']].to_string(float_format=lambda x: f'{x:.4g}'))
sub.to_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/mesothelioma/tead_axis.csv')

# TAZ-vs-YAP split among TEAD1-dependent mesothelioma lines
def col(gene):
    return [x for x in eff.columns if x.startswith(gene+' (')][0]
t1 = g[col('TEAD1')]; yap = g[col('YAP1')]; taz = g[col('WWTR1')]
t1d = t1 < -0.5
print('\nTEAD1-dependent mesothelioma lines:', int(t1d.sum()), '/', len(t1))
print('  TAZ-dependent among them:', int(((taz < -0.5) & t1d).sum()))
print('  YAP-dependent among them:', int(((yap < -0.5) & t1d).sum()))
print('  both TAZ and YAP:', int(((taz < -0.5) & (yap < -0.5) & t1d).sum()))
print('  neither:', int(((taz >= -0.5) & (yap >= -0.5) & t1d).sum()))
print('  TAZ-only:', int(((taz < -0.5) & (yap >= -0.5) & t1d).sum()))
print('  YAP-only:', int(((yap < -0.5) & (taz >= -0.5) & t1d).sum()))
print('Overall TAZ-dep:', int((taz < -0.5).sum()), 'YAP-dep:', int((yap < -0.5).sum()), 'of', len(t1))
print('\nTEAD1-dependent lines:')
for mid in t1[t1d].index:
    print(' ', mid, 'TEAD1', round(t1[mid],3), 'YAP1', round(yap[mid],3), 'TAZ', round(taz[mid],3))
# pan-essential check for the 8 genes (genome-wide dep frac)
gw = (eff.iloc[:,1:] < -0.5).mean()
for gene in genes:
    cc = col(gene)
    print(gene, 'genome-wide dep frac:', round(float(gw[cc]),4))
# scored-line counts for the 8 genes in group
for gene in genes:
    cc = col(gene)
    print(gene, 'scored in group:', int(g[cc].notna().sum()), '/', len(g))
