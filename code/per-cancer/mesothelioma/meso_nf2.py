import pandas as pd
mut_path = '/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv'
meso = ['ACH-000086','ACH-000092','ACH-000153','ACH-000319','ACH-000331','ACH-000335','ACH-000367','ACH-000569','ACH-000645','ACH-000648','ACH-001547','ACH-001555','ACH-001556','ACH-001557','ACH-001558','ACH-001559','ACH-001560','ACH-001562','ACH-001977','ACH-001992','ACH-002033']
# stream through the 642MB file, keep only NF2 rows for the 21 lines
chunks = pd.read_csv(mut_path, usecols=['HugoSymbol','ModelID','VariantType','VepImpact','ProteinChange','Hotspot','AMClass','LikelyLoF'], chunksize=200000)
nf2 = []
for ch in chunks:
    sub = ch[(ch['HugoSymbol']=='NF2') & (ch['ModelID'].isin(meso))]
    if len(sub): nf2.append(sub)
nf2 = pd.concat(nf2) if nf2 else pd.DataFrame()
print('NF2 rows for the 21 mesothelioma lines:', len(nf2))
print(nf2.to_string())
print('\nVariantType values:', nf2['VariantType'].value_counts().to_dict() if len(nf2) else {})
print('LikelyLoF values:', nf2['LikelyLoF'].value_counts().to_dict() if len(nf2) else {})
