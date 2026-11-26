"""GDC fetch for ESCC tumour-vs-normal (design (b), pre-registered before seeing results).

Primary TCGA leg: ALL ESCA squamous-cell primary tumours (95 STAR-count files)
vs ALL ESCA adjacent normals (13 files), GDC STAR counts, TPM, log2(TPM+1).
TCGA-ESCA mixes squamous (95 tumours) and adenocarcinoma (89); tumours are
filtered to squamous via the GDC cases API (primary_diagnosis containing
'squamous'), reusing the esophagogastric run's file + case dumps so no query
is repeated. No subsampling: the study is small (95 tumours), so all squamous
tumours are used. The 13 normals are fewer than HNSC's 44 but far more than
PDAC's 4 — the axis can discriminate, weakly; the bar is applied exactly as
pre-defined and the normal count is stated in the report.
"""
import json, requests

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophageal-squamous/downloads'
ESGA_DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/downloads'
API = 'https://api.gdc.cancer.gov'

cases = json.load(open(f'{ESGA_DL}/gdc_esca_cases.json'))
sq = set()
for c in cases:
    diags = [d.get('primary_diagnosis', '') for d in c.get('diagnoses', [])]
    if any('squamous' in d.lower() for d in diags):
        sq.add(c['case_id'])
print('ESCA squamous cases:', len(sq))

files = json.load(open(f'{ESGA_DL}/gdc_tcga-esca_files.json'))
tum_sq, norm, tum_other = [], [], []
for h in files:
    st = h['cases'][0]['samples'][0]['sample_type']
    cid = h['cases'][0]['case_id']
    if st == 'Primary Tumor':
        (tum_sq if cid in sq else tum_other).append(h)
    elif st == 'Solid Tissue Normal':
        norm.append(h)
print('ESCA squamous tumours:', len(tum_sq), 'normals:', len(norm), 'excluded non-squamous tumours:', len(tum_other))
subset = tum_sq + norm
ids = [h['file_id'] for h in subset]
json.dump(ids, open(f'{DL}/gdc_subset_ids.json', 'w'))
json.dump(subset, open(f'{DL}/gdc_subset_files.json', 'w'))
print('subset files:', len(subset), 'total size MB:', round(sum(h['file_size'] for h in subset) / 1e6, 1))

# download via /data POST with {"ids": [...]} wrapping
r = requests.post(API + '/data', json={'ids': ids}, timeout=600)
r.raise_for_status()
open(f'{DL}/gdc_subset.tar.gz', 'wb').write(r.content)
print('downloaded bytes:', len(r.content))
