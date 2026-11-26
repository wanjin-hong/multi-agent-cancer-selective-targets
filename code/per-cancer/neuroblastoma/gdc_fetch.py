"""GDC fetch for neuroblastoma patient RNA — pre-registered design (c),
decided before seeing any patient result.

There is no TCGA neuroblastoma study and TARGET-NBL has no adjacent
normals, so design (b) tumour-vs-normal is impossible. Design (c):
TARGET-NBL primary tumours (ALL 153 with STAR counts) vs pooled primary
tumours of 6 other cancers (BRCA, LUAD, LUSC, COADREAD, STAD, LIHC —
25 randomly sampled per study, seed 20261121, total 150), GDC STAR
counts, TPM, log2(TPM+1). Both sides are the same GDC-harmonised
pipeline, so the units are identical — the SCLC ucologne scale trap
does not apply here; a scale check (cohort median log2(TPM+1) across
fetched genes) is run before the bar is applied, and the comparison
is rejected if the medians differ by >1.0 log2 unit.

No protein leg exists: no NBL RPPA or CPTAC profile in cBioPortal
(nbl_target_2018_pub and nbl_target_gdc carry mRNA/CNA/mutations
only, checked 21 Nov 2026) — recorded as unavailable, not retried.
"""
import json, random, requests

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/downloads'
API = 'https://api.gdc.cancer.gov'

def query_files(project, sample_types):
    q = {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "files.data_format", "value": ["tsv"]}},
            {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["RNA-Seq"]}},
            {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": sample_types}},
        ]},
        "fields": "file_id,file_name,file_size,cases.case_id,cases.submitter_id,cases.samples.sample_type",
        "format": "JSON", "size": 2000,
    }
    r = requests.post(API + '/files', json=q, timeout=120)
    r.raise_for_status()
    return r.json()['data']['hits']

nbl = query_files('TARGET-NBL', ['Primary Tumor'])
print('TARGET-NBL primary-tumour STAR files:', len(nbl),
      'total MB:', round(sum(h['file_size'] for h in nbl) / 1e6, 1))
json.dump(nbl, open(f'{DL}/gdc_nbl_files.json', 'w'))

random.seed(20261121)
pooled = []
for proj in ['TCGA-BRCA', 'TCGA-LUAD', 'TCGA-LUSC', 'TCGA-COAD', 'TCGA-STAD', 'TCGA-LIHC']:
    hits = query_files(proj, ['Primary Tumor'])
    tum = [h for h in hits if h['cases'][0]['samples'][0]['sample_type'] == 'Primary Tumor']
    sample = random.sample(tum, 25)
    for h in sample:
        h = dict(h)
        h['project_id'] = proj
        pooled.append(h)
    print(proj, 'tumours available:', len(tum), 'sampled 25')
json.dump(pooled, open(f'{DL}/gdc_pooled_files.json', 'w'))
json.dump({'nbl': nbl, 'pooled': pooled}, open(f'{DL}/gdc_subset_files.json', 'w'))
print('subset total files:', len(nbl) + len(pooled),
      'total MB:', round((sum(h['file_size'] for h in nbl) + sum(h['file_size'] for h in pooled)) / 1e6, 1))
