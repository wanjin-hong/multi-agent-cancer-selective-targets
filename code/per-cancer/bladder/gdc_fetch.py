"""GDC fetch for bladder tumour-vs-normal (design (b), pre-registered before
seeing results). TCGA-BLCA: 150 random primary tumours (seed 20261121) vs ALL
solid tissue normals, GDC STAR counts, TPM, log2(TPM+1). Same design as
HNSC/HCC/RCC.
"""
import json, random, requests

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/bladder/downloads'
API = 'https://api.gdc.cancer.gov'

q = {
    "filters": {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-BLCA"]}},
        {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "files.data_format", "value": ["tsv"]}},
        {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["RNA-Seq"]}},
        {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor", "Solid Tissue Normal"]}},
    ]},
    "fields": "file_id,file_name,file_size,cases.case_id,cases.submitter_id,cases.samples.sample_type",
    "format": "JSON", "size": 1000,
}
r = requests.post(API + '/files', json=q, timeout=120)
r.raise_for_status()
hits = r.json()['data']['hits']
print('TCGA-BLCA STAR-count files:', len(hits))
json.dump(hits, open(f'{DL}/gdc_blca_files.json', 'w'))
tum = [h for h in hits if h['cases'][0]['samples'][0]['sample_type'] == 'Primary Tumor']
norm = [h for h in hits if h['cases'][0]['samples'][0]['sample_type'] == 'Solid Tissue Normal']
print('tumours:', len(tum), 'normals:', len(norm))
random.seed(20261121)
sample = random.sample(tum, min(150, len(tum)))
subset = sample + norm
print('subset:', len(sample), 'tumours +', len(norm), 'normals =', len(subset))
ids = [h['file_id'] for h in subset]
json.dump(ids, open(f'{DL}/gdc_subset_ids.json', 'w'))
json.dump(subset, open(f'{DL}/gdc_subset_files.json', 'w'))
print('total size MB:', round(sum(h['file_size'] for h in subset) / 1e6, 1))
