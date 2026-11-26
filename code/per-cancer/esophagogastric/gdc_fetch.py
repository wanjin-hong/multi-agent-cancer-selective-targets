"""GDC fetch for ESGA tumour-vs-normal (design (b), pre-registered before seeing results).

Primary TCGA leg: pooled STAD primary tumours + ESCA adenocarcinoma primary
tumours vs pooled adjacent normals (STAD + ESCA), GDC STAR counts, TPM,
log2(TPM+1). STAD is stomach adenocarcinoma by definition; ESCA is mixed
squamous + adenocarcinoma, so ESCA tumours are filtered to adenocarcinoma via
the GDC cases API (primary_diagnosis containing 'adenocarcinoma').
Subset: 150 random STAD primary tumours (seed 20261121) + ALL ESCA-adeno
primary tumours + ALL normals from both studies.
"""
import json, random, requests

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/downloads'
API = 'https://api.gdc.cancer.gov'

def files_query(project):
    q = {
        "filters": {"op": "and", "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
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
    return r.json()['data']['hits']

all_files = []
for proj in ['TCGA-STAD', 'TCGA-ESCA']:
    hits = files_query(proj)
    print(proj, 'files:', len(hits))
    json.dump(hits, open(f'{DL}/gdc_{proj.lower()}_files.json', 'w'))
    all_files.extend(hits)

# ESCA adenocarcinoma case ids via cases API
cq = {
    "filters": {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-ESCA"]}},
    "fields": "case_id,submitter_id,diagnoses.primary_diagnosis",
    "format": "JSON", "size": 1000,
}
r = requests.post(API + '/cases', json=cq, timeout=120)
r.raise_for_status()
cases = r.json()['data']['hits']
json.dump(cases, open(f'{DL}/gdc_esca_cases.json', 'w'))
adeno = set()
other_diag = {}
for c in cases:
    diags = [d.get('primary_diagnosis', '') for d in c.get('diagnoses', [])]
    if any('adenocarcinoma' in d.lower() for d in diags):
        adeno.add(c['case_id'])
    else:
        other_diag[c['case_id']] = diags
print('ESCA cases:', len(cases), 'adenocarcinoma cases:', len(adeno))
print('ESCA non-adeno diagnoses sample:', list(other_diag.values())[:5])

# classify files
stad_tum, stad_norm, esca_tum_adeno, esca_tum_other, esca_norm = [], [], [], [], []
for h in all_files:
    proj = h['cases'][0]['project']['project_id'] if 'project' in h['cases'][0] else None
    # project not in fields; infer from file list origin via case submitter prefix is unreliable — use per-project lists instead
    pass
# redo classification from the per-project dumps (project known by file)
def classify(hits, project):
    tum, norm, tum_other = [], [], []
    for h in hits:
        st = h['cases'][0]['samples'][0]['sample_type']
        cid = h['cases'][0]['case_id']
        if st == 'Primary Tumor':
            if project == 'TCGA-ESCA' and cid not in adeno:
                tum_other.append(h)
            else:
                tum.append(h)
        elif st == 'Solid Tissue Normal':
            norm.append(h)
    return tum, norm, tum_other

for proj, fn in [('TCGA-STAD', 'gdc_tcga-stad_files.json'), ('TCGA-ESCA', 'gdc_tcga-esca_files.json')]:
    hits = json.load(open(f'{DL}/{fn}'))
    tum, norm, other = classify(hits, proj)
    print(proj, 'tumours:', len(tum), 'normals:', len(norm), 'excluded non-adeno tumours:', len(other))
    if proj == 'TCGA-STAD':
        stad_tum, stad_norm = tum, norm
    else:
        esca_tum_adeno, esca_norm, esca_tum_other = tum, norm, other

random.seed(20261121)
stad_sample = random.sample(stad_tum, min(150, len(stad_tum)))
subset = stad_sample + esca_tum_adeno + stad_norm + esca_norm
print('subset: STAD tumours', len(stad_sample), '+ ESCA-adeno tumours', len(esca_tum_adeno),
      '+ STAD normals', len(stad_norm), '+ ESCA normals', len(esca_norm), '=', len(subset))
ids = [h['file_id'] for h in subset]
json.dump(ids, open(f'{DL}/gdc_subset_ids.json', 'w'))
json.dump(subset, open(f'{DL}/gdc_subset_files.json', 'w'))
print('total size MB:', round(sum(h['file_size'] for h in subset) / 1e6, 1))
