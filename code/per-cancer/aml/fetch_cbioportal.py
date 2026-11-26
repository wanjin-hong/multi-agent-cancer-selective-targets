"""Fetch cBioPortal data for the AML selective-targets run.

Profiles (study-scoped POST API):
- laml_tcga_rna_seq_v2_mrna (+ same profile in 6 other TCGA studies, cross-cancer
  RNA — THE TCGA leg of the over-expression bar, design (c): TCGA-LAML has no
  adjacent-normal RNA samples — GDC STAR-count RNA files for TCGA-LAML are 151,
  all from primary blood-derived cancer samples — so tumour-vs-normal is
  impossible, same design as melanoma/colon/glioma/ovarian).
- NO RPPA profile exists for any AML study in cBioPortal (laml_tcga,
  aml_ohsu_2018/2022, aml_tcga_gdc, laml_tcga_pan_can_atlas_2018, aml_target,
  aml_stjude_2024 all checked 21 Nov 2026) — no protein bar at all.
- NO CPTAC AML study exists in cBioPortal (no laml_cptac_* study). A CPTAC AML
  proteomics cohort exists in the literature (Pino et al., PMID 38232702, 210
  Beat AML patients) but is not available through cBioPortal, and the CPTAC
  portal API is a dead route (HTTP 500, playbook) — recorded, not fetched.
Fetch set: 463 relaxed-bar genes (union of pan-cancer and haematology
comparators) + check genes = 515 genes.
Uses curl (Python urllib fails with IncompleteRead on chunked responses).
"""
import json, subprocess, os

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/aml/downloads'
BASE = 'https://www.cbioportal.org/api'
entrez_map = json.load(open(os.path.join(DL, 'entrez_map.json')))
eids = sorted(entrez_map.values())
print('fetch genes:', len(eids))

def curl_post(url, payload, outpath):
    with open('/tmp/cb_payload.json', 'w') as f:
        json.dump(payload, f)
    r = subprocess.run(['curl', '-s', '-w', '\nHTTP %{http_code}', '-X', 'POST', url,
                        '-H', 'Content-Type: application/json', '-d', '@/tmp/cb_payload.json',
                        '-o', outpath, '--max-time', '300'], capture_output=True, text=True)
    code = r.stdout.strip().split()[-1]
    print(url.split('/api/')[-1], '->', outpath, 'HTTP', code, 'size', os.path.getsize(outpath))
    return code == '200'

for c in ['laml', 'brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']:
    profile = f'{c}_tcga_rna_seq_v2_mrna'
    ok = curl_post(f'{BASE}/molecular-profiles/{profile}/molecular-data/fetch',
                   {'entrezGeneIds': eids, 'sampleListId': profile},
                   os.path.join(DL, f'tcga_{c}_mrna.json'))
    if ok:
        try:
            d = json.load(open(os.path.join(DL, f'tcga_{c}_mrna.json')))
            print('   rows:', len(d), 'genes:', len(set(x['entrezGeneId'] for x in d)),
                  'samples:', len(set(x['sampleId'] for x in d)))
        except Exception as e:
            print('   parse failed:', e)
print('done')
