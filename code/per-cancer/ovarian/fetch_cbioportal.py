"""Fetch cBioPortal data for the ovarian selective-targets run.

Profiles (study-scoped POST API):
- ov_tcga_rna_seq_v2_mrna (+ same profile in 6 other TCGA studies, cross-cancer
  RNA — the TCGA leg of the over-expression bar, design (c): TCGA-OV has ZERO
  adjacent normals in GDC, so tumour-vs-normal is impossible and the pooled
  other-cancer comparison is the pre-registered TCGA design, as in
  melanoma/colon/glioma)
- ov_tcga_rppa (+ same 5 other studies; LIHC has no RPPA profile — RPPA is
  the protein bar)
- ov_tcga_protein_quantification (CPTAC mass spec, within-study; sample
  composition checked in process_tcga.py — if no normals, descriptive only)
Fetch set: 94 relaxed-bar genes + 28 headline/check genes = 116 genes
(overlap removed). Uses curl (Python urllib fails with IncompleteRead).
"""
import json, subprocess, os

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/ovarian/downloads'
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

profiles = []
for c in ['ov', 'brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']:
    profiles.append((f'{c}_tcga_rna_seq_v2_mrna', f'{c}_tcga_rna_seq_v2_mrna', f'tcga_{c}_mrna.json'))
for c in ['ov', 'brca', 'luad', 'lusc', 'coadread', 'stad']:
    profiles.append((f'{c}_tcga_rppa', f'{c}_tcga_rppa', f'tcga_{c}_rppa.json'))
profiles.append(('ov_tcga_protein_quantification', 'ov_tcga_protein_quantification', 'tcga_ov_cptac.json'))
for profile, sample_list, fname in profiles:
    ok = curl_post(f'{BASE}/molecular-profiles/{profile}/molecular-data/fetch',
                   {'entrezGeneIds': eids, 'sampleListId': sample_list},
                   os.path.join(DL, fname))
    if ok:
        try:
            d = json.load(open(os.path.join(DL, fname)))
            print('   rows:', len(d), 'genes:', len(set(x['entrezGeneId'] for x in d)),
                  'samples:', len(set(x['sampleId'] for x in d)))
        except Exception as e:
            print('   parse failed:', e)
print('done')
