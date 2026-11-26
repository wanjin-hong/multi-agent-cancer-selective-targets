"""Fetch cBioPortal data for the ESGA selective-targets run.

Profiles (study-scoped POST API):
- stad_tcga_rna_seq_v2_mrna + esca_tcga_rna_seq_v2_mrna (+ same profile in 6
  other TCGA studies, cross-cancer RNA context only — NOT part of the
  over-expression bar; the bar's TCGA leg is the GDC tumour-vs-normal axis)
- stad_tcga_rppa + esca_tcga_rppa (+ same 5 other studies; LIHC has no RPPA
  profile — RPPA is the protein bar; pooled STAD+ESCA-adeno primary tumours
  vs pooled other-cancer primary tumours)
- No CPTAC protein profile exists for STAD or ESCA in cBioPortal (no
  stad_cptac_* / esca_cptac_* study; the *_cptac_gdc studies carry mRNA only)
  — protein rests on RPPA alone.
Fetch set: 241 relaxed-bar genes + 39 headline/check genes = 280 genes.
Uses curl (Python urllib fails with IncompleteRead on chunked responses).
"""
import json, subprocess, os

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/esophagogastric/downloads'
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
for c in ['stad', 'esca', 'brca', 'luad', 'lusc', 'coadread', 'lihc', 'hnsc']:
    profiles.append((f'{c}_tcga_rna_seq_v2_mrna', f'{c}_tcga_rna_seq_v2_mrna', f'tcga_{c}_mrna.json'))
for c in ['stad', 'esca', 'brca', 'luad', 'lusc', 'coadread', 'hnsc']:
    profiles.append((f'{c}_tcga_rppa', f'{c}_tcga_rppa', f'tcga_{c}_rppa.json'))
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
