"""Fetch cBioPortal data for the PDAC selective-targets run.

Profiles (study-scoped API, as in the HCC run):
- paad_tcga_rna_seq_v2_mrna (+ same profile in 6 other TCGA studies, cross-cancer context)
- paad_tcga_rppa (+ same 6 other studies)
- paad_cptac_2021_protein_quantification (CPTAC mass-spec protein; tumour vs normal
  within the same study — never compared across CPTAC studies, per the playbook trap)
- sample types for paad_cptac_2021 to split tumour vs normal

Fetch set: 165 relaxed-bar genes + 16 headline/check genes (TEAD1-4, YAP1, WWTR1,
PTK2, PKN2, KRAS, GRB2, ITGAV, FERMT2, RAC1, BRAF, MAPK1, CCND1).
Uses curl (Python urllib fails with IncompleteRead on chunked responses).
"""
import json, subprocess, os, sys

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/pancreatic/downloads'
BASE = 'https://www.cbioportal.org/api'
fetch_genes = json.load(open(os.path.join(DL, 'fetch_genes.json')))

# --- Entrez map from DepMap column names ---
import pandas as pd
scan = pd.read_csv('/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/pancreatic/dependency_scan.csv')
entrez = {}
for col in scan['Gene']:
    sym, eid = col.rsplit(' (', 1)
    entrez[sym] = int(eid.rstrip(')'))
missing = [g for g in fetch_genes if g not in entrez]
print('fetch genes:', len(fetch_genes), 'missing Entrez:', missing)
entrez_map = {g: entrez[g] for g in fetch_genes if g in entrez}
json.dump(entrez_map, open(os.path.join(DL, 'entrez_map.json'), 'w'))
print('entrez map size:', len(entrez_map))

def curl_post(url, payload, outpath):
    with open('/tmp/cb_payload.json', 'w') as f:
        json.dump(payload, f)
    r = subprocess.run(['curl', '-s', '-w', '\nHTTP %{http_code}', '-X', 'POST', url,
                        '-H', 'Content-Type: application/json', '-d', '@/tmp/cb_payload.json',
                        '-o', outpath], capture_output=True, text=True)
    code = r.stdout.strip().split()[-1]
    print(url.split('/api/')[-1], '->', outpath, 'HTTP', code, 'size', os.path.getsize(outpath))
    if code != '200':
        print('FAILED', r.stdout, r.stderr)
    return code == '200'

def curl_get(url, outpath):
    r = subprocess.run(['curl', '-s', '-w', '\nHTTP %{http_code}', url, '-o', outpath],
                       capture_output=True, text=True)
    code = r.stdout.strip().split()[-1]
    print(url.split('/api/')[-1], '->', outpath, 'HTTP', code, 'size', os.path.getsize(outpath))
    return code == '200'

eids = sorted(entrez_map.values())

# --- CPTAC sample types (tumour vs normal) ---
curl_post(f'{BASE}/samples/fetch',
          {'studyIds': ['paad_cptac_2021'], 'sampleListIds': ['paad_cptac_2021_protein_quantification']},
          os.path.join(DL, 'cptac_paad_samples.json'))
try:
    samps = json.load(open(os.path.join(DL, 'cptac_paad_samples.json')))
    from collections import Counter
    print('CPTAC PAAD protein samples:', Counter(s.get('sampleType') for s in samps))
    print('sample ids example:', samps[0])
except Exception as e:
    print('sample parse failed:', e)

# --- molecular data fetches ---
profiles = [
    ('paad_tcga_rna_seq_v2_mrna', 'paad_tcga_rna_seq_v2_mrna', 'tcga_paad_mrna.json'),
    ('brca_tcga_rna_seq_v2_mrna', 'brca_tcga_rna_seq_v2_mrna', 'tcga_brca_mrna.json'),
    ('luad_tcga_rna_seq_v2_mrna', 'luad_tcga_rna_seq_v2_mrna', 'tcga_luad_mrna.json'),
    ('lusc_tcga_rna_seq_v2_mrna', 'lusc_tcga_rna_seq_v2_mrna', 'tcga_lusc_mrna.json'),
    ('coadread_tcga_rna_seq_v2_mrna', 'coadread_tcga_rna_seq_v2_mrna', 'tcga_coadread_mrna.json'),
    ('stad_tcga_rna_seq_v2_mrna', 'stad_tcga_rna_seq_v2_mrna', 'tcga_stad_mrna.json'),
    ('lihc_tcga_rna_seq_v2_mrna', 'lihc_tcga_rna_seq_v2_mrna', 'tcga_lihc_mrna.json'),
    ('paad_tcga_rppa', 'paad_tcga_rppa', 'tcga_paad_rppa.json'),
    ('brca_tcga_rppa', 'brca_tcga_rppa', 'tcga_brca_rppa.json'),
    ('luad_tcga_rppa', 'luad_tcga_rppa', 'tcga_luad_rppa.json'),
    ('lusc_tcga_rppa', 'lusc_tcga_rppa', 'tcga_lusc_rppa.json'),
    ('coadread_tcga_rppa', 'coadread_tcga_rppa', 'tcga_coadread_rppa.json'),
    ('stad_tcga_rppa', 'stad_tcga_rppa', 'tcga_stad_rppa.json'),
    ('lihc_tcga_rppa', 'lihc_tcga_rppa', 'tcga_lihc_rppa.json'),
    ('paad_cptac_2021_protein_quantification', 'paad_cptac_2021_protein_quantification', 'cptac_paad_protein.json'),
]
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
