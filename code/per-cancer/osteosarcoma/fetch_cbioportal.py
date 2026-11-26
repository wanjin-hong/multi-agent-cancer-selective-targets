"""Fetch cBioPortal data for the osteosarcoma selective-targets run.

Patient-RNA leg (design (c), cross-cancer — pre-registered): there is NO
TCGA osteosarcoma study and NO TCGA-OS adjacent normals, so tumour-vs-
normal is impossible. The tumour set is TARGET osteosarcoma
(os_target_gdc, 159 primary solid tumours) and the comparator is pooled
primary tumours of 6 other cancers (BRCA, LUAD, LUSC, COADREAD, STAD,
LIHC), all on the GDC-harmonised TPM scale (os_target_gdc_mrna_seq_tpm
vs *_tcga_gdc_mrna_seq_tpm) — the same TPM unit on both sides, unlike the
legacy TCGA RSEM profiles, so no cross-scale trap. A scale check (median
log2(TPM+1) of the two sets) is run in process_tcga.py; if it fails, the
leg is rejected as SCLC's was.
Protein: NO protein leg exists. There is no CPTAC osteosarcoma study in
cBioPortal (checked 21 Nov 2026 — CPTAC studies cover brain, breast,
colon, GBM, pancreas, endometrium, lung, H&N, ovary, kidney, uterus
only) and no RPPA profile for os_target_gdc (its profile list is CNA /
mRNA / mutations only). The +1.0 protein bonus therefore cannot fire for
any gene; recorded, not imputed.
Fetch set: 198 relaxed-bar genes + 49 headline/check genes = 247 genes.
Uses curl (Python urllib fails with IncompleteRead).
"""
import json, subprocess, os

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/osteosarcoma/downloads'
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

profiles = [('os_target_gdc_mrna_seq_tpm', 'os_target_gdc_tpm', 'tcga_os_tpm.json')]
for c in ['brca', 'luad', 'lusc', 'coad', 'stad', 'hcc']:
    profiles.append((f'{c}_tcga_gdc_mrna_seq_tpm', f'{c}_tcga_gdc_tpm', f'tcga_{c}_tpm.json'))
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
