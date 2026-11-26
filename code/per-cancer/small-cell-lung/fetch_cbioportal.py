"""Fetch cBioPortal patient RNA for the SCLC selective-targets run.

Patient RNA leg (design (c), PRE-REGISTERED before seeing results): SCLC
primary tumours from sclc_ucologne_2015 (George et al., Nature 2015; 81
tumours, RNA-seq RSEM — log2-transformed here) vs pooled primary tumours of
6 other TCGA cancers (BRCA, LUAD, LUSC, COADREAD, STAD, LIHC), cBioPortal
RNA Seq V2 RSEM, log2-transformed here. There is NO TCGA SCLC study and no
tumour-vs-normal design is possible (ucologne has tumours only, no normals),
so design (c) is the only patient leg available — same design as
melanoma/colon/glioma/ovarian.

Protein: NO protein leg exists for SCLC in cBioPortal — no SCLC RPPA profile
and no SCLC CPTAC study (checked 21 Nov 2026: the sclc studies carry
mutations/RNA/SV only). The +1.0 protein bonus therefore cannot fire for any
gene, as in the glioma run. Recorded as an unavailable route, not retried.

Fetch set: 202 relaxed-bar genes + check/headline genes.
Uses curl (Python urllib fails with IncompleteRead).
"""
import json, subprocess, os
import pandas as pd

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/small-cell-lung/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/small-cell-lung/'
BASE = 'https://www.cbioportal.org/api'

scan = pd.read_csv(OUT + 'dependency_scan.csv')
entrez_map = {}
for g in scan['Gene']:
    sym, eid = g.rsplit(' (', 1)
    entrez_map[sym] = int(eid.rstrip(')'))
json.dump(entrez_map, open(os.path.join(DL, 'entrez_map.json'), 'w'))
print('entrez map genes:', len(entrez_map))

relaxed = pd.read_csv(OUT + 'relaxed_bar.csv')
relaxed['Symbol'] = relaxed['Gene'].str.split(' (', regex=False).str[0]
check = ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'GRB2',
         'ITGAV', 'FERMT2', 'ITGB5', 'ITGB1', 'RAC1', 'BRAF', 'MAPK1', 'CCND1', 'TP63', 'SOX2',
         'EGFR', 'PIK3CA', 'BRCA1', 'BRCA2', 'PARP1', 'MYC', 'MYCL', 'MYCN', 'ASCL1', 'NEUROD1',
         'POU2F3', 'INSM1', 'DLL3', 'BCL2', 'EZH2', 'KDM1A', 'RB1', 'TP53', 'SOX10', 'KLF5',
         'FOSL1', 'NEUROG2', 'CHGA', 'SYP', 'FOXM1', 'CCNE1', 'CCNE2', 'E2F3', 'SKP2', 'SLC2A1',
         'NAMPT', 'SMARCA4', 'EP300', 'DOT1L', 'TIMM50', 'RNASEK', 'CXADR', 'CTRB2', 'FOXA2',
         'MIB1', 'TFDP1', 'SDHB', 'SDHA', 'FDPS', 'MVD', 'DHFR', 'PARP2', 'AURKA', 'AURKB',
         'CDK2', 'CDK4', 'CDK6', 'CHEK1', 'WEE1', 'PLK1', 'TOP2A', 'SLFN11', 'BCL2L1', 'MCL1']
fetch_genes = sorted(set(relaxed['Symbol'].tolist()) | set(check))
print('fetch genes:', len(fetch_genes), '(relaxed', len(relaxed), '+ check)')
json.dump(fetch_genes, open(os.path.join(DL, 'fetch_genes.json'), 'w'))
eids = sorted(entrez_map[g] for g in fetch_genes if g in entrez_map)
print('entrez ids:', len(eids), 'missing:', [g for g in fetch_genes if g not in entrez_map])

def curl_post(url, payload, outpath):
    with open('/tmp/cb_payload.json', 'w') as f:
        json.dump(payload, f)
    r = subprocess.run(['curl', '-s', '-w', '\nHTTP %{http_code}', '-X', 'POST', url,
                        '-H', 'Content-Type: application/json', '-d', '@/tmp/cb_payload.json',
                        '-o', outpath, '--max-time', '300'], capture_output=True, text=True)
    code = r.stdout.strip().split()[-1]
    print(url.split('/api/')[-1], '->', outpath, 'HTTP', code, 'size', os.path.getsize(outpath))
    return code == '200'

profiles = [('sclc_ucologne_2015_rna_seq_mrna', 'sclc_ucologne_2015_all', 'tcga_sclc_mrna.json')]
for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']:
    profiles.append((f'{c}_tcga_rna_seq_v2_mrna', f'{c}_tcga_rna_seq_v2_mrna', f'tcga_{c}_mrna.json'))
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
