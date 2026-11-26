"""Process cBioPortal TCGA data for AML.

- Cross-cancer mRNA (THE TCGA leg of the over-expression bar, design (c)):
  LAML primary tumours vs pooled primary tumours of 6 other cancers
  (BRCA, LUAD, LUSC, COADREAD, STAD, LIHC), RSEM log2-transformed here.
  This is the pre-registered TCGA design: TCGA-LAML has ZERO adjacent-normal
  RNA samples (GDC STAR-count RNA files for TCGA-LAML: 151, all from primary
  blood-derived cancer samples; the 'Solid Tissue Normal' samples in the case
  facet are DNA germline controls, not RNA) — so tumour-vs-normal is
  impossible, same design as melanoma/colon/glioma/ovarian.
  Scale check included (the SCLC leg was rejected on a scale mismatch):
  median raw RSEM per study is printed before any comparison.
- NO protein leg: no RPPA profile exists for any AML study in cBioPortal and
  no CPTAC AML study exists there (checked 21 Nov 2026). A CPTAC AML
  proteomics cohort exists in the literature (Pino et al., PMID 38232702,
  210 Beat AML patients) but is not available through cBioPortal; the CPTAC
  portal API is a dead route (HTTP 500, playbook). The +1.0 protein bonus
  therefore never fires in this run.
Primary-tumour samples only, by TCGA barcode sample-type 01.
"""
import json, os, numpy as np, pandas as pd
from scipy import stats

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/aml/downloads'
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/aml/'
entrez = json.load(open(os.path.join(DL, 'entrez_map.json')))
inv = {v: k for k, v in entrez.items()}

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p); out = np.full(n, np.nan)
    mask = ~np.isnan(p); pp = p[mask]; order = np.argsort(pp); ranked = pp[order]
    q = ranked * n / np.arange(1, len(pp) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out[mask] = np.minimum(q, 1.0)[np.argsort(order)]
    return out

def load_profile(fname):
    return pd.DataFrame(json.load(open(os.path.join(DL, fname))))

def is_primary_tumour(sid):
    # TCGA barcode sample types: 01 = primary solid tumour, 03 = primary
    # blood derived cancer. LAML tumours are type 03 (peripheral blood /
    # bone marrow), NOT 01 — a 01-only filter returns zero LAML samples.
    return isinstance(sid, str) and len(sid) >= 15 and sid[13:15] in ('01', '03')

# ---------- scale check (SCLC trap) ----------
for c in ['laml', 'brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']:
    d = load_profile(f'tcga_{c}_mrna.json')
    v = d['value'].astype(float)
    print(f'{c}: rows {len(d)} samples {d["sampleId"].nunique()} median raw RSEM {v.median():.1f} mean {v.mean():.1f}')

# ---------- cross-cancer mRNA (the TCGA leg) ----------
aml = load_profile('tcga_laml_mrna.json')
print('LAML mRNA rows:', len(aml), 'samples:', aml['sampleId'].nunique(),
      'primary tumours:', aml[aml['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
others = pd.concat([load_profile(f'tcga_{c}_mrna.json') for c in ['brca', 'luad', 'lusc', 'coadread', 'stad', 'lihc']])
print('other-cancer mRNA rows:', len(others), 'samples:', others['sampleId'].nunique(),
      'primary tumours:', others[others['sampleId'].apply(is_primary_tumour)]['sampleId'].nunique())
aml = aml[aml['sampleId'].apply(is_primary_tumour)].copy()
others = others[others['sampleId'].apply(is_primary_tumour)].copy()
aml['value'] = np.log2(aml['value'].astype(float) + 1)
others['value'] = np.log2(others['value'].astype(float) + 1)
rows = []
for eid, g in inv.items():
    a = aml.loc[aml['entrezGeneId'] == eid, 'value'].values
    b = others.loc[others['entrezGeneId'] == eid, 'value'].values
    am = a.mean() if len(a) else np.nan; bm = b.mean() if len(b) else np.nan
    p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
    rows.append(dict(Gene=g, Entrez=eid, TCGA_RNA_n_AML=len(a), TCGA_RNA_n_other=len(b),
                     TCGA_RNA_mean_AML=am, TCGA_RNA_mean_other=bm,
                     TCGA_RNA_log2FC=am - bm, TCGA_RNA_p=p))
xc = pd.DataFrame(rows)
xc['TCGA_RNA_q'] = bh(xc['TCGA_RNA_p'].values)
xc = xc.sort_values('TCGA_RNA_log2FC', ascending=False)
xc.to_csv(OUT + 'tcga_rna.csv', index=False)
bar = xc[(xc.TCGA_RNA_log2FC >= 0.5) & (xc.TCGA_RNA_p < 0.05)]
print('cross-cancer TCGA bar pass:', len(bar))
print('genes with no LAML TCGA data:', xc[xc.TCGA_RNA_n_AML == 0]['Gene'].tolist())
print(xc.head(15)[['Gene', 'TCGA_RNA_n_AML', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_q']].to_string(index=False))
for g in ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'YAP1', 'WWTR1', 'PTK2', 'PKN2', 'KRAS', 'MYB', 'CBFB', 'SPI1', 'RUNX1', 'CEBPA', 'CDK6', 'FLT3', 'IDH1', 'IDH2', 'NPM1', 'KMT2A', 'BCL2', 'CD33', 'IL3RA', 'KLF5', 'FOSL1', 'HNF1B', 'PARD6B', 'CXADR', 'JAK2', 'MCL1', 'BCL2L1', 'KIT', 'DOT1L', 'MEN1', 'MEIS1', 'HOXA9', 'GATA2', 'ZEB2', 'LMO2', 'CSF2RB']:
    r = xc[xc.Gene == g]
    print(g, r[['TCGA_RNA_log2FC', 'TCGA_RNA_p', 'TCGA_RNA_q']].to_string(index=False) if len(r) else 'NOT FOUND')
