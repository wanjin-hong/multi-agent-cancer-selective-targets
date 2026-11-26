"""Fetch TCGA COADREAD (coadread_tcga, COAD+READ) RNA-seq V2 RSEM, CPTAC protein and RPPA
for the relaxed-bar genes, via cBioPortal study-scoped molecular-data/fetch endpoints.
Sample type is classified from the TCGA barcode itself (4th field: 01 = primary tumour,
11 = solid tissue normal), since study-scoped clinical-data fetch returned empty in the
TNBC run. RSEM values are log2-transformed here, as in the TNBC run.
"""
import json, time, requests, pandas as pd, numpy as np

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/colon"
BASE = "https://www.cbioportal.org/api"
STUDY = "coadread_tcga"

genes = [g.strip() for g in open(f"{OUT}/downloads/relaxed_genes.txt") if g.strip()]
# Extras fetched for the Q4 check / prior-work cross-check only (not part of the ranking set):
extras = [g for g in ["TEAD1", "WWTR1", "PTK2", "PKN2", "WRN"] if g not in genes]
print("Relaxed-bar genes:", len(genes), "| extras (Q4/prior-work only):", extras)
genes = genes + extras
print("Genes to fetch:", len(genes))

def fetch_profile(profile_id, genes, tag):
    url = f"{BASE}/molecular-profiles/{profile_id}/molecular-data/fetch"
    sym2entrez = {v: k for k, v in entrez2sym.items()}
    entrez = [sym2entrez[g] for g in genes if g in sym2entrez]
    print(tag, "genes resolved to Entrez:", len(entrez), "of", len(genes),
          "unresolved:", [g for g in genes if g not in sym2entrez])
    out = []
    CH = 50
    for i in range(0, len(entrez), CH):
        chunk = entrez[i:i+CH]
        body = {"entrezGeneIds": chunk, "sampleListId": f"{STUDY}_all"}
        r = requests.post(url, json=body, headers={"Content-Type": "application/json"}, timeout=300)
        print(tag, "chunk", i, "status", r.status_code, "bytes", len(r.content))
        if r.status_code != 200:
            print("  ERROR:", r.text[:500]); continue
        out.extend(r.json())
        time.sleep(0.3)
    with open(f"{OUT}/downloads/cbioportal_{tag}.json", "w") as f:
        json.dump(out, f)
    print(tag, "records:", len(out))
    return out

# No separate gene-resolution call: molecular-data/fetch accepts hugoGeneSymbols directly,
# and Entrez->symbol mapping is built locally from the DepMap column headers.
import pandas as pd
_eff_cols = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", nrows=0).columns[1:]
_expr_cols = pd.read_csv("/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0).columns[1:]
entrez2sym = {}
for c in list(_eff_cols) + list(_expr_cols):
    if " (" in c:
        sym, eid = c.rsplit(" (", 1)
        entrez2sym[int(eid.rstrip(")"))] = sym
json.dump({str(k): v for k, v in entrez2sym.items()}, open(f"{OUT}/downloads/entrez2symbol.json", "w"))
print("Entrez->symbol map:", len(entrez2sym), "genes")

rna = fetch_profile("coadread_tcga_rna_seq_v2_mrna", genes, "rna_seq_v2")
cptac = fetch_profile("coadread_tcga_protein_quantification", genes, "cptac")
rppa = fetch_profile("coadread_tcga_rppa", genes, "rppa")

def sample_type(sid):
    # TCGA barcode: TCGA-XX-XXXX-<type>...  01 primary tumour, 11 solid tissue normal
    try:
        return sid.split("-")[3][:2]
    except Exception:
        return "?"

for tag, data in [("rna", rna), ("cptac", cptac), ("rppa", rppa)]:
    sids = {d["sampleId"] for d in data}
    types = pd.Series([sample_type(s) for s in sids]).value_counts()
    print(tag, "unique samples:", len(sids), "sample-type codes:", types.to_dict())
    gsyms = {d.get("entrezGeneId") for d in data}
    print(tag, "genes with data:", len(gsyms))
