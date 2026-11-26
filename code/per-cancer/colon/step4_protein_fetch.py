"""Step 4: CPTAC mass-spec protein + RPPA for the 175 fetched genes, PanCancer Atlas.
Only 3 PanCancer Atlas studies have a protein_quantification profile (surveyed in step 3):
brca, coadread, ov. Protein comparison is therefore COADREAD tumours vs pooled
BRCA + OV tumours — a thin comparator, flagged as supporting-only in the report.
RPPA is fetched for the same 3 studies as the secondary protein source.
"""
import json, time, requests, numpy as np, pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/colon"
BASE = "https://www.cbioportal.org/api"
STUDIES = ["coadread_tcga_pan_can_atlas_2018", "brca_tcga_pan_can_atlas_2018", "ov_tcga_pan_can_atlas_2018"]

genes = [g.strip() for g in open(f"{OUT}/downloads/relaxed_genes.txt") if g.strip()]
extras = [g for g in ["TEAD1", "WWTR1", "PTK2", "PKN2", "WRN"] if g not in genes]
genes = genes + extras
entrez2sym = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/entrez2symbol.json")).items()}
sym2entrez = {v: k for k, v in entrez2sym.items()}
entrez = [sym2entrez[g] for g in genes]

def sample_type(sid):
    try: return sid.split("-")[3][:2]
    except Exception: return "?"

def fetch(study, profile_suffix):
    url = f"{BASE}/molecular-profiles/{study}_{profile_suffix}/molecular-data/fetch"
    out = []
    for i in range(0, len(entrez), 50):
        body = {"entrezGeneIds": entrez[i:i+50], "sampleListId": f"{study}_all"}
        r = requests.post(url, json=body, headers={"Content-Type": "application/json"}, timeout=300)
        if r.status_code != 200:
            print(study, profile_suffix, "chunk", i, "ERROR", r.status_code, r.text[:200]); continue
        out.extend(r.json()); time.sleep(0.2)
    return out

for suffix, tag in [("protein_quantification", "cptac"), ("rppa", "rppa")]:
    all_data = {}
    for st in STUDIES:
        d = fetch(st, suffix)
        all_data[st] = d
        sids = {x["sampleId"] for x in d}
        print(tag, st, "records:", len(d), "samples:", len(sids),
              "types:", pd.Series([sample_type(s) for s in sids]).value_counts().to_dict(),
              "genes:", len({x["entrezGeneId"] for x in d}))
    json.dump(all_data, open(f"{OUT}/downloads/cbioportal_pca_{tag}.json", "w"))
