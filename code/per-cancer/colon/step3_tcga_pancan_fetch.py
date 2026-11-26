"""Step 3: TCGA PanCancer Atlas RNA — colorectal (coadread) vs all other TCGA cancers.
Profiles: {study}_rna_seq_v2_mrna (RSEM), harmonised PanCancer Atlas pipeline, so the
colorectal group and the comparator come from the same pipeline (unlike mixing Firehose
Legacy coadread with PanCancer Atlas comparators).
Comparator samples: primary tumours only (barcode type 01). Per-gene running sums are
kept for the comparator (sum, sumsq, n) so the full record set need not be held in memory.
Also surveys which PanCancer Atlas studies have a protein_quantification profile.
"""
import json, time, requests, numpy as np, pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/colon"
BASE = "https://www.cbioportal.org/api"
STUDIES = [l.strip() for l in open("/tmp/pca_studies.txt") if l.strip() and "_tcga_pan_can_atlas_2018" in l]
STUDIES = sorted(set(STUDIES))
print("PanCancer Atlas studies:", len(STUDIES))

genes = [g.strip() for g in open(f"{OUT}/downloads/relaxed_genes.txt") if g.strip()]
extras = [g for g in ["TEAD1", "WWTR1", "PTK2", "PKN2", "WRN"] if g not in genes]
genes = genes + extras
entrez2sym = {int(k): v for k, v in json.load(open(f"{OUT}/downloads/entrez2symbol.json")).items()}
sym2entrez = {v: k for k, v in entrez2sym.items()}
entrez = [sym2entrez[g] for g in genes]
print("genes:", len(genes))

def sample_type(sid):
    try: return sid.split("-")[3][:2]
    except Exception: return "?"

# --- survey protein profiles ---
prot_studies = []
for st in STUDIES:
    r = requests.get(f"{BASE}/studies/{st}/molecular-profiles", timeout=60)
    profs = [p["molecularProfileId"] for p in r.json()] if r.status_code == 200 else []
    if f"{st}_protein_quantification" in profs:
        prot_studies.append(st)
    time.sleep(0.1)
print("Studies with protein_quantification profile:", prot_studies)
json.dump(prot_studies, open(f"{OUT}/downloads/pca_protein_studies.json", "w"))

# --- colorectal group: full raw fetch, saved ---
COAD = "coadread_tcga_pan_can_atlas_2018"
def fetch_study(study, entrez_list, tag):
    url = f"{BASE}/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data/fetch"
    out = []
    for i in range(0, len(entrez_list), 50):
        body = {"entrezGeneIds": entrez_list[i:i+50], "sampleListId": f"{study}_all"}
        r = requests.post(url, json=body, headers={"Content-Type": "application/json"}, timeout=300)
        if r.status_code != 200:
            print(tag, study, "chunk", i, "ERROR", r.status_code, r.text[:200]); continue
        out.extend(r.json())
        time.sleep(0.2)
    return out

coad = fetch_study(COAD, entrez, "coad")
json.dump(coad, open(f"{OUT}/downloads/cbioportal_pca_coadread_rna.json", "w"))
coad_s = {d["sampleId"] for d in coad}
print("COADREAD PanCan RNA records:", len(coad), "samples:", len(coad_s),
      "types:", pd.Series([sample_type(s) for s in coad_s]).value_counts().to_dict())

# --- comparators: aggregate on the fly ---
agg = {}   # entrez -> [sum, sumsq, n]
study_counts = {}
for st in STUDIES:
    if st == COAD: continue
    data = fetch_study(st, entrez, "comp")
    sids = {d["sampleId"] for d in data}
    prim = {s for s in sids if sample_type(s) == "01"}
    study_counts[st] = {"samples_with_data": len(sids), "primary_tumour": len(prim),
                        "genes": len({d["entrezGeneId"] for d in data})}
    for d in data:
        if sample_type(d["sampleId"]) != "01": continue
        v = np.log2(d["value"] + 1)
        a = agg.setdefault(d["entrezGeneId"], [0.0, 0.0, 0])
        a[0] += v; a[1] += v * v; a[2] += 1
    print(st, "primary tumours:", len(prim), "genes:", study_counts[st]["genes"])
json.dump({str(k): v for k, v in agg.items()}, open(f"{OUT}/downloads/pca_comparator_agg.json", "w"))
json.dump(study_counts, open(f"{OUT}/downloads/pca_comparator_study_counts.json", "w"), indent=1)
tot_n = sum(v[2] for v in agg.values())
print("Comparator gene-sample observations:", tot_n, "genes aggregated:", len(agg))
print("Comparator primary-tumour samples total:", sum(c["primary_tumour"] for c in study_counts.values()))
