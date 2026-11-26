"""Fetch TCGA LUAD data from cBioPortal (study-scoped endpoints, as in TNBC run).
1. EGFR mutations in luad_tcga (proteinChange per sample) -> activating vs WT.
2. RNA-seq V2 RSEM for relaxed-bar genes (+ TEAD-programme genes) for those samples.
3. RPPA protein for the same genes.
4. CPTAC LUAD (luad_cptac_2020) profiles check + mutation + protein fetch if available.
Saves raw JSON to downloads/.
"""
import json, requests, pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/egfr-mutant-lung/downloads"
BASE = "https://www.cbioportal.org/api"
H = {"Content-Type": "application/json", "Accept": "application/json"}

def post(path, payload, name):
    r = requests.post(BASE + path, headers=H, json=payload, timeout=300)
    print(name, r.status_code, "bytes", len(r.content))
    if r.status_code != 200:
        print(r.text[:500]); return None
    d = r.json()
    with open(f"{OUT}/{name}.json", "w") as f: json.dump(d, f)
    return d

def get(path, name):
    r = requests.get(BASE + path, headers=H, timeout=120)
    print(name, r.status_code, "bytes", len(r.content))
    if r.status_code != 200:
        print(r.text[:500]); return None
    d = r.json()
    with open(f"{OUT}/{name}.json", "w") as f: json.dump(d, f)
    return d

# 1. EGFR mutations, sequenced sample list
mut = post("/molecular-profiles/luad_tcga_mutations/mutations/fetch",
           {"entrezGeneIds": [1956], "sampleListId": "luad_tcga_sequenced"},
           "luad_tcga_egfr_mutations")
if mut is not None:
    print("EGFR mutation rows:", len(mut))
    df = pd.DataFrame(mut)
    print(df[["sampleId", "proteinChange", "mutationType", "mutationStatus"]].to_string() if len(df) else "empty")
    print("distinct samples with EGFR mutation:", df["sampleId"].nunique() if len(df) else 0)

# sample lists for RNA
sl = get("/studies/luad_tcga/sample-lists", "luad_tcga_sample_lists")
# 2. RNA for relaxed-bar genes + programme genes
scan = pd.read_csv("/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/egfr-mutant-lung/dependency_scan.csv")
relaxed = scan[scan["RelaxedBar"]]["Gene"].tolist()
extra = ["TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "PTK2", "PKN2", "EGFR", "ERBB2", "MET", "KRAS", "GRB2"]
genes = sorted(set(relaxed) | set(extra))
print("genes to fetch (RNA/RPPA):", len(genes))
with open(f"{OUT}/tcga_fetch_genes.txt", "w") as f: f.write("\n".join(genes))
# entrez mapping from the DepMap matrix column headers ("GENE (ENTREZ)"),
# no API call needed (cBioPortal /genes/fetch rejects this payload shape here)
import re
ge_head = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", nrows=0)
entrez = {}
for c in ge_head.columns[1:]:
    m2 = re.match(r"(.+) \((\d+)\)$", c)
    if m2:
        entrez[m2.group(1)] = int(m2.group(2))
print("mapped entrez from matrix:", len(entrez), "unmapped:", sorted(set(genes) - set(entrez)))
with open(f"{OUT}/tcga_entrez.json", "w") as f: json.dump(entrez, f)
eids = [entrez[g] for g in genes if g in entrez]
rna = post("/molecular-profiles/luad_tcga_rna_seq_v2_mrna/molecular-data/fetch",
           {"entrezGeneIds": eids, "sampleListId": "luad_tcga_rna_seq_v2_mrna"},
           "luad_tcga_rna")
if rna is not None:
    print("RNA rows:", len(rna), "samples:", len({x['sampleId'] for x in rna}), "genes:", len({x['entrezGeneId'] for x in rna}))

rppa = post("/molecular-profiles/luad_tcga_rppa/molecular-data/fetch",
            {"entrezGeneIds": eids, "sampleListId": "luad_tcga_rppa"},
            "luad_tcga_rppa")
if rppa is not None:
    print("RPPA rows:", len(rppa), "samples:", len({x['sampleId'] for x in rppa}), "genes:", len({x['entrezGeneId'] for x in rppa}))

# 4. CPTAC LUAD profiles
prof = get("/studies/luad_cptac_2020/molecular-profiles", "luad_cptac_2020_profiles")
if prof:
    for x in prof:
        print(x["molecularProfileId"], x.get("datatype"), x.get("molecularAlterationType"))
