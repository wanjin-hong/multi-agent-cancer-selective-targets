"""Fetch TCGA LUAD data from cBioPortal for the KRAS-mutant lung repeat.

1. KRAS mutations in luad_tcga (proteinChange per sample) -> activating vs WT.
   Activating rule (cBioPortal proteinChange, no 'p.' prefix): missense at
   codon G12, G13 or Q61 (regex ^(G12|G13|Q61)[A-Z]); same strict rule as the
   cell-line analysis (excludes G138V and all non-activating variants).
2. RNA-seq V2 RSEM for relaxed-bar genes (+ programme genes).
3. RPPA protein for the same genes.
4. CPTAC LUAD (luad_cptac_2020): KRAS mutations + protein quantification.
Saves raw JSON to downloads/.
"""
import json, re
import requests, pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung/downloads"
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

# 1. KRAS mutations in luad_tcga, sequenced sample list
mut = post("/molecular-profiles/luad_tcga_mutations/mutations/fetch",
           {"entrezGeneIds": [3845], "sampleListId": "luad_tcga_sequenced"},
           "luad_tcga_kras_mutations")
if mut is not None:
    df = pd.DataFrame(mut)
    print("KRAS mutation rows:", len(df))
    print(df[["sampleId", "proteinChange", "mutationType"]].to_string() if len(df) else "empty")
    print("distinct samples with KRAS mutation:", df["sampleId"].nunique() if len(df) else 0)

seq = post("/samples/fetch", {"sampleListIds": ["luad_tcga_sequenced"]},
           "luad_tcga_sequenced_samples")
if seq is not None:
    print("sequenced samples:", len(seq))

# 2. RNA for relaxed-bar genes + programme genes
scan = pd.read_csv("/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung/dependency_scan.csv")
relaxed = scan[scan["RelaxedBar"]]["Gene"].tolist()
extra = ["TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "PTK2", "PKN2", "EGFR", "ERBB2", "MET", "KRAS", "GRB2",
         "BRAF", "RAF1", "MAPK1", "MAP2K1", "DOCK5", "ILK"]
genes = sorted(set(relaxed) | set(extra))
print("genes to fetch (RNA/RPPA/CPTAC):", len(genes))
with open(f"{OUT}/tcga_fetch_genes.txt", "w") as f: f.write("\n".join(genes))
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

# 4. CPTAC LUAD
prof = get("/studies/luad_cptac_2020/molecular-profiles", "luad_cptac_2020_profiles")
if prof:
    for x in prof:
        print(x["molecularProfileId"], x.get("datatype"), x.get("molecularAlterationType"))
sl = get("/studies/luad_cptac_2020/sample-lists", "luad_cptac_2020_sample_lists")
if sl:
    for x in sl:
        print(x.get("sampleListId"), x.get("sampleCount"))

cmut = post("/molecular-profiles/luad_cptac_2020_mutations/mutations/fetch",
            {"entrezGeneIds": [3845], "sampleListId": "luad_cptac_2020_sequenced"},
            "luad_cptac_kras_mutations")
if cmut is None:
    cmut = post("/molecular-profiles/luad_cptac_2020_mutations/mutations/fetch",
                {"entrezGeneIds": [3845]},
                "luad_cptac_kras_mutations")
if cmut is not None:
    df = pd.DataFrame(cmut)
    print("CPTAC KRAS mutation rows:", len(df))
    if len(df):
        print(df[["sampleId", "proteinChange", "mutationType"]].to_string())

cprot = post("/molecular-profiles/luad_cptac_2020_protein_quantification/molecular-data/fetch",
             {"entrezGeneIds": eids, "sampleListId": "luad_cptac_2020_protein_quantification"},
             "luad_cptac_protein")
if cprot is None:
    cprot = post("/molecular-profiles/luad_cptac_2020_protein_quantification/molecular-data/fetch",
                 {"entrezGeneIds": eids, "sampleListId": "luad_cptac_2020_all"},
                 "luad_cptac_protein")
if cprot is not None:
    print("CPTAC protein rows:", len(cprot), "samples:", len({x['sampleId'] for x in cprot}),
          "genes:", len({x['entrezGeneId'] for x in cprot}))
print("done")
