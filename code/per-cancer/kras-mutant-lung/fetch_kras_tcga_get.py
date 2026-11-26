"""Fetch TCGA LUAD data via cBioPortal GET endpoints only (POST is denied
in this preparation stage). Reuses the EGFR run's sequenced-samples list
(same study, same sample list luad_tcga_sequenced, 230 samples).
"""
import json, re, subprocess
import pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung/downloads"
BASE = "https://www.cbioportal.org/api"

def curl_get(url, outname):
    out = f"{OUT}/{outname}"
    r = subprocess.run(["curl", "-s", "-w", "\n%{http_code}", url, "-H", "Accept: application/json",
                        "--max-time", "300", "-o", out],
                       capture_output=True, text=True)
    code = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "?"
    print(outname, "HTTP", code)
    if code != "200":
        print(open(out).read()[:300])
        return None
    d = json.load(open(out))
    print("  rows:", len(d))
    return d

# mutations already fetched by GET into test_mut_get.json — rename to canonical name
import os
os.rename(f"{OUT}/test_mut_get.json", f"{OUT}/luad_tcga_kras_mutations.json")
try:
    os.remove(f"{OUT}/test_samp_get.json")
except FileNotFoundError:
    pass

scan = pd.read_csv("/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung/dependency_scan.csv")
relaxed = scan[scan["RelaxedBar"]]["Gene"].tolist()
extra = ["TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "PTK2", "PKN2", "EGFR", "ERBB2", "MET", "KRAS", "GRB2",
         "BRAF", "RAF1", "MAPK1", "MAP2K1", "DOCK5", "ILK"]
genes = sorted(set(relaxed) | set(extra))
print("genes to fetch:", len(genes))
with open(f"{OUT}/tcga_fetch_genes.txt", "w") as f: f.write("\n".join(genes))
ge_head = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", nrows=0)
entrez = {}
for c in ge_head.columns[1:]:
    m2 = re.match(r"(.+) \((\d+)\)$", c)
    if m2:
        entrez[m2.group(1)] = int(m2.group(2))
print("unmapped:", sorted(set(genes) - set(entrez)))
with open(f"{OUT}/tcga_entrez.json", "w") as f: json.dump(entrez, f)
eids = [entrez[g] for g in genes if g in entrez]
qs = "&".join(f"entrezGeneId={e}" for e in eids)

rna = curl_get(f"{BASE}/molecular-profiles/luad_tcga_rna_seq_v2_mrna/molecular-data?sampleListId=luad_tcga_rna_seq_v2_mrna&{qs}",
               "luad_tcga_rna.json")
if rna is not None:
    print("  RNA samples:", len({x['sampleId'] for x in rna}), "genes:", len({x['entrezGeneId'] for x in rna}))

rppa = curl_get(f"{BASE}/molecular-profiles/luad_tcga_rppa/molecular-data?sampleListId=luad_tcga_rppa&{qs}",
                "luad_tcga_rppa.json")
if rppa is not None:
    print("  RPPA samples:", len({x['sampleId'] for x in rppa}), "genes:", len({x['entrezGeneId'] for x in rppa}))

cmut = curl_get(f"{BASE}/molecular-profiles/luad_cptac_2020_mutations/mutations?sampleListId=luad_cptac_2020_sequenced&entrezGeneId=3845",
                "luad_cptac_kras_mutations.json")
if cmut is not None:
    df = pd.DataFrame(cmut)
    print("  CPTAC KRAS rows:", len(df))
    if len(df):
        print(df[["sampleId", "proteinChange", "mutationType"]].to_string())

cprot = curl_get(f"{BASE}/molecular-profiles/luad_cptac_2020_protein_quantification/molecular-data?sampleListId=luad_cptac_2020_protein_quantification&{qs}",
                 "luad_cptac_protein.json")
if cprot is not None:
    print("  CPTAC protein samples:", len({x['sampleId'] for x in cprot}), "genes:", len({x['entrezGeneId'] for x in cprot}))
print("done")
