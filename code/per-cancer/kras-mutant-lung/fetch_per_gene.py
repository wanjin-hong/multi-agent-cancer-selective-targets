"""Per-gene GET fetches for TCGA LUAD RNA / RPPA / CPTAC protein.

cBioPortal's GET molecular-data endpoint honours only ONE entrezGeneId per
call, so fetch gene by gene (fast, small responses) into per-gene files,
then concatenate. Uses curl (GET only).
"""
import json, os, re, subprocess
import pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung/downloads"
BASE = "https://www.cbioportal.org/api"

genes = open(f"{OUT}/tcga_fetch_genes.txt").read().split()
entrez = json.load(open(f"{OUT}/tcga_entrez.json"))
pairs = [(g, entrez[g]) for g in genes if g in entrez]
print("genes:", len(pairs))

profiles = [
    ("luad_tcga_rna_seq_v2_mrna", "luad_tcga_rna_seq_v2_mrna", "rna"),
    ("luad_tcga_rppa", "luad_tcga_rppa", "rppa"),
    ("luad_cptac_2020_protein_quantification", "luad_cptac_2020_protein_quantification", "cptac"),
]
for profile, sample_list, tag in profiles:
    d = os.path.join(OUT, f"per_gene_{tag}")
    os.makedirs(d, exist_ok=True)
    cmds = []
    for g, eid in pairs:
        url = (f"{BASE}/molecular-profiles/{profile}/molecular-data"
               f"?sampleListId={sample_list}&entrezGeneId={eid}")
        cmds.append(f"curl -s --max-time 60 '{url}' -H 'Accept: application/json' -o {d}/{eid}.json")
    # run 8 at a time
    import concurrent.futures
    def run(cmd):
        subprocess.run(cmd, shell=True, capture_output=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(run, cmds))
    rows = []
    for g, eid in pairs:
        p = f"{d}/{eid}.json"
        try:
            data = json.load(open(p))
        except Exception as e:
            print(tag, g, "BAD JSON", e); continue
        if isinstance(data, list):
            rows.extend(data)
        elif isinstance(data, dict):
            print(tag, g, "dict response:", str(data)[:200])
    with open(f"{OUT}/luad_tcga_{tag}.json" if tag != "cptac" else f"{OUT}/luad_cptac_protein.json", "w") as f:
        json.dump(rows, f)
    print(tag, "total rows:", len(rows), "samples:", len({x['sampleId'] for x in rows}) if rows else 0,
          "genes:", len({x['entrezGeneId'] for x in rows}) if rows else 0)
print("done")
