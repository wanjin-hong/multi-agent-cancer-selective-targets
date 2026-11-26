"""Step 2: TCGA fetch via cBioPortal study-scoped molecular-data/fetch.

TCGA comparison (pre-defined, before any result was seen):
- SKCM tumours (skcm_tcga) vs pooled tumours of the 10 largest other TCGA studies
  by sample count: brca, coadread, gbm, ov, luad, ucec, kirc, hnsc, lgg, thca.
  Reason: SKCM has ZERO adjacent-normal samples in TCGA (480 samples: 111 primary
  + 369 metastatic, 0 normal), so tumour-vs-normal within SKCM is impossible.
  The pooled other-cancer comparator mirrors the cell-line comparator
  (melanoma vs all non-melanoma lines).
- RNA: rna_seq_v2_mrna (RSEM, log2-transformed here as in the TNBC run).
- Protein: RPPA (skcm_tcga_rppa + same 10 studies). CPTAC: no SKCM profile
  exists — attempted and recorded.
- Genes: 197 relaxed-bar genes + YAP1, PKN2 (Q4 extras; TEAD1/WWTR1/PTK2
  already in the relaxed set).
Sample type from TCGA barcode 4th field (01 primary, 06 metastatic, 11 normal).
"""
import json, time, requests, pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/melanoma"
BASE = "https://www.cbioportal.org/api"
SKCM = "skcm_tcga"
OTHERS = ["brca_tcga", "coadread_tcga", "gbm_tcga", "ov_tcga", "luad_tcga",
          "ucec_tcga", "kirc_tcga", "hnsc_tcga", "lgg_tcga", "thca_tcga"]

relax = pd.read_csv(f"{OUT}/_relax_genes.csv")["Gene"].tolist()
extras = [g for g in ["YAP1", "PKN2", "TEAD1", "WWTR1", "PTK2"] if g not in relax]
genes = relax + extras
print("relaxed:", len(relax), "extras:", extras, "total fetch:", len(genes))
with open(f"{OUT}/downloads/fetch_genes.txt", "w") as f:
    f.write("\n".join(genes))

_eff = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", nrows=0).columns[1:]
_expr = pd.read_csv("/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nrows=0).columns[1:]
entrez2sym = {}
for c in list(_eff) + list(_expr):
    if " (" in c:
        sym, eid = c.rsplit(" (", 1)
        entrez2sym[int(eid.rstrip(")"))] = sym
sym2entrez = {v: k for k, v in entrez2sym.items()}
json.dump({str(k): v for k, v in entrez2sym.items()}, open(f"{OUT}/downloads/entrez2symbol.json", "w"))
unres = [g for g in genes if g not in sym2entrez]
print("unresolved to Entrez:", unres)
entrez = [sym2entrez[g] for g in genes if g in sym2entrez]

def profiles(study):
    r = requests.get(f"{BASE}/studies/{study}/molecular-profiles", timeout=60)
    return {p["molecularProfileId"]: p for p in r.json()} if r.status_code == 200 else {}

def fetch_profile(study, profile_id, tag):
    url = f"{BASE}/molecular-profiles/{profile_id}/molecular-data/fetch"
    out = []
    for i in range(0, len(entrez), 50):
        body = {"entrezGeneIds": entrez[i:i+50], "sampleListId": f"{study}_all"}
        r = requests.post(url, json=body, headers={"Content-Type": "application/json"}, timeout=300)
        if r.status_code != 200:
            print(f"  {tag} chunk {i}: HTTP {r.status_code} {r.text[:200]}")
            continue
        out.extend(r.json()); time.sleep(0.2)
    with open(f"{OUT}/downloads/cbioportal_{tag}.json", "w") as f:
        json.dump(out, f)
    sids = {d["sampleId"] for d in out}
    print(f"  {tag}: {len(out)} records, {len(sids)} samples, {len({d['entrezGeneId'] for d in out})} genes")
    return out

def sample_type(sid):
    try:
        return sid.split("-")[3][:2]
    except Exception:
        return "?"

# --- SKCM RNA + RPPA ---
profs = profiles(SKCM)
print("SKCM profiles:", sorted(profs))
assert "skcm_tcga_rna_seq_v2_mrna" in profs
fetch_profile(SKCM, "skcm_tcga_rna_seq_v2_mrna", "skcm_rna_seq_v2")
if "skcm_tcga_rppa" in profs:
    fetch_profile(SKCM, "skcm_tcga_rppa", "skcm_rppa")
# CPTAC attempt (expected: no such profile for SKCM)
if "skcm_tcga_protein_quantification" in profs:
    fetch_profile(SKCM, "skcm_tcga_protein_quantification", "skcm_cptac")
else:
    print("CPTAC: no skcm_tcga_protein_quantification profile exists — CPTAC melanoma unavailable via cBioPortal (recorded, not worked around).")
    with open(f"{OUT}/downloads/cptac_skcm_attempt.txt", "w") as f:
        f.write("No skcm_tcga_protein_quantification molecular profile in cBioPortal study skcm_tcga (profiles listed above). CPTAC protein data for SKCM unavailable.\n")

# --- Comparator studies: RNA + RPPA ---
for study in OTHERS:
    profs = profiles(study)
    rna_id = f"{study}_rna_seq_v2_mrna"
    if rna_id in profs:
        fetch_profile(study, rna_id, f"{study}_rna_seq_v2")
    else:
        print(f"  {study}: NO rna_seq_v2_mrna profile; available: {sorted(profs)}")
    rppa_id = f"{study}_rppa"
    if rppa_id in profs:
        fetch_profile(study, rppa_id, f"{study}_rppa")
    else:
        print(f"  {study}: no RPPA profile")

print("STEP2 FETCH DONE")
