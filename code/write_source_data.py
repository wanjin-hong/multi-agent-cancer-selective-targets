#!/usr/bin/env python3
"""write_source_data.py — per-figure Source Data files for the manuscript
(26 Nov 2026). Nothing here edits manuscript.md, figures or tables.

Reads the same verified source CSVs that build_figs_tables.py uses and writes
one Source Data CSV per data figure into manuscript/source_data/:

  Source_Data_Fig2.csv  — 20-cancer x 269-gene top-20 matrix (combined scores)
  Source_Data_Fig3.csv  — recurrence counts (25 genes >=3 cancers + singletons)
  Source_Data_Fig4.csv  — TEAD-axis enrichment ranks (9 genes x 21 settings)
  Source_Data_Fig5.csv  — tau + switch-like values, 269 selective + 300 background
                          (CORRECTED true-CPM summaries — NOT the inflated run)

Fig 1 is a schematic and has no source data (Data Availability statement).
"""
import csv, os, statistics
from collections import Counter

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
TNBC = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-selective-targets-20261120")
MS = os.path.join(BASE, "manuscript")
OUT = os.path.join(MS, "source_data")
os.makedirs(OUT, exist_ok=True)

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def write_csv(name, fieldnames, rows):
    path = os.path.join(OUT, name)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {name} ({len(rows)} rows)")
    return path

CANCERS = [
    ("tnbc", "Triple-negative breast cancer", "TNBC"),
    ("egfr-mutant-lung", "EGFR-mutant lung cancer", "EGFR lung"),
    ("hcc", "Hepatocellular carcinoma", "HCC"),
    ("melanoma", "Melanoma", "Melanoma"),
    ("colon", "Colon cancer", "Colon"),
    ("glioma", "Glioma", "Glioma"),
    ("pancreatic", "Pancreatic cancer", "PDAC"),
    ("kras-mutant-lung", "KRAS-mutant lung cancer", "KRAS lung"),
    ("head-and-neck", "Head and neck squamous cell carcinoma", "HNSC"),
    ("ovarian", "Ovarian cancer", "Ovarian"),
    ("small-cell-lung", "Small-cell lung cancer", "SCLC"),
    ("esophagogastric", "Esophagogastric adenocarcinoma", "ESGA"),
    ("esophageal-squamous", "Oesophageal squamous cell carcinoma", "ESCC"),
    ("bladder", "Bladder cancer", "Bladder"),
    ("renal-cell", "Renal cell carcinoma", "RCC"),
    ("ewing", "Ewing sarcoma", "Ewing"),
    ("endometrial", "Endometrial cancer", "Endometrial"),
    ("aml", "Acute myeloid leukaemia", "AML"),
    ("neuroblastoma", "Neuroblastoma", "Neuroblastoma"),
    ("osteosarcoma", "Osteosarcoma", "Osteosarcoma"),
]
CKEYS = [c[0] for c in CANCERS]
CLABEL = {c[0]: c[2] for c in CANCERS}
CFULL = {c[0]: c[1] for c in CANCERS}

def cancer_dir(key):
    return TNBC if key == "tnbc" else os.path.join(BASE, key)

def gene_col(rows):
    return "Symbol" if "Symbol" in rows[0] else "Gene"

def score_col(rows):
    return "CombinedScore" if "CombinedScore" in rows[0] else "Score"

# ---------------------------------------------------------------- top-20 lists
top20 = {}
for key in CKEYS:
    rows = read_csv(os.path.join(cancer_dir(key), "top20.csv"))
    gc, sc = gene_col(rows), score_col(rows)
    top20[key] = [(r[gc].split(" (")[0].strip(), float(r[sc])) for r in rows]
rec = Counter(g for v in top20.values() for g, _ in v)
assert sum(len(v) for v in top20.values()) == 385, "slots != 385"
assert len(rec) == 269, "distinct genes != 269"

# druggability strip: calibration target list only (20 targets)
tsum = read_csv(os.path.join(BASE, "drug-target-calibration", "target_summary.csv"))
DRUG = {}
for r in tsum:
    st = r["Drug_status"]
    DRUG[r["Target"]] = "approved" if st.startswith("Approved") else ("clinical" if "Clinical" in st else "none")

ADHESION = {"ITGAV", "ITGB5", "FERMT2", "PTK2"}
LINEAGE = {"KLF5", "SOX9", "HNF1B", "TP63", "PAX8"}
POLARITY = {"PARD3", "PARD6B", "ARHGAP29"}
def cluster(g):
    if g in ADHESION: return "adhesion axis"
    if g in LINEAGE: return "lineage factor"
    if g in POLARITY: return "polarity module"
    return "none"

# ---------------------------------------------------------------- Fig 2
score_of = {}
for k in CKEYS:
    for g, s in top20[k]:
        score_of[(g, k)] = s
genes_sorted = sorted(rec.keys(), key=lambda g: (-rec[g],
    -max(s for kk in CKEYS for gg, s in top20[kk] if gg == g), g))
fields = (["gene", "cluster", "druggability_at_gene_calibration_list",
           "n_cancers_in_top20"]
          + [CFULL[k] for k in CKEYS])
rows = []
for g in genes_sorted:
    r = {"gene": g, "cluster": cluster(g),
         "druggability_at_gene_calibration_list": DRUG.get(g, "none"),
         "n_cancers_in_top20": rec[g]}
    for k in CKEYS:
        r[CFULL[k]] = score_of.get((g, k), "")
    rows.append(r)
write_csv("Source_Data_Fig2.csv", fields, rows)

# ---------------------------------------------------------------- Fig 3
rec3 = {g: c for g, c in rec.items() if c >= 3}
assert len(rec3) == 25, "recurring >=3 != 25"
rows = []
for g, c in sorted(rec3.items(), key=lambda kv: (-kv[1], kv[0])):
    rows.append({"gene": g, "n_cancers_in_top20": c,
                 "druggability_at_gene_calibration_list": DRUG.get(g, "none"),
                 "n_genes_represented": 1})
rows.append({"gene": "200 singleton genes (aggregate)",
             "n_cancers_in_top20": 1,
             "druggability_at_gene_calibration_list": "not applicable (aggregate)",
             "n_genes_represented": 200})
write_csv("Source_Data_Fig3.csv",
          ["gene", "n_cancers_in_top20", "druggability_at_gene_calibration_list",
           "n_genes_represented"], rows)

# ---------------------------------------------------------------- Fig 4
TEAD_GENES = ["TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "PTK2", "PKN2"]
TEAD_LABEL = {"TEAD1": "TEAD1", "TEAD2": "TEAD2", "TEAD3": "TEAD3", "TEAD4": "TEAD4",
              "YAP1": "YAP1", "WWTR1": "WWTR1 (TAZ)", "PTK2": "PTK2 (FAK)", "PKN2": "PKN2"}
def enrichment_ranks(key):
    rows = read_csv(os.path.join(cancer_dir(key), "dependency_scan.csv"))
    gc = gene_col(rows)
    if "EnrichRank" in rows[0]:
        return {r[gc].split(" (")[0].strip(): int(float(r["EnrichRank"])) for r in rows}
    srt = sorted(rows, key=lambda r: float(r["Enrichment"]))
    return {r[gc].split(" (")[0].strip(): i + 1 for i, r in enumerate(srt)}
tead_rank = {key: enrichment_ranks(key) for key in CKEYS}
meso = {r["Gene"]: int(r["EnrichRank"]) for r in read_csv(os.path.join(BASE, "mesothelioma", "tead_axis.csv"))}
tead_rank["mesothelioma"] = meso
SETTINGS = CKEYS + ["mesothelioma"]
SLABEL = dict(CLABEL); SLABEL["mesothelioma"] = "Mesothelioma"
CLASS = {}
for k in ["tnbc", "melanoma", "glioma", "ovarian", "renal-cell", "osteosarcoma"]:
    CLASS[k] = "TAZ cancer"
for k in ["colon", "esophagogastric", "head-and-neck", "mesothelioma"]:
    CLASS[k] = "YAP1 cancer / leaning"
for k in ["small-cell-lung", "esophageal-squamous", "ewing", "endometrial", "aml", "neuroblastoma"]:
    CLASS[k] = "null / anti-TEAD pole"
for k in ["egfr-mutant-lung", "hcc", "pancreatic", "kras-mutant-lung", "bladder"]:
    CLASS[k] = "other"
assert tead_rank["mesothelioma"]["TEAD1"] == 43
assert tead_rank["aml"]["TEAD1"] == 17802
assert tead_rank["bladder"]["TEAD3"] == 53
fields = ["gene"] + [CFULL[k] for k in CKEYS] + ["Mesothelioma"] + ["top20_in"]
rows = []
for g in TEAD_GENES:
    r = {"gene": TEAD_LABEL[g]}
    for k in CKEYS:
        r[CFULL[k]] = tead_rank[k][g]
    r["Mesothelioma"] = tead_rank["mesothelioma"][g]
    r["top20_in"] = "; ".join(CLABEL[k] for k in CKEYS if g in [x for x, _ in top20[k]]) or "-"
    rows.append(r)
# cancer-class strip as a labelled annotation row (documented in the README)
rows.append({"gene": "CANCER CLASS (annotation strip, from manuscript text)",
             **{CFULL[k]: CLASS[k] for k in CKEYS},
             "Mesothelioma": CLASS["mesothelioma"], "top20_in": "-"})
write_csv("Source_Data_Fig4.csv", fields, rows)

# ---------------------------------------------------------------- Fig 5
tg = read_csv(os.path.join(BASE, "tabula-sapiens", "per_gene_summary_truecpm.csv"))
q269 = [r for r in tg if r["in_269"] == "True"]
bg = read_csv(os.path.join(BASE, "tabula-sapiens", "background_summary_truecpm.csv"))
assert len(tg) == 290 and len(q269) == 269 and len(bg) == 300
itgav = next(r for r in tg if r["gene"] == "ITGAV")
assert abs(float(itgav["worst_mean_cpm10k"]) - 27.2) < 1.0, \
    f"ITGAV worst CPM is {itgav['worst_mean_cpm10k']}, expected ~27.2 (true CPM)"
assert "retinal pigment" in itgav["worst_celltype"]
fields = ["gene", "set", "tau", "switch_like",
          "n_types_on_frac_ge_0.50", "n_types_off_frac_le_0.01",
          "worst_celltype", "worst_mean_cpm10k", "worst_frac",
          "max_tpm_fc"]
rows = []
for r in q269:
    rows.append({"gene": r["gene"], "set": "269 selective genes", "tau": r["tau"],
                 "switch_like": r["switch_like"],
                 "n_types_on_frac_ge_0.50": r["n_on_frac_ge_0.50"],
                 "n_types_off_frac_le_0.01": r["n_off_frac_le_0.01"],
                 "worst_celltype": r["worst_celltype"],
                 "worst_mean_cpm10k": r["worst_mean_cpm10k"],
                 "worst_frac": r["worst_frac"],
                 "max_tpm_fc": r["max_tpm_fc"]})
for r in bg:
    rows.append({"gene": r["gene"], "set": "300 background genes", "tau": r["tau"],
                 "switch_like": r["switch_like"],
                 "n_types_on_frac_ge_0.50": r["n_on_frac_ge_0.50"],
                 "n_types_off_frac_le_0.01": r["n_off_frac_le_0.01"],
                 "worst_celltype": r["worst_celltype"],
                 "worst_mean_cpm10k": r["worst_mean_cpm10k"],
                 "worst_frac": r["worst_frac"],
                 "max_tpm_fc": ""})
write_csv("Source_Data_Fig5.csv", fields, rows)

print("SOURCE DATA DONE")
print(f"  Fig 2: {len(genes_sorted)} genes x {len(CKEYS)} cancers")
print(f"  Fig 3: {len(rec3)} recurring genes + singletons aggregate")
print(f"  Fig 4: {len(TEAD_GENES)} genes x {len(SETTINGS)} settings")
print(f"  Fig 5: {len(q269)} selective + {len(bg)} background genes")
print(f"  ITGAV worst: {itgav['worst_celltype']} "
      f"{float(itgav['worst_mean_cpm10k']):.1f} CPM10k (true CPM), frac {itgav['worst_frac']}")
