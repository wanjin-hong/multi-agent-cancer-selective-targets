#!/usr/bin/env python3
"""verify_source_data.py — checks the Source Data files in manuscript/source_data/
against the manuscript's headline values and, for a sample of cells, against the
original source CSVs. Independent of write_source_data.py's own asserts: it reads
only the written files (plus the originals for spot-checks).
"""
import csv, math, os, statistics
from collections import Counter

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
TNBC = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-selective-targets-20261120")
SD = os.path.join(BASE, "manuscript", "source_data")
PASS, FAIL = 0, 0
FAILURES = []

def check(name, got, want, tol=None):
    global PASS, FAIL
    ok = abs(got - want) <= tol if tol is not None else got == want
    if ok:
        PASS += 1
        print(f"  PASS {name}: {got}")
    else:
        FAIL += 1
        FAILURES.append((name, got, want))
        print(f"  FAIL {name}: got {got}, want {want}")

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

# ---------------------------------------------------------------- Fig 2
print("== Source_Data_Fig2.csv ==")
f2 = read_csv(os.path.join(SD, "Source_Data_Fig2.csv"))
check("Fig2 rows (genes)", len(f2), 269)
cancer_cols = [c for c in f2[0] if c not in
               ("gene", "cluster", "druggability_at_gene_calibration_list", "n_cancers_in_top20")]
check("Fig2 cancer columns", len(cancer_cols), 20)
filled = sum(1 for r in f2 for c in cancer_cols if r[c] != "")
check("Fig2 filled cells (top-20 slots)", filled, 385)
check("Fig2 singletons", sum(1 for r in f2 if int(r["n_cancers_in_top20"]) == 1), 200)
for g, w in [("KLF5", 8), ("HNF1B", 6), ("WWTR1", 6), ("FERMT2", 6),
             ("EGFR", 5), ("ITGAV", 5), ("CDK6", 4)]:
    check(f"Fig2 {g} recurrence", int(next(r for r in f2 if r["gene"] == g)["n_cancers_in_top20"]), w)
check("Fig2 TEAD1 absent", "TEAD1" in [r["gene"] for r in f2], False)
check("Fig2 KLF5 in 8 cancers",
      sum(1 for c in cancer_cols if next(r for r in f2 if r["gene"] == "KLF5")[c] != ""), 8)
check("Fig2 EGFR drug strip", next(r for r in f2 if r["gene"] == "EGFR")["druggability_at_gene_calibration_list"], "approved")
check("Fig2 KLF5 cluster", next(r for r in f2 if r["gene"] == "KLF5")["cluster"], "lineage factor")
check("Fig2 ITGAV cluster", next(r for r in f2 if r["gene"] == "ITGAV")["cluster"], "adhesion axis")
# spot-checks against the original top20.csv files
def top20_value(key, gene):
    d = TNBC if key == "tnbc" else os.path.join(BASE, key)
    rows = read_csv(os.path.join(d, "top20.csv"))
    gc = "Symbol" if "Symbol" in rows[0] else "Gene"
    sc = "CombinedScore" if "CombinedScore" in rows[0] else "Score"
    for r in rows:
        if r[gc].split(" (")[0].strip() == gene:
            return float(r[sc])
    return None
CFULL = {"tnbc": "Triple-negative breast cancer", "egfr-mutant-lung": "EGFR-mutant lung cancer",
         "hcc": "Hepatocellular carcinoma", "melanoma": "Melanoma", "colon": "Colon cancer",
         "glioma": "Glioma", "pancreatic": "Pancreatic cancer", "kras-mutant-lung": "KRAS-mutant lung cancer",
         "head-and-neck": "Head and neck squamous cell carcinoma", "ovarian": "Ovarian cancer",
         "small-cell-lung": "Small-cell lung cancer", "esophagogastric": "Esophagogastric adenocarcinoma",
         "esophageal-squamous": "Oesophageal squamous cell carcinoma", "bladder": "Bladder cancer",
         "renal-cell": "Renal cell carcinoma", "ewing": "Ewing sarcoma", "endometrial": "Endometrial cancer",
         "aml": "Acute myeloid leukaemia", "neuroblastoma": "Neuroblastoma", "osteosarcoma": "Osteosarcoma"}
for key, gene in [("tnbc", "MT1E"), ("egfr-mutant-lung", "EGFR"), ("hcc", "ITGAV"),
                  ("glioma", "FERMT2"), ("pancreatic", "KLF5"), ("kras-mutant-lung", "KRAS"),
                  ("head-and-neck", "TP63"), ("ovarian", "PAX8"), ("osteosarcoma", "RUNX2"),
                  ("esophagogastric", "KLF5")]:
    got = float(next(r for r in f2 if r["gene"] == gene)[CFULL[key]])
    check(f"Fig2 spot-check {gene} in {key}", got, top20_value(key, gene), 0.005)
# every non-empty cell matches its original top20.csv value
mismatch = 0
for r in f2:
    for k, col in CFULL.items():
        if r[col] != "":
            if abs(float(r[col]) - top20_value(k, r["gene"])) > 0.005:
                mismatch += 1
check("Fig2 all 385 cells match original top20.csv", mismatch, 0)

# ---------------------------------------------------------------- Fig 3
print("== Source_Data_Fig3.csv ==")
f3 = read_csv(os.path.join(SD, "Source_Data_Fig3.csv"))
check("Fig3 rows (25 genes + aggregate)", len(f3), 26)
check("Fig3 KLF5", int(f3[0]["n_cancers_in_top20"]), 8)
check("Fig3 first gene is KLF5", f3[0]["gene"], "KLF5")
check("Fig3 genes >= 3 cancers", sum(1 for r in f3[:-1] if int(r["n_cancers_in_top20"]) >= 3), 25)
agg = f3[-1]
check("Fig3 aggregate row genes", int(agg["n_genes_represented"]), 200)
check("Fig3 aggregate n_cancers", int(agg["n_cancers_in_top20"]), 1)
check("Fig3 gene set matches Fig2 >=3",
      sorted(r["gene"] for r in f3[:-1]),
      sorted(r["gene"] for r in f2 if int(r["n_cancers_in_top20"]) >= 3))

# ---------------------------------------------------------------- Fig 4
print("== Source_Data_Fig4.csv ==")
f4 = read_csv(os.path.join(SD, "Source_Data_Fig4.csv"))
genes4 = [r for r in f4 if not r["gene"].startswith("CANCER CLASS")]
check("Fig4 gene rows", len(genes4), 8)
check("Fig4 annotation row present",
      sum(1 for r in f4 if r["gene"].startswith("CANCER CLASS")), 1)
def cell4(gene, cancer_col):
    return int(next(r for r in genes4 if r["gene"].startswith(gene))[cancer_col])
check("Fig4 TEAD1 mesothelioma", cell4("TEAD1", "Mesothelioma"), 43)
check("Fig4 TEAD1 AML", cell4("TEAD1", "Acute myeloid leukaemia"), 17802)
check("Fig4 TEAD1 KRAS lung", cell4("TEAD1", "KRAS-mutant lung cancer"), 8248)
check("Fig4 TEAD3 bladder", cell4("TEAD3", "Bladder cancer"), 53)
check("Fig4 TEAD1 glioma", cell4("TEAD1", "Glioma"), 291)
check("Fig4 YAP1 mesothelioma", cell4("YAP1", "Mesothelioma"), 29)
check("Fig4 WWTR1 mesothelioma", cell4("WWTR1", "Mesothelioma"), 360)
check("Fig4 PTK2 HCC", cell4("PTK2", "Hepatocellular carcinoma"), 3)
check("Fig4 PTK2 HNSC", cell4("PTK2", "Head and neck squamous cell carcinoma"), 17845)
check("Fig4 WWTR1 top-20 in 6 cancers",
      next(r for r in genes4 if r["gene"].startswith("WWTR1"))["top20_in"],
      "TNBC; Melanoma; Glioma; Ovarian; RCC; Osteosarcoma")
check("Fig4 TEAD1 top-20 in none",
      next(r for r in genes4 if r["gene"] == "TEAD1")["top20_in"], "-")
check("Fig4 TEAD3 top-20 bladder only",
      next(r for r in genes4 if r["gene"] == "TEAD3")["top20_in"], "Bladder")
# spot-check ranks against the original dependency_scan.csv files
def rank_orig(key, gene):
    d = TNBC if key == "tnbc" else os.path.join(BASE, key)
    rows = read_csv(os.path.join(d, "dependency_scan.csv"))
    gc = "Symbol" if "Symbol" in rows[0] else "Gene"
    if "EnrichRank" in rows[0]:
        for r in rows:
            if r[gc].split(" (")[0].strip() == gene:
                return int(float(r["EnrichRank"]))
    srt = sorted(rows, key=lambda r: float(r["Enrichment"]))
    return [r[gc].split(" (")[0].strip() for r in srt].index(gene) + 1
for key, col, gene in [("glioma", "Glioma", "TEAD1"), ("melanoma", "Melanoma", "TEAD1"),
                       ("egfr-mutant-lung", "EGFR-mutant lung cancer", "TEAD1"),
                       ("pancreatic", "Pancreatic cancer", "PTK2"),
                       ("head-and-neck", "Head and neck squamous cell carcinoma", "TEAD3")]:
    check(f"Fig4 spot-check {gene} in {key}", cell4(gene, col), rank_orig(key, gene))
check("Fig4 cancer-class strip: TNBC is TAZ cancer",
      next(r for r in f4 if r["gene"].startswith("CANCER CLASS"))["Triple-negative breast cancer"],
      "TAZ cancer")

# ---------------------------------------------------------------- Fig 5
print("== Source_Data_Fig5.csv ==")
f5 = read_csv(os.path.join(SD, "Source_Data_Fig5.csv"))
sel = [r for r in f5 if r["set"] == "269 selective genes"]
bgr = [r for r in f5 if r["set"] == "300 background genes"]
check("Fig5 selective rows", len(sel), 269)
check("Fig5 background rows", len(bgr), 300)
check("Fig5 switch-like selective", sum(1 for r in sel if r["switch_like"] == "True"), 25)
check("Fig5 switch-like background", sum(1 for r in bgr if r["switch_like"] == "True"), 63)
tau_s = [float(r["tau"]) for r in sel]
tau_b = [float(r["tau"]) for r in bgr]
check("Fig5 median tau selective", round(statistics.median(tau_s), 3), 0.891, 0.001)
check("Fig5 median tau background", round(statistics.median(tau_b), 3), 0.904, 0.001)
check("Fig5 tau>=0.9 selective", sum(1 for v in tau_s if v >= 0.9), 121)
check("Fig5 tau>=0.9 background", sum(1 for v in tau_b if v >= 0.9), 152)
check("Fig5 median worst normal",
      round(statistics.median(float(r["worst_mean_cpm10k"]) for r in sel), 1), 2.2, 0.5)
check("Fig5 worst <100", sum(1 for r in sel if float(r["worst_mean_cpm10k"]) < 100), 267)
check("Fig5 worst >=500", sum(1 for r in sel if float(r["worst_mean_cpm10k"]) >= 500), 0)
def fc(r):
    return float(r["max_tpm_fc"]) if r["max_tpm_fc"] not in ("", None) else None
wide = [r for r in sel if fc(r) is not None and fc(r) >= 3 and float(r["worst_mean_cpm10k"]) < 100]
narrow = [r for r in sel if fc(r) is not None and fc(r) >= 3 and float(r["worst_mean_cpm10k"]) >= 500]
check("Fig5 wide-window genes", len(wide), 90)
check("Fig5 narrow-window genes", len(narrow), 0)
switch_genes = sorted(r["gene"] for r in sel if r["switch_like"] == "True")
check("Fig5 25 switch-like genes", switch_genes,
      ["ASCL1", "CDX2", "CEP55", "CKAP2L", "CSF2RB", "FOXA2", "GFI1", "GPR61",
       "HAND2", "HNF1A", "HNF1B", "INSM1", "ISL1", "KIF18B", "LGALS9B", "MYBL2",
       "MYCN", "PHOX2A", "RCOR2", "SOX10", "SOX11", "SOX2", "SPRR1B", "TFAP2A", "TUBA3D"])
itgav = next(r for r in f5 if r["gene"] == "ITGAV")
check("Fig5 ITGAV worst cell type is RPE", "retinal pigment" in itgav["worst_celltype"], True)
check("Fig5 ITGAV worst CPM10k (true CPM)", float(itgav["worst_mean_cpm10k"]), 27.2, 1.0)
check("Fig5 ITGAV value is NOT the inflated 1227.5",
      abs(float(itgav["worst_mean_cpm10k"]) - 1227.5) > 100, True)
check("Fig5 ITGAV worst frac", round(float(itgav["worst_frac"]), 2), 0.93, 0.01)
# cross-check every Fig5 row against the corrected Tabula summaries
tg = read_csv(os.path.join(BASE, "tabula-sapiens", "per_gene_summary_truecpm.csv"))
bg = read_csv(os.path.join(BASE, "tabula-sapiens", "background_summary_truecpm.csv"))
tgd = {r["gene"]: r for r in tg}
bgd = {r["gene"]: r for r in bg}
mismatch = 0
for r in sel:
    o = tgd[r["gene"]]
    if (abs(float(r["tau"]) - float(o["tau"])) > 1e-9 or r["switch_like"] != o["switch_like"]
            or abs(float(r["worst_mean_cpm10k"]) - float(o["worst_mean_cpm10k"])) > 1e-6):
        mismatch += 1
for r in bgr:
    o = bgd[r["gene"]]
    if (abs(float(r["tau"]) - float(o["tau"])) > 1e-9 or r["switch_like"] != o["switch_like"]
            or abs(float(r["worst_mean_cpm10k"]) - float(o["worst_mean_cpm10k"])) > 1e-6):
        mismatch += 1
check("Fig5 all 569 rows match corrected Tabula summaries", mismatch, 0)

print(f"\nSOURCE DATA VERIFICATION: {PASS} passed, {FAIL} failed")
if FAILURES:
    print("FAILURES:")
    for n, g, w in FAILURES:
        print(f"  - {n}: got {g}, want {w}")
    raise SystemExit(1)
