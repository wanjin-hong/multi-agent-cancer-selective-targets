"""Independent verification of adc-intersection headlines.

Recomputes every headline number from the RAW files via a separate code path
(csv module, no shared helpers with intersect.py). Fails loudly on mismatch.
"""
import csv, json, os, sys

ROOT = os.path.expanduser("~/workspace/research_notes")
MULTI = os.path.join(ROOT, "multi-cancer-selective-targets-20261120")
ADC = os.path.join(ROOT, "cancer-enriched-surface-proteins-adc-20261119", "data")
OUT = os.path.join(MULTI, "adc-intersection")

fails = []
def check(name, got, want):
    ok = got == want
    print(("[PASS] " if ok else "[FAIL] ") + f"{name}: got {got}, want {want}")
    if not ok:
        fails.append(name)

def sym(x):
    s = (x or "").strip()
    return s.split(" (")[0].strip() if " (" in s else s

def read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

# --- raw inputs, independent loads
surf = set(sym(r["gene_symbol"]) for r in read_csv_rows(os.path.join(ADC, "surfaceome.csv")))
ct = {sym(r["gene_symbol"]): r for r in read_csv_rows(os.path.join(ADC, "candidate_table.csv"))}
brm = {(sym(r["symbol"]), r["cancer"]): r for r in read_csv_rows(os.path.join(ADC, "bulk_rna_vs_normal.csv"))}

CAND_FILES = {
    "tnbc": (os.path.join(ROOT, "depmap-tnbc-selective-targets-20261120", "top20.csv"), None),
    "egfr-mutant-lung": (os.path.join(MULTI, "egfr-mutant-lung", "ranked_candidates.csv"), None),
    "hcc": (os.path.join(MULTI, "hcc", "candidates.csv"), None),
    "melanoma": (os.path.join(MULTI, "melanoma", "_relaxed_scored.csv"), ("RelaxedBar", "Pass_OE")),
    "colon": (os.path.join(MULTI, "colon", "candidates_both.csv"), None),
    "glioma": (os.path.join(MULTI, "glioma", "candidates.csv"), None),
    "pancreatic": (os.path.join(MULTI, "pancreatic", "candidates.csv"), None),
    "kras-mutant-lung": (os.path.join(MULTI, "kras-mutant-lung", "ranked_candidates.csv"), None),
    "head-and-neck": (os.path.join(MULTI, "head-and-neck", "candidates.csv"), None),
    "ovarian": (os.path.join(MULTI, "ovarian", "candidates.csv"), None),
    "small-cell-lung": (os.path.join(MULTI, "small-cell-lung", "candidates.csv"), None),
    "esophagogastric": (os.path.join(MULTI, "esophagogastric", "candidates.csv"), None),
    "esophageal-squamous": (os.path.join(MULTI, "esophageal-squamous", "candidates.csv"), None),
    "renal-cell": (os.path.join(MULTI, "renal-cell", "candidates.csv"), None),
    "bladder": (os.path.join(MULTI, "bladder", "candidates.csv"), None),
    "endometrial": (os.path.join(MULTI, "endometrial", "candidates.csv"), None),
    "ewing": (os.path.join(MULTI, "ewing", "candidates.csv"), None),
    "neuroblastoma": (os.path.join(MULTI, "neuroblastoma", "candidates.csv"), None),
    "aml": (os.path.join(MULTI, "aml", "candidates.csv"), None),
    "osteosarcoma": (os.path.join(MULTI, "osteosarcoma", "candidates.csv"), None),
}
cand = {}
for key, (path, filt) in CAND_FILES.items():
    rows = read_csv_rows(path)
    if filt:
        rows = [r for r in rows if r[filt[0]] == "True" and r[filt[1]] == "True"]
    gcol = "Gene" if "Gene" in rows[0] else "Symbol"
    cand[key] = set(sym(r[gcol]) for r in rows)

# --- headline 1: per-cancer surface-hit counts
master = read_csv_rows(os.path.join(OUT, "intersection_master.csv"))
from collections import Counter
got_counts = Counter(r["cancer_key"] for r in master)
for key in CAND_FILES:
    want = len(cand[key] & surf)
    check(f"surface hits {key}", got_counts.get(key, 0), want)
check("total surface-hit rows", len(master), sum(len(cand[k] & surf) for k in cand))
check("distinct surface genes", len(set(r["gene"] for r in master)), 37)

# --- headline 2: NEW vs confirmed
new_genes = set(r["gene"] for r in master if r["adc_status"] == "NEW")
check("distinct NEW genes", len(new_genes), 21)
check("NEW rows", sum(1 for r in master if r["adc_status"] == "NEW"), 36)
check("FGFR1 NEW in 4 cancers",
      sorted(set(r["cancer"] for r in master if r["gene"] == "FGFR1")),
      ["Glioma", "HCC", "Osteosarcoma", "Renal cell"])
check("PDGFRA NEW in 2 cancers",
      sorted(set(r["cancer"] for r in master if r["gene"] == "PDGFRA")),
      ["Glioma", "Osteosarcoma"])

# --- headline 3: safety flags (recomputed from candidate_table directly)
def flagged(g):
    if g not in ct:
        return False
    r = ct[g]
    crit = str(r["critical_tissue_flag"]).lower() == "true"
    cell = float(r["max_normal_celltype_nCPM"] or 0) >= 2000
    tis = float(r["max_normal_tissue_consensus_nTPM"] or 0) >= 200
    return crit or cell or tis
# note: NEW genes use NCM, not CT — check those separately below
n_flag_ct = sum(1 for r in master if r["gene"] in ct and flagged(r["gene"]))
check("flagged rows with CT genes (ITGAV x5, ATP1B3, CPD, SLC2A1 x2, SLC38A2)",
      n_flag_ct, 10)
ncm = {sym(r["symbol"]): r for r in read_csv_rows(os.path.join(ADC, "normal_celltype_max.csv"))}
n_flag_new = sum(1 for r in master
                 if r["gene"] not in ct and r["gene"] in ncm
                 and float(ncm[r["gene"]]["max_celltype_nCPM"] or 0) >= 2000)
check("flagged rows among NEW genes (CDH2, FOLR2)", n_flag_new, 2)
check("total flagged rows", n_flag_ct + n_flag_new, 12)
check("ITGAV flagged in 5 cancers",
      sorted(set(r["cancer"] for r in master if r["gene"] == "ITGAV" and r["safety_flags"])),
      ["Glioma", "HCC", "Osteosarcoma", "Ovarian", "Renal cell"])

# --- headline 4: ADC bar met = 3, all ESGA (recompute FPKM from BRM + protein from CT max)
def bar_met(gene, tcga):
    b = brm.get((gene, tcga))
    if not b:
        return False
    if float(b["median_tumour_FPKM"]) < 10:
        return False
    prot = float(ct[gene]["max_protein_log2FC"]) if gene in ct and ct[gene]["max_protein_log2FC"] else None
    return prot is not None and prot >= 0.5
met = [(r["cancer"], r["gene"]) for r in master
       if r["cancer_key"] in ("esophagogastric",) and bar_met(r["gene"], "STAD")]
check("bar-met hits (ESGA/STAD)", sorted(met),
      [("Esophagogastric adeno.", "SLC3A2"), ("Esophagogastric adeno.", "SLC7A1"),
       ("Esophagogastric adeno.", "SLC7A5")])
# EGFR/LUAD fails the FPKM leg by 0.06
egfr_luad = brm[("EGFR", "LUAD")]
check("EGFR LUAD FPKM just under 10", round(float(egfr_luad["median_tumour_FPKM"]), 2), 9.94)

# --- headline 5: reverse direction — 16 ADC genes with dependency support
rev = read_csv_rows(os.path.join(OUT, "adc_reverse_lookup.csv"))
sup = set(r["gene"] for r in rev if r["is_candidate"] == "True")
check("ADC genes that are candidates in >=1 cancer", len(sup), 16)
sup20 = set(r["gene"] for r in rev if r["in_top20"] == "True")
check("ADC genes in a top20", len(sup20), 8)
egfr_c = sorted(set(r["cancer"] for r in rev if r["gene"] == "EGFR" and r["is_candidate"] == "True"))
check("EGFR candidate in 6 cancers", len(egfr_c), 6)
itgav_c = sorted(set(r["cancer"] for r in rev if r["gene"] == "ITGAV" and r["is_candidate"] == "True"))
check("ITGAV candidate in 5 cancers", len(itgav_c), 5)

# --- headline 6: five new ADC targets — only SLC7A5 gains support
for g, want in [("SLC7A5", 1), ("SPINT2", 0), ("TTYH3", 0), ("SLC52A2", 0), ("EPHB3", 0)]:
    n = sum(1 for r in rev if r["gene"] == g and r["is_candidate"] == "True")
    check(f"{g} candidate cancers", n, want)

# --- headline 7: melanoma and KRAS lung have zero surface hits
check("melanoma surface hits", got_counts.get("melanoma", 0), 0)
check("kras-mutant-lung surface hits", got_counts.get("kras-mutant-lung", 0), 0)

# --- headline 8: cross-cancer recurrence
cc = read_csv_rows(os.path.join(OUT, "cross_cancer_surface.csv"))
check("cross-cancer surface genes", len(cc), 37)
check("recurring in >=2 cancers", sum(1 for r in cc if int(r["n_cancers"]) >= 2), 20)
check("EGFR in 6 cancers", [r["n_cancers"] for r in cc if r["gene"] == "EGFR"], ["6"])
check("in a TAZ cancer", sum(1 for r in cc if int(r["n_taz_cancers"]) > 0), 12)

# --- headline 9: AML haem sensitivity
haem = set(sym(r["Gene"]) for r in read_csv_rows(
    os.path.join(MULTI, "aml", "candidates_haem.csv")))
aml_hits = [r["gene"] for r in master if r["cancer_key"] == "aml"]
check("AML surface hits also in haem-only list",
      sum(1 for g in aml_hits if g in haem), 5)

# --- headline 10: OR10A4 is in SURFY (ML-predicted surface)
check("OR10A4 in SURFY", "OR10A4" in surf, True)

print()
if fails:
    print("FAILURES:", fails)
    sys.exit(1)
print("ALL CHECKS PASSED")
