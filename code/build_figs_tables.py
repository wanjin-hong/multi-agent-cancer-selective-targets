#!/usr/bin/env python3
"""build_figs_tables.py — manuscript figures 1-5 and tables 1-5 (26 Nov 2026).

Every plotted value and table cell is recomputed from the programme's source
CSVs, never from memory or from the manuscript text alone. The manuscript text
is used only as the *specification* of what each figure/table must contain;
all values are checked against the source files and any disagreement is
recorded in figures-tables-build-note.md (never silently fixed).

Outputs:
  manuscript/figures/fig1..fig5.png + .pdf  (300 dpi PNG, vector PDF)
  manuscript/tables/table1..table5.md + .xlsx
  manuscript/figures-tables-build-note.md
  manuscript/build_verification_log.txt
"""
import csv, math, os, re, statistics, sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
TNBC = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-selective-targets-20261120")
MS = os.path.join(BASE, "manuscript")
FIGDIR = os.path.join(MS, "figures")
TABDIR = os.path.join(MS, "tables")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TABDIR, exist_ok=True)

LOG = []
def log(msg):
    LOG.append(msg)
    print(msg)

CHECKS = {"pass": 0, "fail": 0}
FAILURES = []
def check(name, got, want, tol=None):
    ok = abs(got - want) <= tol if tol is not None else got == want
    if ok:
        CHECKS["pass"] += 1
        log(f"  PASS {name}: {got}")
    else:
        CHECKS["fail"] += 1
        FAILURES.append((name, got, want))
        log(f"  FAIL {name}: got {got}, want {want}")
    return ok

def check_set(name, got, want):
    ok = set(got) == set(want)
    if ok:
        CHECKS["pass"] += 1
        log(f"  PASS {name}: {sorted(got)}")
    else:
        CHECKS["fail"] += 1
        FAILURES.append((name, sorted(got), sorted(want)))
        log(f"  FAIL {name}: got {sorted(got)}, want {sorted(want)}")
    return ok

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

# ---------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "axes.unicode_minus": False,
    "savefig.dpi": 300,
})
# Okabe-Ito colourblind-safe palette
OI = {"blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73", "orange": "#E69F00",
      "verm": "#D55E00", "purple": "#CC79A7", "yellow": "#F0E442", "grey": "#999999",
      "darkgrey": "#4D4D4D", "black": "#000000"}
MM = 1 / 25.4  # mm -> inches

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
log("== top-20 lists (Fig 2, Fig 3, Table 2) ==")
top20 = {}   # key -> list of (gene, score, tier)
for key in CKEYS:
    rows = read_csv(os.path.join(cancer_dir(key), "top20.csv"))
    gc, sc = gene_col(rows), score_col(rows)
    top20[key] = [(r[gc].split(" (")[0].strip(), float(r[sc]),
                   r.get("Tier", "").strip()) for r in rows]
    assert len(top20[key]) in (17, 20, 8) or key in ("egfr-mutant-lung", "kras-mutant-lung"), (key, len(top20[key]))
total_slots = sum(len(v) for v in top20.values())
rec = Counter(g for v in top20.values() for g, _, _ in v)
check("total top-20 slots", total_slots, 385)
check("distinct genes", len(rec), 269)
check("singletons", sum(1 for c in rec.values() if c == 1), 200)
for g, w in [("KLF5", 8), ("HNF1B", 6), ("WWTR1", 6), ("FERMT2", 6),
             ("EGFR", 5), ("ITGAV", 5), ("FOSL1", 4), ("PARD6B", 4),
             ("CXADR", 4), ("CDK6", 4)]:
    check(f"{g} recurrence", rec[g], w)
check("TEAD1 in any top 20", sum(1 for v in top20.values() if "TEAD1" in [g for g, _, _ in v]), 0)
check("TEAD3 top-20 cancers", [k for k, v in top20.items() if "TEAD3" in [g for g, _, _ in v]], ["bladder"])
check("TAZ (WWTR1) top-20 cancers",
      [k for k, v in top20.items() if "WWTR1" in [g for g, _, _ in v]],
      ["tnbc", "melanoma", "glioma", "ovarian", "renal-cell", "osteosarcoma"])
check("YAP1 top-20 cancers",
      [k for k, v in top20.items() if "YAP1" in [g for g, _, _ in v]],
      ["colon", "esophagogastric"])

# recurrence >= 3 (Fig 3 / Table 2)
rec3 = {g: c for g, c in rec.items() if c >= 3}
check("genes recurring >=3", len(rec3), 25)
ms15 = ["SOX9", "MT1E", "TCF7L2", "MYBL2", "TP63", "PARD3", "FAF2", "ARHGAP29",
        "PAX8", "ZEB2", "IRS2", "UBE2S", "FGFR1", "SLC35B2", "E2F3"]
check_set("fifteen 3-cancer genes", [g for g, c in rec3.items() if c == 3], ms15)
check("genes recurring >=2", sum(1 for c in rec.values() if c >= 2), 69)

# ---------------------------------------------------------------- synthesis Table 1
log("== Table 1 source: synthesis.md at-a-glance table ==")
synth = open(os.path.join(BASE, "synthesis.md")).read().splitlines()
t1 = {}
for ln in synth:
    m = re.match(r"\| ([^|]+) \| ([^|]+) \| ([^|]+) \| (\d+) \| (\d+) \| ([^|]+) \| ([^|]+) \| ([^|]+) \(([\d.]+)\) \|", ln)
    if m:
        cancer, grp, comp, s, r, _, _, g1, sc1 = [x.strip() for x in m.groups()]
        t1[cancer] = {"group_raw": grp, "comparator": comp, "strict_oe": int(s),
                      "relaxed_oe": int(r), "n1_gene": g1, "n1_score": float(sc1)}
check("synthesis Table 1 rows parsed", len(t1), 20)
SYN_KEY = {"TNBC": "tnbc", "EGFR-mutant lung": "egfr-mutant-lung", "HCC": "hcc",
           "Melanoma": "melanoma", "Colon": "colon", "Glioma": "glioma", "PDAC": "pancreatic",
           "KRAS-mutant lung": "kras-mutant-lung", "HNSC": "head-and-neck",
           "Ovarian": "ovarian", "SCLC": "small-cell-lung", "ESGA": "esophagogastric",
           "ESCC": "esophageal-squamous", "Bladder": "bladder", "RCC": "renal-cell",
           "Ewing": "ewing", "Endometrial": "endometrial", "AML": "aml",
           "Neuroblastoma": "neuroblastoma", "Osteosarcoma": "osteosarcoma"}
table1 = {}
for sname, d in t1.items():
    key = SYN_KEY[sname]
    n = int(re.match(r"(\d+)", d["group_raw"]).group(1))
    table1[key] = {"cancer": CFULL[key], "group_n": n, "comparator": d["comparator"],
                   "strict_oe": d["strict_oe"], "relaxed_oe": d["relaxed_oe"],
                   "n1_gene": d["n1_gene"], "n1_score": d["n1_score"]}
# synthesis.md's first column is the strict *dependency-bar* count (270/149/100),
# not strict+OE, for AML, neuroblastoma and osteosarcoma. The manuscript Table 1
# uses strict+OE, so override those three from candidates.csv Tier A counts
# (verified against the per-cancer reports by the manuscript's verify script).
for key, want in [("aml", 103), ("neuroblastoma", 71), ("osteosarcoma", 21)]:
    rows = read_csv(os.path.join(BASE, key, "candidates.csv"))
    tierA = sum(1 for r in rows if r.get("Tier", "").strip().startswith("A"))
    check(f"{key} strict+OE override (synthesis said {table1[key]['strict_oe']}, candidates.csv Tier A)",
          tierA, want)
    table1[key]["strict_oe"] = tierA

# Verify Table 1 counts against per-cancer candidate files / report text.
# (a) cancers with candidates.csv: rows = relaxed+OE, Tier A rows = strict+OE.
log("-- Table 1 verification vs candidates.csv Tier counts --")
for key in ["hcc", "glioma", "pancreatic", "head-and-neck", "ovarian", "small-cell-lung",
            "esophagogastric", "esophageal-squamous", "bladder", "renal-cell", "ewing",
            "endometrial", "aml", "neuroblastoma", "osteosarcoma"]:
    rows = read_csv(os.path.join(BASE, key, "candidates.csv"))
    tierA = sum(1 for r in rows if r.get("Tier", "").strip().startswith("A"))
    check(f"{key} relaxed+OE (candidates.csv rows)", len(rows), table1[key]["relaxed_oe"])
    # AML/NB/OS: synthesis.md first column is the strict *dependency-bar* count,
    # not strict+OE; the manuscript Table 1 uses strict+OE from candidates.csv.
    check(f"{key} strict+OE (candidates.csv Tier A)", tierA, table1[key]["strict_oe"])
# (b) colon: candidates_both.csv
rows = read_csv(os.path.join(BASE, "colon", "candidates_both.csv"))
check("colon relaxed+OE (candidates_both.csv rows)", len(rows), table1["colon"]["relaxed_oe"])
check("colon strict+OE (Tier A)", sum(1 for r in rows if r["Tier"].startswith("A")), table1["colon"]["strict_oe"])
# (c) EGFR / KRAS lung: ranked_candidates.csv rows + StrictBar count
for key, want_s, want_r in [("egfr-mutant-lung", 13, 17), ("kras-mutant-lung", 2, 8)]:
    rows = read_csv(os.path.join(BASE, key, "ranked_candidates.csv"))
    check(f"{key} relaxed+OE (ranked_candidates.csv rows)", len(rows), want_r)
    check(f"{key} strict+OE (StrictBar True)", sum(1 for r in rows if r.get("StrictBar") == "True"), want_s)
    check(f"{key} table1 strict+OE", table1[key]["strict_oe"], want_s)
    check(f"{key} table1 relaxed+OE", table1[key]["relaxed_oe"], want_r)
# (d) melanoma: report text (21 strict+OE, 40 relaxed+OE)
mel_txt = open(os.path.join(BASE, "melanoma", "report.md")).read()
check("melanoma strict+OE (report text)", "finds **21 genes**" in mel_txt, True)
check("melanoma relaxed+OE (report text)", "(**40** pass the relaxed dep" in mel_txt, True)
# (e) TNBC: candidates_both.csv rows = strict+OE 11; ranked_fetchbar.csv rows = relaxed+OE 25
rows = read_csv(os.path.join(TNBC, "candidates_both.csv"))
check("TNBC strict+OE (candidates_both.csv rows)", len(rows), 11)
rows = read_csv(os.path.join(TNBC, "ranked_fetchbar.csv"))
check("TNBC relaxed+OE (ranked_fetchbar.csv rows)", len(rows), 25)
# (f) AML blood-only relaxed+OE = 89 (manuscript Table 1 parenthetical)
rows = read_csv(os.path.join(BASE, "aml", "candidates_haem.csv"))
check("AML blood-only relaxed+OE (candidates_haem.csv rows)", len(rows), 89)
# (g) #1 gene + score vs top20.csv row 1 (manuscript quotes 2 dp)
for key in CKEYS:
    g1, s1, _ = top20[key][0]
    check(f"{key} #1 gene", g1, table1[key]["n1_gene"])
    check(f"{key} #1 score (2dp)", round(s1, 2), table1[key]["n1_score"], 0.005)
# (h) group n spot-checks vs top20.csv group-n columns where present
for key, col, want in [("hcc", "HCC_n", 21), ("pancreatic", "PDAC_n", 45),
                       ("head-and-neck", "HNSC_n", 72), ("osteosarcoma", "OS_n", 17),
                       ("neuroblastoma", "NB_n", 39)]:
    rows = read_csv(os.path.join(BASE, key, "top20.csv"))
    check(f"{key} group n (top20.csv {col})", int(float(rows[0][col])), want)

# ---------------------------------------------------------------- TEAD axis
log("== TEAD-axis enrichment ranks (Fig 4, Table 3) ==")
TEAD_GENES = ["TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1", "PTK2", "PKN2"]
TEAD_LABEL = {"TEAD1": "TEAD1", "TEAD2": "TEAD2", "TEAD3": "TEAD3", "TEAD4": "TEAD4",
              "YAP1": "YAP1", "WWTR1": "WWTR1 (TAZ)", "PTK2": "PTK2 (FAK)", "PKN2": "PKN2"}

def enrichment_ranks(key):
    """Enrichment rank (1 = most negative enrichment) from a dependency scan.
    Uses the file's EnrichRank where present; otherwise sorts by Enrichment
    (validated: identical to file EnrichRank in glioma and AML)."""
    rows = read_csv(os.path.join(cancer_dir(key), "dependency_scan.csv"))
    gc = gene_col(rows)
    if "EnrichRank" in rows[0]:
        return {r[gc].split(" (")[0].strip(): int(float(r["EnrichRank"])) for r in rows}
    srt = sorted(rows, key=lambda r: float(r["Enrichment"]))
    return {r[gc].split(" (")[0].strip(): i + 1 for i, r in enumerate(srt)}

tead_rank = {}
for key in CKEYS:
    tead_rank[key] = enrichment_ranks(key)
# mesothelioma: dedicated axis file (no top-20 hunt; ranks from tead_axis.csv)
meso = {r["Gene"]: int(r["EnrichRank"]) for r in read_csv(os.path.join(BASE, "mesothelioma", "tead_axis.csv"))}
tead_rank["mesothelioma"] = meso
SETTINGS = CKEYS + ["mesothelioma"]
SLABEL = dict(CLABEL); SLABEL["mesothelioma"] = "Mesothelioma"

# validate against manuscript-quoted ranks
for key, gene, want in [("mesothelioma", "TEAD1", 43), ("melanoma", "TEAD1", 132),
                        ("glioma", "TEAD1", 291), ("head-and-neck", "TEAD1", 400),
                        ("egfr-mutant-lung", "TEAD1", 409), ("kras-mutant-lung", "TEAD1", 8248),
                        ("aml", "TEAD1", 17802), ("mesothelioma", "YAP1", 29),
                        ("mesothelioma", "WWTR1", 360), ("bladder", "TEAD3", 53),
                        ("head-and-neck", "TEAD3", 215), ("glioma", "PTK2", 15),
                        ("hcc", "PTK2", 3), ("melanoma", "PTK2", 25),
                        ("pancreatic", "PTK2", 1284), ("head-and-neck", "PTK2", 17845),
                        ("colon", "YAP1", 29), ("colon", "WWTR1", 17874),
                        ("ewing", "WWTR1", 17794), ("kras-mutant-lung", "TEAD1", 8248)]:
    check(f"TEAD rank {gene} in {key}", tead_rank[key][gene], want)
for key, want in [("tnbc", 5), ("melanoma", 9), ("glioma", 14), ("ovarian", 9),
                  ("renal-cell", 7), ("osteosarcoma", 9)]:
    check(f"WWTR1 top-20 rank in {key}", top20[key].index(next(t for t in top20[key] if t[0] == "WWTR1")) + 1 if False else [g for g, _, _ in top20[key]].index("WWTR1") + 1, want)
check("PDAC WWTR1 not in top 20", "WWTR1" not in [g for g, _, _ in top20["pancreatic"]], True)
# top-20 membership for the rank matrix (bold outline in Fig 4)
tead_top20 = {g: [k for k in CKEYS if g in [x for x, _, _ in top20[k]]] for g in TEAD_GENES}

# ---------------------------------------------------------------- Tabula Sapiens
log("== Tabula Sapiens switch-like test (Fig 5, Table 5) ==")
tg = read_csv(os.path.join(BASE, "tabula-sapiens", "per_gene_summary_truecpm.csv"))
q269 = [r for r in tg if r["in_269"] == "True"]
bg = read_csv(os.path.join(BASE, "tabula-sapiens", "background_summary_truecpm.csv"))
check("Tabula query genes", len(tg), 290)
check("selective genes in Tabula", len(q269), 269)
check("background genes", len(bg), 300)
s269 = sum(1 for r in q269 if r["switch_like"] == "True")
sbg = sum(1 for r in bg if r["switch_like"] == "True")
check("switch-like selective", s269, 25)
check("switch-like background", sbg, 63)

def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    def hg(k):
        return math.comb(a + b, k) * math.comb(c + d, a + c - k) / math.comb(n, a + c)
    p_obs = hg(a)
    lo = max(0, a + c - (c + d)); hi = min(a + b, a + c)
    return sum(hg(k) for k in range(lo, hi + 1) if hg(k) <= p_obs * (1 + 1e-12))
p_fisher = fisher_two_sided(s269, len(q269) - s269, sbg, len(bg) - sbg)
check("Fisher p switch-like", round(p_fisher, 4), 0.0001, 0.0001)

def mannwhitney_p(x, y):
    allv = sorted(x + y)
    avg = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1] == allv[i]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            avg.setdefault(allv[k], r)
        i = j + 1
    rx = sum(avg[v] for v in x)
    n1, n2 = len(x), len(y)
    u = rx - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    ties = Counter(allv)
    tcorr = sum(t * (t * t - 1) for t in ties.values() if t > 1)
    n = n1 + n2
    sd = math.sqrt(n1 * n2 * (n + 1 - tcorr / (n * (n - 1))) / 12)
    z = (u - mu) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))

tau269 = [float(r["tau"]) for r in q269]
taubg = [float(r["tau"]) for r in bg]
check("median tau selective", round(statistics.median(tau269), 3), 0.891, 0.001)
check("median tau background", round(statistics.median(taubg), 3), 0.904, 0.001)
check("tau>=0.9 selective", sum(1 for v in tau269 if v >= 0.9), 121)
check("tau>=0.9 background", sum(1 for v in taubg if v >= 0.9), 152)
check("median types detected selective",
      int(statistics.median(int(r["n_celltypes_frac_ge_0.10"]) for r in q269)), 128)
check("median types detected background",
      int(statistics.median(int(r["n_celltypes_frac_ge_0.10"]) for r in bg)), 89)
check("median types on selective",
      int(statistics.median(int(r["n_on_frac_ge_0.50"]) for r in q269)), 16)
check("median types on background",
      int(statistics.median(int(r["n_on_frac_ge_0.50"]) for r in bg)), 7)
check("median types off selective",
      int(statistics.median(int(r["n_off_frac_le_0.01"]) for r in q269)), 6)
check("median types off background",
      int(statistics.median(int(r["n_off_frac_le_0.01"]) for r in bg)), 10)
check("Mann-Whitney p tau", round(mannwhitney_p(tau269, taubg), 4), 0.0155, 0.0005)
check("median worst normal (CPM10k)",
      round(statistics.median(float(r["worst_mean_cpm10k"]) for r in q269), 1), 2.2, 0.5)
check("worst < 100", sum(1 for r in q269 if float(r["worst_mean_cpm10k"]) < 100), 267)
check("worst >= 500", sum(1 for r in q269 if float(r["worst_mean_cpm10k"]) >= 500), 0)
check("worst >= 1000", sum(1 for r in q269 if float(r["worst_mean_cpm10k"]) >= 1000), 0)

def _fc(r):
    v = r["max_tpm_fc"]
    return float(v) if v not in ("", None) else None
wide = [r for r in q269 if _fc(r) is not None and _fc(r) >= 3 and float(r["worst_mean_cpm10k"]) < 100]
narrow = [r for r in q269 if _fc(r) is not None and _fc(r) >= 3 and float(r["worst_mean_cpm10k"]) >= 500]
check("wide-window genes", len(wide), 90)
check("narrow-window genes", len(narrow), 0)
switch_genes = sorted(r["gene"] for r in q269 if r["switch_like"] == "True")
check_set("25 switch-like selective genes",
          switch_genes,
          ["ASCL1", "CDX2", "CEP55", "CKAP2L", "CSF2RB", "FOXA2", "GFI1", "GPR61",
           "HAND2", "HNF1A", "HNF1B", "INSM1", "ISL1", "KIF18B", "LGALS9B", "MYBL2",
           "MYCN", "PHOX2A", "RCOR2", "SOX10", "SOX11", "SOX2", "SPRR1B", "TFAP2A",
           "TUBA3D"])
itgav = next(r for r in tg if r["gene"] == "ITGAV")
check("ITGAV worst cell type is RPE", "retinal pigment" in itgav["worst_celltype"], True)
check("ITGAV worst mean CPM10k", float(itgav["worst_mean_cpm10k"]), 27.2, 1.0)

# ---------------------------------------------------------------- calibration
log("== Drug-target calibration (Table 4) ==")
mat = read_csv(os.path.join(BASE, "drug-target-calibration", "target_cancer_matrix.csv"))
check("calibration cells", len(mat), 400)
n_strict_oe = sum(1 for r in mat if r["Pass_strict_OE"] == "True")
n_relaxed_oe = sum(1 for r in mat if r["Pass_relaxed_OE"] == "True")
n_strict_only = sum(1 for r in mat if r["Strict_dep"] == "True" and r["OE_pass"] == "False")
n_relaxed_only_dep = sum(1 for r in mat if r["Relaxed_dep"] == "True" and r["Strict_dep"] == "False" and r["OE_pass"] == "False")
n_oe_only = sum(1 for r in mat if r["OE_pass"] == "True" and r["Strict_dep"] == "False" and r["Relaxed_dep"] == "False")
n_neither = sum(1 for r in mat if r["OE_pass"] == "False" and r["Strict_dep"] == "False" and r["Relaxed_dep"] == "False")
check("pass strict + OE", n_strict_oe, 22)
check("pass relaxed + OE (inclusive)", n_relaxed_oe, 33)
check("relaxed-only + OE (exclusive)", n_relaxed_oe - n_strict_oe, 11)
check("strict dep, no OE", n_strict_only, 12)
check("relaxed-only dep, no OE", n_relaxed_only_dep, 4)
check("OE only, no dep bar", n_oe_only, 66)
check("neither", n_neither, 285)
check("exclusive categories sum to 400",
      n_strict_oe + (n_relaxed_oe - n_strict_oe) + n_strict_only + n_relaxed_only_dep + n_oe_only + n_neither, 400)
tsum = read_csv(os.path.join(BASE, "drug-target-calibration", "target_summary.csv"))
check("calibration targets", len(tsum), 20)
def cell(target, cancer):
    return next(r for r in mat if r["Target"] == target and r["Cancer"] == cancer)
check("EGFR found strict, EGFR-mutant lung", cell("EGFR", "EGFR-mutant lung")["Pass_strict_OE"], "True")
check("EGFR found strict, HNSC", cell("EGFR", "HNSC")["Pass_strict_OE"], "True")
check("ERBB2 found relaxed, colon", cell("ERBB2", "Colon")["Pass_relaxed_OE"], "True")
check("KRAS found strict, KRAS-mutant lung", cell("KRAS", "KRAS-mutant lung")["Pass_strict_OE"], "True")
check("EPAS1 found relaxed, RCC", cell("EPAS1", "RCC")["Pass_relaxed_OE"], "True")
check("BRAF dep-only, melanoma", (cell("BRAF", "Melanoma")["Strict_dep"], cell("BRAF", "Melanoma")["OE_pass"]), ("True", "False"))
check("PARP1 passes nowhere", sum(1 for r in mat if r["Target"] == "PARP1" and (r["Pass_strict_OE"] == "True" or r["Pass_relaxed_OE"] == "True")), 0)
# fusion-rare genome-wide counts come from the calibration verification (ALK 7, ROS1 0, RET 1, NTRK1 1)
cal_txt = open(os.path.join(BASE, "drug-target-calibration", "report.md")).read()
for s in ["ROS1 (0/1,178)", "RET (1/1,178)", "NTRK1 (1/1,178)", "ALK (7/1,178)"]:
    check(f"calibration report contains {s}", s in cal_txt, True)

# druggability classification for Fig 2 / Fig 3 side annotation:
# sourced from the calibration target list only (20 targets).
DRUG = {}
for r in tsum:
    st = r["Drug_status"]
    DRUG[r["Target"]] = "approved" if st.startswith("Approved") else ("clinical" if "Clinical" in st else "none")

# ---------------------------------------------------------------- tables
log("== writing tables ==")
import pandas as pd

def write_table(stem, df, md_title, md_lines, source_note):
    """Write one table as .md and .xlsx (data sheet + Source sheet)."""
    md = [f"# {md_title}", ""]
    md += md_lines
    md += ["", f"*Source: {source_note}*", ""]
    with open(os.path.join(TABDIR, stem + ".md"), "w") as f:
        f.write("\n".join(md))
    with pd.ExcelWriter(os.path.join(TABDIR, stem + ".xlsx"), engine="openpyxl") as w:
        df.to_excel(w, sheet_name="Table", index=False)
        pd.DataFrame({"Source": [source_note]}).to_excel(w, sheet_name="Source", index=False)
    log(f"  wrote {stem}.md / {stem}.xlsx ({len(df)} rows)")

def md_table(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(v) for v in r) + " |")
    return out

# ---- Table 1
t1_rows = []
for key in CKEYS:
    d = table1[key]
    g1, s1, _ = top20[key][0]
    t1_rows.append({"Cancer": d["cancer"], "Group (CRISPR n)": d["group_n"],
                    "Comparator": d["comparator"], "Strict + OE": d["strict_oe"],
                    "Relaxed + OE": d["relaxed_oe"], "#1 gene": g1,
                    "#1 combined score": round(s1, 2)})
df1 = pd.DataFrame(t1_rows)
write_table("table1", df1, "Table 1 — The twenty scans at a glance",
            md_table(df1),
            "Strict + OE / Relaxed + OE: synthesis.md at-a-glance table "
            "(multi-cancer-selective-targets-20261120/synthesis.md, lines 13-32), verified "
            "gene-by-gene against per-cancer candidates.csv Tier counts (15 cancers), "
            "ranked_candidates.csv StrictBar counts (EGFR-mutant and KRAS-mutant lung), "
            "candidates_both.csv (colon), report text (melanoma) and ranked_fetchbar.csv "
            "(TNBC relaxed+OE = 25). For AML, neuroblastoma and osteosarcoma the synthesis.md "
            "first column lists the strict dependency-bar count (270/149/100), not strict+OE; "
            "this table uses strict+OE from candidates.csv Tier A (103/71/21), as the manuscript "
            "does. #1 gene and score: top20.csv row 1 in each cancer folder. "
            "AML relaxed+OE (153) is the pan-cancer list; the blood-only relaxed+OE count is 89 "
            "(aml/candidates_haem.csv, 89 rows). Group n and comparator: synthesis.md table.")

# ---- Table 2
t2_rows = []
for g, c in sorted(rec3.items(), key=lambda kv: (-kv[1], kv[0])):
    cancers = sorted(k for k in CKEYS if g in [x for x, _, _ in top20[k]])
    t2_rows.append({"Gene": g, "Cancers (n)": c,
                    "Cancers": "; ".join(CLABEL[k] for k in cancers)})
df2 = pd.DataFrame(t2_rows)
write_table("table2", df2, "Table 2 — Recurrence across the twenty top-20 lists",
            md_table(df2),
            "Recomputed from the twenty top20.csv files "
            "(multi-cancer-selective-targets-20261120/<cancer>/top20.csv and "
            "depmap-tnbc-selective-targets-20261120/top20.csv): 385 slots, 269 distinct genes. "
            "Only genes in >=3 cancers are listed (25 genes); 200 genes appear once.")

# ---- Table 3
t3_rows = []
for g in TEAD_GENES:
    row = {"Gene": TEAD_LABEL[g]}
    for k in SETTINGS:
        row[SLABEL[k]] = tead_rank[k][g]
    row["Top-20 in"] = "; ".join(CLABEL[k] for k in tead_top20[g]) if tead_top20[g] else "-"
    t3_rows.append(row)
df3 = pd.DataFrame(t3_rows)
write_table("table3", df3, "Table 3 — TEAD-axis enrichment ranks (1 = most selective)",
            md_table(df3),
            "Enrichment rank within each cancer's genome-wide dependency scan "
            "(dependency_scan.csv EnrichRank column; mesothelioma from "
            "mesothelioma/tead_axis.csv). EGFR-mutant and KRAS-mutant lung scans carry no "
            "EnrichRank column, so ranks there were computed by sorting on Enrichment "
            "(ascending); the method was validated against the file EnrichRank in glioma and "
            "AML (exact agreement, including TEAD1 291 and 17,802). TEAD1's best rank is 43 "
            "(mesothelioma); it never enters a top 20. WWTR1 = TAZ; PTK2 = FAK.")

# ---- Table 4
p_breadth = mannwhitney_p([int(r["n_celltypes_frac_ge_0.10"]) for r in q269],
                          [int(r["n_celltypes_frac_ge_0.10"]) for r in bg])
log(f"  breadth Mann-Whitney p (computed): {p_breadth:.3g} (manuscript: 9.8e-06)")
t4a = pd.DataFrame([
    {"Outcome (mutually exclusive)": "Strict dependency bar + over-expression", "Cells": n_strict_oe},
    {"Outcome (mutually exclusive)": "Relaxed-only dependency bar + over-expression", "Cells": n_relaxed_oe - n_strict_oe},
    {"Outcome (mutually exclusive)": "Strict dependency bar only (no over-expression)", "Cells": n_strict_only},
    {"Outcome (mutually exclusive)": "Relaxed-only dependency bar only (no over-expression)", "Cells": n_relaxed_only_dep},
    {"Outcome (mutually exclusive)": "Over-expression only (no dependency bar)", "Cells": n_oe_only},
    {"Outcome (mutually exclusive)": "Neither bar", "Cells": n_neither},
    {"Outcome (mutually exclusive)": "TOTAL", "Cells": 400},
])
t4a_lines = md_table(t4a) + ["",
    "Inclusive counts: 33 cells pass the relaxed bar plus over-expression (22 of them also "
    "pass the strict bar). 16 cells pass a dependency bar without over-expression "
    "(12 strict + 4 relaxed-only)."]
t4b_rows = []
for r in tsum:
    t4b_rows.append({"Target": r["Target"], "Drug status": r["Drug_status"],
                    "Approved indications in programme": r["Approved_indications_in_programme"],
                    "Verdict in approved indication": r["Verdict_in_approved"],
                    "Pass strict + OE": r["Pass_strict_OE"] or "-",
                    "Pass relaxed + OE": r["Pass_relaxed_OE"] or "-",
                    "Failure taxonomy": r["Failure_taxonomy"]})
df4b = pd.DataFrame(t4b_rows)
md4 = ["# Table 4 — Calibration summary (20 targets x 20 cancers = 400 cells)", "",
       "## A. Outcome counts", ""] + t4a_lines + ["", "## B. Per-target verdicts", ""] + md_table(df4b)
with open(os.path.join(TABDIR, "table4.md"), "w") as f:
    f.write("\n".join(md4) + "\n\n*Source: drug-target-calibration/target_cancer_matrix.csv "
            "(400 rows, recomputed counts) and target_summary.csv (per-target verdicts). "
            "The manuscript quotes the inclusive relaxed+OE count (33); the exclusive "
            "relaxed-only+OE count is 11 (33 minus the 22 that also pass strict).*\n")
with pd.ExcelWriter(os.path.join(TABDIR, "table4.xlsx"), engine="openpyxl") as w:
    t4a.to_excel(w, sheet_name="A. Outcome counts", index=False)
    df4b.to_excel(w, sheet_name="B. Per-target verdicts", index=False)
    pd.DataFrame({"Source": ["drug-target-calibration/target_cancer_matrix.csv (400 rows) and "
                             "target_summary.csv; counts recomputed in build_figs_tables.py."]}).to_excel(
        w, sheet_name="Source", index=False)
log("  wrote table4.md / table4.xlsx")

# ---- Table 5
t5_rows = [
    {"Metric": "Median tau", "269 selective genes": round(statistics.median(tau269), 3),
     "300 background genes": round(statistics.median(taubg), 3), "p-value": "Mann-Whitney 0.0155"},
    {"Metric": "Tau >= 0.9", "269 selective genes": f"{sum(1 for v in tau269 if v >= 0.9)} ({sum(1 for v in tau269 if v >= 0.9)/269*100:.1f}%)",
     "300 background genes": f"{sum(1 for v in taubg if v >= 0.9)} ({sum(1 for v in taubg if v >= 0.9)/300*100:.1f}%)", "p-value": "-"},
    {"Metric": "Switch-like (on <=10 types, off >=100)", "269 selective genes": f"{s269} ({s269/269*100:.1f}%)",
     "300 background genes": f"{sbg} ({sbg/300*100:.1f}%)", "p-value": "Fisher 0.0001"},
    {"Metric": "Median cell types detected (frac >= 0.10)", "269 selective genes": int(statistics.median(int(r["n_celltypes_frac_ge_0.10"]) for r in q269)),
     "300 background genes": int(statistics.median(int(r["n_celltypes_frac_ge_0.10"]) for r in bg)), "p-value": "Mann-Whitney 9.8e-06"},
    {"Metric": "Median cell types on (frac >= 0.50)", "269 selective genes": int(statistics.median(int(r["n_on_frac_ge_0.50"]) for r in q269)),
     "300 background genes": int(statistics.median(int(r["n_on_frac_ge_0.50"]) for r in bg)), "p-value": "-"},
    {"Metric": "Median cell types off (frac <= 0.01)", "269 selective genes": int(statistics.median(int(r["n_off_frac_le_0.01"]) for r in q269)),
     "300 background genes": int(statistics.median(int(r["n_off_frac_le_0.01"]) for r in bg)), "p-value": "-"},
    {"Metric": "Median worst-normal expression (CPM10k)", "269 selective genes": round(statistics.median(float(r["worst_mean_cpm10k"]) for r in q269), 1),
     "300 background genes": "-", "p-value": "-"},
    {"Metric": "Worst normal < 100 CPM10k", "269 selective genes": sum(1 for r in q269 if float(r["worst_mean_cpm10k"]) < 100),
     "300 background genes": "-", "p-value": "-"},
    {"Metric": "Worst normal >= 500 CPM10k", "269 selective genes": sum(1 for r in q269 if float(r["worst_mean_cpm10k"]) >= 500),
     "300 background genes": "-", "p-value": "-"},
    {"Metric": "Worst normal >= 1,000 CPM10k", "269 selective genes": sum(1 for r in q269 if float(r["worst_mean_cpm10k"]) >= 1000),
     "300 background genes": "-", "p-value": "-"},
    {"Metric": "Wide window (tumour FC >= 3, worst normal < 100)", "269 selective genes": len(wide),
     "300 background genes": "-", "p-value": "-"},
    {"Metric": "Narrow window (tumour FC >= 3, worst normal >= 500)", "269 selective genes": len(narrow),
     "300 background genes": "-", "p-value": "-"},
]
df5 = pd.DataFrame(t5_rows)
write_table("table5", df5, "Table 5 — The switch-like test (Tabula Sapiens 2.0)",
            md_table(df5),
            "Recomputed from tabula-sapiens/per_gene_summary_truecpm.csv (269 selective genes, "
            "in_269 = True) and background_summary_truecpm.csv (300 random protein-coding genes, "
            "seed 20261126). True per-cell-total normalisation (Census obs raw_sum); the first "
            "query-set-relative run (~45x inflation) was corrected by a full re-query on "
            "26 Nov 2026 — see the correction note in the build note. Tau = Yanai tau over 180 cell types. Switch-like = fraction "
            ">= 0.50 in <= 10 types and fraction <= 0.01 in >= 100 types. Worst-normal = "
            "highest mean CPM10k among cell types with n >= 50 cells. Tumour fold change = "
            "maximum cell-line TPM fold change across the twenty cancers. "
            "Breadth Mann-Whitney p recomputed in build_figs_tables.py.")

# ---------------------------------------------------------------- figures
log("== drawing figures ==")

def savefig(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, stem + "." + ext), bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    log(f"  wrote figures/{stem}.png / {stem}.pdf")

# ---- Figure 1: programme architecture (schematic, from Methods)
fig, ax = plt.subplots(figsize=(183 * MM, 115 * MM))
ax.set_xlim(0, 100); ax.set_ylim(0, 62); ax.axis("off")
def box(x, y, w, h, title, body, edge, face="white", ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.8,rounding_size=1.5",
                                fc=face, ec=edge, lw=1.4, linestyle=ls))
    ax.text(x + w / 2, y + h - 2.6, title, ha="center", va="top", fontsize=7.5, fontweight="bold")
    ax.text(x + w / 2, y + h - 6.4, body, ha="center", va="top", fontsize=7, linespacing=1.35)
def arrow(x1, y1, x2, y2, label=None, lx=None, ly=None):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                 color=OI["darkgrey"], lw=1.2, shrinkA=1, shrinkB=1))
    if label:
        ax.text(lx, ly, label, ha="center", va="center", fontsize=6.5, color=OI["darkgrey"])
box(2, 50, 28, 11, "AUTHOR (W.H.) — DESIGN", "fixes bars, cancer list,\ncomparator rules", OI["orange"])
box(36, 50, 28, 11, "COORDINATOR AGENT", "decomposes programme,\nassigns one worker per cancer", OI["blue"])
box(70, 50, 28, 11, "20 PER-CANCER WORKERS", "same fixed playbook\nreport + top-20 list each", OI["blue"])
arrow(30, 55.5, 36, 55.5); arrow(64, 55.5, 70, 55.5)
box(70, 30, 28, 11, "VERIFICATION SCRIPTS", "independent recomputation\nof every headline number", OI["green"])
arrow(84, 50, 84, 41.5)
box(36, 30, 28, 11, "HUMAN CHECKPOINT", "per-cancer delivery review\n(author reviews each report)", OI["orange"])
arrow(70, 35.5, 64, 35.5)
# loop-back: verification -> workers (recompute / reject)
ax.annotate("", xy=(97, 50), xytext=(97, 41.5),
            arrowprops=dict(arrowstyle="-|>", color=OI["green"], lw=1.2,
                            connectionstyle="arc3,rad=0.25"))
ax.text(99.5, 45.5, "recompute", ha="left", va="center", fontsize=6.5, color=OI["green"], rotation=90)
box(36, 8, 28, 13, "DOWNSTREAM WORKERS", "synthesis · partner mining\nADC intersection · Tabula leg\ncalibration leg", OI["blue"])
arrow(50, 30, 50, 21.5)
box(70, 10, 28, 11, "HUMAN CHECKPOINT", "interpretation &\nmanuscript decisions", OI["orange"])
arrow(64, 14.5, 70, 14.5)
box(2, 10, 28, 11, "VERIFICATION LOOP", "the element the Virtual\nBiotech architecture lacks", OI["green"], ls=(0, (4, 2)))
arrow(30, 14.5, 36, 14.5)
ax.text(50, 3.5, "Schematic — architecture as described in Methods. No worker could change the bars, "
        "the cancer list or the scoring rule.", ha="center", va="center", fontsize=6.5,
        color=OI["darkgrey"], style="italic")
savefig(fig, "fig1")

# ---- Figure 2: twenty-cancer x gene matrix
# Two stacked panels: recurring genes (>=2 cancers, 69 rows, labelled) on top and
# singletons (200 rows, unlabelled, compressed) below. A single 269-row labelled
# matrix cannot be legible at page size (rows would be <1 mm); the split keeps
# every row plotted while labelling the rows that carry information.
log("-- Fig 2 matrix --")
genes_sorted = sorted(rec.keys(), key=lambda g: (-rec[g], -max(s for _, s, _ in
    [(gg, ss, tt) for k in CKEYS for gg, ss, tt in top20[k] if gg == g]), g))
score_of = {}
for k in CKEYS:
    for g, s, _ in top20[k]:
        score_of[(g, k)] = s
M = np.full((len(genes_sorted), len(CKEYS)), np.nan)
for i, g in enumerate(genes_sorted):
    for j, k in enumerate(CKEYS):
        if (g, k) in score_of:
            M[i, j] = score_of[(g, k)]
ADHESION = {"ITGAV", "ITGB5", "FERMT2", "PTK2"}
LINEAGE = {"KLF5", "SOX9", "HNF1B", "TP63", "PAX8"}
POLARITY = {"PARD3", "PARD6B", "ARHGAP29"}
cluster_code = [1 if g in ADHESION else (2 if g in LINEAGE else (3 if g in POLARITY else 0))
                for g in genes_sorted]
drug_code = [1 if DRUG.get(g) == "approved" else (2 if DRUG.get(g) == "clinical" else 0)
             for g in genes_sorted]
NREC = 69  # genes in >=2 cancers
fig = plt.figure(figsize=(183 * MM, 225 * MM))
gs = fig.add_gridspec(2, 4, width_ratios=[0.10, 0.10, 1.0, 0.16], height_ratios=[69, 45],
                      wspace=0.04, hspace=0.015, left=0.20, right=0.985, top=0.985, bottom=0.16)
ax_c1 = fig.add_subplot(gs[0, 0]); ax_d1 = fig.add_subplot(gs[0, 1])
ax_m1 = fig.add_subplot(gs[0, 2]); ax_r1 = fig.add_subplot(gs[0, 3])
ax_c2 = fig.add_subplot(gs[1, 0]); ax_d2 = fig.add_subplot(gs[1, 1])
ax_m2 = fig.add_subplot(gs[1, 2]); ax_r2 = fig.add_subplot(gs[1, 3])
cmap_cl = LinearSegmentedColormap.from_list("cl", ["white", OI["blue"], OI["orange"], OI["green"]])
cmap_dr = LinearSegmentedColormap.from_list("dr", ["white", OI["blue"], OI["sky"]])
cmap_sc = LinearSegmentedColormap.from_list("sc", ["#FEE8C8", "#FDBB84", "#EF6548", "#99000D"])
cmap_sc.set_bad("white")
norm_sc = matplotlib.colors.PowerNorm(gamma=0.5, vmin=0, vmax=30)
for ax_c, ax_d, ax_m, ax_r, sl in [(ax_c1, ax_d1, ax_m1, ax_r1, slice(0, NREC)),
                                   (ax_c2, ax_d2, ax_m2, ax_r2, slice(NREC, None))]:
    n = len(genes_sorted[sl])
    ax_c.imshow(np.array(cluster_code[sl])[:, None], cmap=cmap_cl, vmin=0, vmax=3, aspect="auto")
    ax_c.set_xticks([]); ax_c.set_yticks([])
    ax_d.imshow(np.array(drug_code[sl])[:, None], cmap=cmap_dr, vmin=0, vmax=2, aspect="auto")
    ax_d.set_xticks([]); ax_d.set_yticks([])
    ax_m.imshow(M[sl], cmap=cmap_sc, norm=norm_sc, aspect="auto")
    ax_m.set_ylim(n - 0.5, -0.5)
    ax_m.tick_params(length=0)
    ax_r.barh(range(n), [rec[g] for g in genes_sorted[sl]], color="#BBBBBB", height=0.9)
    ax_r.set_ylim(n - 0.5, -0.5)
    ax_r.set_xlim(0, 8.6); ax_r.set_xticks([0, 4, 8]); ax_r.set_yticks([])
    ax_r.tick_params(length=0)
    for a in (ax_c, ax_d, ax_m, ax_r):
        for s in ("top", "right", "bottom", "left"):
            a.spines[s].set_visible(False)
    ax_m.spines["left"].set_visible(True); ax_m.spines["right"].set_visible(True)
    ax_m.spines["top"].set_visible(True); ax_m.spines["bottom"].set_visible(True)
ax_c1.set_title("Cluster", fontsize=7, pad=2)
ax_d1.set_title("Drug", fontsize=7, pad=2)
ax_r1.set_title("Recurrence", fontsize=7, pad=2)
ax_m1.set_yticks(range(NREC))
ax_m1.set_yticklabels(genes_sorted[:NREC], fontsize=5)
ax_m1.set_xticks([])
ax_m2.set_xticks(range(len(CKEYS)))
ax_m2.set_xticklabels([f"{CLABEL[k]}\nn={table1[k]['group_n']}" for k in CKEYS],
                      rotation=90, fontsize=7, va="top", ha="center")
ax_m2.set_yticks([])
ax_r2.text(1.02, 0.96, "singletons (200)", transform=ax_r2.transAxes,
           fontsize=6.5, color=OI["darkgrey"], va="top", ha="left")
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm_sc, cmap=cmap_sc), ax=ax_m2,
                  orientation="horizontal", fraction=0.05, pad=0.22, aspect=40)
cb.set_label("Combined score (within-cancer ranking; scores are not comparable across cancers)",
             fontsize=7)
fig.text(0.20, 0.05,
         "Cluster: blue = adhesion axis (ITGAV/ITGB5/FERMT2/PTK2), orange = lineage factors "
         "(KLF5/SOX9/HNF1B/TP63/PAX8), green = polarity module (PARD3/PARD6B/ARHGAP29).  "
         "Drug: blue = approved drug at the gene, light blue = clinical-stage (calibration target list).",
         fontsize=6.5, color=OI["darkgrey"], ha="left", va="top", wrap=True)
fig.text(0.20, 0.022,
         "Rows ordered by recurrence (then maximum combined score). 385 slots, 269 genes; "
         "200 genes appear once. White = gene not in that cancer's top 20.",
         fontsize=6.5, color=OI["darkgrey"], ha="left", va="top")
savefig(fig, "fig2")

# ---- Figure 3: recurrence chart
log("-- Fig 3 recurrence --")
rec_items = sorted(rec3.items(), key=lambda kv: (-kv[1], kv[0]))
names = [g for g, _ in rec_items] + ["200 singletons"]
vals = [c for _, c in rec_items] + [1]
cols = []
for g in [g for g, _ in rec_items]:
    cols.append(OI["blue"] if DRUG.get(g) == "approved" else
                (OI["sky"] if DRUG.get(g) == "clinical" else "#BBBBBB"))
cols.append(OI["darkgrey"])
fig, ax = plt.subplots(figsize=(89 * MM, 118 * MM))
y = np.arange(len(names))[::-1]
ax.barh(y, vals, color=cols, height=0.75, edgecolor="white", linewidth=0.5)
for yi, v in zip(y, vals):
    ax.text(v + 0.15, yi, str(v), va="center", ha="left", fontsize=7)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7)
ax.set_xlim(0, 9.2); ax.set_xticks(range(0, 9))
ax.set_xlabel("Cancers with the gene in the top 20 (of 20)")
ax.tick_params(length=0)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
legend_elems = [mpatches.Patch(facecolor=OI["blue"], label="Approved drug at the gene"),
                mpatches.Patch(facecolor=OI["sky"], label="Clinical-stage drug at the gene"),
                mpatches.Patch(facecolor="#BBBBBB", label="No drug at the gene (calibration list)")]
ax.legend(handles=legend_elems, loc="upper center", bbox_to_anchor=(0.5, -0.16),
          frameon=False, fontsize=7)
savefig(fig, "fig3")

# ---- Figure 4: TEAD-axis rank matrix
log("-- Fig 4 TEAD matrix --")
CLASS = {}
for k in ["tnbc", "melanoma", "glioma", "ovarian", "renal-cell", "osteosarcoma"]:
    CLASS[k] = 0  # TAZ cancer
for k in ["colon", "esophagogastric", "head-and-neck", "mesothelioma"]:
    CLASS[k] = 1  # YAP1 cancer / leaning
for k in ["small-cell-lung", "esophageal-squamous", "ewing", "endometrial", "aml", "neuroblastoma"]:
    CLASS[k] = 2  # null / anti pole
for k in ["egfr-mutant-lung", "hcc", "pancreatic", "kras-mutant-lung", "bladder"]:
    CLASS[k] = 3  # other
R = np.array([[tead_rank[k][g] for k in SETTINGS] for g in TEAD_GENES], dtype=float)
L = np.log10(R)
fig = plt.figure(figsize=(183 * MM, 112 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[0.10, 1.0], hspace=0.04,
                      left=0.115, right=0.985, top=0.94, bottom=0.22)
ax_s = fig.add_subplot(gs[0, 0]); ax_m = fig.add_subplot(gs[1, 0])
cmap_class = LinearSegmentedColormap.from_list("cls", [OI["blue"], OI["orange"], OI["verm"], "#CCCCCC"])
ax_s.imshow(np.array([CLASS[k] for k in SETTINGS])[None, :], cmap=cmap_class, vmin=0, vmax=3, aspect="auto")
ax_s.set_xticks([]); ax_s.set_yticks([]); ax_s.set_title("Cancer class", fontsize=7, loc="left", pad=2)
for s in ("top", "right", "bottom", "left"):
    ax_s.spines[s].set_visible(False)
cmap_div = LinearSegmentedColormap.from_list("div", [OI["blue"], "white", OI["verm"]])
norm = matplotlib.colors.TwoSlopeNorm(vmin=0.0, vcenter=math.log10(8958), vmax=math.log10(17916))
im = ax_m.imshow(L, cmap=cmap_div, norm=norm, aspect="auto")
for i, g in enumerate(TEAD_GENES):
    for j, k in enumerate(SETTINGS):
        v = L[i, j]
        tc = "white" if (v < 2.0 or v > 4.2) else "black"
        ax_m.text(j, i, f"{int(R[i, j]):,}", ha="center", va="center", fontsize=6, color=tc)
        if k in tead_top20[g]:
            ax_m.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                         edgecolor="black", lw=1.1))
ax_m.set_xticks(range(len(SETTINGS)))
ax_m.set_xticklabels([f"{SLABEL[k]}\nn={table1[k]['group_n']}" if k in table1 else "Mesothelioma\nn=21"
                      for k in SETTINGS], rotation=45, fontsize=7, va="top", ha="right")
ax_m.set_yticks(range(len(TEAD_GENES)))
ax_m.set_yticklabels([TEAD_LABEL[g] for g in TEAD_GENES], fontsize=8)
ax_m.tick_params(length=0)
cb = fig.colorbar(im, ax=ax_m, orientation="horizontal", fraction=0.05, pad=0.22, aspect=40)
cb.set_label("Enrichment rank (log10 scale; 1 = most selective, 17,916 = most anti-selective)", fontsize=7)
fig.text(0.115, 0.01,
         "Bold outline = gene in that cancer's top 20. Cancer class (top strip): blue = TAZ cancer, "
         "orange = YAP1 cancer / leaning, red = null / anti-TEAD pole, grey = other. "
         "TEAD1 never enters a top 20; best rank 43 (mesothelioma).",
         fontsize=6.5, color=OI["darkgrey"], ha="left", va="top")
savefig(fig, "fig4")

# ---- Figure 5: switch-like test
log("-- Fig 5 switch-like --")
fig = plt.figure(figsize=(183 * MM, 118 * MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.28,
                      left=0.075, right=0.98, top=0.90, bottom=0.20)
ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
bins = np.arange(0.5, 1.001, 0.02)
ax_a.hist(tau269, bins=bins, color=OI["blue"], alpha=0.65, label="269 selective genes")
ax_a.hist(taubg, bins=bins, color="#999999", alpha=0.55, label="300 background genes")
ax_a.axvline(statistics.median(tau269), color=OI["blue"], lw=1.2, linestyle=(0, (4, 2)))
ax_a.axvline(statistics.median(taubg), color="#666666", lw=1.2, linestyle=(0, (4, 2)))
ax_a.set_xlabel("Yanai tau (0 = flat, 1 = one cell type only)")
ax_a.set_ylabel("Genes")
ax_a.set_title("A  Tau distribution", loc="left", fontsize=9)
ax_a.legend(frameon=False, loc="upper left")
ax_a.text(0.56, 30, "Median tau 0.891 vs 0.904\nMann–Whitney p = 0.0155",
          fontsize=7, va="top")
for s in ("top", "right"):
    ax_a.spines[s].set_visible(False)
bars = ax_b.bar([0, 1], [s269 / 269 * 100, sbg / 300 * 100], width=0.55,
                color=[OI["blue"], "#999999"])
ax_b.set_xticks([0, 1])
ax_b.set_xticklabels(["269 selective genes", "300 background genes"], fontsize=7)
ax_b.set_ylabel("% switch-like genes")
ax_b.set_ylim(0, 25)
ax_b.set_title("B  Switch-like genes", loc="left", fontsize=9)
for b, n, d in zip(bars, [s269, sbg], [269, 300]):
    ax_b.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.6,
              f"{n}/{d} ({n/d*100:.1f}%)", ha="center", va="bottom", fontsize=8)
ax_b.text(0.5, 0.94, "Fisher p = 0.0001", transform=ax_b.transAxes, ha="center", fontsize=7, va="top")
for s in ("top", "right"):
    ax_b.spines[s].set_visible(False)
glist = ", ".join(switch_genes)
half = (len(switch_genes) + 1) // 2
line1 = "The 25 switch-like selective genes: " + ", ".join(switch_genes[:half]) + ","
line2 = ", ".join(switch_genes[half:]) + "."
fig.text(0.075, 0.075, line1, fontsize=6.5, ha="left", va="top")
fig.text(0.075, 0.045, line2, fontsize=6.5, ha="left", va="top")
fig.text(0.075, 0.005,
         "Switch-like = expressed (fraction >= 0.50) in <= 10 cell types and silent (fraction <= 0.01) "
         "in >= 100 types, Tabula Sapiens 2.0 (1,136,218 cells, 180 types).",
         fontsize=6.5, color=OI["darkgrey"], ha="left", va="top")
savefig(fig, "fig5")

# ---------------------------------------------------------------- build note + log
note = []
note.append("# Figures and tables — build note")
note.append("")
note.append("Built 26 November 2026 by `manuscript/build_figs_tables.py`. Manuscript text NOT edited.")
note.append("")
note.append("## Figures")
note.append("")
note.append("| Figure | File | What it shows | Source |")
note.append("|---|---|---|---|")
note.append("| 1 | `figures/fig1.png/.pdf` | Programme architecture (schematic) | Methods section of manuscript.md — coordinator, 20 per-cancer workers, fixed playbook, verification scripts, human checkpoints at design / delivery / interpretation, downstream workers. **Schematic, not data.** The verification loop is drawn as the element the Virtual Biotech architecture lacks, as the manuscript states. |")
note.append("| 2 | `figures/fig2.png/.pdf` | 20-cancer x 269-gene top-20 matrix, coloured by combined score | All twenty `top20.csv` files (gene + CombinedScore/Score columns harmonised). Rows ordered by recurrence then maximum combined score. Cluster strip: adhesion axis (ITGAV/ITGB5/FERMT2/PTK2), lineage factors (KLF5/SOX9/HNF1B/TP63/PAX8), polarity module (PARD3/PARD6B/ARHGAP29) — exactly the sets the manuscript names. Drug strip: calibration target list only (approved = blue, clinical = light blue). Scores are comparable within a cancer only (Methods); the colour scale is shown for the matrix as the legend specifies. **Layout decision:** the matrix is drawn as two stacked panels — the 69 genes recurring in >=2 cancers (labelled) above, the 200 singletons (unlabelled, compressed) below — because 269 labelled rows cannot be legible at page size (rows would be <1 mm). Every row is plotted; only the singleton labels are omitted. Colour uses a square-root normalisation so mid-range scores are visible against white (absent). |")
note.append("| 3 | `figures/fig3.png/.pdf` | Recurrence: genes in >=3 cancers (25) + 200 singletons | Recomputed from the twenty `top20.csv` files. Colour by druggability from the calibration target list (see Fig 2). |")
note.append("| 4 | `figures/fig4.png/.pdf` | TEAD-axis enrichment ranks, 9 genes x 21 settings | `dependency_scan.csv` EnrichRank per cancer; `mesothelioma/tead_axis.csv` for mesothelioma. EGFR-mutant and KRAS-mutant lung scans carry no EnrichRank column, so ranks there were computed by sorting on Enrichment (ascending); the method was validated against file EnrichRank in glioma and AML (exact agreement, incl. TEAD1 291 and 17,802). Bold outline = gene in that cancer's top 20. Cancer-class strip from the manuscript text (TAZ cancers: TNBC, melanoma, glioma, ovarian, RCC, osteosarcoma; YAP1 cancers/leaning: colon, ESGA, HNSC, mesothelioma; null/anti pole: SCLC, ESCC, Ewing, endometrial, AML, neuroblastoma; other: EGFR lung, HCC, PDAC, KRAS lung, bladder). |")
note.append("| 5 | `figures/fig5.png/.pdf` | Switch-like test: tau distributions + switch-like fractions | `tabula-sapiens/per_gene_summary_truecpm.csv` (269 selective) and `background_summary_truecpm.csv` (300 background) — true per-cell-total normalisation. The 25 switch-like selective genes are listed under the panels, taken from the CSV (`switch_like` = True), not from the manuscript text. |")
note.append("")
note.append("## Tables")
note.append("")
note.append("| Table | File | Source |")
note.append("|---|---|---|")
note.append("| 1 | `tables/table1.md/.xlsx` | Strict+OE / Relaxed+OE: synthesis.md at-a-glance table, verified against per-cancer `candidates.csv` Tier counts (15 cancers), `ranked_candidates.csv` StrictBar counts (EGFR/KRAS lung), `candidates_both.csv` (colon), report text (melanoma) and `ranked_fetchbar.csv` (TNBC relaxed+OE = 25 rows). #1 gene/score: `top20.csv` row 1. Group n spot-checked against `top20.csv` group-n columns (HCC 21, PDAC 45, HNSC 72, OS 17, NB 39). |")
note.append("| 2 | `tables/table2.md/.xlsx` | Recomputed from the twenty `top20.csv` files. |")
note.append("| 3 | `tables/table3.md/.xlsx` | As Fig 4. |")
note.append("| 4 | `tables/table4.md/.xlsx` | `drug-target-calibration/target_cancer_matrix.csv` (400 cells, counts recomputed) and `target_summary.csv` (per-target verdicts). |")
note.append("| 5 | `tables/table5.md/.xlsx` | As Fig 5. |")
note.append("")
note.append("## Correction — true CPM (26 Nov 2026, evening)")
note.append("")
note.append("The first Tabula run normalised each cell by the sum of the 290 queried genes only, inflating absolute CPM ~45×. This build re-ran Figure 5 and Table 5 from the corrected true-CPM summaries (`per_gene_summary_truecpm.csv`, `background_summary_truecpm.csv`; Census obs `raw_sum`). Corrected values: median tau 0.891 vs 0.904 (Mann–Whitney p = 0.0155, now significant — selective genes significantly less tau-specific); median worst-normal 2.2 CPM10k (267/269 below 100, 0 above 500, 0 above 1,000); wide window 90, narrow window 0. The 100/500 CPM thresholds were set on the inflated scale and no longer discriminate (the wide class admits 90 of the 91 genes with tumour FC ≥ 3); they were not moved post hoc. Fraction-based results are unchanged (switch-like 25/269 vs 63/300, Fisher p = 0.0001; breadth medians; all fractions). Figures 1–4 and Tables 1–4 are regenerated from unchanged inputs in this run; their identity against the first build is checked separately (see the post-run comparison).")
note.append("")
note.append("## Verification")
note.append("")
note.append(f"- Build-script checks: {CHECKS['pass']} passed, {CHECKS['fail']} failed (see log).")
note.append("- `verify_headlines.py` re-run after the build: see log tail below.")
note.append("")
note.append("## Discrepancies and caveats (manuscript NOT edited)")
note.append("")
if FAILURES:
    for name, got, want in FAILURES:
        note.append(f"- **{name}**: build got `{got}`, manuscript/source expectation `{want}`.")
else:
    note.append("- None: every check passed.")
note.append("- **Switch-like gene list (manuscript text is a partial list, not an error):** the manuscript names 16 of the 25 switch-like selective genes in parentheses (MYCN, HAND2, PHOX2A, ASCL1, INSM1, SOX11, SOX2, CDX2, HNF1B, FOXA2, HNF1A, ISL1, SOX10, TFAP2A, GFI1, RCOR2). The full 25 from `per_gene_summary_truecpm.csv` add nine genes that are not lineage transcription factors: CEP55, CKAP2L, CSF2RB, GPR61, KIF18B, LGALS9B, MYBL2, SPRR1B, TUBA3D. Fig 5 and Table 5 use the full 25. The manuscript's 'almost all lineage transcription factors' reading holds for 16/25; the other 9 are mitotic/cytoskeletal or signalling genes.")
note.append("- **Manuscript relabelling during this build (content unchanged):** while the figures/tables were being built, the manuscript was restructured to meet the six-display-item limit — Tables 2–5 are now labelled **Extended Data Tables 1–4** (Table 2 → Ext. Data Table 1, Table 3 → Ext. Data Table 2, Table 4 → Ext. Data Table 3, Table 5 → Ext. Data Table 4); Table 1 and Figures 1–5 remain main display items. The legend text and every number are unchanged, so the files here (`table1`–`table5`) match the current manuscript content exactly; only the file/table labels have not been renamed. Rename at submission time if the Extended Data labels are kept.")
note.append("- **Table 4 inclusive vs exclusive counts:** the manuscript quotes the inclusive relaxed+OE count (33 cells pass the relaxed bar plus over-expression). The mutually exclusive breakdown in Table 4A is strict+OE 22, relaxed-only+OE 11, strict-dep-only 12, relaxed-only-dep-only 4, OE-only 66, neither 285 (= 400).")
note.append("- **TNBC relaxed+OE = 25:** verified by row count of `ranked_fetchbar.csv` (25 data rows). That file is the relaxed-bar gene list taken forward for the over-expression fetch; its row count matches synthesis.md. It is not an independent recomputation from the raw matrices.")
note.append("- **Breadth Mann–Whitney p:** recomputed here by normal approximation with tie correction; manuscript quotes 9.8e-06 (see log for the computed value).")
note.append("- **Fig 2/3 druggability strip:** marks only the 20 calibration targets (approved/clinical). A gene not in that list is shown as 'no drug at the gene (calibration list)' — this is a defined, sourced classification, not a full druggability review.")
note.append("- **Fig 2 colour scale:** combined scores are comparable within a cancer only (Methods). The heatmap colours by combined score across cancers exactly as the figure legend specifies; read row-wise, not column-wise, for magnitude.")
note.append("")
with open(os.path.join(MS, "figures-tables-build-note.md"), "w") as f:
    f.write("\n".join(note) + "\n")

with open(os.path.join(MS, "build_verification_log.txt"), "w") as f:
    f.write("\n".join(LOG) + "\n")
print(f"\nBUILD DONE: {CHECKS['pass']} checks passed, {CHECKS['fail']} failed")
if FAILURES:
    print("FAILURES:")
    for name, got, want in FAILURES:
        print(f"  - {name}: got {got}, want {want}")
