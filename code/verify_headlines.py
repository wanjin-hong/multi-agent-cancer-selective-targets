#!/usr/bin/env python3
"""verify_headlines.py — manuscript (26 Nov 2026).

Recomputes every number cited in manuscript.md from the source CSVs, through a
separate code path from the pipeline and from the per-cancer verification scripts.
All checks must pass. Exit 0 = pass, 1 = fail.
"""
import csv, math, os, re, statistics, sys
from collections import Counter

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
TNBC = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-selective-targets-20261120")
HERE = os.path.dirname(os.path.abspath(__file__))

passed, failed = 0, 0
def check(name, got, want, tol=None):
    global passed, failed
    ok = abs(got - want) <= tol if tol is not None else got == want
    if ok:
        passed += 1
        print(f"  PASS {name}: {got}")
    else:
        failed += 1
        print(f"  FAIL {name}: got {got}, want {want}")

def check_set(name, got, want):
    global passed, failed
    if set(got) == set(want):
        passed += 1
        print(f"  PASS {name}: {sorted(got)}")
    else:
        failed += 1
        print(f"  FAIL {name}: got {sorted(got)}, want {sorted(want)}")

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

# ---------- top-20 lists ----------
print("== programme overview (top20.csv files + TNBC report table) ==")
cancers = ["egfr-mutant-lung", "hcc", "melanoma", "colon", "glioma", "pancreatic",
           "kras-mutant-lung", "head-and-neck", "ovarian", "small-cell-lung",
           "esophagogastric", "esophageal-squamous", "bladder", "renal-cell",
           "ewing", "endometrial", "aml", "neuroblastoma", "osteosarcoma"]
lists = {}
txt = open(os.path.join(TNBC, "report.md")).read()
m = re.search(r"\| Rank \| Gene \| Tier \| TNBC mean.*?((?:\n\| \d+ \|[^\n]+){20})", txt, re.S)
genes = []
for row in m.group(1).strip().split("\n"):
    parts = [p.strip() for p in row.split("|")]
    genes.append(parts[2].split(" (")[0])
lists["tnbc"] = genes
check("TNBC top-20 length", len(genes), 20)
check("TNBC #1 gene", genes[0], "MT1E")
for key in cancers:
    rows = read_csv(os.path.join(BASE, key, "top20.csv"))
    col = "Symbol" if "Symbol" in rows[0] else "Gene"
    lists[key] = [r[col].split(" (")[0].strip() for r in rows]
total = sum(len(v) for v in lists.values())
check("total top-20 slots", total, 385)
cnt = Counter(g for v in lists.values() for g in v)
check("distinct genes", len(cnt), 269)
check("singletons", sum(1 for c in cnt.values() if c == 1), 200)
for g, w in [("KLF5", 8), ("HNF1B", 6), ("WWTR1", 6), ("FERMT2", 6),
             ("EGFR", 5), ("ITGAV", 5), ("FOSL1", 4), ("PARD6B", 4),
             ("CXADR", 4), ("CDK6", 4)]:
    check(f"{g} recurrence", cnt[g], w)
check("TEAD1 in any top 20", sum(1 for v in lists.values() if "TEAD1" in v), 0)
check("TEAD2 in any top 20", sum(1 for v in lists.values() if "TEAD2" in v), 0)
check("TEAD4 in any top 20", sum(1 for v in lists.values() if "TEAD4" in v), 0)
check_set("TEAD3 top-20 cancers", [k for k, v in lists.items() if "TEAD3" in v], ["bladder"])
check_set("YAP1 top-20 cancers", [k for k, v in lists.items() if "YAP1" in v],
          ["colon", "esophagogastric"])
check_set("TAZ (WWTR1) top-20 cancers", [k for k, v in lists.items() if "WWTR1" in v],
          ["tnbc", "melanoma", "glioma", "ovarian", "renal-cell", "osteosarcoma"])
for k, rank in [("tnbc", 5), ("melanoma", 9), ("glioma", 14), ("ovarian", 9),
                ("renal-cell", 7), ("osteosarcoma", 9)]:
    check(f"TAZ rank in {k} top 20", lists[k].index("WWTR1") + 1, rank)

# FAK enrichment ranks from the dependency scans
print("== FAK enrichment ranks (dependency_scan.csv) ==")
for k, want in [("hcc", 3), ("glioma", 15), ("melanoma", 25),
                ("pancreatic", 1284), ("head-and-neck", 17845)]:
    rows = read_csv(os.path.join(BASE, k, "dependency_scan.csv"))
    col = "Symbol" if "Symbol" in rows[0] else "Gene"
    r = next(x for x in rows if x[col].split(" (")[0] == "PTK2")
    check(f"FAK EnrichRank {k}", int(r["EnrichRank"]), want)

# ---------- mesothelioma ----------
print("== mesothelioma TEAD1 check ==")
rows = read_csv(os.path.join(BASE, "mesothelioma", "tead_axis.csv"))
tead1 = next(r for r in rows if r["Gene"] == "TEAD1")
check("mesothelioma TEAD1 EnrichRank", int(tead1["EnrichRank"]), 43)
check("mesothelioma TEAD1 dep frac", float(tead1["GroupDepFrac"]), 10 / 21, 0.001)
check("mesothelioma TEAD1 q", float(tead1["q_BH"]), 0.144, 0.001)

# ---------- calibration ----------
print("== drug-target calibration (target_cancer_matrix.csv) ==")
mat = read_csv(os.path.join(BASE, "drug-target-calibration", "target_cancer_matrix.csv"))
check("calibration cells", len(mat), 400)
check("pass strict + OE", sum(1 for r in mat if r["Pass_strict_OE"] == "True"), 22)
check("pass relaxed + OE", sum(1 for r in mat if r["Pass_relaxed_OE"] == "True"), 33)
check("strict dep, no OE", sum(1 for r in mat if r["Strict_dep"] == "True" and r["OE_pass"] == "False"), 12)
check("OE only, no dep bar", sum(1 for r in mat if r["OE_pass"] == "True" and r["Strict_dep"] == "False" and r["Relaxed_dep"] == "False"), 66)
check("neither", sum(1 for r in mat if r["OE_pass"] == "False" and r["Strict_dep"] == "False" and r["Relaxed_dep"] == "False"), 285)
check("relaxed-only dep, no OE", sum(1 for r in mat if r["Relaxed_dep"] == "True" and r["Strict_dep"] == "False" and r["OE_pass"] == "False"), 4)
def cell(target, cancer):
    return next(r for r in mat if r["Target"] == target and r["Cancer"] == cancer)
check("EGFR found strict, EGFR-mutant lung", cell("EGFR", "EGFR-mutant lung")["Pass_strict_OE"], "True")
check("EGFR found strict, HNSC", cell("EGFR", "HNSC")["Pass_strict_OE"], "True")
check("ERBB2 found relaxed, colon", cell("ERBB2", "Colon")["Pass_relaxed_OE"], "True")
check("KRAS found strict, KRAS-mutant lung", cell("KRAS", "KRAS-mutant lung")["Pass_strict_OE"], "True")
check("EPAS1 found relaxed, RCC", cell("EPAS1", "RCC")["Pass_relaxed_OE"], "True")
check("BRAF dep-only, melanoma", (cell("BRAF", "Melanoma")["Strict_dep"], cell("BRAF", "Melanoma")["OE_pass"]), ("True", "False"))
check("FLT3 dep-only, AML", (cell("FLT3", "AML")["OE_pass"], cell("FLT3", "AML")["Strict_dep"], cell("FLT3", "AML")["Relaxed_dep"]), ("True", "False", "False"))
check("PARP1 passes nowhere", sum(1 for r in mat if r["Target"] == "PARP1" and (r["Pass_strict_OE"] == "True" or r["Pass_relaxed_OE"] == "True")), 0)

# ---------- Tabula Sapiens (corrected 26 Nov 2026: true per-cell-total CPM) ----------
print("== Tabula Sapiens switch-like test ==")
tg = read_csv(os.path.join(BASE, "tabula-sapiens", "per_gene_summary_truecpm.csv"))
q269 = [r for r in tg if r["in_269"] == "True"]
check("Tabula query genes", len(tg), 290)
check("selective genes in Tabula", len(q269), 269)
bg = read_csv(os.path.join(BASE, "tabula-sapiens", "background_summary_truecpm.csv"))
check("background genes", len(bg), 300)
s269 = sum(1 for r in q269 if r["switch_like"] == "True")
sbg = sum(1 for r in bg if r["switch_like"] == "True")
check("switch-like selective", s269, 25)
check("switch-like background", sbg, 63)
# Fisher exact, two-sided, via hypergeometric sum
def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    def hg(k):
        return math.comb(a + b, k) * math.comb(c + d, a + c - k) / math.comb(n, a + c)
    p_obs = hg(a)
    lo = max(0, a + c - (c + d)); hi = min(a + b, a + c)
    return sum(hg(k) for k in range(lo, hi + 1) if hg(k) <= p_obs * (1 + 1e-12))
p_f = fisher_two_sided(s269, len(q269) - s269, sbg, len(bg) - sbg)
check("Fisher p switch-like", round(p_f, 4), 0.0001, 0.0001)
tau269 = statistics.median(float(r["tau"]) for r in q269)
taubg = statistics.median(float(r["tau"]) for r in bg)
check("median tau selective", round(tau269, 3), 0.891, 0.001)
check("median tau background", round(taubg, 3), 0.904, 0.001)
check("tau>=0.9 selective", sum(1 for r in q269 if float(r["tau"]) >= 0.9), 121)
check("tau>=0.9 background", sum(1 for r in bg if float(r["tau"]) >= 0.9), 152)
check("median types detected selective", int(statistics.median(int(r["n_celltypes_frac_ge_0.10"]) for r in q269)), 128)
check("median types detected background", int(statistics.median(int(r["n_celltypes_frac_ge_0.10"]) for r in bg)), 89)
check("median types on selective", int(statistics.median(int(r["n_on_frac_ge_0.50"]) for r in q269)), 16)
check("median types on background", int(statistics.median(int(r["n_on_frac_ge_0.50"]) for r in bg)), 7)
check("median types off selective", int(statistics.median(int(r["n_off_frac_le_0.01"]) for r in q269)), 6)
check("median types off background", int(statistics.median(int(r["n_off_frac_le_0.01"]) for r in bg)), 10)
# Mann-Whitney U, normal approximation with tie correction
def mannwhitney_p(x, y):
    import itertools
    allv = sorted(x + y)
    ranks = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1] == allv[i]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks.setdefault(allv[k], []).append(r)
        i = j + 1
    # average ranks per value
    avg = {v: sum(rs) / len(rs) for v, rs in ranks.items()}
    rx = sum(avg[v] for v in x)
    n1, n2 = len(x), len(y)
    u = rx - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    # tie correction
    ties = Counter(allv)
    tcorr = sum(t * (t * t - 1) for t in ties.values() if t > 1)
    n = n1 + n2
    sd = math.sqrt(n1 * n2 * (n + 1 - tcorr / (n * (n - 1))) / 12)
    z = (u - mu) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
p_mw = mannwhitney_p([float(r["tau"]) for r in q269], [float(r["tau"]) for r in bg])
check("Mann-Whitney p tau", round(p_mw, 4), 0.0155, 0.0005)
check("median worst normal (CPM10k)", round(statistics.median(float(r["worst_mean_cpm10k"]) for r in q269), 1), 2.2, 0.5)
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
itgav = next(r for r in tg if r["gene"] == "ITGAV")
check("ITGAV worst cell type is RPE", "retinal pigment" in itgav["worst_celltype"], True)
check("ITGAV worst mean CPM10k", float(itgav["worst_mean_cpm10k"]), 27.2, 1.0)
tead3 = next(r for r in tg if r["gene"] == "TEAD3")
check("TEAD3 worst mean CPM10k", float(tead3["worst_mean_cpm10k"]), 0.7, 0.5)
tead1 = next(r for r in tg if r["gene"] == "TEAD1")
check("TEAD1 worst mean CPM10k", float(tead1["worst_mean_cpm10k"]), 5.8, 0.5)

# ---------- ADC intersection ----------
print("== ADC intersection ==")
adc = read_csv(os.path.join(BASE, "adc-intersection", "intersection_master.csv"))
check("ADC surface hits", len(adc), 71)
check("ADC distinct genes", len(set(r["gene"] for r in adc)), 37)
check("ADC bar met", sum(1 for r in adc if r["adc_bar"] == "met"), 3)
check_set("ADC bar-met genes", [r["gene"] for r in adc if r["adc_bar"] == "met"],
          ["SLC7A5", "SLC7A1", "SLC3A2"])
check("ADC new distinct genes", len(set(r["gene"] for r in adc if "NEW" in r.get("adc_status", ""))), 21)
check("ADC TAZ-cancer distinct genes", len(set(r["gene"] for r in adc if r["taz_cancer"] == "True")), 12)
rev = read_csv(os.path.join(BASE, "adc-intersection", "adc_genes_with_dependency_support.csv"))
check("ADC genes gaining dependency support", len(rev), 16)

# ---------- partner mining ----------
print("== partner mining tiers ==")
tazp = read_csv(os.path.join(BASE, "klf5-taz-partners", "taz_partners.csv"))
klf5p = read_csv(os.path.join(BASE, "klf5-taz-partners", "klf5_partners.csv"))
check("TAZ Tier 1 partners", sum(1 for r in tazp if r["tier"] == "1"), 3)
check_set("TAZ Tier 1 genes", [r["partner"] for r in tazp if r["tier"] == "1"],
          ["PTK2", "TEAD1", "CDK6"])
check("KLF5 Tier 1 partners", sum(1 for r in klf5p if r["tier"] == "1"), 9)
check_set("KLF5 Tier 1 genes", [r["partner"] for r in klf5p if r["tier"] == "1"],
          ["EGFR", "CTNNB1", "CDK6", "ERBB2", "CCNE1", "PIK3CA", "MAPK1", "WWTR1", "YAP1"])

# ---------- Table 1 intersection counts (strict+OE / relaxed+OE) ----------
# NB: synthesis.md's table lists the strict *dependency-bar* count (270/149/100)
# in the first column for AML, neuroblastoma and osteosarcoma; the per-cancer
# reports give strict+OE as 103 / 71 / 21. The manuscript Table 1 uses strict+OE.
print("== Table 1 intersection counts (per-cancer reports) ==")
aml_txt = open(os.path.join(BASE, "aml", "report.md")).read()
check("AML strict+OE", "103 pass the over-expression bar" in aml_txt, True)
check("AML relaxed+OE", "153 pass" in aml_txt, True)
nb_txt = open(os.path.join(BASE, "neuroblastoma", "report.md")).read()
check("NB strict+OE", "71 pass the over-expression bar" in nb_txt, True)
check("NB relaxed+OE", "126 pass" in nb_txt, True)
os_txt = open(os.path.join(BASE, "osteosarcoma", "report.md")).read()
check("OS relaxed+OE", "33 genes" in os_txt and "pass both a relaxed" in os_txt, True)
check("OS strict+OE", "21 pass the strict dependency bar" in os_txt, True)

# ---------- draft text spot checks ----------
print("== draft text contains the headline numbers ==")
ms = open(os.path.join(HERE, "manuscript.md")).read()
for s in ["385 top-20 list slots", "269 distinct genes", "200 singletons", "25/269 (9.3%)",
          "63/300 (21.0%)", "Fisher p = 0.0001", "0.891", "0.904", "SLC7A5, SLC7A1 and SLC3A2",
          "4 of the 7", "22 pass the strict bar", "285 pass neither", "1,136,218 cells",
          "27.2 CPM10k", "sole author", "not an author"]:
    check(f"draft contains {s!r}", s in ms, True)

print(f"\n{passed} passed, {failed} failed")
sys.exit(1 if failed else 0)
