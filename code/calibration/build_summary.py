#!/usr/bin/env python3
"""Build target_summary.csv for the drug-target calibration run from the matrix."""
import csv, os
CAL = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120/drug-target-calibration")
rows = list(csv.DictReader(open(os.path.join(CAL, "target_cancer_matrix.csv"))))
from collections import defaultdict
by = defaultdict(list)
for r in rows:
    by[r["Target"]].append(r)

def passes(tr, col):
    return "; ".join(sorted(r["Cancer"] for r in tr if r[col] == "True"))

SUMMARY = {
 "EGFR":  ("Approved", "EGFR-mutant NSCLC; HNSC; CRC RAS-WT; PDAC", "FOUND strict (EGFR-mutant lung, HNSC); FOUND relaxed (PDAC); dep-only in colon",
           "Expression-blind in one approved indication (colon)"),
 "ERBB2": ("Approved", "HER2+ gastric/GEJ (ESGA); HER2+ CRC (colon); HER2 IHC 3+ solid tumours (endometrial, bladder)",
           "FOUND relaxed (colon); NOT FOUND in ESGA (dep fails), endometrial, bladder",
           "Dilution — amplified subset diluted in whole-cancer group (ESGA)"),
 "KRAS":  ("Approved (G12C)", "KRAS-mutant NSCLC; KRAS G12C CRC (colon); KRAS G12C PDAC (no PDAC-specific approval)",
           "FOUND strict (KRAS-mutant lung); dep-only in colon and PDAC",
           "Expression-blind (PDAC, colon)"),
 "BRAF":  ("Approved", "BRAF V600E melanoma; BRAF V600E CRC (colon)",
           "Dep-only in both (melanoma strict dep, colon relaxed dep); OE fails both",
           "Expression-blind"),
 "ALK":   ("Approved", "ALK+ NSCLC — no ALK+ subgroup in programme",
           "NOT FOUND anywhere; nearest signal neuroblastoma (relaxed fails on dep frac 0.15)",
           "Fusion-rare (GW dep frac 0.006, ~7 lines) + unrepresented indication"),
 "ROS1":  ("Approved", "ROS1+ NSCLC — no ROS1+ subgroup in programme",
           "NOT FOUND anywhere; zero lines dependent genome-wide",
           "Fusion-rare, purest case (GW dep frac 0.000) + unrepresented indication"),
 "MET":   ("Approved", "MET exon-14 NSCLC — no subgroup in programme; HCC (cabozantinib, not MET-selected)",
           "NOT FOUND in any approved indication; passes relaxed in glioma (new signal)",
           "Unrepresented indication + amplification-rare (GW 0.031)"),
 "RET":   ("Approved", "RET fusion NSCLC / thyroid — no subgroup in programme",
           "NOT FOUND anywhere; one line dependent genome-wide (K-562, -0.515, not a RET-fusion tumour)",
           "Fusion-rare (GW dep frac 0.001) + unrepresented indication"),
 "NTRK1": ("Approved (tumour-agnostic)", "NTRK fusion solid tumours — no subgroup in programme",
           "NOT FOUND anywhere; one line dependent genome-wide (NCI-H1573, -0.531)",
           "Fusion-rare (GW dep frac 0.001) + unrepresented indication"),
 "FGFR1": ("Approved (FGFR1-rearranged MLN; FGFR2/3 indications are not FGFR1)",
           "None for FGFR1 itself in programme (bladder/cholangiocarcinoma approvals are FGFR2/3)",
           "Passes strict in HCC, glioma, osteosarcoma; relaxed in RCC — all new-indication signals",
           "Unrepresented indication; found elsewhere"),
 "PDGFRA":("Approved", "GIST — not in programme",
           "Passes strict in glioma, osteosarcoma — new-indication signals",
           "Unrepresented indication; found elsewhere"),
 "KIT":   ("Approved", "GIST; systemic mastocytosis — neither in programme (KIT-mutant melanoma is guideline-level only)",
           "NOT FOUND anywhere; GW dep frac 0.003 (~3 lines)",
           "Mutation-rare + unrepresented indication"),
 "FLT3":  ("Approved", "FLT3-mutant AML",
           "Dep-only in AML (rank 521; fails relaxed by one line, dep frac 0.18); OE passes strongly",
           "Mutation-blind — dependency tracks the ITD mutation, not expression"),
 "CDK4":  ("Approved", "HR+/HER2- breast — no HR+ breast group in programme (TNBC only)",
           "NOT FOUND anywhere; GW dep frac 0.667 exceeds the 0.60 cap in every cancer",
           "Pan-essential (near) + unrepresented indication"),
 "CDK6":  ("Approved", "HR+/HER2- breast — no HR+ breast group in programme (TNBC only)",
           "Passes strict in AML, HNSC, neuroblastoma, RCC; relaxed in TNBC — all new-indication signals",
           "Unrepresented indication; found elsewhere (PDAC verdict is TCGA-design-dependent)"),
 "PARP1": ("Approved", "BRCA-mutant/HRD ovarian; germline BRCA breast (incl. TNBC); germline BRCA PDAC — all mutation subsets, not whole-cancer groups",
           "NOT FOUND anywhere; GW dep frac 0.025",
           "Synthetic-lethal (+ knockout-vs-drug gap: even BRCA-LoF lines are only 8% PARP1-dependent)"),
 "TIPARP":("Clinical only (RBN-2397 Phase 1)", "None approved",
           "Passes strict in HNSC, PDAC; relaxed in RCC",
           "Found (no approved indication to test)"),
 "PIK3CA":("Approved", "HR+/HER2- PIK3CA-mutant breast — no HR+ breast group in programme",
           "Passes strict in ESCC, ESGA; relaxed in HNSC — new-indication signals; dep-only strict in bladder, Ewing, neuroblastoma",
           "Unrepresented indication; found elsewhere"),
 "MDM2":  ("Clinical only (no approval)", "None approved",
           "Passes strict in melanoma, RCC — driven by the TP53-WT fraction",
           "Conditional (TP53-WT) but found: 74% of TP53-WT lines vs 11% of TP53-LoF lines MDM2-dependent"),
 "EPAS1": ("Approved", "VHL-associated RCC; advanced RCC (belzutifan)",
           "FOUND relaxed in RCC (fails strict on dep frac by one line: 7/30, strict needs 8; q=0.11)",
           "Conditional (VHL-loss): dependency concentrated in VHL-LoF lines but only 29% of them"),
}

with open(os.path.join(CAL, "target_summary.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Target", "Drug_status", "Approved_indications_in_programme", "Verdict_in_approved",
                "Pass_strict_OE", "Pass_relaxed_OE", "Failure_taxonomy"])
    for t in ["EGFR","ERBB2","KRAS","BRAF","ALK","ROS1","MET","RET","NTRK1","FGFR1",
              "PDGFRA","KIT","FLT3","CDK4","CDK6","PARP1","TIPARP","PIK3CA","MDM2","EPAS1"]:
        status, ind, verdict, tax = SUMMARY[t]
        w.writerow([t, status, ind, verdict, passes(by[t], "Pass_strict_OE"),
                    passes(by[t], "Pass_relaxed_OE"), tax])
print("wrote target_summary.csv")
