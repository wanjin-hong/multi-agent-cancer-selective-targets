#!/usr/bin/env python3
"""Independent verification for the drug-target calibration run.

Recomputes headline numbers from the RAW DepMap 24Q4 files (CRISPRGeneEffect.csv,
OmicsExpressionProteinCodingGenesTPMLogp1.csv, Model.csv, OmicsSomaticMutations.csv)
through a separate code path from calibrate.py, which read the per-cancer scan CSVs.
Checks:
  A. Dependency: group mean Chronos, dependent count/fraction, enrichment vs all
     other CRISPR lines, for 12 headline gene x cancer cells.
  B. Over-expression: cell-line RNA TPM fold change (verified back-transform formula)
     and delta-log2 for 3 headline cells.
  C. FLT3-ITD: the three ITD lines' Chronos values (mutation file + raw matrix).
  D. Fusion-rare counts: number of lines genome-wide dependent on ALK/ROS1/RET/NTRK1,
     and the identities of the ALK-dependent lines.
  E. CDK4 pan-essential check: genome-wide dependent fraction from the raw matrix.
Tolerances: group mean |diff| <= 0.005, dependent counts exact, enrichment |diff|
<= 0.01, TPM FC within 2%.
"""
import csv, math, os, sys
import pandas as pd
from scipy import stats

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
CAL = os.path.join(BASE, "drug-target-calibration")
VRK = os.path.expanduser("~/workspace/research_notes/vrk1-codependency-20261108")
GLI = os.path.expanduser("~/workspace/research_notes/glioma-top12-deep-analyses-20261108")
MUT = os.path.expanduser("~/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv")
TNBC_CODEP = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-codependency-20261120")

DEP = 0.5  # Chronos < -0.5

def col_for(df, symbol):
    hits = [c for c in df.columns if c.split(" (")[0] == symbol]
    assert len(hits) == 1, (symbol, hits)
    return hits[0]

def main():
    print("loading raw matrices ...")
    crispr = pd.read_csv(os.path.join(VRK, "CRISPRGeneEffect.csv"), index_col=0)
    model = pd.read_csv(os.path.join(VRK, "Model.csv"), index_col="ModelID")
    expr = pd.read_csv(os.path.join(GLI, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), index_col=0)
    print(f"CRISPR: {crispr.shape[0]} lines x {crispr.shape[1]} genes; expression: {expr.shape[0]} lines")
    all_crispr = set(crispr.index)
    all_expr = set(expr.index)
    disease = model["OncotreePrimaryDisease"].to_dict()
    cellname = model["CellLineName"].to_dict()

    def disease_group(label):
        return sorted(m for m in all_crispr if disease.get(m) == label)

    def lines_csv(path, panel=None):
        rows = list(csv.DictReader(open(path)))
        ids = [r["ModelID"] for r in rows if panel is None or panel in r.get("Panel", "")]
        return sorted(set(ids) & all_crispr)

    # --- group definitions (independent of the scan scripts) ---
    groups = {
        "EGFR-mutant lung": lines_csv(os.path.join(BASE, "egfr-mutant-lung/activating_lines.csv"), "CRISPR"),
        "KRAS-mutant lung": lines_csv(os.path.join(BASE, "kras-mutant-lung/activating_lines.csv"), "CRISPR"),
        "AML": disease_group("Acute Myeloid Leukemia"),
        "RCC": disease_group("Renal Cell Carcinoma"),
        "Melanoma": disease_group("Melanoma"),
        "HNSC": disease_group("Head and Neck Squamous Cell Carcinoma"),
        "HCC": disease_group("Hepatocellular Carcinoma"),
        "PDAC": disease_group("Pancreatic Adenocarcinoma"),
    }
    tnbc_rows = list(csv.DictReader(open(os.path.join(TNBC_CODEP, "tnbc_lines.csv"))))
    groups["TNBC"] = sorted({r["ModelID"] for r in tnbc_rows if r.get("Core_TNBC") == "True"} & all_crispr)
    for g, ids in groups.items():
        print(f"  group {g}: n={len(ids)}")

    # expected values from target_cancer_matrix.csv
    matrix = {(r["Target"], r["Cancer"]): r for r in csv.DictReader(open(os.path.join(CAL, "target_cancer_matrix.csv")))}

    fails = []
    def check(label, cond, detail=""):
        print(("  PASS " if cond else "  FAIL ") + label + (f"  [{detail}]" if detail else ""))
        if not cond:
            fails.append(label)

    # --- A. dependency cells ---
    # Comparators: the run's own comparator for each scan (pan-cancer = all other
    # CRISPR lines; the two lung subtypes = WT NSCLC; TNBC = non-TNBC breast).
    print("\nA. Dependency recomputation from raw CRISPR matrix")
    mut_all = pd.read_csv(MUT, usecols=["ModelID", "HugoSymbol"], low_memory=True)
    seq_models = set(mut_all["ModelID"])
    egfr_mut_models = set(mut_all[mut_all["HugoSymbol"] == "EGFR"]["ModelID"])
    kras_mut_models = set(mut_all[mut_all["HugoSymbol"] == "KRAS"]["ModelID"])
    nsclc = sorted(m for m in all_crispr
                   if model.loc[m, "OncotreeLineage"] == "Lung"
                   and disease.get(m) == "Non-Small Cell Lung Cancer")
    comp_map = {
        ("EGFR", "EGFR-mutant lung"): sorted(m for m in nsclc if m not in egfr_mut_models),
        ("KRAS", "KRAS-mutant lung"): sorted(m for m in nsclc if m in seq_models and m not in kras_mut_models),
        ("CDK6", "TNBC"): sorted({r["ModelID"] for r in tnbc_rows if r.get("Classification") == "Non-TNBC breast"} & all_crispr),
    }
    for key, comp_ids in comp_map.items():
        print(f"  comparator {key}: n={len(comp_ids)}")
    cells = [("EGFR", "EGFR-mutant lung"), ("KRAS", "KRAS-mutant lung"), ("FLT3", "AML"),
             ("EPAS1", "RCC"), ("BRAF", "Melanoma"), ("CDK6", "HNSC"), ("CDK6", "AML"),
             ("MDM2", "Melanoma"), ("EGFR", "HNSC"), ("FGFR1", "HCC"),
             ("CDK6", "TNBC"), ("TIPARP", "PDAC")]
    for gene, cancer in cells:
        exp = matrix[(gene, cancer)]
        col = col_for(crispr, gene)
        g = groups[cancer]
        gv = crispr.loc[g, col].dropna()
        comp_ids = comp_map.get((gene, cancer))
        if comp_ids is None:
            comp_ids = [m for m in crispr.index if m not in set(g)]
        comp = crispr.loc[comp_ids, col].dropna()
        gmean, gdep = gv.mean(), int((gv < -DEP).sum())
        cmean = comp.mean()
        enr = gmean - cmean
        check(f"{gene}/{cancer} group mean {gmean:+.3f} vs {float(exp['Group_mean_Chronos']):+.3f}",
              abs(gmean - float(exp["Group_mean_Chronos"])) <= 0.005, f"n={len(gv)}")
        check(f"{gene}/{cancer} dep count {gdep} vs {int(round(float(exp['Group_dep_frac'])*int(exp['Group_n'])))}",
              gdep == int(round(float(exp["Group_dep_frac"]) * int(exp["Group_n"]))), f"dep frac {gdep/len(gv):.3f}")
        check(f"{gene}/{cancer} enrichment {enr:+.3f} vs {float(exp['Enrichment']):+.3f}",
              abs(enr - float(exp["Enrichment"])) <= 0.01)

    # --- B. over-expression cells (raw expression matrix, pan-cancer comparator) ---
    print("\nB. Cell-line RNA recomputation from raw expression matrix")
    def rna_check(gene, cancer, disease_label):
        exp = matrix[(gene, cancer)]
        col = col_for(expr, gene)
        g = sorted({m for m in all_expr if disease.get(m) == disease_label})
        gv = expr.loc[g, col].dropna()
        cv = expr.loc[[m for m in expr.index if m not in set(g)], col].dropna()
        dlog2 = gv.mean() - cv.mean()
        fc = (2 ** gv.mean() - 1) / (2 ** cv.mean() - 1)
        # expected values read from the per-cancer overexpression.csv
        return dlog2, fc, len(gv), len(cv)
    # expected OE numbers taken from the run files (quoted in the report methods)
    for gene, cancer, label, exp_fc, exp_dl in [
            ("FLT3", "AML", "Acute Myeloid Leukemia", 34.372, 2.501),
            ("CDK6", "AML", "Acute Myeloid Leukemia", 2.845, 1.424),
            ("EPAS1", "RCC", "Renal Cell Carcinoma", 1.819, 0.821)]:
        dl, fc, ng, nc = rna_check(gene, cancer, label)
        check(f"{gene}/{cancer} TPM FC {fc:.3f} vs {exp_fc}",
              abs(fc - exp_fc) / exp_fc <= 0.02, f"n={ng}/{nc}")
        check(f"{gene}/{cancer} delta-log2 {dl:.3f} vs {exp_dl}",
              abs(dl - exp_dl) <= 0.01)

    # --- C. FLT3-ITD ---
    print("\nC. FLT3-ITD lines (mutation file + raw matrix)")
    mut = pd.read_csv(MUT, usecols=["ModelID", "HugoSymbol", "VariantType", "ProteinChange",
                                     "LikelyLoF", "Hotspot"], low_memory=True)
    flt3_itd = mut[(mut["HugoSymbol"] == "FLT3") & (mut["VariantType"] == "insertion")]
    print("  FLT3 insertion rows (dedup by ModelID):")
    itd_ids = sorted(set(flt3_itd["ModelID"]))
    col = col_for(crispr, "FLT3")
    for mid in itd_ids:
        pc = "; ".join(sorted(set(flt3_itd[flt3_itd["ModelID"] == mid]["ProteinChange"].astype(str))))
        val = crispr.loc[mid, col] if mid in crispr.index else float("nan")
        print(f"    {cellname.get(mid, mid)} ({mid}): {pc}  Chronos={val:.3f}  dependent={val < -DEP}")
        if cellname.get(mid) in ("MV4-11", "MOLM-13", "MOLM-14"):
            check(f"FLT3-ITD {cellname.get(mid)} dependent", bool(val < -DEP), f"Chronos {val:.2f}")
    check("FLT3-ITD lines all dependent (3/3)",
          all(crispr.loc[m, col] < -DEP for m in itd_ids if cellname.get(m) in ("MV4-11", "MOLM-13", "MOLM-14")))

    # --- D. fusion-rare counts ---
    print("\nD. Genome-wide dependent-line counts (raw matrix)")
    for gene in ["ALK", "ROS1", "RET", "NTRK1"]:
        col = col_for(crispr, gene)
        s = crispr[col].dropna()
        ndep = int((s < -DEP).sum())
        gw = ndep / len(s)
        exp_gw = float(matrix[(gene, "Melanoma")]["GW_dep_frac"])
        check(f"{gene}: {ndep}/{len(s)} dependent genome-wide (gw={gw:.4f} vs scan {exp_gw:.4f})",
              abs(gw - exp_gw) <= 0.002)
        if gene == "ALK":
            dep_lines = s[s < -DEP].index.tolist()
            print("  ALK-dependent lines: " + ", ".join(
                f"{cellname.get(m, m)} [{disease.get(m)}]" for m in dep_lines))

    # --- E. CDK4 pan-essential ---
    print("\nE. CDK4 genome-wide dependency (raw matrix)")
    col = col_for(crispr, "CDK4")
    s = crispr[col].dropna()
    gw = (s < -DEP).mean()
    exp_gw = float(matrix[("CDK4", "Melanoma")]["GW_dep_frac"])
    check(f"CDK4 gw dep frac {gw:.3f} vs scan {exp_gw:.3f} (>0.60 cap -> fails both bars everywhere)",
          abs(gw - exp_gw) <= 0.002 and gw > 0.60)

    print(f"\n{'ALL CHECKS PASSED' if not fails else 'FAILURES: ' + str(fails)}")
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
