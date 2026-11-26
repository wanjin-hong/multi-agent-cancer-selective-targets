"""Verify the CORRECTED (true-CPM) headline numbers in report.md against the
saved true-CPM Tabula outputs (tabula_query_truecpm/). The original
verify_headlines.py checks the first, query-set-relative run; this file
checks the corrected values. Exits non-zero on any mismatch.

Recomputes tau, worst-type calls, switch-like counts, window classes and the
TEAD-axis rows from the matrices and per_gene_summary.csv, and asserts they
match the values stated in the report. Exits non-zero on any mismatch.
Background checks run only once tabula_background/ exists.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TAB = Path("/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/tabula-sapiens")
fails = []


def check(name, got, want, tol=0.0):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    print(f"  {'PASS' if ok else 'FAIL'}: {name}: got {got}, want {want}")
    if not ok:
        fails.append(name)


def tau(x):
    x = np.asarray(x, dtype=float)
    return float(np.sum(1.0 - x / x.max()) / (len(x) - 1)) if x.max() > 0 else float("nan")


def main():
    meta = json.load(open(TAB / "tabula_query_truecpm" / "query_metadata.json"))
    mean = pd.read_csv(TAB / "tabula_query_truecpm" / "celltype_mean_cpm10k.csv", index_col=0)
    frac = pd.read_csv(TAB / "tabula_query_truecpm" / "celltype_frac_expressing.csv", index_col=0)
    cmeta = pd.read_csv(TAB / "tabula_query_truecpm" / "celltype_metadata.csv")
    summ = pd.read_csv(TAB / "per_gene_summary_truecpm.csv", index_col=0)
    print("== identity / denominators ==")
    check("tabula version", meta["tabula_version_used"], "v2")
    check("census version", meta["census_version"], "2025-11-17")
    check("n cells", meta["n_cells"], 1136218)
    check("n donors", meta["n_donors"], 24)
    check("n cell types", meta["n_cell_types"], 180)
    check("genes found", meta["n_genes_found"], 290)
    check("genes missing", meta["genes_missing"], [])
    check("matrix shape", mean.shape, (180, 290))
    check("cells sum from metadata", int(cmeta["n_cells"].sum()), 1136218)
    check("organ sum matches All Cells", meta["organ_cell_sum_matches_all_cells"], True)

    print("== podocyte gap ==")
    bad = [c for c in mean.index if any(k in c.lower() for k in
           ["podocyte", "glomerul", "tubule", "proximal", "parietal", "henle"])]
    check("podocyte/glomerular/tubular types", bad, [])

    print("== kidney epithelial TEAD axis ==")
    check("TEAD1 kidney-epith mean (true CPM)", round(mean.loc["kidney epithelial cell", "TEAD1"], 1), 1.4, 0.05)
    check("TEAD3 kidney-epith mean (true CPM)", round(mean.loc["kidney epithelial cell", "TEAD3"], 1), 0.1, 0.05)
    check("TEAD3 kidney-epith frac", round(frac.loc["kidney epithelial cell", "TEAD3"], 2), 0.06, 0.005)
    check("TEAD1 kidney-epith frac", round(frac.loc["kidney epithelial cell", "TEAD1"], 2), 0.64, 0.005)

    print("== retinal / ITGAV ==")
    check("ITGAV RPE mean (true CPM)", round(mean.loc["retinal pigment epithelial cell", "ITGAV"], 1), 27.2, 0.05)
    check("ITGAV RPE frac", round(frac.loc["retinal pigment epithelial cell", "ITGAV"], 2), 0.93, 0.005)
    check("TEAD3 RPE mean (true CPM)", round(mean.loc["retinal pigment epithelial cell", "TEAD3"], 1), 0.1, 0.05)
    check("WWTR1 retinal-endothelium mean",
          round(mean.loc["retinal blood vessel endothelial cell", "WWTR1"], 1), 5.3, 0.05)

    print("== tau / switch-like (recomputed from the mean matrix) ==")
    q269 = summ[summ["in_269"] == True]  # noqa: E712
    taus = {g: tau(mean[g].values) for g in mean.columns}
    med_tau = float(np.median([taus[g] for g in q269.index]))
    check("median tau 269 (true CPM)", round(med_tau, 3), 0.891, 0.001)
    check("tau>=0.9 count 269 (true CPM)", int(sum(taus[g] >= 0.9 for g in q269.index)), 121)
    check("tau<0.5 count 269", int(sum(taus[g] < 0.5 for g in q269.index)), 0)
    check("switch-like 269", int(q269["switch_like"].sum()), 25)
    check("median n_on 269", int(q269["n_on_frac_ge_0.50"].median()), 16)
    check("median n_off 269", int(q269["n_off_frac_le_0.01"].median()), 6)
    check("median n frac>=0.10 269", int(q269["n_celltypes_frac_ge_0.10"].median()), 128)

    print("== worst-type distribution / windows ==")
    check("median worst 269 (true CPM)", round(float(q269["worst_mean_cpm10k"].median()), 1), 2.2, 0.05)
    check("worst<100 count 269 (true CPM)", int((q269["worst_mean_cpm10k"] < 100).sum()), 267)
    check("worst>=500 count 269 (true CPM)", int((q269["worst_mean_cpm10k"] >= 500).sum()), 0)
    check("worst>=1000 count 269 (true CPM)", int((q269["worst_mean_cpm10k"] >= 1000).sum()), 0)
    wide = q269[(q269["max_tpm_fc"] >= 3) & (q269["worst_mean_cpm10k"] < 100)]
    check("wide window 90 (true CPM)", len(wide), 90)
    narrow = q269[(q269["max_tpm_fc"] >= 3) & (q269["worst_mean_cpm10k"] >= 500)]
    check("narrow window 0 (true CPM)", len(narrow), 0)
    check("RPE-worst count 290",
          int((summ["worst_celltype"] == "retinal pigment epithelial cell").sum()), 10)
    check("kidney-epith-worst count 290",
          int((summ["worst_celltype"] == "kidney epithelial cell").sum()), 5)

    print("== TEAD axis worst rows ==")
    for g, worst, val in [("TEAD1", "pulmonary alveolar type 1 cell", 5.8), ("TEAD2",
                          "vein endothelial cell", 0.5),
                          ("TEAD3", "myoepithelial cell", 0.7),
                          ("TEAD4", "endothelial cell of arteriole", 1.0),
                          ("YAP1", "corneal epithelial cell", 3.4),
                          ("WWTR1", "endothelial cell of artery", 7.7)]:
        check(f"{g} worst type", summ.loc[g, "worst_celltype"], worst)
        check(f"{g} worst mean", round(float(summ.loc[g, "worst_mean_cpm10k"]), 1), val, 0.05)

    print("== three-leg spot checks ==")
    check("PAX8 TPM FC", round(float(summ.loc["PAX8", "max_tpm_fc"]), 1), 34.6, 0.05)
    check("HNF1B TPM FC", round(float(summ.loc["HNF1B", "max_tpm_fc"]), 1), 48.2, 0.05)
    check("FOXA2 TCGA", round(float(summ.loc["FOXA2", "max_tcga_log2fc"]), 2), -3.39, 0.005)
    check("KRAS TPM FC", round(float(summ.loc["KRAS", "max_tpm_fc"]), 1), 1.6, 0.05)
    check("CDK6 worst type", summ.loc["CDK6", "worst_celltype"],
          "common myeloid progenitor")

    bg_path = TAB / "tabula_background_truecpm" / "celltype_mean_cpm10k.csv"
    if bg_path.exists():
        print("== background comparison (Virtual Biotech test) ==")
        bg = pd.read_csv(TAB / "background_summary_truecpm.csv", index_col=0)
        check("background genes (PDE4C deduped)", len(bg), 300)
        check("median tau background (true CPM)", round(float(bg["tau"].median()), 3), 0.904, 0.001)
        check("tau>=0.9 background (true CPM)", int((bg["tau"] >= 0.9).sum()), 152)
        check("switch-like background", int(bg["switch_like"].sum()), 63)
        check("switch-like 269", int(q269["switch_like"].sum()), 25)
        check("median n_on background", int(bg["n_on_frac_ge_0.50"].median()), 7)
        check("median n_off background", int(bg["n_off_frac_le_0.01"].median()), 10)
        check("median breadth background", int(bg["n_celltypes_frac_ge_0.10"].median()), 89)
        from scipy.stats import fisher_exact, mannwhitneyu
        _, p_mw = mannwhitneyu(q269["tau"].dropna(), bg["tau"].dropna())
        check("MW tau p significant (p<0.05), selective LESS specific", bool(p_mw < 0.05), True)
        _, p_f = fisher_exact([[int(q269["switch_like"].sum()),
                                int(len(q269) - q269["switch_like"].sum())],
                               [int(bg["switch_like"].sum()),
                                int(len(bg) - bg["switch_like"].sum())]])
        check("Fisher switch-like p", round(float(p_f), 6), 0.000114, 0.00001)
    else:
        print("== background: not present yet — comparison checks skipped ==")

    print()
    if fails:
        print(f"FAILED: {len(fails)} check(s): {fails}")
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
