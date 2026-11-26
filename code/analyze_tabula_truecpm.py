"""Analyse the TRUE-CPM Tabula re-query (tabula_query_truecpm/ +
tabula_background_truecpm/) exactly as analyze_tabula.py analysed the
original query-set-relative run, and print the same headline comparisons.

Writes per_gene_summary_truecpm.csv and background_summary_truecpm.csv
(one row per gene: tau, switch-like, worst type, cancer legs). The original
per_gene_summary.csv / background_summary.csv are NOT overwritten.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_tabula as A

TAB = A.TAB


def main():
    gene_meta = pd.read_csv(TAB / "gene_list.csv")
    q, _, _, _ = A.summarise("tabula_query_truecpm", gene_meta)
    legs = A.load_cancer_legs()
    q = q.join(legs, how="left")
    q.to_csv(TAB / "per_gene_summary_truecpm.csv")
    print(f"\nQuery genes summarised (true CPM): {len(q)}")

    bg = None
    if (TAB / "tabula_background_truecpm" / "celltype_mean_cpm10k.csv").exists():
        bg, _, _, _ = A.summarise("tabula_background_truecpm", None)
        bg.to_csv(TAB / "background_summary_truecpm.csv")
        print(f"Background genes summarised (true CPM): {len(bg)}")
    else:
        print("\nNOTE: tabula_background_truecpm/ not present yet — background "
              "comparison skipped; re-run after the background query lands.")

    q269 = q[q["in_269"] == True]  # noqa: E712
    adc = q[q["in_adc37"] == True]  # noqa: E712
    print("\n--- tau specificity (TRUE CPM; 0=ubiquitous, 1=single-cell-type) ---")
    for name, s in [("269 selective", q269["tau"]), ("37 ADC surface", adc["tau"]),
                    ("290 query union", q["tau"])]:
        print(f"  {name}: median tau {s.median():.3f}, mean {s.mean():.3f}, "
              f"n={s.notna().sum()}, tau>=0.9: {(s >= 0.9).sum()}, tau<0.5: {(s < 0.5).sum()}")
    if bg is not None:
        s = bg["tau"]
        print(f"  300 background: median tau {s.median():.3f}, mean {s.mean():.3f}, "
              f"n={s.notna().sum()}, tau>=0.9: {(s >= 0.9).sum()}, tau<0.5: {(s < 0.5).sum()}")
        u, p = mannwhitneyu(q269["tau"].dropna(), bg["tau"].dropna(), alternative="two-sided")
        print(f"  Mann-Whitney 269 vs background tau: U={u:.0f}, p={p:.3g}")

    print("\n--- breadth: n cell types with frac >= 0.10 (frac-based, should match original) ---")
    print(f"  269 selective: median {q269['n_celltypes_frac_ge_0.10'].median():.0f}")
    if bg is not None:
        print(f"  300 background: median {bg['n_celltypes_frac_ge_0.10'].median():.0f}")
        u, p = mannwhitneyu(q269["n_celltypes_frac_ge_0.10"], bg["n_celltypes_frac_ge_0.10"])
        print(f"  Mann-Whitney breadth: U={u:.0f}, p={p:.3g}")

    print("\n--- switch-like (frac-based, should match original) ---")
    print(f"  269 selective: {q269['switch_like'].sum()}/{len(q269)} "
          f"({100 * q269['switch_like'].mean():.1f}%)")
    print(f"  37 ADC surface: {adc['switch_like'].sum()}/{len(adc)}")
    print(f"  269 selective: median n_on {q269['n_on_frac_ge_0.50'].median():.0f}, "
          f"median n_off {q269['n_off_frac_le_0.01'].median():.0f}")
    if bg is not None:
        print(f"  300 background: {bg['switch_like'].sum()}/{len(bg)} "
              f"({100 * bg['switch_like'].mean():.1f}%)")
        print(f"  300 background: median n_on {bg['n_on_frac_ge_0.50'].median():.0f}, "
              f"median n_off {bg['n_off_frac_le_0.01'].median():.0f}")
        a = int(q269["switch_like"].sum()); b = int(len(q269) - a)
        c = int(bg["switch_like"].sum()); d = int(len(bg) - c)
        _, p = fisher_exact([[a, b], [c, d]])
        print(f"  Fisher exact 269 vs background switch-like: p={p:.3g}")

    print("\n--- worst-type distribution / windows (TRUE CPM) ---")
    print(f"  median worst 269: {q269['worst_mean_cpm10k'].median():.1f}")
    print(f"  worst<100: {(q269['worst_mean_cpm10k'] < 100).sum()}, "
          f"worst>=500: {(q269['worst_mean_cpm10k'] >= 500).sum()}, "
          f"worst>=1000: {(q269['worst_mean_cpm10k'] >= 1000).sum()}")
    wide = q269[(q269["max_tpm_fc"] >= 3) & (q269["worst_mean_cpm10k"] < 100)]
    narrow = q269[(q269["max_tpm_fc"] >= 3) & (q269["worst_mean_cpm10k"] >= 500)]
    print(f"  wide window (TPM FC>=3, worst<100): {len(wide)}")
    print(f"  narrow window (TPM FC>=3, worst>=500): {len(narrow)}")
    print(f"  RPE-worst count 290: {(q['worst_celltype'] == 'retinal pigment epithelial cell').sum()}")
    print(f"  kidney-epith-worst count 290: {(q['worst_celltype'] == 'kidney epithelial cell').sum()}")
    print("  wide-window genes:", ", ".join(wide.index.tolist()))
    print("  narrow-window genes:", ", ".join(narrow.index.tolist()))

    print("\n--- TEAD axis in Tabula v2 (TRUE CPM mean / frac) ---")
    mean_q = pd.read_csv(TAB / "tabula_query_truecpm" / "celltype_mean_cpm10k.csv", index_col=0)
    frac_q = pd.read_csv(TAB / "tabula_query_truecpm" / "celltype_frac_expressing.csv", index_col=0)
    for g in ["TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1", "WWTR1"]:
        r = q.loc[g]
        ke = (f"kidney-epith {mean_q[g].get('kidney epithelial cell', float('nan')):.1f} "
              f"(frac {frac_q[g].get('kidney epithelial cell', float('nan')):.2f})")
        rpe = (f"RPE {mean_q[g].get('retinal pigment epithelial cell', float('nan')):.1f} "
               f"(frac {frac_q[g].get('retinal pigment epithelial cell', float('nan')):.2f})")
        print(f"  {g}: tau={r['tau']:.2f}, worst={r['worst_celltype']} "
              f"({r['worst_mean_cpm10k']:.1f}, frac {r['worst_frac']:.2f}), {ke}, {rpe}, "
              f"retinal {r['retinal_mean_cpm10k']:.1f} (frac {r['retinal_frac']:.2f})")

    print("\n--- ADC 37 RPE / worst values (TRUE CPM) ---")
    for g in ["ITGAV", "SLC3A2", "SLC7A5", "SLC7A1", "SLC2A1", "LRP8", "ERBB2", "CDH2", "ERBB3"]:
        r = q.loc[g]
        print(f"  {g}: worst={r['worst_celltype']} ({r['worst_mean_cpm10k']:.1f}, "
              f"frac {r['worst_frac']:.2f}), switch_like={bool(r['switch_like'])}")
    print(f"  ADC median tau: {adc['tau'].median():.3f}; ADC switch-like: "
          f"{adc['switch_like'].sum()}/37")
    print(f"  RPE-worst among ADC 37: "
          f"{(adc['worst_celltype'] == 'retinal pigment epithelial cell').sum()}")

    print("\n--- spot checks vs original (frac must match; means must NOT) ---")
    old = pd.read_csv(TAB / "per_gene_summary.csv", index_col=0)
    for g in ["ITGAV", "TEAD1", "TEAD3", "WWTR1", "CDK6", "EGFR", "PAX8"]:
        print(f"  {g}: old worst {old.loc[g, 'worst_celltype']} "
              f"({old.loc[g, 'worst_mean_cpm10k']:.1f}) -> new worst "
              f"{q.loc[g, 'worst_celltype']} ({q.loc[g, 'worst_mean_cpm10k']:.1f}); "
              f"worst_frac {old.loc[g, 'worst_frac']:.2f} -> {q.loc[g, 'worst_frac']:.2f}")
    frac_old = pd.read_csv(TAB / "tabula_query" / "celltype_frac_expressing.csv", index_col=0)
    frac_new = pd.read_csv(TAB / "tabula_query_truecpm" / "celltype_frac_expressing.csv", index_col=0)
    print(f"  frac matrices identical: "
          f"{np.allclose(frac_old.values, frac_new.values, equal_nan=True)} "
          f"(max abs diff {np.abs(frac_old.values - frac_new.values).max():.2e})")


if __name__ == "__main__":
    main()
