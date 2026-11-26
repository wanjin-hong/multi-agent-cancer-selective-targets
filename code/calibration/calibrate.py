#!/usr/bin/env python3
"""Drug-target calibration: 20 approved/clinical targets through the programme's bars.

Reads the 20 finished per-cancer dependency_scan.csv + overexpression.csv files
(DepMap Public 24Q4, dependency = Chronos < -0.5) and, for each of the 20 targets,
extracts its dependency statistics in every cancer, recomputes the strict and
relaxed dependency bars from the numeric columns (fixed bars, pipeline-playbook.md
section 3), and takes the over-expression pass from each run's own pass column
(recomputed only for glioma and TNBC, whose overexpression.csv files carry no pass
columns, using the playbook section 4 bar on the same raw columns).

Outputs:
  target_cancer_matrix.csv  (20 targets x 20 cancers)
  target_summary.csv        (per-target approved-indication verdict + passes elsewhere)
"""
import csv, math, os, re, sys
from collections import defaultdict

BASE = os.path.expanduser("~/workspace/research_notes/multi-cancer-selective-targets-20261120")
TNBC = os.path.expanduser("~/workspace/research_notes/depmap-tnbc-selective-targets-20261120")
OUT = os.path.join(BASE, "drug-target-calibration")

TARGETS = ["EGFR","ERBB2","KRAS","BRAF","ALK","ROS1","MET","RET","NTRK1","FGFR1",
           "PDGFRA","KIT","FLT3","CDK4","CDK6","PARP1","TIPARP","PIK3CA","MDM2","EPAS1"]

# cancer -> (folder, dep cols, oe pass column or None->recompute, oe recompute spec)
# dep cols: gene, group_mean, group_dep_frac, group_n, enrichment, fracdiff, gw_dep_frac, panessential, q, rank
CANCERS = {
 "TNBC":            (TNBC, dict(gene="Gene", mean="TNBC_mean", depfrac="TNBC_dep_frac", n="TNBC_N_scored",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), None, "tnbc"),
 "EGFR-mutant lung":(os.path.join(BASE,"egfr-mutant-lung"), dict(gene="Gene", mean="GroupMean", depfrac="GroupDepFrac", n="GroupN",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWideDepFrac",
                                pe="PanEssential", q="BH_Q", rank=None), "Overexp_pass", None),
 "HCC":             (os.path.join(BASE,"hcc"), dict(gene="Gene", mean="HCC_mean", depfrac="HCC_dep_frac", n="HCC_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "Melanoma":        (os.path.join(BASE,"melanoma"), dict(gene="Gene", mean="Mel_mean", depfrac="Mel_dep_frac", n="Mel_N_scored",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "Pass_OE", None),
 "Colon":           (os.path.join(BASE,"colon"), dict(gene="Gene", mean="Colon_mean", depfrac="Colon_dep_frac", n="Colon_N_scored",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "Pass_overexpression_any", None),
 "Glioma":          (os.path.join(BASE,"glioma"), dict(gene="Gene", mean="Glioma_mean", depfrac="Glioma_dep_frac", n="Glioma_N_scored",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), None, "glioma"),
 "PDAC":            (os.path.join(BASE,"pancreatic"), dict(gene="Gene", mean="PDAC_mean", depfrac="PDAC_dep_frac", n="PDAC_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "KRAS-mutant lung":(os.path.join(BASE,"kras-mutant-lung"), dict(gene="Gene", mean="GroupMean", depfrac="GroupDepFrac", n="GroupN",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWideDepFrac",
                                pe="PanEssential", q="BH_Q", rank=None), "Overexp_pass", None),
 "HNSC":            (os.path.join(BASE,"head-and-neck"), dict(gene="Gene", mean="HNSC_mean", depfrac="HNSC_dep_frac", n="HNSC_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "Ovarian":         (os.path.join(BASE,"ovarian"), dict(gene="Gene", mean="OV_mean", depfrac="OV_dep_frac", n="OV_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "SCLC":            (os.path.join(BASE,"small-cell-lung"), dict(gene="Gene", mean="SCLC_mean", depfrac="SCLC_dep_frac", n="SCLC_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "ESGA":            (os.path.join(BASE,"esophagogastric"), dict(gene="Gene", mean="ESGA_mean", depfrac="ESGA_dep_frac", n="ESGA_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "ESCC":            (os.path.join(BASE,"esophageal-squamous"), dict(gene="Gene", mean="ESCC_mean", depfrac="ESCC_dep_frac", n="ESCC_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "Bladder":         (os.path.join(BASE,"bladder"), dict(gene="Gene", mean="Bladder_mean", depfrac="Bladder_dep_frac", n="Bladder_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "RCC":             (os.path.join(BASE,"renal-cell"), dict(gene="Gene", mean="RCC_mean", depfrac="RCC_dep_frac", n="RCC_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "Ewing":           (os.path.join(BASE,"ewing"), dict(gene="Gene", mean="Ewing_mean", depfrac="Ewing_dep_frac", n="Ewing_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q", rank=None), "OverExpr_pass", None),
 "Endometrial":     (os.path.join(BASE,"endometrial"), dict(gene="Gene", mean="ENDO_mean", depfrac="ENDO_dep_frac", n="ENDO_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "AML":             (os.path.join(BASE,"aml"), dict(gene="Gene", mean="AML_mean", depfrac="AML_dep_frac", n="AML_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "Neuroblastoma":   (os.path.join(BASE,"neuroblastoma"), dict(gene="Gene", mean="NB_mean", depfrac="NB_dep_frac", n="NB_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q_BH", rank="EnrichRank"), "OverExpr_pass", None),
 "Osteosarcoma":    (os.path.join(BASE,"osteosarcoma"), dict(gene="Gene", mean="OS_mean", depfrac="OS_dep_frac", n="OS_n",
                                enr="Enrichment", fd="FracDiff", gw="GenomeWide_dep_frac",
                                pe="PanEssential", q="q", rank="EnrichRank"), "OverExpr_pass", None),
}

def sym(s):
    """Strip ' (Entrez)' suffix -> plain symbol."""
    return re.sub(r"\s*\(\d+\)\s*$", "", (s or "").strip())

def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None

def fbool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("true", "1", "t", "yes"): return True
    if s in ("false", "0", "f", "no", "", "nan"): return False
    return None

def strict_bar(enr, depfrac, fd, pe, gw):
    return (enr is not None and depfrac is not None and fd is not None and gw is not None
            and enr <= -0.10 and depfrac >= 0.25 and fd >= 0.10 and pe is False and gw <= 0.60)

def relaxed_bar(enr, depfrac, fd, pe, gw):
    return (enr is not None and depfrac is not None and fd is not None and gw is not None
            and enr <= -0.05 and depfrac >= 0.20 and fd >= 0.10 and pe is False and gw <= 0.60)

def oe_recompute(row, kind):
    """Recompute the playbook over-expression bar from raw OE columns (glioma/TNBC only)."""
    if kind == "glioma":
        rna = (fnum(row.get("RNA_TPM_FC")) or 0) >= 1.5 and (fnum(row.get("RNA_delta_log2")) or 0) >= 0.5 \
              and (fnum(row.get("RNA_p_welch")) is not None and fnum(row.get("RNA_p_welch")) < 0.05)
        tcga = (fnum(row.get("TCGA_RNA_log2FC")) or 0) >= 0.5 \
               and (fnum(row.get("TCGA_RNA_p")) is not None and fnum(row.get("TCGA_RNA_p")) < 0.05)
        return bool(rna or tcga)   # protein leg dropped in the glioma run
    if kind == "tnbc":
        rna = (fnum(row.get("RNA_TPM_FC")) or 0) >= 1.5 and (fnum(row.get("RNA_delta_log2")) or 0) >= 0.5 \
              and (fnum(row.get("RNA_p_welch")) is not None and fnum(row.get("RNA_p_welch")) < 0.05)
        tcga = (fnum(row.get("TCGA_RNA_log2FC")) or 0) >= 0.5 \
               and (fnum(row.get("TCGA_RNA_p")) is not None and fnum(row.get("TCGA_RNA_p")) < 0.05)
        cptac = (fnum(row.get("CPTAC_delta")) or 0) >= 0.5 \
                and (fnum(row.get("CPTAC_p")) is not None and fnum(row.get("CPTAC_p")) < 0.05)
        rppa = (fnum(row.get("RPPA_delta")) or 0) >= 0.5 \
               and (fnum(row.get("RPPA_p")) is not None and fnum(row.get("RPPA_p")) < 0.05)
        return bool(rna or tcga or cptac or rppa)
    raise ValueError(kind)

def load_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def main():
    matrix = []          # one row per target x cancer
    crosscheck = []      # recomputed bars vs each file's own Strict/Relaxed columns
    missing = []
    for cancer, (folder, dc, oe_col, oe_kind) in CANCERS.items():
        dep = load_rows(os.path.join(folder, "dependency_scan.csv"))
        oe = load_rows(os.path.join(folder, "overexpression.csv"))
        # index OE rows by symbol (prefer Symbol column, else Gene stripped)
        oe_by_sym = {}
        for r in oe:
            s = sym(r.get("Symbol") or r.get("Gene"))
            oe_by_sym[s] = r
        # file's own bar columns, for cross-check (names vary)
        dep0 = dep[0]
        own_strict = own_relaxed = None
        for cand in ("StrictBar", "Strict_pass"):
            if cand in dep0: own_strict = cand
        for cand in ("RelaxedBar", "Relaxed_pass"):
            if cand in dep0: own_relaxed = cand
        oe0 = oe[0]
        own_oe_strict = own_oe_relaxed = None
        for cand in ("StrictBar", "Strict_pass"):
            if cand in oe0: own_oe_strict = cand
        for cand in ("RelaxedBar", "Relaxed_pass"):
            if cand in oe0: own_oe_relaxed = cand
        for t in TARGETS:
            drow = None
            for r in dep:
                if sym(r.get(dc["gene"])) == t:
                    drow = r; break
            if drow is None:
                missing.append((t, cancer, "dependency_scan"))
                continue
            orow = oe_by_sym.get(t)
            if orow is None:
                missing.append((t, cancer, "overexpression"))
                continue
            enr = fnum(drow.get(dc["enr"])); depfrac = fnum(drow.get(dc["depfrac"]))
            fd = fnum(drow.get(dc["fd"])); gw = fnum(drow.get(dc["gw"]))
            pe = fbool(drow.get(dc["pe"])); q = fnum(drow.get(dc["q"]))
            n = drow.get(dc["n"]); mean = fnum(drow.get(dc["mean"]))
            rank = drow.get(dc["rank"]) if dc["rank"] else None
            s_bar = strict_bar(enr, depfrac, fd, pe, gw)
            r_bar = relaxed_bar(enr, depfrac, fd, pe, gw)
            if oe_col:
                oe_pass = fbool(orow.get(oe_col))
            else:
                oe_pass = oe_recompute(orow, oe_kind)
            matrix.append(dict(Target=t, Cancer=cancer, Group_n=n, Group_mean_Chronos=mean,
                               Group_dep_frac=depfrac, Enrichment=enr, FracDiff=fd,
                               GW_dep_frac=gw, PanEssential=pe, q_BH=q, EnrichRank=rank,
                               Strict_dep=s_bar, Relaxed_dep=r_bar, OE_pass=oe_pass,
                               Pass_strict_OE=bool(s_bar and oe_pass),
                               Pass_relaxed_OE=bool(r_bar and oe_pass)))
            # cross-check recomputed bars against the file's own columns
            for src, own_col, mine in (("dep", own_strict, s_bar), ("dep", own_relaxed, r_bar),
                                       ("oe", own_oe_strict, s_bar), ("oe", own_oe_relaxed, r_bar)):
                if own_col:
                    src_row = drow if src == "dep" else orow
                    theirs = fbool(src_row.get(own_col))
                    if theirs is not None and theirs != mine:
                        crosscheck.append((t, cancer, src, own_col, theirs, mine))
    # write matrix
    fields = ["Target","Cancer","Group_n","Group_mean_Chronos","Group_dep_frac","Enrichment",
              "FracDiff","GW_dep_frac","PanEssential","q_BH","EnrichRank","Strict_dep",
              "Relaxed_dep","OE_pass","Pass_strict_OE","Pass_relaxed_OE"]
    with open(os.path.join(OUT, "target_cancer_matrix.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in matrix: w.writerow(r)
    print(f"matrix rows: {len(matrix)} (expected {len(TARGETS)*len(CANCERS)})")
    print(f"missing lookups: {missing}")
    print(f"bar cross-check disagreements vs files' own columns: {len(crosscheck)}")
    for c in crosscheck: print("  ", c)
    # per-target summary of passes
    by_target = defaultdict(list)
    for r in matrix: by_target[r["Target"]].append(r)
    print("\nPer-target passes (strict+OE / relaxed+OE):")
    for t in TARGETS:
        s = sorted([r["Cancer"] for r in by_target[t] if r["Pass_strict_OE"]])
        r_ = sorted([r["Cancer"] for r in by_target[t] if r["Pass_relaxed_OE"]])
        print(f"  {t}: strict {s} | relaxed {r_}")

if __name__ == "__main__":
    main()
