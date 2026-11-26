"""ADC x DepMap intersection driver.

Direction 1 (dependency-first): per-cancer selective candidates (published
candidate files) intersected with the SURFY-2886 surfaceome, ADC calibrated bar
and normal-tissue safety check applied.
Direction 2 (ADC-first): each of the 375 ADC-nominated genes looked up in all
20 cancers' dependency_scan.csv files.
Direction 3 (cross-cancer): recurrence of surface dependencies; TAZ-cancer flag.

Definitions (stated atop every table in the report):
- DepMap Public 24Q4; dependency = Chronos < -0.5.
- Candidate = published per-cancer candidate list (strict + relaxed dependency
  bar + over-expression bar, as in synthesis.md). TNBC uses top20.csv only
  (20 genes; the 5 further relaxed+OE candidates have no machine-readable
  candidate file, so TNBC counts are a floor — flagged in the report).
- AML candidates.csv is the pan-cancer-comparator list (153); the 89-gene
  blood-only list (candidates_haem.csv) is used as a sensitivity check.
- Surface = SURFY in-silico human surfaceome, 2,886 rows / 2,799 symbols.
- ADC calibrated bar (19 Nov 2026 programme): surface + tumour >= 10 FPKM +
  protein support. Assessable only where enriched_per_cancer.csv has data
  (BRCA, COAD, LIHC, LUAD, LUSC, READ, STAD).
- Safety flags (explicit): CRITICAL = candidate_table critical_tissue_flag;
  CELLTYPE = max normal cell-type nCPM >= 2000 (the ADC Tier 1 cutoff);
  TISSUE = max normal tissue consensus nTPM >= 200.
"""
import json, os, re
import pandas as pd
import numpy as np

ROOT = os.path.expanduser("~/workspace/research_notes")
MULTI = os.path.join(ROOT, "multi-cancer-selective-targets-20261120")
ADC = os.path.join(ROOT, "cancer-enriched-surface-proteins-adc-20261119", "data")
OUT = os.path.join(MULTI, "adc-intersection")
os.makedirs(os.path.join(OUT, "per_cancer"), exist_ok=True)

def sym(x):
    """Normalise a gene symbol: strip ' (Entrez)' suffixes, whitespace."""
    if pd.isna(x):
        return None
    s = str(x).strip()
    if " (" in s:
        s = s.split(" (")[0].strip()
    return s

def pick_col(df, patterns, default=None):
    for pat in patterns:
        for c in df.columns:
            if re.fullmatch(pat, c, flags=re.I):
                return c
    return default

# ---------------------------------------------------------------- ADC inputs
surf = pd.read_csv(os.path.join(ADC, "surfaceome.csv"))
SURFY = set(sym(s) for s in surf["gene_symbol"] if sym(s))
print("SURFY symbols:", len(SURFY))

ct = pd.read_csv(os.path.join(ADC, "candidate_table.csv"))
ct["symbol"] = ct["gene_symbol"].map(sym)
CT = ct.set_index("symbol")
print("candidate_table rows:", len(ct), "| tiers:", ct["tier"].value_counts().to_dict())

epc = pd.read_csv(os.path.join(ADC, "enriched_per_cancer.csv"))
epc["symbol"] = epc["symbol"].map(sym)
EPC = {(r.symbol, r.cancer): r for r in epc.itertuples()}
print("enriched_per_cancer rows:", len(epc), "| cancers:", sorted(epc.cancer.unique()))

ncm = pd.read_csv(os.path.join(ADC, "normal_celltype_max.csv"))
ncm["symbol"] = ncm["symbol"].map(sym)
NCM = ncm.set_index("symbol")

# Full-coverage tumour RNA (all genes x 7 cancers) — the right source for the
# ADC tumour-FPKM bar; enriched_per_cancer.csv is a filtered subset and must
# NOT be used for the bar (absence there means "not enriched", not "low").
brm = pd.read_csv(os.path.join(ADC, "bulk_rna_vs_normal.csv"))
brm["symbol"] = brm["symbol"].map(sym)
BRM = {(r.symbol, r.cancer): r for r in brm.itertuples()}
print("bulk_rna_vs_normal rows:", len(brm))

# CPTAC protein log2FC (LUAD, LUSC, BRCA, COAD only); best isoform kept.
cpt = pd.read_csv(os.path.join(ADC, "cptac_protein_fc.csv"))
cpt["symbol"] = cpt["gene_symbol"].map(sym)
CPT = {}
for (g, c), sub in cpt.groupby(["symbol", "cancer"]):
    CPT[(g, c)] = float(sub["log2FC"].max())
print("CPTAC gene-cancer pairs:", len(CPT))

ihc = pd.read_csv(os.path.join(ADC, "hpa_ihc_summary.csv"))
ihc["symbol"] = ihc["gene_symbol"].map(sym)
IHC = ihc.set_index("symbol")

# ---------------------------------------------------------------- cancers
# (key, display, folder, candidate-file spec, group n, TCGA code for ADC
#  mapping or None, TAZ cancer?)
CANCERS = [
    ("tnbc", "TNBC", os.path.join(ROOT, "depmap-tnbc-selective-targets-20261120"),
     "top20.csv", 26, "BRCA", True),
    ("egfr-mutant-lung", "EGFR-mutant lung", os.path.join(MULTI, "egfr-mutant-lung"),
     "ranked_candidates.csv", 6, "LUAD", False),
    ("hcc", "HCC", os.path.join(MULTI, "hcc"),
     "candidates.csv", 21, "LIHC", False),
    ("melanoma", "Melanoma", os.path.join(MULTI, "melanoma"),
     "_relaxed_scored.csv", 66, None, True),
    ("colon", "Colon", os.path.join(MULTI, "colon"),
     "candidates_both.csv", 63, "COAD", False),
    ("glioma", "Glioma", os.path.join(MULTI, "glioma"),
     "candidates.csv", 72, None, True),
    ("pancreatic", "PDAC", os.path.join(MULTI, "pancreatic"),
     "candidates.csv", 45, None, False),
    ("kras-mutant-lung", "KRAS-mutant lung", os.path.join(MULTI, "kras-mutant-lung"),
     "ranked_candidates.csv", 36, "LUAD", False),
    ("head-and-neck", "HNSC", os.path.join(MULTI, "head-and-neck"),
     "candidates.csv", 72, None, False),
    ("ovarian", "Ovarian", os.path.join(MULTI, "ovarian"),
     "candidates.csv", 58, None, True),
    ("small-cell-lung", "SCLC", os.path.join(MULTI, "small-cell-lung"),
     "candidates.csv", 26, None, False),
    ("esophagogastric", "Esophagogastric adeno.", os.path.join(MULTI, "esophagogastric"),
     "candidates.csv", 45, "STAD", False),
    ("esophageal-squamous", "Esophageal SCC", os.path.join(MULTI, "esophageal-squamous"),
     "candidates.csv", 24, None, False),
    ("renal-cell", "Renal cell", os.path.join(MULTI, "renal-cell"),
     "candidates.csv", 30, None, True),
    ("bladder", "Bladder", os.path.join(MULTI, "bladder"),
     "candidates.csv", 45, None, False),
    ("endometrial", "Endometrial", os.path.join(MULTI, "endometrial"),
     "candidates.csv", 28, None, False),
    ("ewing", "Ewing sarcoma", os.path.join(MULTI, "ewing"),
     "candidates.csv", 27, None, False),
    ("neuroblastoma", "Neuroblastoma", os.path.join(MULTI, "neuroblastoma"),
     "candidates.csv", 39, None, False),
    ("aml", "AML", os.path.join(MULTI, "aml"),
     "candidates.csv", 28, None, False),
    ("osteosarcoma", "Osteosarcoma", os.path.join(MULTI, "osteosarcoma"),
     "candidates.csv", 17, None, True),
]
CANCER_BY_KEY = {c[0]: c for c in CANCERS}

def load_candidates(key, folder, fname):
    """Return dict symbol -> {enrichment, dep_frac, dep_n, group_n, q, rank, score}."""
    p = os.path.join(folder, fname)
    d = pd.read_csv(p)
    if fname == "_relaxed_scored.csv":
        d = d[(d["RelaxedBar"] == True) & (d["Pass_OE"] == True)]
    gcol = "Gene" if "Gene" in d.columns else "Symbol"
    ecol = "Enrichment"
    fcol = pick_col(d, [r".*_dep_frac", r"GroupDepFrac"])
    ncol = pick_col(d, [r".*_dep$", r"GroupDep", r".*_dep_n"])
    qcol = pick_col(d, [r"q_BH", r"BH_Q"])
    rcol = pick_col(d, [r"EnrichRank"])
    scol = pick_col(d, [r"Score", r"CombinedScore"])
    out = {}
    for r in d.itertuples(index=False):
        g = sym(getattr(r, gcol))
        if not g:
            continue
        dd = d.columns
        def val(col):
            return getattr(r, col) if col else np.nan
        out[g] = {
            "enrichment": val(ecol),
            "dep_frac": val(fcol),
            "dep_n": val(ncol),
            "q": val(qcol),
            "rank": val(rcol),
            "score": val(scol),
        }
    return out

CAND = {}
for key, disp, folder, fname, n, tcga, taz in CANCERS:
    CAND[key] = load_candidates(key, folder, fname)
    print(f"{key}: {len(CAND[key])} candidates")

# top20 sets (for reverse lookup)
TOP20 = {}
for key, disp, folder, fname, n, tcga, taz in CANCERS:
    t = pd.read_csv(os.path.join(folder, "top20.csv"))
    gcol = "Gene" if "Gene" in t.columns else "Symbol"
    TOP20[key] = set(sym(g) for g in t[gcol] if sym(g))
    print(f"{key}: top20 n={len(TOP20[key])}")

# ================================================================ DIRECTION 1
# Dependency-first: candidates ∩ surfaceome, ADC bar + safety per hit.
# ================================================================
def safety_for(symbol):
    """Return dict of safety numbers/flags from ADC tables (or NaN)."""
    s = {"normal_tissue_max": np.nan, "normal_tissue": None,
         "celltype_max": np.nan, "celltype": None,
         "critical": False, "ihc_level": None, "ihc_n_high": np.nan,
         "tissue_source": None}
    if symbol in CT.index:
        r = CT.loc[symbol]
        s["normal_tissue_max"] = r["max_normal_tissue_consensus_nTPM"]
        s["normal_tissue"] = r["max_normal_tissue"]
        s["celltype_max"] = r["max_normal_celltype_nCPM"]
        s["celltype"] = r["max_normal_celltype"]
        s["critical"] = bool(r["critical_tissue_flag"]) if pd.notna(r["critical_tissue_flag"]) else False
        s["ihc_level"] = r["ihc_max_normal_level"]
        s["ihc_n_high"] = r["ihc_n_normal_tissues_high"]
        s["tissue_source"] = "candidate_table (max over 51 tissues)"
    else:
        if symbol in NCM.index:
            r = NCM.loc[symbol]
            s["celltype_max"] = r["max_celltype_nCPM"]
            s["celltype"] = r["max_celltype"]
        if symbol in IHC.index:
            r = IHC.loc[symbol]
            s["ihc_level"] = r["max_normal_level"]
            s["ihc_n_high"] = r["n_tissues_high"]
        s["tissue_source"] = ("celltype/IHC only — no 51-tissue max in "
                               "candidate_table; full normal-tissue check still needed")
    return s

def safety_flags(s):
    flags = []
    if s["critical"]:
        flags.append("CRITICAL")
    if pd.notna(s["celltype_max"]) and s["celltype_max"] >= 2000:
        flags.append("CELLTYPE>=2000")
    if pd.notna(s["normal_tissue_max"]) and s["normal_tissue_max"] >= 200:
        flags.append("TISSUE>=200")
    return ";".join(flags)

master_rows = []
per_cancer = {}
for key, disp, folder, fname, n, tcga, taz in CANCERS:
    rows = []
    for gene, dep in CAND[key].items():
        if gene not in SURFY:
            continue
        adc_tier = int(CT.loc[gene]["tier"]) if gene in CT.index else None
        nominated = "Tier %d" % adc_tier if adc_tier else "NEW"
        # ADC per-cancer evidence (mapped cancers only: BRCA, COAD, LIHC,
        # LUAD, LUSC, READ, STAD). Tumour FPKM + fold from bulk_rna_vs_normal
        # (full coverage); protein support from enriched_per_cancer if present,
        # else CPTAC (LUAD/LUSC/BRCA/COAD), else candidate_table max (labelled).
        fpkm, fold, prot_sup, prot_fc, prot_src, adc_bar = (
            np.nan, np.nan, None, np.nan, None, "not assessable")
        if tcga:
            b = BRM.get((gene, tcga))
            if b is not None:
                fpkm = b.median_tumour_FPKM
                fold = b.fold_vs_normal_approx
            if (gene, tcga) in EPC:
                e = EPC[(gene, tcga)]
                prot_sup = bool(e.protein_support) if pd.notna(e.protein_support) else None
                prot_fc = e.protein_log2FC
                prot_src = "enriched_per_cancer"
            elif (gene, tcga) in CPT:
                prot_fc = CPT[(gene, tcga)]
                prot_sup = bool(prot_fc >= 0.5)
                prot_src = "CPTAC best isoform"
            elif gene in CT.index:
                prot_fc = CT.loc[gene]["max_protein_log2FC"]
                prot_sup = bool(pd.notna(prot_fc) and prot_fc >= 0.5)
                prot_src = "candidate_table max (other cancers)"
            if pd.notna(fpkm) and prot_sup is not None:
                adc_bar = ("met" if (fpkm >= 10 and prot_sup) else "not met")
            elif pd.notna(fpkm):
                adc_bar = "not met (no protein data)"
        elif gene in CT.index:
            # unmapped cancer: record the gene's best ADC evidence elsewhere
            r = CT.loc[gene]
            fpkm = r["max_tumour_FPKM"]
            fold = r["max_fold_floored"]
            prot_fc = r["max_protein_log2FC"]
            adc_bar = "other-cancer evidence only"
        s = safety_for(gene)
        fl = safety_flags(s)
        rows.append({
            "cancer": disp, "cancer_key": key, "gene": gene,
            "adc_status": nominated, "adc_tier": adc_tier,
            "enrichment": dep["enrichment"], "dep_frac": dep["dep_frac"],
            "dep_n": dep["dep_n"], "group_n": n, "q": dep["q"],
            "enrich_rank": dep["rank"], "score": dep["score"],
            "tumour_fpkm_adc": fpkm, "fold_adc": fold,
            "protein_support_adc": prot_sup, "protein_log2fc_adc": prot_fc,
            "protein_support_source": prot_src,
            "adc_bar": adc_bar,
            "normal_tissue_max_ntpm": s["normal_tissue_max"],
            "normal_tissue_max_tissue": s["normal_tissue"],
            "normal_celltype_max_ncpm": s["celltype_max"],
            "normal_celltype_max_type": s["celltype"],
            "critical_tissue_flag": s["critical"],
            "ihc_max_normal_level": s["ihc_level"],
            "ihc_n_tissues_high": s["ihc_n_high"],
            "safety_flags": fl,
            "taz_cancer": taz,
        })
    # rank: unflagged first, then by enrichment (most negative first)
    rows.sort(key=lambda r: (1 if r["safety_flags"] else 0,
                             r["enrichment"] if pd.notna(r["enrichment"]) else 0))
    per_cancer[key] = rows
    master_rows.extend(rows)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "per_cancer", key + "_surface_hits.csv"), index=False)
    n_new = sum(1 for r in rows if r["adc_status"] == "NEW")
    n_flag = sum(1 for r in rows if r["safety_flags"])
    print(f"{disp}: {len(rows)} surface hits ({n_new} NEW, {n_flag} safety-flagged)")

master = pd.DataFrame(master_rows)
master.to_csv(os.path.join(OUT, "intersection_master.csv"), index=False)
print("MASTER rows:", len(master), "| distinct genes:",
      master["gene"].nunique())
print("NEW nominations:", int((master["adc_status"] == "NEW").sum()),
      "| distinct NEW genes:", master[master["adc_status"] == "NEW"]["gene"].nunique())

# ================================================================ DIRECTION 2
# ADC-first: every nominated gene looked up in all 20 dependency scans.
# ================================================================
def load_scan(key, folder):
    p = os.path.join(folder, "dependency_scan.csv")
    d = pd.read_csv(p)
    gcol = d.columns[0]
    ecol = "Enrichment"
    fcol = pick_col(d, [r".*_dep_frac", r"GroupDepFrac"])
    qcol = pick_col(d, [r"q_BH", r"BH_Q"])
    rcol = pick_col(d, [r"EnrichRank"])
    out = {}
    for r in d.itertuples(index=False):
        g = sym(getattr(r, gcol))
        if not g:
            continue
        out[g] = {
            "enrichment": getattr(r, ecol),
            "dep_frac": getattr(r, fcol) if fcol else np.nan,
            "q": getattr(r, qcol) if qcol else np.nan,
            "rank": getattr(r, rcol) if rcol else np.nan,
        }
    return out

SCANS = {}
for key, disp, folder, fname, n, tcga, taz in CANCERS:
    SCANS[key] = load_scan(key, folder)
    print(f"scan {key}: {len(SCANS[key])} genes")

rev_rows = []
for gene in CT.index:
    tier = int(CT.loc[gene]["tier"])
    for key, disp, folder, fname, n, tcga, taz in CANCERS:
        s = SCANS[key].get(gene)
        if s is None:
            rev_rows.append({"gene": gene, "adc_tier": tier, "cancer": disp,
                             "cancer_key": key, "enrichment": np.nan,
                             "dep_frac": np.nan, "q": np.nan, "rank": np.nan,
                             "is_candidate": False, "in_top20": False,
                             "in_scan": False, "taz_cancer": taz})
        else:
            rev_rows.append({"gene": gene, "adc_tier": tier, "cancer": disp,
                             "cancer_key": key, "enrichment": s["enrichment"],
                             "dep_frac": s["dep_frac"], "q": s["q"],
                             "rank": s["rank"],
                             "is_candidate": gene in CAND[key],
                             "in_top20": gene in TOP20[key],
                             "in_scan": True, "taz_cancer": taz})
rev = pd.DataFrame(rev_rows)
rev.to_csv(os.path.join(OUT, "adc_reverse_lookup.csv"), index=False)
print("reverse rows:", len(rev))

# summary per ADC gene: in how many cancers a candidate / top20
sup = rev[rev["is_candidate"]].groupby("gene").agg(
    n_cand_cancers=("cancer", "nunique"),
    cand_cancers=("cancer", lambda x: "; ".join(sorted(set(x)))),
).reset_index()
sup20 = rev[rev["in_top20"]].groupby("gene").agg(
    n_top20_cancers=("cancer", "nunique"),
    top20_cancers=("cancer", lambda x: "; ".join(sorted(set(x)))),
).reset_index()
sup = sup.merge(sup20, on="gene", how="left")
sup["adc_tier"] = sup["gene"].map(lambda g: int(CT.loc[g]["tier"]))
sup = sup.sort_values(["n_cand_cancers", "n_top20_cancers"], ascending=False)
sup.to_csv(os.path.join(OUT, "adc_genes_with_dependency_support.csv"), index=False)
print("ADC genes that are candidates in >=1 cancer:", len(sup))
print("ADC genes in a top20 in >=1 cancer:", int(sup["n_top20_cancers"].fillna(0).gt(0).sum()))

# five new targets specifically
for g in ["SLC7A5", "SPINT2", "TTYH3", "SLC52A2", "EPHB3"]:
    sub = rev[(rev["gene"] == g) & (rev["is_candidate"])]
    sub20 = rev[(rev["gene"] == g) & (rev["in_top20"])]
    print(f"{g}: candidate in {sorted(sub['cancer'])}; top20 in {sorted(sub20['cancer'])}")

# ================================================================ DIRECTION 3
# Cross-cancer view over dependency-first surface hits.
# ================================================================
cc = master.groupby("gene").agg(
    n_cancers=("cancer", "nunique"),
    cancers=("cancer", lambda x: "; ".join(sorted(set(x)))),
    n_taz_cancers=("taz_cancer", lambda x: int(sum(1 for v in x if v))),
    adc_tier=("adc_tier", "first"),
    adc_status=("adc_status", "first"),
    best_enrichment=("enrichment", "min"),
).reset_index()
cc = cc.sort_values(["n_cancers", "best_enrichment"], ascending=[False, True])
cc.to_csv(os.path.join(OUT, "cross_cancer_surface.csv"), index=False)
print("cross-cancer surface genes:", len(cc))
print("recurring in >=2 cancers:", int((cc["n_cancers"] >= 2).sum()))
print("in a TAZ cancer:", int((cc["n_taz_cancers"] > 0).sum()))
print(cc.head(15).to_string())

# save a JSON of headline numbers for the verifier + report
headlines = {
    "n_cancers": 20,
    "total_candidate_slots": sum(len(CAND[k]) for k in CAND),
    "surface_hits_total": int(len(master)),
    "surface_genes_distinct": int(master["gene"].nunique()),
    "new_nominations_rows": int((master["adc_status"] == "NEW").sum()),
    "new_genes_distinct": int(master[master["adc_status"] == "NEW"]["gene"].nunique()),
    "confirmed_rows": int((master["adc_status"] != "NEW").sum()),
    "safety_flagged_rows": int((master["safety_flags"] != "").sum()),
    "adc_bar_met_rows": int((master["adc_bar"] == "met").sum()),
    "adc_genes_with_dependency": int(len(sup)),
    "adc_genes_top20": int(sup["n_top20_cancers"].fillna(0).gt(0).sum()),
    "cross_cancer_genes": int(len(cc)),
    "cross_cancer_ge2": int((cc["n_cancers"] >= 2).sum()),
    "cross_cancer_taz": int((cc["n_taz_cancers"] > 0).sum()),
    "per_cancer_hits": {k: len(per_cancer[k]) for k in per_cancer},
    "per_cancer_new": {k: int(sum(1 for r in per_cancer[k] if r["adc_status"] == "NEW"))
                       for k in per_cancer},
}
with open(os.path.join(OUT, "headlines.json"), "w") as f:
    json.dump(headlines, f, indent=2)
print(json.dumps(headlines, indent=2))
