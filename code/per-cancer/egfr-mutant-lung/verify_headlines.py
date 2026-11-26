"""INDEPENDENT verification: recompute headline numbers from the raw files
with a separate code path (direct row indexing, no shared functions with
compute_egfr_lung.py / merge_tcga.py). Prints CHECK lines to compare against
the report; any mismatch must be resolved before delivery.
"""
import json, re
import numpy as np
import pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/egfr-mutant-lung"
D = OUT + "/downloads"

# --- groups, from Model.csv + mutation file directly ---
model = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/Model.csv")
mut = pd.read_csv("/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv",
                  usecols=["HugoSymbol", "ModelID", "ProteinChange"])
eg = mut[mut.HugoSymbol == "EGFR"]
ACT = {"p.E746_A750del", "p.L747_P753delinsS", "p.L747_E749del", "p.L858R", "p.T790M"}
act = set(eg[eg.ProteinChange.isin(ACT)].ModelID)
anymut = set(eg.ModelID)
ge = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", index_col=0)
nsclc = set(model[(model.OncotreeLineage == "Lung") & (model.OncotreePrimaryDisease == "Non-Small Cell Lung Cancer")].ModelID)
cr = [m for m in ge.index if m in nsclc]
a = [m for m in cr if m in act]; b = [m for m in cr if m in anymut]; w = [m for m in cr if m not in anymut]
print("CHECK groups CRISPR NSCLC:", len(cr), "activating:", len(a), "broad:", len(b), "WT:", len(w))
assert (len(cr), len(a), len(b), len(w)) == (98, 6, 21, 77)

ex = pd.read_csv("/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
seq = set(mut.ModelID)
er = [m for m in ex.index if m in nsclc]
ae = [m for m in er if m in act]; be = [m for m in er if m in anymut]
we = [m for m in er if m in seq and m not in anymut]
print("CHECK groups RNA NSCLC:", len(er), "activating:", len(ae), "broad:", len(be), "WT:", len(we))
assert (len(er), len(ae), len(be), len(we)) == (143, 9, 28, 106)

# --- EGFR + GRB2 + KRAS Chronos, direct ---
def col(gene):
    c = [c for c in ge.columns if c.startswith(gene + " (")][0]
    return ge[c]
for gene in ["EGFR", "GRB2", "KRAS", "TEAD1", "WWTR1", "PTK2"]:
    v = col(gene)
    va, vw = v.loc[a], v.loc[w]
    dep_a = int((va < -0.5).sum()); dep_w = int((vw < -0.5).sum())
    print(f"CHECK {gene}: mean act {va.mean():.3f} ({dep_a}/{va.notna().sum()}) vs WT {vw.mean():.3f} "
          f"({dep_w}/{vw.notna().sum()}), enrich {va.mean()-vw.mean():.3f}, "
          f"genome-wide dep frac {(v < -0.5).mean():.3f}")
# KRAS Welch p recomputed
v = col("KRAS"); t, p = stats.ttest_ind(v.loc[a].dropna(), v.loc[w].dropna(), equal_var=False)
print("CHECK KRAS Welch p (raw):", f"{p:.3g}")

# --- counts from saved scans ---
p = pd.read_csv(f"{OUT}/dependency_scan.csv"); s = pd.read_csv(f"{OUT}/dependency_scan_secondary.csv")
print("CHECK primary strict/relaxed:", int(p.StrictBar.sum()), int(p.RelaxedBar.sum()),
      "| q<0.05:", int((p.BH_Q < 0.05).sum()), "| q<0.05 mutant-direction:",
      int(((p.BH_Q < 0.05) & (p.Enrichment < 0)).sum()))
print("CHECK secondary strict/relaxed:", int(s.StrictBar.sum()), int(s.RelaxedBar.sum()),
      "| q<0.05:", int((s.BH_Q < 0.05).sum()))
print("CHECK missing entries in effect matrix:", int(ge.isna().sum().sum()))
cand = pd.read_csv(f"{OUT}/ranked_candidates.csv")
print("CHECK candidates:", len(cand), "| strict among them:", int(cand.StrictBar.sum()))
print("CHECK EGFR score:", round(float(cand[cand.Gene == 'EGFR'].Score.iloc[0]), 2))

# --- TCGA EGFR RNA, from raw JSON ---
rna = pd.DataFrame(json.load(open(f"{D}/luad_tcga_rna.json")))
rna = rna[rna.entrezGeneId == 1956].copy()
rna["v"] = np.log2(rna["value"] + 1)
rna = rna[rna.sampleId.str.endswith("-01")]
muts = pd.DataFrame(json.load(open(f"{D}/luad_tcga_egfr_mutations.json")))
act_samp = set(muts[muts.proteinChange.isin(
    ["L858R", "L861Q", "G719A", "S768I", "T790M", "H773dup", "E746_A750del", "L747_E749del",
     "L747_T751del", "T751_E758del", "E709_T710delinsD", "E746_S752delinsA", "L747_A750del"])].sampleId)
seq_samp = {x["sampleId"] for x in json.load(open(f"{D}/luad_tcga_sequenced_samples.json"))}
any_samp = set(muts.sampleId)
va = rna[rna.sampleId.isin(act_samp)]["v"]; vw = rna[rna.sampleId.isin(seq_samp - any_samp)]["v"]
t, pp = stats.ttest_ind(va, vw, equal_var=False)
print(f"CHECK TCGA EGFR RNA: log2FC {va.mean()-vw.mean():.3f}, p {pp:.3g}, n {len(va)} vs {len(vw)}")
prot = pd.DataFrame(json.load(open(f"{D}/luad_cptac_protein.json")))
prot = prot[prot.entrezGeneId == 1956]
cm = pd.DataFrame(json.load(open(f"{D}/luad_cptac_egfr_mutations.json")))
cact = set(cm[cm.proteinChange.isin(
    ["L858R", "L861Q", "G719A", "G719C", "S768I", "T790M", "H773dup", "E746_A750del", "L747_E749del",
     "L747_T751del", "T751_E758del", "E709_T710delinsD", "E746_S752delinsA", "L747_A750del"])].sampleId)
cany = set(cm.sampleId)
pa = prot[prot.sampleId.isin(cact)]["value"]; pw = prot[prot.sampleId.isin(set(prot.sampleId) - cany)]["value"]
t, pp = stats.ttest_ind(pa.dropna(), pw.dropna(), equal_var=False)
print(f"CHECK CPTAC EGFR protein: delta {pa.mean()-pw.mean():.3f}, p {pp:.3g}, n {pa.notna().sum()} vs {pw.notna().sum()}")
print("CHECK CPTAC activating samples:", len(cact), "any-mutant:", len(cany))
print("ALL CHECKS DONE")
