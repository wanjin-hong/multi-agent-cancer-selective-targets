"""INDEPENDENT verification: recompute headline numbers from the raw files
with a separate code path (direct row indexing, no shared functions with
compute_kras_lung.py / merge_kras_tcga.py). Prints CHECK lines to compare
against the report; any mismatch must be resolved before delivery.
"""
import json, re
import numpy as np
import pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung"
D = OUT + "/downloads"

# --- groups, from Model.csv + mutation file directly ---
model = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/Model.csv")
mut = pd.read_csv("/home/hatch/workspace/browser_downloads/sess-1145554613/browser-download-20261108T082841.719973602Z-1-OmicsSomaticMutations.csv",
                  usecols=["HugoSymbol", "ModelID", "ProteinChange"])
kr = mut[mut.HugoSymbol == "KRAS"]
PAT = re.compile(r"^p\.(G12[A-Z]|G13[A-Z]|Q61[A-Z])")
act = set(kr[kr.ProteinChange.fillna("").apply(lambda s: bool(PAT.match(s)))].ModelID)
anymut = set(kr.ModelID)
ge = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", index_col=0)
nsclc = set(model[(model.OncotreeLineage == "Lung") & (model.OncotreePrimaryDisease == "Non-Small Cell Lung Cancer")].ModelID)
cr = [m for m in ge.index if m in nsclc]
a = [m for m in cr if m in act]; b = [m for m in cr if m in anymut]; w = [m for m in cr if m not in anymut]
print("CHECK groups CRISPR NSCLC:", len(cr), "activating:", len(a), "broad:", len(b), "WT:", len(w))
assert (len(cr), len(a), len(b), len(w)) == (98, 36, 36, 62)

ex = pd.read_csv("/home/hatch/workspace/research_notes/glioma-top12-deep-analyses-20261108/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
seq = set(mut.ModelID)
er = [m for m in ex.index if m in nsclc]
ae = [m for m in er if m in act]; be = [m for m in er if m in anymut]
we = [m for m in er if m in seq and m not in anymut]
print("CHECK groups RNA NSCLC:", len(er), "activating:", len(ae), "broad:", len(be), "WT:", len(we))
assert (len(er), len(ae), len(be), len(we)) == (143, 42, 43, 91)

def col(df, gene):
    c = [c for c in df.columns if c.startswith(gene + " (")][0]
    return df[c]

for gene in ["KRAS", "DOCK5", "TEAD1", "WWTR1", "PTK2", "PKN2", "FOSL1", "UFM1"]:
    v = col(ge, gene)
    va, vw = v.loc[a], v.loc[w]
    dep_a = int((va < -0.5).sum()); dep_w = int((vw < -0.5).sum())
    print(f"CHECK {gene}: mean act {va.mean():.3f} ({dep_a}/{va.notna().sum()}) vs WT {vw.mean():.3f} "
          f"({dep_w}/{vw.notna().sum()}), enrich {va.mean()-vw.mean():.3f}, "
          f"genome-wide dep frac {(v < -0.5).mean():.3f}")
v = col(ge, "KRAS"); t, p = stats.ttest_ind(v.loc[a].dropna(), v.loc[w].dropna(), equal_var=False)
print("CHECK KRAS Welch p (raw):", f"{p:.3g}")

# cell-line RNA FC, direct
for gene in ["KRAS", "FOSL1", "DOCK5"]:
    e = col(ex, gene)
    ea, ew = e.loc[ae], e.loc[we]
    fc = (2.0 ** ea.mean() - 1) / (2.0 ** ew.mean() - 1)
    print(f"CHECK RNA {gene}: FC {fc:.2f}x dlog2 {ea.mean()-ew.mean():.2f} (act n={len(ea)}, WT n={len(ew)})")

# counts from saved scans
p = pd.read_csv(f"{OUT}/dependency_scan.csv")
print("CHECK primary strict/relaxed:", int(p.StrictBar.sum()), int(p.RelaxedBar.sum()),
      "| q<0.05:", int((p.BH_Q < 0.05).sum()), "| q<0.05 mutant-direction:",
      int(((p.BH_Q < 0.05) & (p.Enrichment < 0)).sum()))
print("CHECK missing entries in effect matrix:", int(ge.isna().sum().sum()))
rel = p[p.RelaxedBar]
print("CHECK relaxed-bar genes with GroupN<36 (low-n trap):", int((rel.GroupN < 36).sum()), "of", len(rel),
      "| strict-bar with GroupN<36:", int((p[p.StrictBar].GroupN < 36).sum()))
cand = pd.read_csv(f"{OUT}/ranked_candidates.csv")
print("CHECK candidates:", len(cand), "| strict among them:", int(cand.StrictBar.sum()))
print("CHECK candidate genes:", ", ".join(cand.Gene.tolist()))

# TCGA DOCK5 / KRAS log2FC from raw JSON, direct
ge_head = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", nrows=0)
id2gene = {}
for c in ge_head.columns[1:]:
    m2 = re.match(r"(.+) \((\d+)\)$", c)
    if m2: id2gene[int(m2.group(2))] = m2.group(1)
rna = pd.DataFrame(json.load(open(f"{D}/luad_tcga_rna.json")))
rna["value"] = np.log2(rna["value"] + 1)
rna = rna[rna["sampleId"].str.endswith("-01")]
rna["Gene"] = rna["entrezGeneId"].map(id2gene)
tmut = pd.DataFrame(json.load(open(f"{D}/luad_tcga_kras_mutations.json")))
ACT = re.compile(r"^(G12|G13|Q61)[A-Z]$")
tmut["activating"] = [bool(ACT.match(str(pc))) and mt == "Missense_Mutation"
                      for pc, mt in zip(tmut["proteinChange"], tmut["mutationType"])]
samp_act = tmut.groupby("sampleId")["activating"].any()
seq_s = {x["sampleId"] for x in json.load(open(f"{D}/luad_tcga_sequenced_samples.json"))}
samp_any = set(tmut["sampleId"])
for gene in ["DOCK5", "KRAS"]:
    sub = rna[rna.Gene == gene]
    piv = sub.pivot_table(index="Gene", columns="sampleId", values="value")
    aa = [s for s in piv.columns if s in samp_act[samp_act].index]
    ww = [s for s in piv.columns if s in seq_s and s not in samp_any]
    d = piv.loc[gene, aa].mean() - piv.loc[gene, ww].mean()
    tt, pp = stats.ttest_ind(piv.loc[gene, aa].dropna(), piv.loc[gene, ww].dropna(), equal_var=False)
    print(f"CHECK TCGA {gene}: log2FC {d:.3f} p={pp:.3g} (act n={len(aa)}, WT n={len(ww)})")

# CPTAC DOCK5 delta from raw JSON, direct
cmut = pd.DataFrame(json.load(open(f"{D}/luad_cptac_kras_mutations.json")))
cmut["activating"] = [bool(ACT.match(str(pc))) and mt == "Missense_Mutation"
                      for pc, mt in zip(cmut["proteinChange"], cmut["mutationType"])]
csamp_act = cmut.groupby("sampleId")["activating"].any()
csamp_any = set(cmut["sampleId"])
prot = pd.DataFrame(json.load(open(f"{D}/luad_cptac_protein.json")))
prot["Gene"] = prot["entrezGeneId"].map(id2gene)
for gene in ["DOCK5", "CEP55"]:
    sub = prot[prot.Gene == gene]
    piv = sub.pivot_table(index="Gene", columns="sampleId", values="value")
    aa = [s for s in piv.columns if s in csamp_act[csamp_act].index]
    ww = [s for s in piv.columns if s not in csamp_any]
    d = piv.loc[gene, aa].mean() - piv.loc[gene, ww].mean()
    tt, pp = stats.ttest_ind(piv.loc[gene, aa].dropna(), piv.loc[gene, ww].dropna(), equal_var=False)
    print(f"CHECK CPTAC {gene}: delta {d:.3f} p={pp:.3g} (act n={len(aa)}, WT n={len(ww)})")
print("done")
