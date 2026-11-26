"""Merge + score: TCGA LUAD RNA / CPTAC / RPPA vs DepMap primary scan (KRAS).

TCGA activating rule (cBioPortal proteinChange, no 'p.' prefix):
  missense at codon G12, G13 or Q61 — regex ^(G12|G13|Q61)[A-Z].
  Same strict rule as the cell-line analysis (excludes G138V-type and all
  non-activating variants, e.g. D33E, K88*). A sample counts as activating
  if ANY of its KRAS changes qualifies.
WT = sequenced sample with no KRAS mutation row at all.
RNA values: luad_tcga_rna_seq_v2_mrna is raw RSEM (checked: range ~0-16k);
log2(RSEM+1) transform here, log2FC = difference of group means.
Protein bar: CPTAC delta >= 0.5 and p < 0.05, or RPPA delta >= 0.5 and p < 0.05.
"""
import json, re
import numpy as np
import pandas as pd
from scipy import stats

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/kras-mutant-lung"
D = OUT + "/downloads"

def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p); order = np.argsort(p); ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1.0)
    return out

ACT_PAT = re.compile(r"^(G12|G13|Q61)[A-Z]$")

def is_activating(pc, mtype):
    if pc is None or (isinstance(pc, float) and np.isnan(pc)):
        return False
    if mtype != "Missense_Mutation":
        return False
    return bool(ACT_PAT.match(str(pc)))

def welch(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, len(a), len(b)
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return float(np.nanmean(a) - np.nanmean(b)), float(p), len(a), len(b)

# ---------- TCGA LUAD mutations ----------
mut = pd.DataFrame(json.load(open(f"{D}/luad_tcga_kras_mutations.json")))
mut["activating"] = [is_activating(pc, mt) for pc, mt in zip(mut["proteinChange"], mut["mutationType"])]
samp_act = mut.groupby("sampleId")["activating"].any()
samp_any = set(mut["sampleId"])
seq_samples = {x["sampleId"] for x in json.load(open(f"{D}/luad_tcga_sequenced_samples.json"))}
print("sequenced samples:", len(seq_samples), "| any KRAS mutation:", len(samp_any),
      "| activating:", int(samp_act.sum()))
print("activating samples:", sorted(samp_act[samp_act].index.tolist()))
print("non-activating mutant samples:", sorted(samp_act[~samp_act].index.tolist()))
print("mutant samples NOT in sequenced list (should be 0):", len(samp_any - seq_samples))
print("samples with >1 KRAS row:", mut["sampleId"].value_counts()[lambda s: s > 1].to_dict())

# ---------- TCGA RNA ----------
rna = pd.DataFrame(json.load(open(f"{D}/luad_tcga_rna.json")))
print("RNA raw value range:", rna["value"].min(), rna["value"].max(), "-> log2(RSEM+1) transform")
rna["value"] = np.log2(rna["value"] + 1)
print("RNA log2 value range:", rna["value"].min(), rna["value"].max())
print("RNA sample suffixes:", rna["sampleId"].str[-3:].value_counts().to_dict())
rna = rna[rna["sampleId"].str.endswith("-01")]  # primary tumour only
print("RNA primary-tumour samples:", rna["sampleId"].nunique())
ge_head = pd.read_csv("/home/hatch/workspace/research_notes/vrk1-codependency-20261108/CRISPRGeneEffect.csv", nrows=0)
id2gene = {}
for c in ge_head.columns[1:]:
    m2 = re.match(r"(.+) \((\d+)\)$", c)
    if m2:
        id2gene[int(m2.group(2))] = m2.group(1)
rna["Gene"] = rna["entrezGeneId"].map(id2gene)
print("RNA genes unmapped:", rna["Gene"].isna().sum(), "| distinct genes:", rna["Gene"].nunique())
fetch_genes = open(f"{D}/tcga_fetch_genes.txt").read().split()
print("fetched genes missing from RNA:", sorted(set(fetch_genes) - set(rna["Gene"].dropna())))
piv = rna.pivot_table(index="Gene", columns="sampleId", values="value")
act_samples = [s for s in piv.columns if s in samp_act[samp_act].index]
wt_samples = [s for s in piv.columns if s in seq_samples and s not in samp_any]
print("RNA: activating n =", len(act_samples), "| WT (sequenced, no KRAS mut) n =", len(wt_samples),
      "| mutant-but-not-activating with RNA:", len([s for s in piv.columns if s in samp_act[~samp_act].index]))
rows = []
for g in piv.index:
    d, p, na, nw = welch(piv.loc[g, act_samples], piv.loc[g, wt_samples])
    rows.append({"Gene": g, "TCGA_log2FC": d, "TCGA_p": p, "TCGA_n_act": na, "TCGA_n_wt": nw})
tcga = pd.DataFrame(rows)
tcga["TCGA_q"] = bh(np.where(tcga["TCGA_p"].isna(), 1.0, tcga["TCGA_p"]))
tcga.to_csv(f"{OUT}/tcga_rna.csv", index=False)
print("TCGA RNA genes passing bar (log2FC>=0.5, p<0.05):",
      int(((tcga["TCGA_log2FC"] >= 0.5) & (tcga["TCGA_p"] < 0.05)).sum()), "of", len(tcga))
print(tcga.sort_values("TCGA_log2FC", ascending=False).head(12).to_string())

# ---------- RPPA ----------
rppa = pd.DataFrame(json.load(open(f"{D}/luad_tcga_rppa.json")))
print("RPPA value range:", rppa["value"].min(), rppa["value"].max())
rppa = rppa[rppa["sampleId"].str.endswith("-01")]
rppa["Gene"] = rppa["entrezGeneId"].map(id2gene)
print("RPPA genes:", sorted(rppa["Gene"].dropna().unique()), "| samples:", rppa["sampleId"].nunique())
rp = rppa.pivot_table(index="Gene", columns="sampleId", values="value")
rrows = []
for g in rp.index:
    a = [s for s in rp.columns if s in samp_act[samp_act].index]
    w = [s for s in rp.columns if s in seq_samples and s not in samp_any]
    d, p, na, nw = welch(rp.loc[g, a], rp.loc[g, w])
    rrows.append({"Gene": g, "RPPA_delta": d, "RPPA_p": p, "RPPA_n_act": na, "RPPA_n_wt": nw})
rppa_df = pd.DataFrame(rrows)
rppa_df.to_csv(f"{OUT}/tcga_rppa.csv", index=False)
print(rppa_df.to_string())

# ---------- CPTAC protein ----------
cmut = pd.DataFrame(json.load(open(f"{D}/luad_cptac_kras_mutations.json")))
cmut["activating"] = [is_activating(pc, mt) for pc, mt in zip(cmut["proteinChange"], cmut["mutationType"])]
csamp_act = cmut.groupby("sampleId")["activating"].any()
csamp_any = set(cmut["sampleId"])
print("CPTAC samples with any KRAS mutation:", len(csamp_any), "| activating:", int(csamp_act.sum()))
print("CPTAC activating:", sorted(csamp_act[csamp_act].index.tolist()))
print("CPTAC non-activating mutant:", sorted(csamp_act[~csamp_act].index.tolist()))
prot = pd.DataFrame(json.load(open(f"{D}/luad_cptac_protein.json")))
prot["Gene"] = prot["entrezGeneId"].map(id2gene)
print("CPTAC genes unmapped:", prot["Gene"].isna().sum(), "| distinct:", prot["Gene"].nunique(),
      "| samples:", prot["sampleId"].nunique())
print("CPTAC fetched genes missing:", sorted(set(fetch_genes) - set(prot["Gene"].dropna())))
pp = prot.pivot_table(index="Gene", columns="sampleId", values="value")
prows = []
for g in pp.index:
    a = [s for s in pp.columns if s in csamp_act[csamp_act].index]
    w = [s for s in pp.columns if s not in csamp_any]
    d, p, na, nw = welch(pp.loc[g, a], pp.loc[g, w])
    prows.append({"Gene": g, "CPTAC_delta": d, "CPTAC_p": p, "CPTAC_n_act": na, "CPTAC_n_wt": nw})
cptac = pd.DataFrame(prows)
cptac.to_csv(f"{OUT}/tcga_cptac.csv", index=False)
print("CPTAC genes passing bar (delta>=0.5, p<0.05):",
      int(((cptac["CPTAC_delta"] >= 0.5) & (cptac["CPTAC_p"] < 0.05)).sum()), "of", len(cptac))
print(cptac.sort_values("CPTAC_delta", ascending=False).head(10).to_string())

# ---------- merge with dependency + cell-line RNA ----------
dep = pd.read_csv(f"{OUT}/dependency_scan.csv")
rna_cl = pd.read_csv(f"{OUT}/celline_rna.csv")
m = dep.merge(rna_cl[["Gene", "DeltaLog2", "TPM_FC", "WelchP", "BH_Q"]].rename(
    columns={"DeltaLog2": "CellLine_DeltaLog2", "TPM_FC": "CellLine_TPM_FC",
             "WelchP": "CellLine_RNA_p", "BH_Q": "CellLine_RNA_q"}), on="Gene", how="left")
m = m.merge(tcga[["Gene", "TCGA_log2FC", "TCGA_p", "TCGA_q", "TCGA_n_act", "TCGA_n_wt"]], on="Gene", how="left")
m = m.merge(cptac[["Gene", "CPTAC_delta", "CPTAC_p", "CPTAC_n_act", "CPTAC_n_wt"]], on="Gene", how="left")
m = m.merge(rppa_df[["Gene", "RPPA_delta", "RPPA_p", "RPPA_n_act", "RPPA_n_wt"]], on="Gene", how="left")
m["CellLine_RNA_pass"] = (m["CellLine_TPM_FC"] >= 1.5) & (m["CellLine_DeltaLog2"] >= 0.5) & (m["CellLine_RNA_p"] < 0.05)
m["TCGA_RNA_pass"] = (m["TCGA_log2FC"] >= 0.5) & (m["TCGA_p"] < 0.05)
m["Protein_pass"] = ((m["CPTAC_delta"] >= 0.5) & (m["CPTAC_p"] < 0.05)) | \
                    ((m["RPPA_delta"] >= 0.5) & (m["RPPA_p"] < 0.05))
m["Overexp_pass"] = m["CellLine_RNA_pass"] | m["TCGA_RNA_pass"] | m["Protein_pass"]
m.to_csv(f"{OUT}/overexpression.csv", index=False)
print("strict-bar genes passing over-expression:", int((m["StrictBar"] & m["Overexp_pass"]).sum()),
      "| relaxed-bar passing:", int((m["RelaxedBar"] & m["Overexp_pass"]).sum()))

# ---------- combined score across relaxed-bar genes ----------
rel = m[m["RelaxedBar"]].copy()
def z(s):
    s = pd.Series(s, dtype=float)
    return (s - s.mean()) / s.std(ddof=1)
rel["Score"] = z(-rel["Enrichment"]) + z(rel["FracDiff"]) + \
    z(rel["TCGA_log2FC"].fillna(0)) + z(rel["CellLine_DeltaLog2"].fillna(0)) + \
    rel["Protein_pass"].fillna(False).astype(float) * 1.0
rel = rel.sort_values("Score", ascending=False)
cand = rel[rel["Overexp_pass"]].copy()
print("relaxed-bar + over-expression candidates:", len(cand))
print(cand[["Gene", "StrictBar", "GroupMean", "GroupDep", "GroupN", "CompDep", "CompN",
            "Enrichment", "CellLine_TPM_FC", "TCGA_log2FC", "TCGA_q", "CPTAC_delta",
            "Score"]].head(30).to_string())
cand.to_csv(f"{OUT}/ranked_candidates.csv", index=False)
rel.to_csv(f"{OUT}/ranked_relaxed.csv", index=False)
print("done")
