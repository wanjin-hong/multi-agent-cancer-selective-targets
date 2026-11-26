"""Build top20.csv for melanoma from the verified _top20_raw.csv + druggability notes."""
import pandas as pd

OUT = "/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/melanoma"
top = pd.read_csv(f"{OUT}/_top20_raw.csv")
druggability = {
    "SOX10": "Lineage transcription factor; no direct inhibitor — undruggable as a small-molecule target.",
    "DUSP4": "Dual-specificity MAPK phosphatase; no selective clinical inhibitor. DUSP4 loss is a known BRAF/MEK-inhibitor resistance mechanism.",
    "MITF": "Lineage transcription factor; no direct inhibitor. HDAC inhibitors modulate the MITF programme indirectly (preclinical).",
    "ZEB2": "Transcription factor; no direct inhibitor.",
    "TFAP2A": "Transcription factor; no direct inhibitor.",
    "IRF4": "Transcription factor; IMiDs (lenalidomide/pomalidomide) lower IRF4 indirectly via IKZF1/3 degradation in myeloma — not a melanoma drug.",
    "FERMT2": "Kindlin-2, integrin-adhesion scaffold; no selective drug.",
    "MDM2": "MDM2 inhibitors in clinic (idasanutlin, milademetan, alrizomadlin/APG-115), Phase 1/2, none approved; activity requires TP53 wild-type.",
    "WWTR1": "TAZ; drugged indirectly via pan-TEAD inhibitors in clinic (VT3989, ISM6331, ODM-212).",
    "MAFF": "Small-MAF transcription factor; no drug.",
    "RNASEK": "Ribonuclease; no drug. Scored in only 7 melanoma lines — least secure rank in the top 20.",
    "EGLN1": "HIF prolyl hydroxylase (PHD2); PHD inhibitors approved for anaemia (roxadustat, daprodustat) — inhibiting EGLN1 stabilises HIF, the opposite direction to a dependency.",
    "TIMM50": "Mitochondrial import machinery; no drug. Scored in only 7 melanoma lines — least secure rank in the top 20.",
    "PEA15": "No selective drug.",
    "P3H3": "Prolyl 3-hydroxylase; no drug. Scored in only 7 melanoma lines — least secure rank in the top 20.",
    "DOCK7": "Rac GEF; no selective drug.",
    "DNAJC13": "No drug.",
    "MAD2L1BP": "Spindle-checkpoint regulator; no drug.",
    "FASN": "Denifanstat (TVB-2640), first-in-class FASN inhibitor, Phase 2 in oncology (astrocytoma + bevacizumab; HER2+ breast) and MASH.",
    "PSME3": "Proteasome activator subunit; no selective drug (proteasome inhibitors act on the 20S core, not PSME3).",
}
top["Druggability"] = top.Gene.map(druggability)
cols = ["Gene", "Tier", "Mel_mean", "Mel_dep_frac", "NonMel_mean", "NonMel_dep_frac",
        "Enrichment", "FracDiff", "GenomeWide_dep_frac", "q_BH",
        "RNA_TPM_FC", "RNA_delta_log2", "RNA_q_BH",
        "TCGA_RNA_log2FC", "TCGA_RNA_q", "TCGA_RNA_n_SKCM", "TCGA_RNA_n_Other",
        "RPPA_delta", "Score", "Druggability"]
out = top[cols].rename(columns={"Score": "CombinedScore"})
out.to_csv(f"{OUT}/top20.csv", index=False)
print(out[["Gene", "Tier", "CombinedScore"]].to_string())
print("top20.csv written:", out.shape)
