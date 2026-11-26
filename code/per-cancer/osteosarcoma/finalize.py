"""Build top20.csv for the osteosarcoma repeat, with druggability notes
(public knowledge only, as of 21 Nov 2026)."""
import pandas as pd
OUT = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/osteosarcoma/'
cand = pd.read_csv(OUT + 'candidates.csv')
top20 = cand.head(20).copy()

DRUG = {
    'RUNX2': 'None — Runt-domain transcription factor (osteoblast master regulator); no direct inhibitor.',
    'FGFR1': 'YES — approved pan-FGFR inhibitors exist (erdafitinib, pemigatinib, futibatinib) in FGFR-altered cancers; not osteosarcoma-labelled. FGFR1 amplification occurs in a minority of osteosarcomas.',
    'PDGFRA': 'YES — approved inhibitors (imatinib, sunitinib, avapritinib for PDGFRA D842V GIST, ripretinib); not osteosarcoma-labelled.',
    'FERMT2': 'None — kindlin-2 adhesion scaffold; no drug.',
    'GPX8': 'None — glutathione peroxidase 8 (ER-resident); no drug.',
    'EXT1': 'None — exostosin glycosyltransferase (heparan sulfate synthesis); no drug.',
    'CBFB': 'None — RUNX co-factor transcription subunit; no drug.',
    'IRS1': 'None — insulin receptor substrate 1 adaptor; no direct inhibitor (IGF1R-axis upstream).',
    'WWTR1': 'INDIRECT only — TAZ is drugged via pan-TEAD inhibitors in clinic (VT3989, ISM6331, ODM-212); no TAZ-selective drug.',
    'SURF4': 'None — ER cargo receptor; no drug.',
    'SEC23A': 'None — COPII vesicle coat component; no drug.',
    'SLC35B2': 'None — nucleotide-sugar transporter; no drug.',
    'ITGAV': 'None selective — cilengitide (αv integrin) already failed Phase 3 in GBM; no approved αv drug.',
    'E2F3': 'None — E2F transcription factor; no direct inhibitor.',
    'SCNM1': 'None — splicing factor; no drug.',
    'UBE2S': 'None — E2 ubiquitin-conjugating enzyme; no drug.',
    'EXT2': 'None — exostosin glycosyltransferase, EXT1 partner; no drug.',
    'FBXO42': 'None — F-box protein (SCF ubiquitin ligase subunit); no drug.',
    'MT2A': 'None — metallothionein 2A; undruggable, culture-sensitive (same family as TNBC/glioma/PDAC MT1E).',
    'H2AC11': 'None — histone H2A family member; no drug.',
}
top20['Rank'] = range(1, 21)
top20['Druggability'] = top20['Symbol'].map(DRUG)
cols = ['Rank', 'Symbol', 'Tier', 'EnrichRank', 'OS_mean', 'OS_n', 'OS_dep', 'OS_dep_frac',
        'Comp_mean', 'Comp_dep_frac', 'Enrichment', 'FracDiff', 'GenomeWide_dep_frac',
        'Dep_p', 'Dep_q', 'OS_mean_log2TPM', 'Comp_mean_log2TPM', 'DeltaLog2', 'TPM_FC',
        'CellRNA_p', 'TCGA_RNA_n_OS', 'TCGA_RNA_n_other', 'TCGA_RNA_log2FC', 'TCGA_RNA_p',
        'TCGA_RNA_q', 'Score', 'Druggability']
top20 = top20[cols]
top20.to_csv(OUT + 'top20.csv', index=False)
print(top20[['Rank', 'Symbol', 'Tier', 'OS_dep', 'Enrichment', 'TPM_FC', 'TCGA_RNA_log2FC', 'Score']].to_string(index=False))
print('druggable at the gene in top 20:', ['FGFR1', 'PDGFRA'])
print('ranks 21-25:', pd.read_csv(OUT + 'candidates.csv').iloc[20:25]['Symbol'].tolist())
