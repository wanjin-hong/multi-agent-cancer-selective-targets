# Selective-targets pipeline playbook — genome-wide dependency × over-expression

The reusable method for future repeats of the TNBC selective-targets analysis. Template: `~/workspace/research_notes/depmap-tnbc-selective-targets-20261120/report.md`. Skill: `~/workspace/custom-skills/depmap_dependency_analysis/SKILL.md` (release/definition/denominator atop every table; counts from bulk files, never the portal; caveats kept in; organise by question; one-table summary alongside the full report).

## 0. What the pipeline is

For one cancer (or subtype), find genes that are BOTH selectively dependent in that cancer's cell lines AND selectively over-expressed in that cancer's samples (RNA and/or protein), ranked by a stated combined score. It is a hypothesis generator ranked by effect size, not a significance-tested discovery — state that in every report.

## 1. Inputs and data sources

| Data | Source (working route) | Notes |
|---|---|---|
| DepMap Public 24Q4 `CRISPRGeneEffect.csv` + `Model.csv` | On disk: `~/workspace/research_notes/vrk1-codependency-20261108/` | Effect matrix 1,178 lines; ~17,916 genes; ~180k missing entries (0.85%) — use NaN-aware denominators everywhere. |
| Expression matrix `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | On disk: `~/workspace/research_notes/glioma-top12-deep-analyses-20261108/` | 1,673 lines, values log₂(TPM+1), 19,193 genes. |
| Somatic mutations `OmicsSomaticMutations.csv` (642 MB) | On disk: `~/workspace/browser_downloads/sess-1145554613/` (find the exact name) | Columns: HugoSymbol, ModelID, VariantType, VepImpact, ProteinChange, Hotspot, AMClass, LikelyLoF. No Variant_Classification. Needed only for mutation-defined groups (e.g. EGFR-mutant lung). |
| TCGA RNA + clinical | cBioPortal study-scoped API (`/api/molecular-profiles/{id}/molecular-data/fetch`, `/api/studies/{id}/...`) + GDC legacy clinical file (`api.gdc.cancer.gov/data/`) | See §4 for the three comparison designs. |
| TCGA tumour-vs-normal counts | GDC STAR counts via `api.gdc.cancer.gov/data/` with `{"ids": [...]}` wrapping | Used for HCC (LIHC). 120 random primary tumours vs all 50 normals, seed recorded. |
| Protein: CPTAC, RPPA | cBioPortal protein profiles (`{study}_protein_quantification`, `{study}_rppa`) | Often unavailable or tiny — see §5 and traps. |
| DepMap proteomics | — | **No bulk file obtainable** (download API returns HTML). Do not waste time on it; record as a failed route. |

## 2. Define the group and comparator — state the exact rule

- **Subtype split (best):** within one tissue, group vs the rest of that tissue (TNBC vs non-TNBC breast; EGFR-mutant vs WT NSCLC). Avoids the lineage tautology. Requires a defensible subtype rule (TNBC: ESR1 < 3.0, PGR < 1.0, ERBB2 < 7.5 log₂(TPM+1), validated against annotations; EGFR: activating mutations on ProteinChange, exon-19 del / L858R / T790M, with exclusions stated).
- **Whole-cancer vs pan-cancer:** group = all lines of the cancer, comparator = all non-cancer lines (HCC, melanoma, colon, glioma). State the tissue-confound explicitly: "selective" means selective vs the pan-cancer average, and lineage markers pass both bars partly by construction. Exclude non-cancerous lines from the group; decide explicitly whether they stay in the comparator and say so.
- **Mutation-defined groups:** state exactly what is included and excluded (AMClass/Hotspot are insufficient on their own — check ProteinChange; e.g. R309G is "likely_pathogenic" and was excluded in EGFR lung). Run a secondary broad-set sensitivity analysis if the primary set is small, and report both.
- **Small-group rules:** record parent/derivative pairs (HCC-827/HCC-827-GR5 — effective n ≈ 5, not 6; PA-TU-8988S/T in PDAC — effective n ≈ 44, harmless there); record lines in the RNA group but not the CRISPR group (PC-9, NCI-H3255, HCC4006 in EGFR lung; KP-1N, PL4 in PDAC); flag genes scored in few lines (CXADR 12/10 in TNBC; RNASEK/TIMM50/P3H3 in 7 melanoma lines) as least secure. **New trap (PDAC, 20 Nov 2026): a block of 127 genes was scored in only 24 of 45 lines — a fixed subset of established lines** (AsPC-1, BxPC-3, Capan-1, …). Any enrichment on those genes rests on a different, smaller line set than the rest of the scan; two of them landed in the PDAC top 20 (SPRR1B #3, CXADR #20) and one was candidate #25 (GGT1). Always report the per-gene scored-line count for the top 20, and check whether the low-n genes cluster in a fixed subset before trusting their ranks.

## 3. Genome-wide dependency scan

For every gene in `CRISPRGeneEffect.csv`, compute: group mean Chronos, dependent fraction (Chronos < −0.5), comparator mean and fraction, **enrichment = group mean minus comparator mean**, dependent-fraction difference, Welch's t-test with Benjamini–Hochberg correction across all ~17,916 genes, genome-wide dependent fraction (across all 1,178 lines). Pan-essential = dependent in ≥ 70% of lines in *both* groups (excluded from selective lists, reported separately).

**Bars (fixed — do not tune them to the results):**
- Strict: enrichment ≤ −0.10, group dep frac ≥ 0.25, frac diff ≥ 0.10, not pan-essential, genome-wide dep frac ≤ 0.60.
- Relaxed (pre-defined before any over-expression result is seen): enrichment ≤ −0.05, dep frac ≥ 0.20, frac diff ≥ 0.10, not pan-essential, genome-wide dep frac ≤ 0.60.

**Statistics handling:** NaN-aware denominators (dependent fractions use scored lines as denominator; Welch tests omit missing values); BH correction with NaN p-values handled (genes never scored get no p-value, not a zero). Report how many genes reach q < 0.05 and in which direction — if none do, say so; the ranking is then effect-size only.

## 4. Over-expression — three axes, and the TCGA comparison design

**Cell-line RNA (DepMap 24Q4 expression matrix):** group vs comparator lines; Δlog₂ = mean log₂(TPM+1) difference; TPM fold change by the verified formula below; Welch p < 0.05.

**TPM fold-change formula (verified — do not use the ratio of mean TPMs):**
`FC = (2^mean_log2TPM1_group − 1) / (2^mean_log2TPM1_comparator − 1)` — the back-transform of the group means. The ratio of arithmetic-mean TPMs gives different numbers (e.g. MT1E 3.52× vs 6.24×).

**TCGA RNA — pick ONE design before seeing results, and state it:**
- (a) **Within-cancer split** (TNBC: IHC subtype; EGFR lung: mutation status) — the cleanest; use the GDC legacy clinical file for subtype/status (HER2-equivocal resolved by FISH) and cBioPortal RNA Seq V2 RSEM, log₂(RSEM+1), primary tumours only.
- (b) **Tumour vs adjacent normal** (HCC: GDC STAR counts, TPM, log₂(TPM+1)) — use when the study has normal samples; it is the check on lineage tautology. **New caveat (PDAC, 20 Nov 2026): TCGA-PAAD has only 4 adjacent normals, so the axis cannot discriminate** — the largest fold changes in the fetched set (SPRR1B 2.34, KLF5 2.23, FOSL1 1.85) all carry p ≈ 0.13–0.20 and only 2 of 177 genes pass the bar. When normals number ≤ ~5, apply the bar exactly as pre-defined (do not relax it to rescue the axis), use all tumours (no subsampling — the study is small), and state plainly that the intersection rests on cell-line RNA and the score's TCGA component is noise-dominated. Check the normal count in the GDC file query *before* promising the check.
- (c) **Tumour vs pooled other-cancer tumours** (melanoma SKCM, colon, glioma) — use only when (a) and (b) are impossible (SKCM has zero adjacent normals in cBioPortal; verify and record that). Choose the pooled set before seeing results. **Scores are comparable within a cancer only** — the TCGA component differs by design across cancers.

**Protein:** CPTAC via cBioPortal where a profile exists; RPPA as a fallback (covers a handful of genes). If no protein comparator exists, drop the protein bar and the +1.0 bonus and say so (glioma did this). **New caveat (PDAC, 20 Nov 2026): a CPTAC study existing does not mean a protein comparison exists** — `paad_cptac_2021` has 140 tumour samples and zero normals, so no within-study tumour-vs-normal protein comparison is possible; check sample types via `/samples/fetch` before assuming one. The GDC-harmonised `pancreas_cptac_gdc` study has no protein profile at all. In that case report CPTAC tumour medians descriptively only (no bar, no bonus) and let RPPA carry the protein leg.

**Over-expression bar (any one suffices):** cell-line RNA TPM FC ≥ 1.5 AND Δlog₂ ≥ 0.5 AND p < 0.05; or TCGA RNA log₂FC ≥ 0.5 AND p < 0.05; or protein Δ ≥ 0.5 AND p < 0.05.

## 5. Intersection, ranking, and the top 20

- Intersect: genes passing a dependency bar AND the over-expression bar. Report both counts (strict + OE; relaxed + OE). If fewer than 20 pass both relaxed bars, the list stops there — do not pad it (EGFR lung stopped at 17).
- **Combined score (fixed):** for each relaxed-bar gene, score = z(−enrichment) + z(frac diff) + z(TCGA log₂FC) + z(cell-line Δlog₂), each z computed across the relaxed-bar genes, **+1.0 if the protein bar passes**. Missing TCGA or cell-line values are z-scored as 0 (no evidence, neither for nor against). Protein gets a fixed bonus, not a z, because its n is too small for a z to be meaningful.
- Rank the top 20 by combined score; label Tier A (passes the strict dependency bar) and Tier B (relaxed bar only). Report ranks 21–25 as well.
- **Druggability:** one line per top-20 gene, public knowledge only (trial registries, literature), as of the run date; do not over-claim. State how many have a drug at the gene itself.
- **Check section (Q4):** where do the programme's prior headline genes rank in this scan (TEAD1/YAP1/WWTR1/PTK2/PKN2 at minimum; plus any earlier top genes for that cancer)? Report ranks and over-expression pass/fail for each — this is the cross-check, not the ranking.
- **Decisive experiment (Q5):** propose one experiment that adjudicates the list (usually drug dose–response across the group + comparator lines, testing whether sensitivity tracks the CRISPR score for the same gene in the same line; plus a mechanism test for the dominant theme, e.g. proliferation adjustment).

## 6. Verification (mandatory before delivery)

Recompute headline numbers from the raw files with an **independent script** (separate code path from the pipeline): dependency means/fractions for the top genes, cell-line RNA Δlog₂ and TPM FC (using the §4 formula), TCGA log₂FC from the raw downloads, protein deltas where applicable. Fix and record any error found (the TNBC run caught a CXADR denominator error; the EGFR run caught an untransformed-RSEM merge error and a wild-type group containing unsequenced samples). Report the verification in the methods appendix.

## 7. Deliverables (per cancer)

`~/workspace/research_notes/multi-cancer-selective-targets-20261120/<cancer>/` (or the equivalent root for a new programme): `report.md` (Q1 dependency scan, Q2 over-expression, Q3 intersection + top 20, Q4 check, Q5 decisive experiment, caveats, methods appendix with failed routes), `one-table-summary.md`, `dependency_scan.csv` (all ~17,916 genes), `overexpression.csv` (all genes, all evidence columns), `top20.csv`. Keep scripts and raw downloads in the folder. No voice-note summary.

## 8. Working download routes (and the dead ones)

- **Works:** cBioPortal study-scoped endpoints — `/api/molecular-profiles/{profileId}/molecular-data/fetch` (POST, `entrezGeneIds`); study-scoped clinical-data fetch works for some studies but returned empty for others — fall back to the GDC legacy clinical file. GDC `/data` endpoint needs `{"ids": [...]}` wrapping; GDC STAR counts for tumour-vs-normal. `curl` for chunked cBioPortal responses (Python `urllib` fails with `IncompleteRead`). **Path-based study endpoints only:** `/api/molecular-profiles?studyId=...` ignores the query parameter in this deployment and returns the first 100 profiles alphabetically regardless of study — use `/api/studies/{studyId}/molecular-profiles` and `/api/studies/{studyId}/sample-lists` instead (PDAC, 20 Nov 2026).
- **Dead — do not retry:** UCSC Xena bulk downloads (S3 403 on all host variants); cBioPortal global `/api/clinical-data/fetch` and `/api/molecular-data/fetch` (empty or error); cBioPortal `/api/genes/fetch` (400 — map via Entrez IDs in the DepMap column names instead); `/api/sample-lists/{id}/samples` (404 — work around via matrix headers and POST `/samples/fetch`); DepMap 24Q4 proteomics (download API returns HTML); CPTAC portal API (HTTP 500) and PDC GraphQL (empty for LIHC).
- **Traps:** CPTAC cross-study values are on different scales (coadread uncentred intensity mean ~19.2 vs BRCA/OV log-ratio mean ~0) — never compare CPTAC deltas across studies; exclude and use RPPA instead. TCGA RSEM must be log₂-transformed before fold changes (untransformed RSEM gives absurd fold changes up to ~27,000). The wild-type group must exclude unsequenced samples. Some genes are absent from the 24Q4 scan (DNA2) or have no expression data (CRIPTO, HAPSTR1) — record, don't impute. **A block of genes may be scored in only a fixed subset of lines** (PDAC: 127 genes in exactly 24 of 45 lines) — check per-gene scored-line counts for the top 20. **A CPTAC study may contain no normal samples** (paad_cptac_2021: 140 tumours, 0 normals) — check sample types before assuming a protein comparison exists. **A ref-normal z-score profile may not exist for the study** (none for PAAD) — the ref-normal z column is optional, not guaranteed. **cBioPortal's `/api/molecular-profiles?studyId=` query parameter is ignored** — use path-based study endpoints.

## 9. Known limitations to carry into every future repeat

- Power tracks group n (n = 6 → nothing significant; n = 66–72 → thousands). At n < 15 the list is a hypothesis generator only.
- The intersection is blind to mutation-driven dependencies that are not over-expressed (BRAF/MAPK1/CCND1, KRAS, EGFR-broad).
- Pan-cancer comparators make lineage markers pass both bars partly by construction; the tumour-vs-normal axis is the check.
- TCGA components are not comparable across cancers (different designs).
- Protein is the weakest leg; where RNA and protein disagree, the disagreement is the finding.
- A knockout dependency is not a drug response (paralog compensation, incomplete knockdown, removal vs inhibition).
