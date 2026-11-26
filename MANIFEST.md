# MANIFEST — staged GitHub repository contents

Staged 26 November 2026 from the working programme folders. All `code/` and `data/` files are **byte-identical copies** of their originals (md5-verified: the original 9 `code/` files spot-checked, and all 147 files added in the second staging pass verified file-by-file, 0 mismatches). Nothing in `code/` or `data/` was edited. Total staged size: **2.4 MB** (201 files: 156 in `code/`, 42 in `data/`, plus README, LICENSE and this manifest). The folder on disk is 3.1 MB because it also holds a local `.git` directory from the first push; `.git` is not staged content and is not uploaded.

Programme folder (P): `/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120`
TNBC folder (T): `/home/hatch/workspace/research_notes/depmap-tnbc-selective-targets-20261120`

## Included files

### code/ (9 files)

| File in repo | Original source path | Notes |
|---|---|---|
| `code/pipeline-playbook.md` | `P/pipeline-playbook.md` | The fixed playbook every worker ran. Contains absolute workspace paths in its text (documentation only). |
| `code/build_figs_tables.py` | `P/manuscript/build_figs_tables.py` | Rebuilds manuscript Figures 1–5 and Tables 1–5. Hard-codes `BASE`/`TNBC` — see README. |
| `code/write_source_data.py` | `P/manuscript/write_source_data.py` | Writes per-figure Source Data files. Hard-codes `BASE`/`TNBC`; also reads `drug-target-calibration/target_summary.csv` and `mesothelioma/tead_axis.csv` (not staged). |
| `code/verify_source_data.py` | `P/manuscript/verify_source_data.py` | 73 checks on the Source Data files. Hard-codes `BASE`/`TNBC`; also reads `tabula-sapiens/per_gene_summary_truecpm.csv` (not staged). |
| `code/verify_headlines.py` | `P/manuscript/verify_headlines.py` | 110 checks on the manuscript headlines. Hard-codes `BASE`/`TNBC`; also reads `mesothelioma/tead_axis.csv`, `drug-target-calibration/target_cancer_matrix.csv` and per-cancer `dependency_scan.csv` files (not staged). |
| `code/census_query.py` | `P/tabula-sapiens/census_query.py` | **Audit only.** Normalises by query-set counts; its absolute CPM values are inflated ~45×. Use `census_query_truecpm.py` instead. |
| `code/census_query_truecpm.py` | `P/tabula-sapiens/census_query_truecpm.py` | Tabula query with true per-cell-total normalisation (Census obs `raw_sum`). Takes paths as CLI arguments. |
| `code/analyze_tabula_truecpm.py` | `P/tabula-sapiens/analyze_tabula_truecpm.py` | Analysis of the true-CPM query output. Path-relative; no editing needed. |
| `code/verify_headlines_truecpm.py` | `P/tabula-sapiens/verify_headlines_truecpm.py` | 62 checks on the corrected Tabula numbers. Hard-codes `TAB` to the original Tabula folder — one line to edit. |

### code/calibration/ (4 files)

| File in repo | Original source path |
|---|---|
| `code/calibration/calibrate.py` | `P/drug-target-calibration/calibrate.py` |
| `code/calibration/conditional_checks.py` | `P/drug-target-calibration/conditional_checks.py` |
| `code/calibration/build_summary.py` | `P/drug-target-calibration/build_summary.py` |
| `code/calibration/verify_headlines.py` | `P/drug-target-calibration/verify_headlines.py` |

### code/partner-mining/ (2 files)

| File in repo | Original source path |
|---|---|
| `code/partner-mining/codepend.py` | `P/klf5-taz-partners/codepend.py` |
| `code/partner-mining/verify_headlines.py` | `P/klf5-taz-partners/verify_headlines.py` |

### code/adc-intersection/ (2 files)

| File in repo | Original source path |
|---|---|
| `code/adc-intersection/intersect.py` | `P/adc-intersection/intersect.py` |
| `code/adc-intersection/verify_headlines.py` | `P/adc-intersection/verify_headlines.py` |

### code/per-cancer/ (139 files, 19 cancers)

One subfolder per cancer, named after the source folder. All files are byte-identical copies of the `.py` scripts in that cancer's programme folder. Glioma has no `.py` files in its programme folder (data and reports only), and TNBC's scan folder (`depmap-tnbc-selective-targets-20261120`) also contains no `.py` files, so neither has a subfolder here.

**`code/per-cancer/aml/` (8 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `merge.py`, `mutation_check.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/bladder/` (11 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `fgfr3_check.py`, `finalize.py`, `gdc_fetch.py`, `gdc_tumour_normal.py`, `luminal_basal_check.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/colon/` (6 files)** — `step1_depmap_scan.py`, `step2_cbioportal_fetch.py`, `step3_tcga_pancan_fetch.py`, `step4_protein_fetch.py`, `step5_merge_score.py`, `verify_headlines.py`

**`code/per-cancer/egfr-mutant-lung/` (4 files)** — `compute_egfr_lung.py`, `fetch_tcga.py`, `merge_tcga.py`, `verify_headlines.py`

**`code/per-cancer/endometrial/` (8 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `gdc_tumour_normal.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/esophageal-squamous/` (9 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `gdc_fetch.py`, `gdc_tumour_normal.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/esophagogastric/` (9 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `gdc_fetch.py`, `gdc_tumour_normal.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/ewing/` (5 files)** — `celline_rna.py`, `fusion_check.py`, `merge.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/hcc/` (6 files)** — `celline_rna.py`, `finalize.py`, `gdc_tumour_normal.py`, `merge.py`, `scan.py`, `verify.py`

**`code/per-cancer/head-and-neck/` (8 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `gdc_tumour_normal.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/kras-mutant-lung/` (6 files)** — `compute_kras_lung.py`, `fetch_kras_tcga.py`, `fetch_kras_tcga_get.py`, `fetch_per_gene.py`, `merge_kras_tcga.py`, `verify_headlines.py`

**`code/per-cancer/melanoma/` (5 files)** — `build_top20.py`, `step1_scan.py`, `step2_fetch_tcga.py`, `step3_merge.py`, `verify_headlines.py`

**`code/per-cancer/mesothelioma/` (2 files)** — `meso_axis.py`, `meso_nf2.py`

**`code/per-cancer/neuroblastoma/` (12 files)** — `celline_rna.py`, `finalize.py`, `gdc_cross_cancer.py`, `gdc_fetch.py`, `gdc_perfile.py`, `gdc_retry_curl.py`, `gdc_threaded.py`, `merge.py`, `mutations.py`, `mycn_split.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/osteosarcoma/` (8 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `merge.py`, `mutations.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/ovarian/` (7 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/pancreatic/` (8 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `gdc_tumour_normal.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify.py`

**`code/per-cancer/renal-cell/` (10 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `gdc_fetch.py`, `gdc_perfile.py`, `gdc_tumour_normal.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

**`code/per-cancer/small-cell-lung/` (7 files)** — `celline_rna.py`, `fetch_cbioportal.py`, `finalize.py`, `merge.py`, `process_tcga.py`, `scan.py`, `verify_headlines.py`

### data/ — top-20 lists (20 files)

| File in repo | Original source path |
|---|---|
| `data/aml-top20.csv` | `P/aml/top20.csv` |
| `data/bladder-top20.csv` | `P/bladder/top20.csv` |
| `data/colon-top20.csv` | `P/colon/top20.csv` |
| `data/egfr-mutant-lung-top20.csv` | `P/egfr-mutant-lung/top20.csv` |
| `data/endometrial-top20.csv` | `P/endometrial/top20.csv` |
| `data/esophageal-squamous-top20.csv` | `P/esophageal-squamous/top20.csv` |
| `data/esophagogastric-top20.csv` | `P/esophagogastric/top20.csv` |
| `data/ewing-top20.csv` | `P/ewing/top20.csv` |
| `data/glioma-top20.csv` | `P/glioma/top20.csv` |
| `data/hcc-top20.csv` | `P/hcc/top20.csv` |
| `data/head-and-neck-top20.csv` | `P/head-and-neck/top20.csv` |
| `data/kras-mutant-lung-top20.csv` | `P/kras-mutant-lung/top20.csv` |
| `data/melanoma-top20.csv` | `P/melanoma/top20.csv` |
| `data/neuroblastoma-top20.csv` | `P/neuroblastoma/top20.csv` |
| `data/osteosarcoma-top20.csv` | `P/osteosarcoma/top20.csv` |
| `data/ovarian-top20.csv` | `P/ovarian/top20.csv` |
| `data/pancreatic-top20.csv` | `P/pancreatic/top20.csv` |
| `data/renal-cell-top20.csv` | `P/renal-cell/top20.csv` |
| `data/small-cell-lung-top20.csv` | `P/small-cell-lung/top20.csv` |
| `data/tnbc-top20.csv` | `T/top20.csv` |

### data/ — candidate tables (20 files + 2 extras)

Genes passing both the dependency bar and the over-expression bar, per cancer. The first five scans and the TNBC scan predate the standardised `candidates.csv` schema, so their equivalents are included under their original filenames (see table).

| File in repo | Original source path | Notes |
|---|---|---|
| `data/aml-candidates.csv` | `P/aml/candidates.csv` | |
| `data/bladder-candidates.csv` | `P/bladder/candidates.csv` | |
| `data/colon-candidates.csv` | `P/colon/candidates_both.csv` | Same concept (both bars), earlier filename. |
| `data/egfr-mutant-lung-candidates.csv` | `P/egfr-mutant-lung/ranked_candidates.csv` | Earlier filename. |
| `data/endometrial-candidates.csv` | `P/endometrial/candidates.csv` | |
| `data/esophageal-squamous-candidates.csv` | `P/esophageal-squamous/candidates.csv` | |
| `data/esophagogastric-candidates.csv` | `P/esophagogastric/candidates.csv` | |
| `data/ewing-candidates.csv` | `P/ewing/candidates.csv` | |
| `data/glioma-candidates.csv` | `P/glioma/candidates.csv` | |
| `data/hcc-candidates.csv` | `P/hcc/candidates.csv` | |
| `data/head-and-neck-candidates.csv` | `P/head-and-neck/candidates.csv` | |
| `data/kras-mutant-lung-candidates.csv` | `P/kras-mutant-lung/ranked_candidates.csv` | Earlier filename. |
| `data/melanoma-candidates.csv` | `P/melanoma/_relaxed_scored.csv` | Relaxed-bar scored list. |
| `data/melanoma-candidates-strict.csv` | `P/melanoma/_strict_genes.csv` | Strict-bar list only (607 bytes). Extra, for completeness. |
| `data/neuroblastoma-candidates.csv` | `P/neuroblastoma/candidates.csv` | |
| `data/osteosarcoma-candidates.csv` | `P/osteosarcoma/candidates.csv` | |
| `data/ovarian-candidates.csv` | `P/ovarian/candidates.csv` | |
| `data/pancreatic-candidates.csv` | `P/pancreatic/candidates.csv` | |
| `data/renal-cell-candidates.csv` | `P/renal-cell/candidates.csv` | |
| `data/small-cell-lung-candidates.csv` | `P/small-cell-lung/candidates.csv` | |
| `data/tnbc-candidates.csv` | `T/candidates_both.csv` | Same concept (both bars), earlier filename. |
| `data/aml-top20-blood-only-comparator.csv` | `P/aml/top20_haem.csv` | **Extra, beyond the brief:** the AML top 20 against a blood-only comparator. Documents the comparator artefact (the TEAD-axis "reversal" seen against the pan-cancer comparator vanishes against the blood-only comparator). |

## Deliberately excluded — and why

| Excluded | Why |
|---|---|
| Raw DepMap matrices (`CRISPRGeneEffect.csv`, `Model.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `OmicsSomaticMutations.csv`) | Public downloads from the DepMap portal (depmap.org/portal); large (the mutation file alone is 642 MB). Named in the README for the user to download. |
| Per-cancer genome-wide intermediates (`dependency_scan.csv`, `overexpression.csv`, `celline_rna.csv`, `tcga_rna.csv`, `tcga_cptac.csv`, `tcga_rppa.csv`) | 17,916-row intermediate outputs per cancer; regenerable from the raw matrices by following the playbook. The derived candidate tables in `data/` are the lists that passed both bars. |
| Per-cancer reports, one-table summaries, `synthesis.md` | Narrative documents, not code or data. The manuscript is the citable record. |
| Manuscript text (`manuscript.md`), figures (PNG/PDF), tables (.md/.xlsx), Source Data files | The submission package. Figures and tables are rebuilt from `data/` by `code/build_figs_tables.py`; Source Data files by `code/write_source_data.py`. Not staged here. |
| Tabula per-gene summaries (`per_gene_summary.csv`, `per_gene_summary_truecpm.csv`, `background_summary*.csv`) and per-organ query outputs | Regenerable with `code/census_query_truecpm.py`. The original `per_gene_summary.csv` (query-set-relative normalisation) is retained in the programme folder for audit. |
| Drug-target calibration outputs (`target_cancer_matrix.csv`, `target_summary.csv`), ADC intersection outputs, KLF5/TAZ partner-mining outputs, `mesothelioma/tead_axis.csv` | The leg **scripts** are now staged (`code/calibration/`, `code/partner-mining/`, `code/adc-intersection/`); only their output data files remain unstaged. Note: `verify_headlines.py`, `write_source_data.py` and `verify_source_data.py` read some of these outputs, so they still cannot run end-to-end from this checkout alone (see README). |
| Melanoma `_top20_raw.csv` | Superseded intermediate; `melanoma-top20.csv` is the final list. |
| `__pycache__/`, run logs (`.log`), build logs | Build artefacts. |
| Personal information | None present in any included file. The only personal name is the author's, in the README citation line and the LICENSE copyright line, as intended. |

## Decisions needed from the author (nothing here is published or uploaded)

1. **Repository name and visibility** (public vs private) — to be decided when creating the GitHub repository.
2. **Whether to include the downstream-leg output data files** (calibration `target_cancer_matrix.csv`/`target_summary.csv`, ADC intersection outputs, partner-mining outputs, `mesothelioma/tead_axis.csv`) — the leg scripts were staged in the second pass (26 Nov 2026 evening); the outputs are still not staged, and they are needed for `verify_headlines.py` to run end-to-end.
3. **Whether the AML blood-only-comparator top 20 and the melanoma strict list should stay in** — both were included as extras beyond the brief; they can be removed with a single delete each.
4. **Zenodo archiving** — a frozen release DOI is still to be created (currently a placeholder in the manuscript); the GitHub release is normally the trigger.
