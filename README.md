# A human-guided multi-agent programme for systematic cancer selective-dependency analysis across twenty cancers

Code and derived data accompanying a manuscript in preparation for *Nature Computational Science* (sole author WanJin Hong). **The manuscript has not been submitted, is not under review, and is not published.** Nothing in this repository should be cited as a published work.

## What this programme is

One coordinating agent and twenty worker agents ran the **same fixed playbook** (`code/pipeline-playbook.md`) in twenty cancers and subtypes: a genome-wide CRISPR dependency screen intersected with an over-expression bar, to find genes a cancer selectively depends on *and* over-expresses.

- **Dependency bar:** group mean Chronos gene-effect score, dependent fraction (Chronos < −0.5), enrichment versus a comparator, Welch's t-test with Benjamini–Hochberg correction across all ~17,916 genes. Strict bar: enrichment ≤ −0.10, group dependent fraction ≥ 0.25, fraction difference ≥ 0.10, not pan-essential, genome-wide dependent fraction ≤ 0.60. A pre-defined relaxed bar is also used. The bars were fixed before any scan ran, and no worker could change them.
- **Over-expression bar (any one suffices):** cell-line RNA TPM fold change ≥ 1.5 with Δlog₂ ≥ 0.5 and p < 0.05; or TCGA RNA log₂FC ≥ 0.5 and p < 0.05; or protein Δ ≥ 0.5 and p < 0.05.
- **Result:** 385 top-20 list slots, 269 distinct genes, 200 singletons. There is **no pan-cancer selective target**. The most recurrent gene, KLF5, appears in eight cancers and is undruggable. WWTR1/TAZ is the best-supported TEAD-axis target (six cancers); TEAD1 never enters a top 20.
- **Legs beyond the twenty scans:** a normal-tissue leg in Tabula Sapiens 2.0, a calibration leg putting 20 approved or clinical drug targets through the identical bars, partner mining, and an antibody–drug conjugate (ADC) surfaceome intersection. The scripts for all of these legs are staged here (`code/calibration/`, `code/partner-mining/`, `code/adc-intersection/`), as are the per-cancer pipeline scripts (`code/per-cancer/` — see the section below). The per-cancer reports, the synthesis, and the leg *output data files* (calibration matrices, intersection tables, partner tables) are not included (see `MANIFEST.md`).

## How to reproduce

### Requirements

- Python 3.12 with `pandas`, `numpy`, `scipy`, `matplotlib`.
- For the Tabula Sapiens leg only: `cellxgene_census` (and its `tiledbsoma` backend) plus network access to the CZ CELLxGENE Census.

### Input data (download; not included)

The raw DepMap matrices are public downloads from the DepMap portal (https://depmap.org/portal) — download them, do not expect them here:

- `CRISPRGeneEffect.csv` — DepMap Public **24Q4**, 1,178 cell lines × 17,916 genes, Chronos scores
- `Model.csv` — lineage and subtype annotation
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` — log₂(TPM+1) expression, 1,673 lines
- `OmicsSomaticMutations.csv` — somatic mutations (loss of function defined by `LikelyLoF = True`; the file repeats each variant 2–3 times, so deduplicate)

TCGA RNA, CPTAC/RPPA protein and normal-tissue data were fetched per cancer from cBioPortal and the GDC (see Data sources). The per-cancer `dependency_scan.csv` and `overexpression.csv` files (17,916 rows each) are also not staged here — they are the intermediate outputs of the per-cancer scans, and the derived `data/*-candidates.csv` tables in this repository are the lists that passed both bars.

### Running the code

1. Place the downloaded DepMap files where the scripts expect them (see the `BASE` / `TNBC` paths at the top of the `manuscript/*.py` scripts — **these contain absolute paths from the original run and must be adjusted to your machine**; see the note below).
2. Follow `code/pipeline-playbook.md` for one cancer: scan, fetch the TCGA/protein legs, merge and score, take the top 20.
3. `code/census_query_truecpm.py <gene_list_csv> <out_prefix>` runs the Tabula Sapiens query with **true per-cell-total normalisation** (see the important note below). `code/analyze_tabula_truecpm.py` analyses the output.
4. `code/build_figs_tables.py` rebuilds the manuscript figures and tables from the `data/` tables; `code/verify_headlines.py` (110 checks) and `code/verify_source_data.py` (73 checks) recompute the headline numbers from the raw files.

### Important notes on the scripts (read before running)

- **Absolute paths.** The manuscript scripts (`build_figs_tables.py`, `write_source_data.py`, `verify_source_data.py`, `verify_headlines.py`) hard-code `BASE` and `TNBC` to the original workspace (`~/workspace/research_notes/multi-cancer-selective-targets-20261120` and `~/workspace/research_notes/depmap-tnbc-selective-targets-20261120`). Edit those two lines to point at your own checkout. `code/verify_headlines_truecpm.py` hard-codes `TAB` to the original Tabula folder — one line to edit. `code/census_query.py`, `code/census_query_truecpm.py` and `code/analyze_tabula_truecpm.py` take paths as arguments or from their own location and need no editing.
- **Verification scripts read inputs not staged here.** `verify_headlines.py` also reads `mesothelioma/tead_axis.csv`, `drug-target-calibration/target_cancer_matrix.csv` and per-cancer `dependency_scan.csv` files; `write_source_data.py` reads `drug-target-calibration/target_summary.csv` and `mesothelioma/tead_axis.csv`; `verify_source_data.py` reads `tabula-sapiens/per_gene_summary_truecpm.csv`. Those output files are not in this repository (only the leg scripts are), so the verification scripts cannot run end-to-end from this checkout alone. See `MANIFEST.md` for the full list.
- **Tabula normalisation — use the `truecpm` script.** `code/census_query.py` is retained for audit only: it normalised each cell by the sum of the *queried genes only*, so its "mean CPM10k" is counts per 10k query-set counts, not true CPM — absolute values are inflated ~45× in the 290-gene run (~200–600× in a 30-gene run). `code/census_query_truecpm.py` fixes this by normalising by total raw counts (Census obs `raw_sum`, validated exact against full-matrix row sums). Fraction-expressing values are unaffected by the bug. All corrected figures in the manuscript use the `truecpm` run (median worst-normal 2.2 CPM10k, not 107.5).
- **Per-cancer file schemas differ.** The first five scans (EGFR-mutant lung, HCC, melanoma, colon, glioma) and the TNBC scan predate the standardised schema; their candidate tables are included under the original filenames recorded in `MANIFEST.md`, and column names vary between cancers. The `*-top20.csv` files are the final ranked lists in all cases.

## Leg and per-cancer scripts

- `code/calibration/` — the drug-target calibration leg: `calibrate.py` puts 20 approved or clinical drug targets through the identical dependency and over-expression bars; `conditional_checks.py` and `build_summary.py` build the summaries; `verify_headlines.py` checks the calibration numbers.
- `code/partner-mining/` — `codepend.py` mines TAZ and KLF5 codependency partners genome-wide (DepMap correlations plus STRING interaction data); `verify_headlines.py` checks the partner numbers.
- `code/adc-intersection/` — `intersect.py` intersects the 269 selective genes with the SURFY surfaceome and applies the calibrated ADC bar; `verify_headlines.py` checks the intersection numbers.
- `code/per-cancer/<cancer>/` — the pipeline and analysis scripts each cancer's worker agent wrote while running the playbook in that cancer (scan, TCGA/protein fetch, merge and score, plus cancer-specific checks and a `verify_headlines.py` / `verify.py`). These are **heterogeneous worker scripts, not a single shared codebase**: each worker implemented the same fixed playbook (`code/pipeline-playbook.md`) in its own way, so file names, step order and column names vary between cancers — the playbook is the canonical pipeline, and these scripts are the record of how each worker ran it. Like the other scripts, they contain absolute workspace paths from the original run and need adjustment to run elsewhere. Two cancers have no scripts to stage: glioma's programme folder holds data and reports only (its scan predates the programme folder structure), and TNBC's scan folder likewise holds data and reports only — their derived tables are in `data/` all the same.

## Repository layout

```
.
├── README.md
├── LICENSE
├── MANIFEST.md
├── code/
│   ├── pipeline-playbook.md          # the fixed playbook every worker ran
│   ├── build_figs_tables.py          # rebuilds manuscript figures + tables
│   ├── write_source_data.py          # writes per-figure Source Data files
│   ├── verify_source_data.py         # 73 checks on the Source Data files
│   ├── verify_headlines.py           # 110 checks on the manuscript headlines
│   ├── census_query.py               # Tabula query (AUDIT ONLY — see note above)
│   ├── census_query_truecpm.py       # Tabula query, true-CPM normalisation
│   ├── analyze_tabula_truecpm.py     # analysis of the true-CPM query output
│   ├── verify_headlines_truecpm.py   # 62 checks on the corrected Tabula numbers
│   ├── calibration/                  # drug-target calibration leg (4 scripts)
│   ├── partner-mining/               # TAZ/KLF5 partner mining (2 scripts)
│   ├── adc-intersection/             # ADC × DepMap intersection (2 scripts)
│   └── per-cancer/                   # per-cancer pipeline scripts (139 files, 19 cancers)
└── data/
    ├── <cancer>-top20.csv            # 20 files: the ranked top-20 list per cancer
    ├── <cancer>-candidates.csv       # 20 files: genes passing both bars per cancer
    ├── melanoma-candidates-strict.csv
    └── aml-top20-blood-only-comparator.csv
```

## Data sources

| Source | Accession / version | Used for |
|---|---|---|
| DepMap | Public **24Q4** — `CRISPRGeneEffect.csv`, `Model.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, `OmicsSomaticMutations.csv` (depmap.org/portal) | Dependency and expression legs, all 20 cancers |
| TCGA | Via cBioPortal study endpoints and GDC legacy clinical / STAR counts | Tumour-versus-normal RNA leg |
| CPTAC / RPPA | Via cBioPortal, where profiles existed | Protein leg |
| Tabula Sapiens 2.0 | CZ CELLxGENE Census release **2025-11-17**, collection `e5f58829-1a66-40b5-a624-9046778e74f5`, dataset "Tabula Sapiens – All Cells" `53d208b0-2cfd-4366-9866-c3c6114081bc`; 1,136,218 cells, 180 cell types, 28 organs, 24 donors | Normal-tissue leg |
| SURFY | SURFY-2886 surfaceome | ADC surfaceome intersection |
| STRING | v12.5 | Protein–protein interaction in partner mining |

## Citation

If you use this code or data, please cite the manuscript once published. Until then:

> Hong W. *A human-guided multi-agent programme for systematic cancer selective-dependency analysis across twenty cancers.* Manuscript in preparation for *Nature Computational Science* (sole author). Code and derived data: this repository.

## Licence

MIT — see `LICENSE`. Copyright (c) 2026 WanJin Hong.
