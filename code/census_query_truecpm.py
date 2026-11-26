"""Query CZ CELLxGENE Census for Tabula Sapiens v2 expression of a gene list,
with TRUE CPM10k normalisation (denominator = per-cell total raw counts).

CORRECTION (26 Nov 2026): the original census_query.py normalised each cell
by the sum of the QUERIED genes only (X was already var-filtered to the gene
list, so X.sum(axis=1) is the query-set total, not the cell total). Its
"mean CPM10k" was counts per 10k query-set counts — inflated ~45x in the
290-gene run and ~200-600x in the 30-gene HCC run. This script instead takes
the denominator from the Census obs column `raw_sum` (total raw counts per
cell), validated 26 Nov 2026 against full-X row sums on retinal pigment
epithelial cells (exact match). Fraction-expressing values are unaffected
(they only use raw count > 0) and were identical between the two runs.

Same per-organ route as the original (28 organ datasets partition All Cells;
per-dataset aggregates saved as they land, resumable). Writes to a NEW
output prefix — never overwrites the original tabula_query/ or
tabula_background/ files.

Usage:
  python -u census_query_truecpm.py <gene_list_csv> <out_prefix>
"""
import gc
import json
import os
import sys
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import scipy.sparse as sp

import cellxgene_census
import tiledbsoma

CENSUS_VERSION = "2025-11-17"
TS_V2_COLLECTION_ID = "e5f58829-1a66-40b5-a624-9046778e74f5"
TS_V2_ALL_CELLS_ID = "53d208b0-2cfd-4366-9866-c3c6114081bc"
ORGAN_DATASETS = [
    ("d45aca3a-ddb8-4185-a462-eff1e0c7af2e", "Ear"),
    ("2423ce2c-3149-4cca-a2ff-cf682ea29b5f", "Kidney"),
    ("0041b9c3-6a49-4bf7-8514-9bc7190067a7", "Skin"),
    ("ff45e623-7f5f-46e3-b47d-56be0341f66b", "Pancreas"),
    ("6d41668c-168c-4500-b06a-4674ccf3e19d", "Liver"),
    ("6ec405bb-4727-4c6d-ab4e-01fe489af7ea", "Uterus"),
    ("d77ec7d6-ef2e-49d6-9e79-05b7f8881484", "Prostate"),
    ("d8732da6-8d1d-42d9-b625-f2416c30054b", "Trachea"),
    ("7357cee7-9f7f-4ab0-8cec-90de8f047e38", "Large_Intestine"),
    ("4f1555bc-4664-46c3-a606-78d34dd10d92", "Bone_Marrow"),
    ("9ba03780-4b13-44bc-a7d3-ce532ea0a856", "Stomach"),
    ("e6a11140-2545-46bc-929e-da243eed2cae", "Heart"),
    ("2ba40233-8576-4dec-a5f1-2adfa115e2dc", "Mammary"),
    ("a0754256-f44b-4c4a-962c-a552e47d3fdc", "Eye"),
    ("0ced5e76-6040-47ff-8a72-93847965afc0", "Thymus"),
    ("f01bdd17-4902-40f5-86e3-240d66dd2587", "Salivary_Gland"),
    ("a357414d-2042-4eb5-95f0-c58604a18bdd", "Small_Intestine"),
    ("bb78adbb-1428-4468-bad0-63602bb974e9", "Ovary"),
    ("55cf0ea3-9d2b-4294-871e-bb4b49a79fc7", "Tongue"),
    ("1c9eb291-6d31-47e1-96b2-129b5e1ae64f", "Muscle"),
    ("a2d4d33e-4c62-4361-b80a-9be53d2e50e8", "Vasculature"),
    ("983d5ec9-40e8-4512-9e65-a572a9c486cb", "Blood"),
    ("cee11228-9f0b-4e57-afe2-cfe15ee56312", "Spleen"),
    ("0d2ee4ac-05ee-40b2-afb6-ebb584caa867", "Lung"),
    ("e5c63d94-593c-4338-a489-e1048599e751", "Bladder"),
    ("5e5e7a2f-8f1c-42ac-90dc-b4f80f38e84c", "Fat"),
    ("18eb630b-a754-4111-8cd4-c24ec80aa5ec", "Lymph_Node"),
    ("ec6ef004-5f30-4a40-b282-980cca20a561", "Testis"),
]


def proxy_context():
    px = os.environ.get("https_proxy") or os.environ.get("HTTPS_PROXY") or ""
    u = urlparse(px)
    cfg = {"vfs.s3.region": "us-west-2",
           "vfs.s3.connect_timeout_ms": "120000",
           "vfs.s3.request_timeout_ms": "600000",
           "vfs.s3.max_attempts": "5"}
    if u.hostname:
        cfg["vfs.s3.proxy_scheme"] = u.scheme or "http"
        cfg["vfs.s3.proxy_host"] = u.hostname
        cfg["vfs.s3.proxy_port"] = str(u.port or 3128)
        if u.username:
            cfg["vfs.s3.proxy_username"] = u.username
        if u.password:
            cfg["vfs.s3.proxy_password"] = u.password
    return tiledbsoma.SOMATileDBContext(tiledb_config=cfg)


def main():
    gene_csv, out_prefix = sys.argv[1], sys.argv[2]
    genes = pd.read_csv(gene_csv)
    gene_col = "gene" if "gene" in genes.columns else genes.columns[0]
    genes = [g.strip() for g in genes[gene_col].astype(str)]
    print(f"{len(genes)} query genes", flush=True)

    out = Path(out_prefix)
    per_ds = out / "per_dataset"
    per_ds.mkdir(parents=True, exist_ok=True)

    organ_filter = os.environ.get("TABULA_ORGANS")
    datasets = ORGAN_DATASETS
    if organ_filter:
        want = set(organ_filter.split(","))
        datasets = [(d, o) for d, o in ORGAN_DATASETS if o in want]
        print(f"Organ subset ({len(datasets)}): {[o for _, o in datasets]}", flush=True)

    census = cellxgene_census.open_soma(
        census_version=CENSUS_VERSION, context=proxy_context())
    print(f"Census URI: {census.uri}", flush=True)

    info = census["census_info"]["datasets"].read().concat().to_pandas()
    v2 = info[info["collection_id"] == TS_V2_COLLECTION_ID].set_index("dataset_id")
    allcells_n = int(v2.loc[TS_V2_ALL_CELLS_ID, "dataset_total_cell_count"])
    organ_sum = 0
    for ds_id, organ in ORGAN_DATASETS:
        row = v2.loc[ds_id]
        assert "Tabula Sapiens" in row["dataset_title"], (ds_id, row["dataset_title"])
        organ_sum += int(row["dataset_total_cell_count"])
    print(f"Organ datasets sum: {organ_sum}; All Cells: {allcells_n}; "
          f"match: {organ_sum == allcells_n}", flush=True)
    assert organ_sum == allcells_n

    var_filter = "feature_name in [" + ",".join(f"'{g}'" for g in genes) + "]"
    found_genes = None
    missing = None
    total_cells = 0
    donors_all = set()
    for ds_id, organ in datasets:
        done = per_ds / f"{organ}.npz"
        if done.exists():
            print(f"  {organ}: already done, skipping", flush=True)
            continue
        adata = cellxgene_census.get_anndata(
            census, organism="Homo sapiens", X_name="raw",
            obs_value_filter=f"dataset_id == '{ds_id}'",
            var_value_filter=var_filter,
            obs_column_names=["cell_type", "cell_type_ontology_term_id",
                              "tissue", "donor_id", "raw_sum"],
            var_column_names=["feature_name"])
        found = list(adata.var["feature_name"])
        if found_genes is None:
            found_genes = found
            missing = sorted(set(genes) - set(found))
            print(f"Genes found: {len(found_genes)}; missing: {missing}", flush=True)
        X = adata.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        X = X.tocsr().astype(np.float32)
        # TRUE CPM10k: denominator is the per-cell TOTAL raw count (obs
        # raw_sum), NOT the sum over the queried genes. The original script
        # used X.sum(axis=1), which is the query-set total only.
        totals = adata.obs["raw_sum"].values.astype(np.float64)
        assert len(totals) == X.shape[0], (organ, len(totals), X.shape)
        n_zero = int((totals <= 0).sum())
        if n_zero:
            print(f"  {organ}: WARNING {n_zero} cells with raw_sum<=0; "
                  f"denominator set to 1 for those cells", flush=True)
            totals = totals.copy()
            totals[totals <= 0] = 1.0
        Xn = X.multiply((1e4 / totals)[:, None]).tocsr()
        ct = adata.obs["cell_type"].astype(str).values
        uniq, inv = np.unique(ct, return_inverse=True)
        onehot = sp.csr_matrix(
            (np.ones(len(inv), dtype=np.float32), (np.arange(len(inv)), inv)),
            shape=(len(inv), len(uniq)))
        sum_here = (onehot.T @ Xn).toarray()
        pos_here = (onehot.T @ (X > 0).astype(np.float32)).toarray()
        n_here = np.asarray(onehot.sum(axis=0)).ravel()
        donors = adata.obs["donor_id"].astype(str).values
        donor_sets = [";".join(sorted(set(donors[inv == i]))) for i in range(len(uniq))]
        tissues = [";".join(sorted(set(adata.obs["tissue"].astype(str).values[inv == i])))
                 for i in range(len(uniq))]
        onto = [str(adata.obs["cell_type_ontology_term_id"].values[inv == i][0])
                for i in range(len(uniq))]
        np.savez_compressed(
            done, cell_types=np.array(uniq), genes=np.array(found),
            sum_cpm=sum_here, n_pos=pos_here, n_cells=n_here,
            donor_sets=np.array(donor_sets), tissues=np.array(tissues),
            ontology=np.array(onto))
        print(f"  {organ}: {len(ct)} cells, {len(uniq)} cell types, "
              f"donors {len(set(donors))}, median raw_sum "
              f"{float(np.median(totals)):.0f}", flush=True)
        total_cells += len(ct)
        donors_all |= set(donors)
        del adata, X, Xn, onehot
        gc.collect()

    present = [o for _, o in ORGAN_DATASETS if (per_ds / f"{o}.npz").exists()]
    if len(present) != len(ORGAN_DATASETS):
        print(f"PARTIAL: {len(present)}/{len(ORGAN_DATASETS)} organ files present; "
              f"missing: {sorted(set(o for _, o in ORGAN_DATASETS) - set(present))}. "
              f"Re-run for the missing organs, then combine.", flush=True)
        return
    ct_sum, ct_pos, ct_n, ct_donors, ct_tissues, ct_onto = {}, {}, {}, {}, {}, {}
    t_sum, t_pos, t_n = {}, {}, {}
    genes_found = None
    for ds_id, organ in ORGAN_DATASETS:
        z = np.load(per_ds / f"{organ}.npz", allow_pickle=True)
        if genes_found is None:
            genes_found = list(z["genes"])
        assert list(z["genes"]) == genes_found, organ
        for i, c in enumerate(z["cell_types"]):
            c = str(c)
            if c not in ct_sum:
                ct_sum[c] = z["sum_cpm"][i].copy()
                ct_pos[c] = z["n_pos"][i].copy()
                ct_n[c] = float(z["n_cells"][i])
                ct_donors[c] = set(str(z["donor_sets"][i]).split(";")) if str(z["donor_sets"][i]) else set()
                ct_tissues[c] = set(str(z["tissues"][i]).split(";")) if str(z["tissues"][i]) else set()
                ct_onto[c] = str(z["ontology"][i])
            else:
                ct_sum[c] += z["sum_cpm"][i]
                ct_pos[c] += z["n_pos"][i]
                ct_n[c] += float(z["n_cells"][i])
                ct_donors[c] |= set(str(z["donor_sets"][i]).split(";"))
                ct_tissues[c] |= set(str(z["tissues"][i]).split(";"))
        t_sum[organ] = z["sum_cpm"].sum(axis=0)
        t_pos[organ] = z["n_pos"].sum(axis=0)
        t_n[organ] = float(z["n_cells"].sum())
    total_cells = int(sum(t_n.values()))
    donors_all = set().union(*ct_donors.values())
    if missing is None:
        missing = sorted(set(genes) - set(genes_found))

    ct_names = sorted(ct_sum.keys())
    mean_mat = np.vstack([ct_sum[c] / ct_n[c] for c in ct_names])
    frac_mat = np.vstack([ct_pos[c] / ct_n[c] for c in ct_names])
    tg_names = sorted(t_sum.keys())
    mean_t = np.vstack([t_sum[t] / t_n[t] for t in tg_names])
    frac_t = np.vstack([t_pos[t] / t_n[t] for t in tg_names])

    pd.DataFrame(mean_mat, index=ct_names, columns=genes_found).to_csv(
        out / "celltype_mean_cpm10k.csv")
    pd.DataFrame(frac_mat, index=ct_names, columns=genes_found).to_csv(
        out / "celltype_frac_expressing.csv")
    pd.DataFrame(mean_t, index=tg_names, columns=genes_found).to_csv(
        out / "tissue_mean_cpm10k.csv")
    pd.DataFrame(frac_t, index=tg_names, columns=genes_found).to_csv(
        out / "tissue_frac_expressing.csv")
    ct_meta = pd.DataFrame([
        {"cell_type": c, "n_cells": int(ct_n[c]), "n_donors": len(ct_donors[c]),
         "tissues": ";".join(sorted(ct_tissues[c])), "ontology_term": ct_onto[c]}
        for c in ct_names])
    ct_meta.to_csv(out / "celltype_metadata.csv", index=False)

    meta = {
        "census_version": CENSUS_VERSION,
        "census_uri": str(census.uri),
        "tabula_version_used": "v2",
        "tabula_collection_id": TS_V2_COLLECTION_ID,
        "tabula_dataset_id": TS_V2_ALL_CELLS_ID,
        "tabula_dataset_title": "Tabula Sapiens - All Cells",
        "tabula_source": "Tabula Sapiens 2.0, Cell 2026, DOI 10.1016/j.cell.2026.08.010",
        "n_cells": total_cells,
        "n_donors": len(donors_all),
        "n_cell_types": int(len(ct_names)),
        "n_organs_queried": len(ORGAN_DATASETS),
        "organ_cell_sum_matches_all_cells": organ_sum == allcells_n,
        "n_query_genes": len(genes),
        "n_genes_found": len(genes_found),
        "genes_missing": missing,
        "chunking": "one get_anndata call per organ dataset (28); full-slice query OOM-killed at 6 GB, tissue_general value-filter slice too slow",
        "normalisation": "TRUE CPM10k: counts per 10k TOTAL raw counts per cell (Census obs raw_sum; validated 26 Nov 2026 against full-X row sums on RPE cells). Corrects census_query.py, which normalised by the queried-genes-only sum (~45x inflation in the 290-gene run).",
    }
    with open(out / "query_metadata.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
