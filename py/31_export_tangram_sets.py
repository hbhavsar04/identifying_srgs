"""Part 1 -> Part 2 handoff (prompts/ob_lineage_plan_v1.md): Tangram reference sets as raw-count h5ad.

Format follows scrna_matching/iter0_data_prep (raw counts in .X, cells x genes; obs with timepoint,
embryo, gene_target, cell_type_sub/broad, nCount_RNA, nFeature_RNA), plus the lineage annotations:
state, state_fine, ob_annot, pool_label, dpt_pseudotime and the fate_* probabilities.

Sets (written to <scrna_matching SHARE_OUT>/iter4_lineage_reference/reference/<set>.h5ad):
  lineage_z72          ZSCAPE 72 hpf refined pool
  lineage_z48_96       ZSCAPE 48-96 hpf refined pool
  lineage_f3           Fabian 3 dpf refined pool
  lineage_z72_sub      all annotated OBs + up to 150 cells per state_fine (others), fixed sample
  lineage_z48_96_sub   same rule on the 48-96 hpf pool (for arm A4)
  lineage_f3_sub       same rule on Fabian 3 dpf (no OBs; for arm A5)

Every set carries all panel genes and the final new HCR genes: genes dropped by the R export
(< 3 expressing cells in that pool) are added as all-zero columns and logged.
Checks: panel detection equals the R table (export/panel_detection_R.csv); the 186 ZSCAPE 72 hpf
OBs are exactly the Iteration 1 zscape_ob72 cells.

Run: uv run --project env python py/31_export_tangram_sets.py
"""
import os
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import common as C

SUB_PER_STATE = 150
OBS_KEEP = ["dataset", "timepoint", "embryo", "gene_target", "pool_label", "pool_cluster", "ob_annot",
            "state", "state_fine", "dpt_pseudotime", "fate_max", "flag_muscle_ambient"]


def share_dir(cfg):
    """SHARE_OUT from scrna_matching/paths.py (read-only), resolved in the scrna_matching env because
    paths.py imports opercle_subclustering modules whose deps (tifffile) are not in this env."""
    sm = cfg["paths"]["scrna_matching_repo"]
    share = subprocess.run(["uv", "run", "--project", os.path.join(sm, "env"), "python", "-c",
                            "import paths; print(paths.SHARE_OUT)"], cwd=sm, capture_output=True,
                           text=True, check=True).stdout.strip().splitlines()[-1]
    d = os.path.join(share, "iter4_lineage_reference", "reference")
    os.makedirs(d, exist_ok=True)
    return d, share


def to_reference(a, must_have):
    """states AnnData -> iter0-format AnnData (raw counts in X); adds missing must-have genes as zeros."""
    X = sp.csr_matrix(a.layers["counts"], dtype=np.float32)
    obs = a.obs[[c for c in OBS_KEEP if c in a.obs] + [c for c in a.obs if c.startswith("fate_") and c != "fate_max"]].copy()
    obs["cell_type_sub"] = a.obs["pool_label"].astype(str)
    obs["cell_type_broad"] = np.where(a.obs["ob_annot"].astype(bool), "osteoblast", "non-OB")
    obs["nCount_RNA"] = np.asarray(X.sum(1)).ravel()
    obs["nFeature_RNA"] = np.asarray((X > 0).sum(1)).ravel()
    for c in obs:
        if obs[c].dtype == object or str(obs[c].dtype) == "category":
            obs[c] = obs[c].astype(str)
    genes = list(a.var_names)
    missing = [g for g in must_have if g not in set(genes)]
    if missing:
        X = sp.hstack([X, sp.csr_matrix((X.shape[0], len(missing)), dtype=np.float32)]).tocsr()
        genes += missing
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes)), missing


def subsample(a, rng):
    """All annotated OBs + up to SUB_PER_STATE cells per state_fine among the rest."""
    ob = a.obs["ob_annot"].astype(bool).values
    keep = list(a.obs_names[ob])
    rest = a.obs[~ob]
    for s, idx in rest.groupby("state_fine", observed=True).groups.items():
        idx = list(idx)
        keep += list(rng.choice(idx, size=min(SUB_PER_STATE, len(idx)), replace=False))
    return a[a.obs_names.isin(keep)].copy()


def main():
    cfg = C.load_config()
    out, P = share_dir(cfg)
    res = cfg["paths"]["results"]
    panel = cfg["panel_genes"] + [g for g in cfg.get("panel_validation_only", []) if g not in cfg["panel_genes"]]
    final_new = list(pd.read_csv(os.path.join(res, "panel", "hcr_panel_v1.csv"))["gene"])
    must = list(dict.fromkeys(panel + final_new))
    r_det = pd.read_csv(os.path.join(res, "export", "panel_detection_R.csv"))
    z = ad.read_h5ad(os.path.join(res, "zscape_states.h5ad"))
    f = ad.read_h5ad(os.path.join(res, "fabian_states.h5ad"))
    tp = z.obs["timepoint"].astype(str)
    full = {"lineage_z72": (z[tp == "72"].copy(), "zscape"),
            "lineage_z48_96": (z.copy(), "zscape"),
            "lineage_f3": (f[f.obs["timepoint"].astype(str) == "3dpf"].copy(), "fabian3")}
    rng = np.random.default_rng(C.SEED)
    sets = dict(full)
    for k in ("lineage_z72", "lineage_z48_96", "lineage_f3"):
        sets[k + "_sub"] = (subsample(full[k][0], rng), full[k][1])

    rows = []
    for name, (a, r_ds) in sets.items():
        ref, missing = to_reference(a, must)
        # detection check vs R (full pools only; the R table is per export dataset)
        if not name.endswith("_sub") and name != "lineage_z72":
            rd = r_det[r_det.dataset == r_ds].set_index("gene")["pct_detected"]
            g = [x for x in rd.index if x in ref.var_names]
            py = pd.Series(100 * np.asarray((ref[:, g].X > 0).mean(0)).ravel(), index=g)
            diff = (py - rd.loc[g]).abs().max()
            assert diff < 0.1, f"{name}: panel detection differs from R by {diff:.3f} points"
        if name == "lineage_z72":
            iter0 = ad.read_h5ad(os.path.join(P, "iter0_data_prep", "reference", "zscape_ob72.h5ad"),
                                 backed="r")
            obs_ob = set(ref.obs_names[ref.obs["ob_annot"].astype(str) == "True"])
            assert obs_ob == set(iter0.obs_names), "lineage_z72 OBs differ from iter0 zscape_ob72"
        assert ref.obs["state"].notna().all() and (ref.obs["state"] != "nan").all(), f"{name}: NA state"
        ref.write_h5ad(os.path.join(out, f"{name}.h5ad"), compression="gzip")
        n_ob = int((ref.obs["ob_annot"].astype(str) == "True").sum())
        rows.append({"set": name, "n_cells": ref.n_obs, "n_genes": ref.n_vars, "n_ob": n_ob,
                     "n_state_fine": ref.obs["state_fine"].nunique(),
                     "timepoints": ";".join(f"{t}:{n}" for t, n in ref.obs["timepoint"].value_counts().sort_index().items()),
                     "zero_filled_genes": ";".join(missing), "path": os.path.join(out, f"{name}.h5ad")})
        C.log(f"{name}: {ref.n_obs} cells x {ref.n_vars} genes, {n_ob} OBs; zero-filled {missing}")
        ref.obs.groupby(["state", "state_fine"], observed=True).size().rename("n_cells").to_csv(
            os.path.join(C.out_dir(cfg, "tangram_sets"), f"{name}_composition.csv"))
    pd.DataFrame(rows).to_csv(os.path.join(C.out_dir(cfg, "tangram_sets"), "manifest.csv"), index=False)


if __name__ == "__main__":
    main()
