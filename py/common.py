"""Shared config, IO and scoring for the OB-lineage Python steps (prompts/ob_lineage_plan_v1.md).

Inputs are the R exports in results/lineage/export/<ds>/ (counts.mtx genes x cells, genes.csv, cells.csv).
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import yaml

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 321


def load_config():
    with open(os.path.join(PROJ, "config", "lineage.yml")) as f:
        cfg = yaml.safe_load(f)
    cfg["paths"]["results"] = os.path.join(PROJ, cfg["paths"]["results"])
    return cfg


def out_dir(cfg, *parts):
    d = os.path.join(cfg["paths"]["results"], *parts)
    os.makedirs(d, exist_ok=True)
    return d


def log(msg):
    print(msg, flush=True)


def lognorm(a):
    """counts in .layers['counts'] -> log1p(1e4-normalized) in .X (same as scrna_matching common.load_cells)."""
    X = sp.csr_matrix(a.layers["counts"], dtype=np.float32)
    lib = np.asarray(X.sum(1)).ravel()
    X = sp.diags(1e4 / lib) @ X
    X.data = np.log1p(X.data)
    a.X = X.tocsr()
    a.obs["n_umi"] = lib
    return a


def load_export(cfg, ds):
    """One R export -> AnnData (cells x genes) with .layers['counts'] and log-normalized .X."""
    d = os.path.join(cfg["paths"]["results"], "export", ds)
    X = scipy.io.mmread(os.path.join(d, "counts.mtx")).T.tocsr().astype(np.float32)
    genes = pd.read_csv(os.path.join(d, "genes.csv"))["gene"].astype(str).values
    obs = pd.read_csv(os.path.join(d, "cells.csv"), index_col=0, low_memory=False)
    obs.index = obs.index.astype(str)
    for c in ["timepoint", "embryo", "gene_target", "pool_cluster", "pool_label", "fabian_cluster"]:
        if c in obs:
            obs[c] = obs[c].astype(str)
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    a.layers["counts"] = a.X.copy()
    return lognorm(a)


def load_dataset(cfg, ds):
    """'zscape' or 'fabian' (= Fabian 3 + 5 dpf merged within Fabian; Gate 1 decision)."""
    if ds == "zscape":
        return load_export(cfg, "zscape")
    if ds == "fabian":
        parts = [load_export(cfg, p) for p in cfg["trajectory"]["fabian_merge"]]
        # Genes absent from one export were dropped there for < 3 expressing cells: fill with 0.
        a = ad.concat(parts, join="outer", fill_value=0, merge="same")
        a.layers["counts"] = sp.csr_matrix(a.layers["counts"])
        return lognorm(a)
    raise ValueError(ds)


def modules(cfg, groups=("lineage", "contaminant", "positional")):
    out = {}
    for g in groups:
        out.update(cfg["modules"][g])
    return out


def score_modules(a, cfg, groups=("lineage", "contaminant", "positional"), tag=""):
    """scanpy score_genes per module into obs['mod_<m>'] and z across cells into obs['modz_<m>'].
    Missing genes are logged; modules with no genes present are skipped (logged)."""
    missing = []
    for m, genes in modules(cfg, groups).items():
        g = [x for x in genes if x in a.var_names]
        missing += [(m, x) for x in genes if x not in a.var_names]
        if not g:
            log(f"[{tag}] module {m}: no genes present, skipped")
            continue
        sc.tl.score_genes(a, g, score_name=f"mod_{m}", random_state=SEED, use_raw=False)
        s = a.obs[f"mod_{m}"].values
        a.obs[f"modz_{m}"] = (s - s.mean()) / (s.std() if s.std() > 0 else 1)
    if missing:
        log(f"[{tag}] module genes missing: " + ", ".join(f"{m}:{g}" for m, g in missing))
    return pd.DataFrame(missing, columns=["module", "gene"])


def knn_smooth(a, X, k=15, use_rep="X_pca"):
    """Mean over each cell's k nearest neighbours (itself included) in `use_rep`. X: cells x features."""
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k).fit(a.obsm[use_rep])
    idx = nn.kneighbors(return_distance=False)
    idx = np.hstack([np.arange(a.n_obs)[:, None], idx[:, : k - 1]])
    rows = np.repeat(np.arange(a.n_obs), k)
    W = sp.csr_matrix((np.full(rows.size, 1.0 / k), (rows, idx.ravel())), shape=(a.n_obs, a.n_obs))
    return W @ X


def detected(a, genes):
    """Boolean cells x genes: >= 1 UMI."""
    g = [x for x in genes if x in a.var_names]
    return np.asarray((a[:, g].layers["counts"] > 0).toarray()), g
