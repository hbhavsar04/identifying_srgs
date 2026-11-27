"""Phase 2: states and OB-lineage trajectory (prompts/ob_lineage_plan_v1.md).

Same code for both datasets, never integrated:
  zscape  ZSCAPE refined pool, 48-96 hpf (k = 30)
  fabian  Fabian 3 dpf + 5 dpf refined pools merged within Fabian (k = 15; Gate 1 decision)

Root = cell with highest kNN-smoothed (nc_mesenchyme - ob_mature) z among the earliest timepoint -> DPT.
CellRank 2: primary 0.8 Pseudotime + 0.2 Connectivity; sensitivity 0.8 CytoTRACE + 0.2 Connectivity
(run with the primary's terminal states, so only the kernel differs). GPCCA n_states chosen by the
largest eigengap within 4..12. No RNA velocity (no spliced/unspliced layers).
Proliferation-named macrostates ('cycling') are never used as terminal states: cycling is not a fate
(first Fabian run called all 8 macrostates terminal, including 2 cycling ones).

Run: uv run --project env python py/20_trajectory.py --dataset zscape|fabian
     uv run --project env python py/20_trajectory.py --concordance   (after both)
Outputs: results/lineage/<ds>_states.h5ad and results/lineage/trajectory/<ds>_*.csv/png;
         results/lineage/trajectory/state_concordance.csv once both datasets are done.
"""
import argparse
import os
import time
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

import common as C  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

K = {"zscape": 30, "fabian": 15}
EARLIEST = {"zscape": "48", "fabian": "3dpf"}
NAME_OF = {"ob_mature": "OB", "ob_early": "OB_early", "chondro": "chondro", "ligament_tendon": "ligament",
           "fibro": "fibro", "nc_mesenchyme": "mesenchyme", "proliferation": "cycling"}
LINEAGE_MODS = list(NAME_OF)
THRESHOLDS = (0.3, 0.4, 0.5)
MAIN_THR = 0.4
MIN_FINE = 15
N_TOP_OB = 30          # cells used if the OB terminal state has to be set by hand
OB_STATES = ("OB_mature", "OB_committed", "OB_progenitor")


# ---- preprocessing --------------------------------------------------------------------
def preprocess(a, ds, cfg):
    C.score_modules(a, cfg, tag=ds)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    sc.pp.pca(a, n_comps=50, mask_var="highly_variable", random_state=C.SEED)
    sc.pp.neighbors(a, n_neighbors=K[ds], n_pcs=50, random_state=C.SEED)
    sc.tl.leiden(a, resolution=1.0, random_state=C.SEED, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(a, random_state=C.SEED)
    sc.tl.diffmap(a, n_comps=15)
    return a


def pick_root(a, ds):
    diff = (a.obs["modz_nc_mesenchyme"] - a.obs["modz_ob_mature"]).values[:, None]
    sm = np.asarray(C.knn_smooth(a, diff, k=15)).ravel()
    early = np.where(a.obs["timepoint"].values == EARLIEST[ds])[0]
    root = early[np.argmax(sm[early])]
    a.uns["iroot"] = int(root)
    sc.tl.dpt(a, n_dcs=10)
    pt = a.obs["dpt_pseudotime"].values
    n_inf = int((~np.isfinite(pt)).sum())
    if n_inf:
        C.log(f"[{ds}] {n_inf} cells unreachable from root (inf pseudotime) -> set to 1.0")
        pt[~np.isfinite(pt)] = 1.0
        a.obs["dpt_pseudotime"] = pt
    return root, n_inf


# ---- CellRank ---------------------------------------------------------------------------
def kernel(a, which):
    from cellrank.kernels import ConnectivityKernel, CytoTRACEKernel, PseudotimeKernel
    ck = ConnectivityKernel(a).compute_transition_matrix()
    if which == "pseudotime":
        k = PseudotimeKernel(a, time_key="dpt_pseudotime").compute_transition_matrix(
            show_progress_bar=False, n_jobs=8)
    else:
        # CytoTRACE needs imputed expression ('Ms' = scVelo moments); kNN mean in PCA space instead.
        a.layers["Ms"] = sp.csr_matrix(C.knn_smooth(a, a.X, k=K[a.uns["ds"]], use_rep="X_pca"))
        k = CytoTRACEKernel(a).compute_cytotrace(layer="Ms").compute_transition_matrix(
            show_progress_bar=False, n_jobs=8)
    return 0.8 * k + 0.2 * ck


def choose_n_states(g, lo=4, hi=12):
    g.compute_eigendecomposition(k=hi + 2, only_evals=True)
    ev = np.sort(np.real(g.eigendecomposition["D"]))[::-1]
    gaps = pd.DataFrame({"n_states": np.arange(1, len(ev)), "eigval_n": ev[:-1], "eigval_n+1": ev[1:],
                         "gap": ev[:-1] - ev[1:]})
    ok = gaps[(gaps.n_states >= lo) & (gaps.n_states <= hi)]
    return int(ok.loc[ok.gap.idxmax(), "n_states"]), gaps


def fit_gpcca(g, n):
    """GPCCA fit; complex-conjugate splits make some n invalid, so try n, n+1, n-1."""
    for m in (n, n + 1, n - 1):
        try:
            g.compute_schur(n_components=max(20, m + 2))
            g.fit(n_states=m, n_cells=30)
            return m
        except Exception as e:  # noqa: BLE001
            C.log(f"  GPCCA n_states={m} failed: {e}")
    raise RuntimeError("GPCCA failed")


def name_states(a, series):
    """Name each category by the lineage module with the highest mean z among its cells."""
    rows, names, used = [], {}, {}
    for cat in series.cat.categories:
        cells = series.index[series == cat]
        mz = a.obs.loc[cells, [f"modz_{m}" for m in LINEAGE_MODS]].mean()
        dom = mz.idxmax()[5:]
        base = NAME_OF[dom]
        used[base] = used.get(base, 0) + 1
        nm = base if used[base] == 1 else f"{base}_{used[base]}"
        names[cat] = nm
        lab = a.obs.loc[cells, "pool_label"].value_counts()
        tp = a.obs.loc[cells, "timepoint"].value_counts()
        rows.append({"gpcca_state": cat, "name": nm, "n_cells": len(cells), "dominant_module": dom,
                     **{f"z_{m}": round(mz[f'modz_{m}'], 2) for m in LINEAGE_MODS},
                     "top_labels": "; ".join(f"{k} ({v})" for k, v in lab.head(3).items()),
                     "timepoints": "; ".join(f"{k} ({v})" for k, v in tp.items())})
    return names, pd.DataFrame(rows)


def fates_df(g, a):
    fp = g.fate_probabilities
    return pd.DataFrame(np.asarray(fp), index=a.obs_names, columns=list(fp.names))


# ---- states -----------------------------------------------------------------------------
def assign_states(fates, terminal, ob_early_z, thr):
    maxf = fates.max(1).values
    arg = fates.idxmax(1).values
    fob = fates["OB"].values
    st = np.where(maxf < thr, "uncommitted_mesenchyme", pd.Series(arg).astype(str) + "_lineage")
    st = np.where((arg == "OB") & (fob >= thr), "OB_progenitor", st)
    st = np.where((fob >= 0.5) & (ob_early_z > 0), "OB_committed", st)
    st = np.where(terminal.astype(str).values == "OB", "OB_mature", st)
    return pd.Series(st, index=fates.index)


def state_fine(a, res=0.5):
    fine = pd.Series("", index=a.obs_names, dtype=object)
    for s in sorted(a.obs["state"].unique()):
        cells = a.obs_names[a.obs["state"] == s]
        if len(cells) < 2 * MIN_FINE:
            fine[cells] = f"{s}.0"
            continue
        sc.tl.leiden(a, resolution=res, restrict_to=("state", [s]), key_added="_sub",
                     random_state=C.SEED, flavor="igraph", n_iterations=2, directed=False)
        sub = a.obs.loc[cells, "_sub"].astype(str).str.split(",").str[-1]
        # merge subclusters < MIN_FINE into the sibling with the nearest PCA centroid
        while True:
            sizes = sub.value_counts()
            small = sizes[sizes < MIN_FINE]
            if small.empty or len(sizes) == 1:
                break
            x = small.index[0]
            pcs = pd.DataFrame(a[cells].obsm["X_pca"][:, :20], index=cells)
            cent = pcs.groupby(sub.values).mean()
            d = ((cent.drop(index=x) - cent.loc[x]) ** 2).sum(1)
            sub[sub == x] = d.idxmin()
        codes = {c: i for i, c in enumerate(sub.value_counts().index)}
        fine[cells] = [f"{s}.{codes[c]}" for c in sub]
        del a.obs["_sub"]
    return fine


# ---- checks -------------------------------------------------------------------------------
def ob_order(a, ds, outd):
    ob = a.obs["state"].isin(OB_STATES).values
    sub = a[ob]
    pt = sub.obs["dpt_pseudotime"].values
    bins = pd.qcut(pd.Series(pt).rank(method="first"), 20, labels=False).values
    genes = {"runx2": ["runx2a", "runx2b"], "sp7": ["sp7"], "col10a1a/ifitm5": ["col10a1a", "ifitm5"]}
    prof = {}
    for k, gs in genes.items():
        gs = [g for g in gs if g in sub.var_names]
        x = np.asarray(sub[:, gs].X.mean(1)).ravel()
        prof[k] = pd.Series(x).groupby(bins).mean().rolling(3, center=True, min_periods=1).mean()
    prof = pd.DataFrame(prof)
    peaks = prof.idxmax()
    ok = bool(peaks["runx2"] <= peaks["sp7"] <= peaks["col10a1a/ifitm5"])
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for k in prof:
        ax.plot(prof.index, prof[k] / prof[k].max(), label=f"{k} (peak bin {peaks[k]})")
    ax.set_xlabel("OB-lineage pseudotime bin (20 quantiles)")
    ax.set_ylabel("smoothed mean log expr (max = 1)")
    ax.set_title(f"{ds}: OB-lineage marker order ({'PASS' if ok else 'FAIL'})")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(outd, f"{ds}_ob_marker_order.png"), dpi=150)
    plt.close(fig)
    return ok, peaks.to_dict(), int(ob.sum())


def positional(a, cfg, ds):
    pos = list(cfg["modules"]["positional"])
    cols = [f"mod_{m}" for m in pos if f"mod_{m}" in a.obs]
    sm = np.asarray(C.knn_smooth(a, a.obs[cols].values, k=15))
    for i, c in enumerate(cols):
        a.obs[f"smooth_{c[4:]}"] = sm[:, i]
    rows = []
    for f, idx in a.obs.groupby("state_fine", observed=True).indices.items():
        r = {"state_fine": f, "n": len(idx)}
        for m in ("hox2", "hox3", "ventral"):
            det, g = C.detected(a[idx], cfg["modules"]["positional"][m])
            r[f"det_{m}_pct"] = round(100 * det.any(1).mean(), 1) if g else np.nan
        for c in cols:
            r[f"smooth_{c[4:]}"] = a.obs[f"smooth_{c[4:]}"].values[idx].mean()
        rows.append(r)
    df = pd.DataFrame(rows)
    for c in [c for c in df if c.startswith("smooth_")]:
        df[c + "_z"] = ((df[c] - df[c].mean()) / df[c].std()).round(2)
    # Opercle = PA2-derived, dorsal: hox2-high, hox3-low, hand2 (ventral)-low. Descriptive flag only.
    df["opercle_candidate"] = (df.get("smooth_hox2_z", 0) > 1) & (df.get("smooth_hox3_z", 0) < 0.5) & \
        (df.get("smooth_ventral_z", 0) <= 0)
    return df


# ---- concordance -----------------------------------------------------------------------------
def concordance(cfg, outd):
    paths = {d: os.path.join(cfg["paths"]["results"], f"{d}_states.h5ad") for d in ("zscape", "fabian")}
    if not all(os.path.exists(p) for p in paths.values()):
        return None
    A = {d: sc.read_h5ad(p) for d, p in paths.items()}
    shared = A["zscape"].var_names.intersection(A["fabian"].var_names)
    lfc = {}
    for d, a in A.items():
        a = a[:, shared]
        X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
        det = a.layers["counts"] > 0
        st = a.obs["state"].values
        means, dets = {}, {}
        for s in np.unique(st):
            m = st == s
            means[s] = np.asarray(X[m].mean(0)).ravel()
            dets[s] = np.asarray(det[m].mean(0)).ravel()
        tot = np.asarray(X.mean(0)).ravel()
        n = len(st)
        lfc[d] = pd.DataFrame({s: means[s] - (tot * n - means[s] * (st == s).sum()) / (n - (st == s).sum())
                               for s in means}, index=shared)
        lfc[d + "_det"] = pd.DataFrame(dets, index=shared).max(1)
    keep = (lfc["zscape_det"] >= 0.05) & (lfc["fabian_det"] >= 0.05)
    z, f = lfc["zscape"][keep], lfc["fabian"][keep]
    cor = pd.DataFrame({fs: [np.corrcoef(z[zs], f[fs])[0, 1] for zs in z] for fs in f}, index=z.columns)
    cor.index.name = "zscape_state"
    cor.round(3).to_csv(os.path.join(outd, "state_concordance.csv"))
    C.log(f"state concordance over {keep.sum()} shared genes (>= 5% in a state in both):\n{cor.round(2)}")
    return cor


# ---- main ---------------------------------------------------------------------------------------
def main():
    import cellrank as cr
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["zscape", "fabian"])
    ap.add_argument("--concordance", action="store_true")
    args = ap.parse_args()
    cfg = C.load_config()
    outd = C.out_dir(cfg, "trajectory")
    if args.concordance:
        concordance(cfg, outd)
        return
    ds = args.dataset
    t0 = time.time()
    tick = lambda msg: C.log(f"[{ds}] {time.time() - t0:7.0f}s {msg}")  # noqa: E731
    np.random.seed(C.SEED)

    a = C.load_dataset(cfg, ds)
    a.uns["ds"] = ds
    C.log(f"[{ds}] {a.n_obs} cells x {a.n_vars} genes; timepoints {a.obs.timepoint.value_counts().to_dict()}")
    preprocess(a, ds, cfg)
    tick("preprocessed")
    root, n_inf = pick_root(a, ds)
    C.log(f"[{ds}] root {a.obs_names[root]} ({a.obs.pool_label.iloc[root]}, {a.obs.timepoint.iloc[root]})")

    # primary kernel
    g = cr.estimators.GPCCA(kernel(a, "pseudotime"))
    n, gaps = choose_n_states(g)
    gaps.to_csv(os.path.join(outd, f"{ds}_eigengap.csv"), index=False)
    tick(f"kernel + eigengap (n = {n})")
    n = fit_gpcca(g, n)
    tick("GPCCA fit")
    names, macro = name_states(a, g.macrostates)
    macro.to_csv(os.path.join(outd, f"{ds}_macrostates.csv"), index=False)
    g.predict_terminal_states()
    term = g.terminal_states
    term_cats = list(term.cat.categories)
    C.log(f"[{ds}] n_states {n}; macrostates {names}; terminal (predicted) {[names[c] for c in term_cats]}")
    states = {names[c]: list(term.index[term == c]) for c in term_cats if not names[c].startswith("cycling")}
    dropped = [names[c] for c in term_cats if names[c].startswith("cycling")]
    if dropped:
        C.log(f"[{ds}] not used as terminal (proliferation): {dropped}")
    ob_set_by_hand = "OB" not in states
    if ob_set_by_hand:
        taken = set().union(*states.values()) if states else set()
        top = a.obs["modz_ob_mature"].drop(index=list(taken)).nlargest(N_TOP_OB).index
        states["OB"] = list(top)
        C.log(f"[{ds}] no OB terminal state predicted -> set top {N_TOP_OB} ob_mature cells as OB")
    g.set_terminal_states(states)
    g.compute_fate_probabilities(show_progress_bar=False)
    fates = fates_df(g, a)
    tick("fate probabilities")
    terminal = g.terminal_states.reindex(a.obs_names)

    term_rows = []
    for s, cells in states.items():
        lab = a.obs.loc[cells, "pool_label"].value_counts()
        term_rows.append({"terminal_state": s, "n_cells": len(cells), "set_by_hand": s == "OB" and ob_set_by_hand,
                          **{f"z_{m}": round(a.obs.loc[cells, f"modz_{m}"].mean(), 2) for m in LINEAGE_MODS},
                          "ob_annot_frac": round(a.obs.loc[cells, "ob_annot"].astype(bool).mean(), 3),
                          "top_labels": "; ".join(f"{k} ({v})" for k, v in lab.head(3).items())})
    pd.DataFrame(term_rows).to_csv(os.path.join(outd, f"{ds}_terminal_states.csv"), index=False)

    drivers = g.compute_lineage_drivers(lineages=["OB"])
    drivers = drivers.sort_values("OB_corr", ascending=False).head(200)
    drivers.index.name = "gene"
    drivers.to_csv(os.path.join(outd, f"{ds}_ob_drivers.csv"))
    tick("lineage drivers")

    # states (main threshold + sensitivity)
    obz = a.obs["modz_ob_early"].values
    st = {t: assign_states(fates, terminal, obz, t) for t in THRESHOLDS}
    a.obs["state"] = pd.Categorical(st[MAIN_THR])
    tab = pd.DataFrame({f"thr_{t}": st[t].value_counts() for t in THRESHOLDS}).fillna(0).astype(int)
    tab.index.name = "state"
    tab.to_csv(os.path.join(outd, f"{ds}_state_table.csv"))
    for c in fates:
        a.obs[f"fate_{c}"] = fates[c].values
    a.obs["fate_max"] = fates.max(1).values

    # sensitivity kernel with the same terminal states
    g2 = cr.estimators.GPCCA(kernel(a, "cytotrace"))
    g2.set_terminal_states(states)
    g2.compute_fate_probabilities(show_progress_bar=False)
    st2 = assign_states(fates_df(g2, a), g2.terminal_states.reindex(a.obs_names), obz, MAIN_THR)
    ob1, ob2 = set(a.obs_names[st[MAIN_THR].isin(OB_STATES)]), set(a.obs_names[st2.isin(OB_STATES)])
    jac = len(ob1 & ob2) / max(1, len(ob1 | ob2))
    a.obs["state_cytotrace"] = pd.Categorical(st2)
    tick("CytoTRACE sensitivity")

    a.obs["state_fine"] = pd.Categorical(state_fine(a))
    fine = a.obs.groupby("state_fine", observed=True).agg(
        state=("state", "first"), n=("state", "size"), ob_annot=("ob_annot", "mean"),
        dpt=("dpt_pseudotime", "median"), fate_OB=("fate_OB", "mean"))
    fine["top_labels"] = [
        "; ".join(f"{k} ({v})" for k, v in a.obs.loc[a.obs.state_fine == f, "pool_label"].value_counts().head(3).items())
        for f in fine.index]
    fine["timepoints"] = [
        "; ".join(f"{k} ({v})" for k, v in a.obs.loc[a.obs.state_fine == f, "timepoint"].value_counts().items())
        for f in fine.index]
    if ds == "zscape":
        flag = str(cfg["flags"]["zscape_muscle_ambient_cluster"])
        a.obs["flag_muscle_ambient"] = a.obs["pool_cluster"].astype(str) == flag
        fine["frac_flag_muscle_ambient"] = a.obs.groupby("state_fine", observed=True)["flag_muscle_ambient"].mean()
    fine.round(3).to_csv(os.path.join(outd, f"{ds}_state_fine.csv"))

    pos = positional(a, cfg, ds)
    pos.round(3).to_csv(os.path.join(outd, f"{ds}_positional.csv"), index=False)

    # checks
    order_ok, peaks, n_ob = ob_order(a, ds, outd)
    checks = {"dataset": ds, "cells": a.n_obs, "n_states": n, "root": a.obs_names[root],
              "root_label": a.obs.pool_label.iloc[root], "unreachable_cells": n_inf,
              "ob_terminal_predicted": not ob_set_by_hand, "ob_lineage_cells": n_ob,
              "terminal_states": list(states), "not_terminal_cycling": dropped,
              "ob_lineage_ob_annot_frac": round(a.obs.loc[a.obs.state.isin(OB_STATES), "ob_annot"].astype(bool).mean(), 3),
              "marker_order_pass": order_ok, "marker_peaks": peaks,
              "kernel_jaccard_ob_lineage": round(jac, 3), "jaccard_target_0.6_met": jac >= 0.6,
              "na_state": int(a.obs["state"].isna().sum())}
    if ds == "zscape":
        p48 = a.obs.dpt_pseudotime[a.obs.timepoint == "48"]
        p72 = a.obs.dpt_pseudotime[a.obs.timepoint == "72"]
        checks.update({"median_dpt_48": round(p48.median(), 3), "median_dpt_72": round(p72.median(), 3),
                       "median_dpt_96": round(a.obs.dpt_pseudotime[a.obs.timepoint == "96"].median(), 3),
                       "mwu_p_48_lt_72": mannwhitneyu(p48, p72, alternative="less").pvalue,
                       "flag_muscle_ambient_states": a.obs.loc[a.obs.flag_muscle_ambient, "state"]
                       .value_counts().to_dict()})
    else:
        checks.update({f"median_dpt_{t}": round(a.obs.dpt_pseudotime[a.obs.timepoint == t].median(), 3)
                       for t in ("3dpf", "5dpf")})
    pd.Series(checks).to_csv(os.path.join(outd, f"{ds}_checks.csv"), header=["value"])
    C.log(pd.Series(checks).to_string())
    assert checks["na_state"] == 0

    # figures
    sc.settings.figdir = outd
    a.obs["root"] = pd.Categorical(np.where(np.arange(a.n_obs) == root, "root", "other"))
    sc.pl.umap(a, color=["state", "timepoint", "pool_label", "dpt_pseudotime", "fate_OB", "fate_max"],
               ncols=3, show=False, size=8, legend_fontsize=6, wspace=0.6)
    plt.savefig(os.path.join(outd, f"{ds}_umap_states.png"), dpi=130, bbox_inches="tight")
    plt.close("all")

    del a.layers["Ms"]
    # CellRank Lineage objects in obsm cannot be written by this anndata version: store plain arrays.
    for k in list(a.obsm):
        a.obsm[k] = np.asarray(a.obsm[k])
    h5 = os.path.join(cfg["paths"]["results"], f"{ds}_states.h5ad")
    try:
        a.write_h5ad(h5, compression="gzip")
    except Exception as e:  # noqa: BLE001  (never lose a long run to a serialization detail)
        C.log(f"[{ds}] full write failed ({e}); retrying without CellRank uns entries")
        for k in [k for k in a.uns if k not in ("ds", "iroot", "neighbors", "pca", "diffmap_evals",
                                                 "leiden", "umap", "hvg")]:
            del a.uns[k]
        a.write_h5ad(h5, compression="gzip")
    tick("wrote states h5ad")
    concordance(cfg, outd)


if __name__ == "__main__":
    main()
