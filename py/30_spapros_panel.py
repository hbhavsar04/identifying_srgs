"""Phase 3: HCR panel (prompts/ob_lineage_plan_v1.md) with spapros, per dataset, then consensus.

Stages (select/evaluate run per dataset and can run in parallel):
  --stage select   --dataset zscape|fabian   prepare cell set + gene universe, run ProbesetSelector
  --stage consensus                          tiered consensus to ~60 new genes
  --stage de       --dataset zscape|fabian   wilcoxon DE ranking per class (for the DE comparators)
  --stage evaluate --dataset zscape|fabian   ProbesetEvaluator on panel-only / +consensus / +random x5 /
                                             +DE_self / +DE_other / +DE_both / +HVG (needs de for both)
  --stage de_panel                           de_both_genes.csv: final panel when panel.final_set = de_both
  --stage table                              hcr_panel_v1.csv (from panel.final_set) + evaluation summary + done-criteria checks

Cell sets: ZSCAPE 72 hpf refined pool, Fabian 3 dpf refined pool; classes = state_fine, plus the
"OB_proximal" class (Gate 2 stopgap; top 1% of non-OB cells by OB fate within the set; WEAK EVIDENCE).
Gene universe: genes shared by the ZSCAPE and Fabian 3 dpf exports, minus mt-/rps/rpl/hemoglobin and
contaminant-module genes, and trainable in ZSCAPE (>= 5% detected in >= 1 class of the 72 hpf set).
spapros `n` counts the preselected genes (installed docstring), so n = n_panel_present + 70.
No prior_genes are passed.

Run: uv run --project env python py/30_spapros_panel.py --stage ... [--dataset ...]
Outputs: results/lineage/panel/
"""
import argparse
import os
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

import common as C

warnings.filterwarnings("ignore")

# Genes kept but flagged for HCR review (not excluded by the HCR-suitability rules).
REVIEW_FLAGS = {"hmgb2a": "proliferation-associated (marks cycling cells)",
                "tmsb4x": "near-ubiquitous actin-sequestering gene (housekeeping-like)"}

OB_STATES = ("OB_mature", "OB_committed", "OB_progenitor")
EVAL_METRICS = ["cluster_similarity", "knn_overlap", "forest_clfs", "gene_corr"]


def pdir(cfg, *p):
    return C.out_dir(cfg, "panel", *p)


# ---- cell set and gene universe ------------------------------------------------------------
def export_genes(cfg, ds):
    return pd.read_csv(os.path.join(cfg["paths"]["results"], "export", ds, "genes.csv"))["gene"].astype(str)


def base_universe(cfg):
    shared = sorted(set(export_genes(cfg, "zscape")) & set(export_genes(cfg, "fabian3")))
    contam = {g for gs in cfg["modules"]["contaminant"].values() for g in gs}
    bad = re.compile(r"^(mt-|rp[sl]\d|rp[sl]p|hb[abz])")
    return [g for g in shared if not bad.match(g) and g not in contam]


def add_ob_proximal(a, top_frac):
    """OB_proximal = top `top_frac` of non-OB cells by fate_OB within this cell set (weak evidence)."""
    a.obs["state"] = a.obs["state"].astype(str)
    a.obs["state_fine"] = a.obs["state_fine"].astype(str)
    non = ~a.obs["state"].isin(OB_STATES)
    n = max(1, int(round(top_frac * non.sum())))
    top = a.obs.loc[non, "fate_OB"].nlargest(n).index
    a.obs["state_pre_ob_proximal"] = a.obs["state"]
    a.obs.loc[top, "state"] = "OB_proximal"
    a.obs.loc[top, "state_fine"] = "OB_proximal.0"
    return top


def class_labels(a, min_cells):
    """spapros classes: state_fine; classes < min_cells merged into '<state>.small'.
    A '.small' class still < min_cells is kept as a label but excluded from the target celltypes."""
    lab = a.obs["state_fine"].astype(str).copy()
    sizes = lab.value_counts()
    small = sizes.index[sizes < min_cells]
    lab[lab.isin(small)] = a.obs.loc[lab.isin(small), "state"].astype(str) + ".small"
    sizes = lab.value_counts()
    a.obs["spapros_ct"] = pd.Categorical(lab)
    return sorted(sizes.index[sizes >= min_cells]), sizes


def class_detection(a, genes, key="spapros_ct"):
    det = (a[:, genes].layers["counts"] > 0)
    out = {}
    for c, idx in a.obs.groupby(key, observed=True).indices.items():
        out[c] = np.asarray(det[idx].mean(0)).ravel()
    return pd.DataFrame(out, index=genes)


def load_set(cfg, ds):
    p = cfg["panel"]
    a = sc.read_h5ad(os.path.join(cfg["paths"]["results"], f"{ds}_states.h5ad"))
    a = a[a.obs["timepoint"].astype(str) == p["spapros_sets"][ds]].copy()
    top = add_ob_proximal(a, cfg["ob_proximal"]["top_frac"])
    cts, sizes = class_labels(a, p["min_class_cells"])
    return a, top, cts, sizes


def universe(cfg):
    """Trainable universe (cached): base universe filtered by ZSCAPE 72 hpf class detection."""
    f = os.path.join(pdir(cfg), "gene_universe.csv")
    if os.path.exists(f):
        return pd.read_csv(f)
    base = base_universe(cfg)
    rows = {}
    for ds in ("zscape", "fabian"):
        a, _, _, _ = load_set(cfg, ds)
        g = [x for x in base if x in a.var_names]
        rows[ds] = class_detection(a, g).max(1).reindex(base)
    u = pd.DataFrame({"gene": base, "max_class_det_zscape72": rows["zscape"].values,
                      "max_class_det_fabian3": rows["fabian"].values})
    th = cfg["panel"]["min_det_trainable"]
    u["trainable"] = u.max_class_det_zscape72 >= th
    u["fabian_only_detectable"] = ~u.trainable & (u.max_class_det_fabian3 >= th)
    u.round(4).to_csv(f, index=False)
    C.log(f"universe: {len(base)} shared/filtered genes; {u.trainable.sum()} trainable in ZSCAPE 72 hpf; "
          f"{u.fabian_only_detectable.sum()} Fabian-only detectable")
    return u


# ---- stage: select ---------------------------------------------------------------------------
def stage_select(cfg, ds):
    import spapros as sp
    u = universe(cfg)
    a, top, cts, sizes = load_set(cfg, ds)
    panel = [g for g in cfg["panel_genes"] if g in a.var_names]
    genes = sorted(set(u.gene[u.trainable]) & set(a.var_names) | set(panel))
    a = a[:, genes].copy()
    sc.pp.highly_variable_genes(a, n_top_genes=min(cfg["panel"]["n_hvg"], a.n_vars))
    a.var["spapros_genes"] = a.var["highly_variable"] | a.var_names.isin(panel)
    n = len(panel) + cfg["panel"]["n_new_select"]
    info = {"dataset": ds, "cells": a.n_obs, "genes_in": a.n_vars, "spapros_genes": int(a.var.spapros_genes.sum()),
            "panel_present": len(panel), "n": n, "classes_targeted": len(cts), "ob_proximal_cells": len(top),
            "ob_proximal_fate_min": round(a.obs.loc[top, "fate_OB"].min(), 4),
            "ob_proximal_fate_median": round(a.obs.loc[top, "fate_OB"].median(), 4)}
    C.log(f"[{ds}] {info}")
    lab = a.obs.loc[top, "state_pre_ob_proximal"].value_counts().to_dict()
    pl = a.obs.loc[top, "pool_label"].value_counts().to_dict()
    pd.DataFrame({"class": sizes.index, "n_cells": sizes.values, "targeted": sizes.index.isin(cts)}).to_csv(
        os.path.join(pdir(cfg), f"{ds}_classes.csv"), index=False)
    pd.Series({**info, "ob_proximal_from_states": lab, "ob_proximal_labels": pl}).to_csv(
        os.path.join(pdir(cfg), f"{ds}_select_info.csv"), header=["value"])
    sel = sp.se.ProbesetSelector(a, celltype_key="spapros_ct", genes_key="spapros_genes", n=n,
                                 preselected_genes=panel, celltypes=cts, seed=C.SEED, verbosity=1,
                                 save_dir=pdir(cfg, f"spapros_{ds}"), n_jobs=12)
    sel.select_probeset()
    ps = sel.probeset.copy()
    ps.index.name = "gene"
    ps.to_csv(os.path.join(pdir(cfg), f"{ds}_probeset.csv"))
    new = ps[ps["selection"] & ~ps.index.isin(panel)]
    C.log(f"[{ds}] selected {int(ps['selection'].sum())} genes, {len(new)} new")


# ---- stage: consensus -----------------------------------------------------------------------
def state_detection(cfg):
    """Detection per coarse state (incl. OB_proximal) per dataset on the spapros cell sets, plus
    full-dataset detection for states absent from the set (e.g. Fabian OB states are 5 dpf only)."""
    out = {}
    for ds in ("zscape", "fabian"):
        a, _, _, _ = load_set(cfg, ds)
        full = sc.read_h5ad(os.path.join(cfg["paths"]["results"], f"{ds}_states.h5ad"))
        g = list(a.var_names)
        d = class_detection(a, g, key="state")
        for s in set(full.obs["state"].astype(str)) - set(d.columns):
            m = (full.obs["state"].astype(str) == s).values
            d[s] = np.asarray((full[m][:, g].layers["counts"] > 0).mean(0)).ravel()
        d["_all_set"] = np.asarray((a.layers["counts"] > 0).mean(0)).ravel()
        out[ds] = d
    return out


def concordant_map(cfg):
    cor = pd.read_csv(os.path.join(cfg["paths"]["results"], "trajectory", "state_concordance.csv"), index_col=0)
    z2f = cor.idxmax(1).to_dict()
    f2z = cor.idxmax(0).to_dict()
    z2f["OB_proximal"] = "OB_proximal"
    f2z["OB_proximal"] = "OB_proximal"
    return z2f, f2z


def excluded_reason(cfg, g):
    explicit = cfg["panel"].get("exclude_genes") or {}
    if g in explicit:
        return f"explicit:{explicit[g]}"
    for reason, rx in cfg["panel"]["exclude_regex"].items():
        if re.match(rx, g):
            return reason
    return ""


def stage_consensus(cfg):
    p = cfg["panel"]
    panel = set(cfg["panel_genes"])
    ps = {ds: pd.read_csv(os.path.join(pdir(cfg), f"{ds}_probeset.csv"), index_col=0) for ds in ("zscape", "fabian")}
    det = state_detection(cfg)
    z2f, f2z = concordant_map(cfg)
    th = p["min_det_trainable"]
    sel = {ds: ps[ds].index[ps[ds]["selection"] & ~ps[ds].index.isin(panel)] for ds in ps}
    nsel = {ds: len(sel[ds]) for ds in sel}

    # spapros 'rank' is identical for every selected gene, so order selected new genes within each
    # dataset by importance_score (descending); normalized rank = position / n selected.
    imp_rank = {ds: {g: (k + 1) / nsel[ds] for k, g in enumerate(
        ps[ds].loc[sel[ds]].sort_values("importance_score", ascending=False).index)} for ds in ps}

    def nrank(ds, g):
        return imp_rank[ds].get(g, np.nan)

    rows = []
    for g in sorted(set(sel["zscape"]) | set(sel["fabian"])):
        inz, inf = g in sel["zscape"], g in sel["fabian"]
        r = {"gene": g, "sel_zscape": inz, "sel_fabian": inf,
             "rank_zscape": round(nrank("zscape", g) * nsel["zscape"]) if inz else np.nan,
             "rank_fabian": round(nrank("fabian", g) * nsel["fabian"]) if inf else np.nan,
             "excluded": excluded_reason(cfg, g),
             "zscape_det": det["zscape"].loc[g, "_all_set"] if g in det["zscape"].index else 0.0}
        src, other = ("zscape", "fabian") if inz else ("fabian", "zscape")
        mp = z2f if src == "zscape" else f2z
        d = det[src].drop(columns="_all_set").loc[g] if g in det[src].index else pd.Series(dtype=float)
        top_state = d.idxmax() if len(d) else None
        conc = mp.get(top_state)
        conc_det = det[other].loc[g, conc] if (conc in det[other].columns and g in det[other].index) else 0.0
        r.update({"top_state_src": top_state, "concordant_state_other": conc, "det_concordant_other": conc_det})
        if r["excluded"]:
            r["tier"], r["source"] = 0, ("both" if inz and inf else src)
            r["mean_nrank"] = np.nanmean([nrank("zscape", g), nrank("fabian", g)])
        elif inz and inf:
            r["tier"], r["source"] = 1, "both"
            r["mean_nrank"] = np.nanmean([nrank("zscape", g), nrank("fabian", g)])
        elif conc_det >= th:
            r["tier"], r["source"] = 2, src
            r["mean_nrank"] = nrank(src, g)
        else:
            r["tier"], r["source"] = 3, src
            r["mean_nrank"] = nrank(src, g)
        rows.append(r)
    cand = pd.DataFrame(rows)
    order = lambda d: d.sort_values(["mean_nrank", "zscape_det"], ascending=[True, False])  # noqa: E731
    picked = list(order(cand[cand.tier == 1]).gene)[: p["n_new_target"]]
    queues = {s: list(order(cand[(cand.tier == 2) & (cand.source == s)]).gene) for s in ("zscape", "fabian")}
    turn = ["zscape", "fabian"]
    i = 0
    while len(picked) < p["n_new_target"] and any(queues.values()):
        q = queues[turn[i % 2]]
        if q:
            picked.append(q.pop(0))
        i += 1
    cand["in_consensus"] = cand.gene.isin(picked)
    cand["consensus_order"] = cand.gene.map({g: k + 1 for k, g in enumerate(picked)})
    cand.round(4).sort_values(["in_consensus", "consensus_order"], ascending=[False, True]).to_csv(
        os.path.join(pdir(cfg), "consensus_candidates.csv"), index=False)
    pd.Series(picked, name="gene").to_csv(os.path.join(pdir(cfg), "consensus_genes.csv"), index=False)
    C.log(f"consensus: {len(picked)} genes (tier1 {int((cand.tier == 1).sum())}, tier2 available "
          f"{int((cand.tier == 2).sum())}, tier3 {int((cand.tier == 3).sum())}); per-dataset new: {nsel}; "
          f"excluded: {cand.loc[cand.tier == 0, ['gene', 'excluded']].values.tolist()}")


# ---- stage: evaluate ------------------------------------------------------------------------
def comparison_sets(cfg, ds, a, panel, cons, u):
    rng = np.random.default_rng(C.SEED)
    pool = [g for g in u.gene[u.trainable] if g in a.var_names and g not in panel and g not in cons
            and not excluded_reason(cfg, g)]
    detz = u.set_index("gene")["max_class_det_zscape72"]
    bins = pd.qcut(detz.loc[pool + [g for g in cons if g in detz.index]], 10, labels=False, duplicates="drop")
    sets = {"panel_only": panel, "panel+consensus": panel + cons}
    by_bin = {b: [g for g in pool if bins.get(g) == b] for b in np.unique(bins.values)}
    for k in range(cfg["panel"]["n_random_draws"]):
        draw = []
        for g in cons:
            cands = [x for x in by_bin.get(bins.get(g), pool) if x not in draw]
            draw.append(rng.choice(cands))
        sets[f"panel+random{k + 1}"] = panel + draw
    n = len(cons)
    # DE / HVG comparators may use consensus genes (only the random draws avoid them).
    elig = pool + [g for g in cons if g in a.var_names]
    both, de, _ = de_both(cfg, set(elig), n)
    other = "fabian" if ds == "zscape" else "zscape"
    sets["panel+DE_self"] = panel + de[ds]
    sets["panel+DE_other"] = panel + de[other]
    sets["panel+DE_both"] = panel + both
    b = a[:, elig].copy()
    sc.pp.highly_variable_genes(b, n_top_genes=n)
    sets["panel+HVG"] = panel + list(b.var_names[b.var["highly_variable"]])
    return sets


def eval_cells(cfg, ds):
    """Evaluation cell set: targeted classes only, genes = trainable universe present + panel."""
    u = universe(cfg)
    a, _, cts, _ = load_set(cfg, ds)
    panel = [g for g in cfg["panel_genes"] if g in a.var_names]
    genes = sorted(set(u.gene[u.trainable]) & set(a.var_names) | set(panel))
    a = a[:, genes].copy()
    return a[a.obs["spapros_ct"].isin(cts)].copy(), panel, u


def stage_de(cfg, ds):
    """Wilcoxon DE per spapros class on the eligible pool (trainable, not panel, not excluded)."""
    a, panel, u = eval_cells(cfg, ds)
    elig = [g for g in a.var_names if g not in panel and not excluded_reason(cfg, g)
            and g in set(u.gene[u.trainable])]
    b = a[:, elig].copy()
    sc.tl.rank_genes_groups(b, "spapros_ct", method="wilcoxon", n_genes=cfg["panel"]["de_n_per_class"])
    names = pd.DataFrame(b.uns["rank_genes_groups"]["names"])
    names.to_csv(os.path.join(pdir(cfg), f"{ds}_de_ranked.csv"), index=False)
    C.log(f"[{ds}] DE ranked: {names.shape[1]} classes x {names.shape[0]} genes from {len(elig)} eligible")


def load_de(cfg, ds):
    f = os.path.join(pdir(cfg), f"{ds}_de_ranked.csv")
    if not os.path.exists(f):
        raise FileNotFoundError(f"{f}: run --stage de --dataset {ds} first")
    return pd.read_csv(f)


def round_robin(names, allowed, n):
    """Take the next-best gene of each class in turn (rank level by rank level) until n genes."""
    out = []
    for row in names.itertuples(index=False):
        for g in row:
            if isinstance(g, str) and g in allowed and g not in out:
                out.append(g)
                if len(out) == n:
                    return out
    return out


def de_both(cfg, allowed, n):
    """Alternate the ZSCAPE and Fabian round-robin DE lists (ZSCAPE first), deduplicated, to n genes."""
    de = {d: round_robin(load_de(cfg, d), allowed, n) for d in ("zscape", "fabian")}
    both, src = [], {}
    for k in range(n):
        for d in ("zscape", "fabian"):
            if k < len(de[d]) and len(both) < n:
                g = de[d][k]
                if g not in both:
                    both.append(g)
                    src[g] = d
    return both, de, src


def de_hit(names, g):
    """Best (rank, class) of gene g in a DE ranking table (rank 1 = top)."""
    best = (np.nan, "")
    for c in names.columns:
        idx = np.flatnonzero(names[c].values == g)
        if idx.size and not (idx[0] + 1 >= best[0]):
            best = (idx[0] + 1, c)
    return best


def stage_de_panel(cfg):
    """Final panel = DE_both (user choice, 2026-11-26) from the trainable, non-excluded, non-panel pool."""
    u = universe(cfg)
    panel = set(cfg["panel_genes"])
    allowed = {g for g in u.gene[u.trainable] if g not in panel and not excluded_reason(cfg, g)}
    n = cfg["panel"]["n_new_target"]
    genes, de, src = de_both(cfg, allowed, n)
    cons = set(pd.read_csv(os.path.join(pdir(cfg), "consensus_genes.csv"))["gene"])
    rows = []
    for k, g in enumerate(genes):
        r = {"gene": g, "order": k + 1, "added_from": src[g], "in_spapros_consensus": g in cons}
        for d in ("zscape", "fabian"):
            rk, c = de_hit(load_de(cfg, d), g)
            r[f"de_rank_{d}"], r[f"de_class_{d}"] = rk, c
            r[f"in_DE_{d}"] = g in de[d]
        rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(pdir(cfg), "de_both_genes.csv"), index=False)
    C.log(f"DE_both panel: {len(genes)} genes, {sum(g in cons for g in genes)} shared with spapros consensus")


def stage_evaluate(cfg, ds):
    import spapros as sp
    a, panel, u = eval_cells(cfg, ds)               # evaluate on targeted classes only
    cons = [g for g in pd.read_csv(os.path.join(pdir(cfg), "consensus_genes.csv"))["gene"] if g in a.var_names]
    sets = comparison_sets(cfg, ds, a, panel, cons, u)
    # Genes with zero variance in this cell set break spapros' gene_corr (NaN correlations): drop them
    # from every set for this dataset and log them.
    x = a.X
    var = np.asarray(x.power(2).mean(0)).ravel() - np.asarray(x.mean(0)).ravel() ** 2
    flat = set(a.var_names[var <= 0])
    dropped = sorted(flat & set(g for v in sets.values() for g in v))
    if dropped:
        C.log(f"[{ds}] zero-variance genes dropped from evaluation sets: {dropped}")
    sets = {k: [g for g in v if g not in flat] for k, v in sets.items()}
    pd.Series(dropped, name="gene", dtype=str).to_csv(os.path.join(pdir(cfg), f"{ds}_eval_zero_var.csv"), index=False)
    pd.DataFrame([{"set": k, "n_genes": len(v), "n_shared_with_consensus": len(set(v) & set(cons)),
                   "genes": ";".join(v)} for k, v in sets.items()]).to_csv(
        os.path.join(pdir(cfg), f"{ds}_eval_sets.csv"), index=False)
    ev = sp.ev.ProbesetEvaluator(a, celltype_key="spapros_ct", scheme="custom", metrics=EVAL_METRICS,
                                 results_dir=pdir(cfg, f"eval_{ds}"), verbosity=0, n_jobs=12)
    for k, v in sets.items():
        C.log(f"[{ds}] evaluating {k} ({len(v)} genes)")
        ev.evaluate_probeset(v, set_id=k)
    ev.summary_results.to_csv(os.path.join(pdir(cfg), f"{ds}_eval_summary.csv"))
    # Per-class recall: spapros forest_clfs stores a row-normalized confusion matrix (true x predicted);
    # its summary 'accuracy' is the mean of the diagonal (balanced accuracy). spapros reports no F1.
    rows = []
    for k in sets:
        cm = ev.results["forest_clfs"][k]
        cm = cm[0] if isinstance(cm, (list, tuple)) else cm
        cm = pd.DataFrame(cm)
        rows.append(pd.Series(np.diag(cm.loc[cm.index, cm.index].values), index=cm.index, name=k))
    pd.concat(rows, axis=1).round(4).to_csv(os.path.join(pdir(cfg), f"{ds}_eval_recall_by_class.csv"))
    C.log(ev.summary_results.round(3).to_string())


# ---- stage: table ------------------------------------------------------------------------------
def positional_regex():
    pats = [r"^hox[abcd][0-9]+[ab]?a?$", r"^dlx[1-6][ab]?$", r"^hand2$", r"^barx1$", r"^pax[0-9]+[ab]?$",
            r"^six[0-9]+[ab]?$", r"^tbx[0-9]+[ab]?$", r"^msx[0-9][a-e]?$", r"^shha?$", r"^en[12][ab]?$",
            r"^pitx[12]$", r"^nr2f[0-9]$", r"^edn1$", r"^jag1b$", r"^irx[0-9][ab]?$", r"^fgf[0-9]+[ab]?$",
            r"^lhx[0-9][ab]?$"]   # same as positional_genes() in R/utils.R
    return re.compile("|".join(pats))


def muscle_detection(cfg, genes):
    import scipy.io
    d = cfg["context72_export"]
    cells = pd.read_csv(os.path.join(d, "cells.csv"), usecols=["cell", "cell_type_sub"])
    g_all = pd.read_csv(os.path.join(d, "genes.csv"))["gene"].astype(str)
    X = scipy.io.mmread(os.path.join(d, "counts.mtx")).tocsr()
    idx = [i for i, g in enumerate(g_all) if g in set(genes)]
    sub = X[idx] > 0
    mus = cells.cell_type_sub.str.startswith("cranial muscle").values
    return pd.DataFrame({"det_cranial_muscle72": np.asarray(sub[:, mus].mean(1)).ravel(),
                         "det_nonmuscle_context72": np.asarray(sub[:, ~mus].mean(1)).ravel()},
                        index=g_all.iloc[idx].values)


def nominee_flags(cfg, genes):
    flags = {}
    for k, f in cfg["nominee_lists"].items():
        try:
            gl = set(pd.read_csv(f)["gene"].astype(str))
        except Exception as e:  # noqa: BLE001
            C.log(f"nominee list {k} unreadable: {e}")
            gl = set()
        flags[f"in_{k}"] = [g in gl for g in genes]
    flags["in_sox10_list"] = [g in set(cfg["sox10_list"]) for g in genes]
    return pd.DataFrame(flags, index=genes)


def stage_table(cfg):
    final = cfg["panel"]["final_set"]
    cand = pd.read_csv(os.path.join(pdir(cfg), "consensus_candidates.csv")).set_index("gene")
    if final == "de_both":
        out = pd.read_csv(os.path.join(pdir(cfg), "de_both_genes.csv")).set_index("gene")
        genes = list(out.index)
        sp_cols = cand.reindex(genes)[["tier", "rank_zscape", "rank_fabian"]].rename(
            columns={"tier": "spapros_tier", "rank_zscape": "spapros_rank_zscape", "rank_fabian": "spapros_rank_fabian"})
        out = out.join(sp_cols)
    else:
        genes = list(pd.read_csv(os.path.join(pdir(cfg), "consensus_genes.csv"))["gene"])
        out = cand.loc[genes, ["tier", "source", "rank_zscape", "rank_fabian", "top_state_src",
                               "concordant_state_other", "det_concordant_other"]].copy()
    det = state_detection(cfg)
    for ds in ("zscape", "fabian"):
        d = det[ds].drop(columns="_all_set").reindex(genes)
        out[f"states_marked_{ds}"] = [";".join(d.columns[(d.loc[g] >= 0.2) & (d.loc[g] >= 2 * d.loc[g].median())])
                                      for g in genes]
        for s in d.columns:
            out[f"pct_{ds}_{s}"] = (100 * d[s]).round(1)
    z = sc.read_h5ad(os.path.join(cfg["paths"]["results"], "zscape_states.h5ad"), backed="r")
    ob72 = ((z.obs["ob_annot"].astype(bool)) & (z.obs["timepoint"].astype(str) == "72")).values
    zo = z[ob72].to_memory()
    gz = [g for g in genes if g in zo.var_names]
    out["pct_zscape72_OB"] = pd.Series(100 * np.asarray((zo[:, gz].layers["counts"] > 0).mean(0)).ravel(),
                                       index=gz).round(1)
    for ds in ("zscape", "fabian"):
        drv = pd.read_csv(os.path.join(cfg["paths"]["results"], "trajectory", f"{ds}_ob_drivers.csv"))
        out[f"ob_driver_rank_{ds}"] = pd.Series(np.arange(1, len(drv) + 1), index=drv.gene).reindex(genes).values
    out = out.join(nominee_flags(cfg, genes))
    rx = positional_regex()
    out["positional_gene"] = [bool(rx.match(g)) for g in genes]
    mus = muscle_detection(cfg, genes).reindex(genes)
    out = out.join(mus)
    zcols = [c for c in out if c.startswith("pct_zscape_")]
    notes = []
    for g in genes:
        n = []
        if (out.loc[g, zcols] >= 30).all():
            n.append("broadly high in every ZSCAPE state")
        m, nm = mus.loc[g, "det_cranial_muscle72"], mus.loc[g, "det_nonmuscle_context72"]
        if pd.notna(m) and m >= 0.2 and m >= 2 * max(nm, 1e-9):
            n.append("also marks cranial muscle")
        if "OB_proximal" in str(out.loc[g, "states_marked_zscape"]) or "OB_proximal" in str(out.loc[g, "states_marked_fabian"]):
            n.append("marks OB_proximal (weak-evidence class)")
        if g in REVIEW_FLAGS:
            n.append(REVIEW_FLAGS[g])
        notes.append("; ".join(n))
    out["notes"] = notes
    out.index.name = "gene"
    out.to_csv(os.path.join(pdir(cfg), "hcr_panel_v1.csv"))

    # done criteria
    panel = set(cfg["panel_genes"])
    u = universe(cfg).set_index("gene")
    checks = {"n_new": len(genes), "in_50_70": 50 <= len(genes) <= 70,
              "none_in_panel": not (set(genes) & panel),
              "all_trainable": bool(u.loc[genes, "trainable"].all())}
    for ds in ("zscape", "fabian"):
        f = os.path.join(pdir(cfg), f"{ds}_eval_summary.csv")
        if os.path.exists(f):
            s = pd.read_csv(f, index_col=0)
            col = [c for c in s.columns if "forest" in c.lower() and "accuracy" in c.lower()]
            if col:
                f1 = s[col[0]]
                rnd = f1[[i for i in f1.index if i.startswith("panel+random")]]
                checks[f"{ds}_balacc_panel_only"] = round(f1["panel_only"], 4)
                checks[f"{ds}_balacc_consensus"] = round(f1["panel+consensus"], 4)
                checks[f"{ds}_balacc_random_max"] = round(rnd.max(), 4)
                fk = "panel+DE_both" if final == "de_both" else "panel+consensus"
                checks[f"{ds}_final_set"] = fk
                checks[f"{ds}_balacc_final"] = round(f1[fk], 4)
                checks[f"{ds}_beats_panel_and_random"] = bool(f1[fk] > max(f1["panel_only"], rnd.max()))
                for k in ("panel+DE_self", "panel+DE_other", "panel+DE_both"):
                    if k in f1.index:
                        checks[f"{ds}_balacc_{k.split('+')[1]}"] = round(f1[k], 4)
    pd.Series(checks).to_csv(os.path.join(pdir(cfg), "panel_checks.csv"), header=["value"])
    C.log(pd.Series(checks).to_string())
    assert checks["in_50_70"] and checks["none_in_panel"] and checks["all_trainable"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["select", "consensus", "de", "de_panel", "evaluate", "table"])
    ap.add_argument("--dataset", choices=["zscape", "fabian"])
    args = ap.parse_args()
    cfg = C.load_config()
    np.random.seed(C.SEED)
    if args.stage == "select":
        stage_select(cfg, args.dataset)
    elif args.stage == "consensus":
        stage_consensus(cfg)
    elif args.stage == "de":
        stage_de(cfg, args.dataset)
    elif args.stage == "de_panel":
        stage_de_panel(cfg)
    elif args.stage == "evaluate":
        stage_evaluate(cfg, args.dataset)
    else:
        stage_table(cfg)


if __name__ == "__main__":
    main()
