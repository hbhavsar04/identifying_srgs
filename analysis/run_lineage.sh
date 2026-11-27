#!/usr/bin/env bash
# Part 1 of prompts/ob_lineage_plan_v1.md: pools -> trajectories -> HCR panel -> Tangram reference sets.
# Run from the identifying_srgs root with the /Volumes/kwongroup share mounted (read only, except the
# Tangram reference sets, which are written as new files under the scrna_matching SHARE_OUT).
# Runtime: pools ~20 min; trajectories ~1.7 h (ZSCAPE) + 40 min (Fabian), dominated by GPCCA's dense
# Schur step; spapros selection and evaluation come on top (not timed).
set -euo pipefail
PY="uv run --project env python"

Rscript analysis/10_render_pools.R zscape fabian3 fabian5      # Phase 1 (Gate 1)
Rscript analysis/11_export_pools.R zscape fabian3 fabian5

$PY py/20_trajectory.py --dataset zscape &                      # Phase 2 (Gate 2)
$PY py/20_trajectory.py --dataset fabian &
wait
$PY py/20_trajectory.py --concordance

$PY py/30_spapros_panel.py --stage select --dataset zscape &    # Phase 3
$PY py/30_spapros_panel.py --stage select --dataset fabian &
wait
$PY py/30_spapros_panel.py --stage consensus
$PY py/30_spapros_panel.py --stage de --dataset zscape
$PY py/30_spapros_panel.py --stage de --dataset fabian
$PY py/30_spapros_panel.py --stage de_panel                     # final panel = DE_both (config panel.final_set)
# spapros caches evaluation results by set name in results/lineage/panel/eval_<ds>/: delete those
# folders if the gene sets changed, or the old results are reused.
$PY py/30_spapros_panel.py --stage evaluate --dataset zscape &
$PY py/30_spapros_panel.py --stage evaluate --dataset fabian &
wait
$PY py/30_spapros_panel.py --stage table                        # hcr_panel_v1.csv + panel_checks.csv

$PY py/31_export_tangram_sets.py                                # Part 1 -> Part 2 handoff
