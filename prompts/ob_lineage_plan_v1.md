# Plan: expand the 3 dpf OB reference to the opercle-forming condensation, choose ~60 new atlas genes, and test the effect on Tangram

> Written for a fresh session (Opus 5.5, medium effort). Everything needed is in this file. Read it top to bottom before touching code.

## Context
**Goal.** Find spatially informative genes in the 3 dpf zebrafish opercle. They should be biologically meaningful and should also improve scRNA → atlas integration (Tangram).

**What exists already.** An HCR atlas of 31 channels:
- 29 genes plus dapi/reference, registered with ANTs SyGN and voxelized;
- the opercle region (33,451 voxels) is subclustered into subdomains A_body, B_pltp-slc8a4b, C_tip, D_posterior, plus a dim runx2a rim E_rim.

Tangram maps ZSCAPE 72 hpf osteoblasts (186 cells) onto that region. This is done in `voxel_measurements_and_analysis/scrna_matching`, Iterations 0–3 plus checks.

**Problem.** Leave-one-gene-out (LOO) works for genes of the committed-OB core and fails toward the periphery:
- median Pearson 0.36 (null 0.15);
- col10a1a 0.85 and col1a1a 0.83;
- runx2a −0.26, inverted. It is the rim gene;
- bcl2a 0.04, sdc2 0.09, mfsd2b 0.08, klf13 0.15.

Rim calls track depth and batch. The gap is two-sided:
- the reference holds only committed OBs;
- the panel was chosen from OB markers, so it holds little information about peripheral or progenitor cells.

**This plan:**
1. builds a broader skeletogenic neural-crest (NC) pool around the OBs in ZSCAPE and in Fabian 2022;
2. uses trajectory and fate analysis to define the OB lineage, including progenitors;
3. picks **~60 new atlas genes** (50–70) with spapros that complement the 29 panel genes;
4. tests in Tangram LOO whether the expanded reference improves peripheral prediction.

**Decisions made with the user:**
- Fabian 2022 is a **co-primary** dataset. Analyse it separately, never integrated: the user has seen little ZSCAPE–Fabian overlap, and Fabian 5 dpf gave poorer Tangram results.
- The scope runs **through Tangram LOO**.
- The panel is **50–70 new genes**.

## Existing work to reuse (read these first; do not rebuild)
**`/Users/msbl/Desktop/Hayden/git_repos/voxel_measurements_and_analysis/`**
- `README.md`: atlas settings and the progress log.
- `scrna_matching/README.md`, `PLAN.md`: the scRNA-matching iterations, gates and decisions.
- `scrna_matching/paths.py`:
  - `RUN_DIR`, `SUB_DIR`, `SHARE_OUT`, `DOMAINS` and `DOMAIN_ORDER`;
  - `ATLAS_TO_SC`: bmp8 → bmp8a, mfsd → mfsd2b, sich211 → si:ch211-106n13.3, znf385 → znf385b;
  - `iter_share_dir()`.
- Underlying paths (from `opercle_subclustering/subdomain_markers.py`):
  - `RUN_DIR = /Volumes/kwongroup/Hayden Bhavsar/voxel_csvs/ants_sygn_pooled_bin1x5x5_rollingball-percentile_mc2_excl-ihhb/` (`voxel_matrix.csv.gz` = whole shared mask, ~162k voxels; `shared_mask.tif`; `voxel_grid_meta.json`)
  - `SUB_DIR = RUN_DIR/seurat_scale_only_rpctFALSE_dims10_res0.3/opercle_subclustering/scale_only/res_0.3`
  - `SHARE_OUT = .../opercle_subclustering/scrna_matching/`
- `scrna_matching/iter0_data_prep/`:
  - `prepare_atlas.py` → `<SHARE_OUT>/iter0_data_prep/atlas_region.csv.gz` (voxel_id, slice, x, y, z_um, y_um, x_um, subcluster, domain, 31 genes);
  - `results/subdomain_profiles_{mean,z}.csv`;
  - `results/iter1_training_genes.csv`: 28 trained genes; ihha is validation-only;
  - `export_reference.R`: the R → mtx/genes.csv/cells.csv export pattern and the context-set rules. The existing `zscape_context72` set has 10,638 cells and includes cranial muscle;
  - `results/zscape72_celltypes.csv`, `results/context72_detection_by_celltype.csv`.
- `scrna_matching/iter1_ob_matching/common.py`:
  - `load_cells`, `load_atlas`, `training_genes`, `dense`;
  - `tangram_map`: cells mode with a uniform density prior, epochs 1000, lr 0.1; it records the effective lambdas;
  - `tangram_cell_probs`, `tangram_scores`, `domain_onehot`.
- `match_tangram.py`: manual LOO (one mapping per held-out gene, prediction = Mᵀx, 200-permutation null, top-bright-subdomain check), with pickle caching in `tg_runs/`.
- `scrna_matching/iter2_stage_checks/run_stage.py`: the generic runner for "map another cell set and score it", with optional inputs.
- `scrna_matching/checks/iter2_pool_impute/{build_sets,run_check,nomination_benchmark}.py`:
  - `run_check.rim_artifact()`: rim call vs UMI/batch;
  - the gate-summary table format;
  - the Level-2 simulated add-back.
- `spatial_stats/common.py`: `rebin()` for coarser voxel grids. Correlation structure holds up to 2.2 µm bins (spatial_stats Iteration 1).
- `scrna_matching/env/`: uv project with python 3.11–3.12, tangram-sc 1.0.4, scanpy 1.11, torch CPU and scvi. **Use it for every Tangram run.**

**`/Users/msbl/Desktop/Hayden/git_repos/Tangram`** is an upstream clone plus `hayden_testing/tutorial.py`. Read only; never modify it.

**Lessons already learned (keep them):**
- tangram-sc 1.0.4 seeds only `if random_state:`, so **use seed 1**. Seed 1 runs are bit-identical, and the noise floor is ≈ 0, so one seed per run is enough.
- kNN imputation and scVI pooling made the ZSCAPE map worse (`checks/iter2_pool_impute`). Don't impute; don't integrate the two datasets.
- In simulated add-back, **detection in the reference predicts a nominee's gain (ρ 0.39); specificity does not (ρ 0.06).** Panel genes have to be detectable in the Tangram reference to help.
- hox PG2 (hoxa2b/hoxb2a) is detected in 0% of ZSCAPE OBs. Opercle OBs can't be told apart from those of other bones by hox code. Positional scores may show only in mesenchyme, if at all.
- **Fabian 3 dpf has no OB cluster:**
  - its best cluster has 1% OB-gene detection;
  - abcc6a and ihha are in 0 cells, runx2a in 0.5%.
  - Fabian 5 dpf has 242 OBs (`sox10_sc_analysis/data/processed/5dpf_RNA_OB.rds`, cluster 14 of `5dpf_scRNA.rds`) at median 2,927 UMI, but they are more mature.
- Existing HCR nominee lists. Report overlap with them; don't use them as selection input.
  - `scrna_matching/iter3_candidates/results/hcr_shortlist.csv`
  - `.../rerank_detection/results/` (18-gene shortlist)
  - `.../exploration/results/`
  - `identifying_srgs/results/ob334/srg_candidates.csv`
  - the sox10-repo list (scxa, lum, mafa, matn4, ifitm5, panx3, col10a1a, slc8a4b, …)

**This repo (`identifying_srgs`).** Reuse from `R/utils.R`:
- `build_clean()` (extend `META_KEEP` with `bio.rep`)
- `embed_features()`, `choose_dims()`, `cluster_sweep()`, `sweep_summary()`, `pick_resolution()`
- `cramers_v()`, `positional_genes()`
- `SEED = 321`

## Data facts (checked)
**ZSCAPE (`saunders.rds`)**
- 636,534 cells × 32,031 genes.
- `gene_target`: ctrl-uninj / ctrl-inj / ctrl-hgfa / ctrl-mafba / ctrl-met / ctrl-noto / tbx16-mut. Drop `tbx16-mut`.
- Late timepoints are only 48, 72 and 96 hpf, with median UMI 391 / 358 / 410.
- 81% of the 72 hpf OBs come from one plate (ctrl-hgfa P17).
- Counts of the relevant `cell_type_sub` labels:

| `cell_type_sub` | 48 hpf | 72 hpf | 96 hpf |
|---|---|---|---|
| pharyngeal arch (contains muscle, early cartilage) | 2124 | 2814 | 329 |
| head mesenchyme (maybe ventral, hand2+) | 2907 | 1687 | 343 |
| jaw chondrocyte | 288 | 1379 | 24 |
| pharyngeal arch (early) | 2800 | 926 | 43 |
| connective tissue-meninges-dermal FB (late) | 6903 | 663 | 18 |
| chondrocranium | 263 | 490 | 141 |
| periocular mesenchyme | 616 | 433 | 76 |
| head mesenchyme (contains sclerotome) | 1727 | 409 | 41 |
| head mesenchyme/PA cartilage | 1125 | 232 | 18 |
| osteoblast | 92 | 186 | 1 |
| cranial neural crest (ntm, foxd3+) | 492 | 176 | 28 |
| head/eye connective tissue | 179 | 52 | 17 |

  Also: `mineralized tissue, bone` 17 cells and `odontoblast` 56 (all timepoints).

**Fabian 3 dpf (`~/Desktop/Hayden/git_repos/sox10_sc_analysis/data/raw/3dpf_scRNA.rds`)**
- Seurat v3 object; needs `UpdateSeuratObject`.
- 6,314 cells × 32,524 genes; median 3,398 UMI / 967 genes.
- Metadata is only `RNA_snn_res.0.8` and `bio.rep` (1–3). There are no cell-type labels.
- 97% of ob334 genes are present. Use the intersected gene universe.

**Fabian 5 dpf:** `.../data/raw/5dpf_scRNA.rds`, with OB cluster 14.

## Where the code goes
**Part 1 (pools, trajectory, panel) goes in this repo**, following its pattern: `R/` helpers, numbered `analysis/` steps, outputs under gitignored `results/`.
```
config/lineage.yml            paths (ZSCAPE, Fabian 3/5 dpf, scrna_matching repo, atlas gene map) + module gene lists
env/pyproject.toml            uv project: python 3.11, scanpy, anndata, cellrank>=2, spapros, pyyaml, matplotlib
                              (if spapros won't resolve, use a second uv project env_spapros/ with python 3.10)
R/lineage.R                   config loader, pool rules, module scoring, contaminant-cluster filter, mtx export
analysis/10_build_pools.Rmd   params: dataset = zscape | fabian3 | fabian5 -> pool + html report
analysis/11_export_pools.R    pools -> results/lineage/export/<set>/{counts.mtx, genes.csv, cells.csv} (same format as export_reference.R)
py/common.py                  config/IO (mtx -> AnnData with counts layer + log-norm X, as common.load_cells), module scores, kNN smoothing
py/20_trajectory.py           --dataset zscape | fabian : states, CellRank fates, OB-lineage drivers -> results/lineage/<ds>_states.h5ad
py/30_spapros_panel.py        per-dataset spapros, consensus, evaluation -> results/lineage/panel/
py/31_export_tangram_sets.py  Tangram reference sets as h5ad (raw counts) for Part 2
analysis/run_lineage.sh       10 -> 11 -> 20 -> 30 -> 31
```

**Part 2 (Tangram LOO) goes in `voxel_measurements_and_analysis/scrna_matching/iter4_lineage_reference/`.** It is a new iteration folder following that repo's conventions:
- scripts plus a `results/` README;
- large files under `P.iter_share_dir("iter4_lineage_reference")`;
- new files only.

It imports `common.py`, `paths.py` and `run_check.py` helpers and does not change them. Add its row to the progress table in `scrna_matching/README.md`.

## Phase 0: Setup (gate: confirm with the user before installing)
1. Write `config/lineage.yml` with:
   - the paths listed above, plus `scrna_matching_repo` and the atlas gene map (`iter0_data_prep/results/atlas_gene_map.csv`);
   - panel genes in scRNA symbols, via `ATLAS_TO_SC`.
2. Set up `env/` with `uv` (the user's convention) and smoke-test `import scanpy, cellrank, spapros`.
3. Write the module gene lists into the YAML, used by both R and Python.
   - Log genes missing from each dataset. Never fail silently.
   - Positional lists are descriptive priors to check against the literature, not filters.

   | Module | Genes |
   |---|---|
   | `ob_mature` | sp7, col10a1a, ifitm5, spp1, bglap, entpd5a, panx3, slc8a4b, cd81b |
   | `ob_early` | runx2a, runx2b, col1a1a, col1a2, sparc, alpl |
   | `nc_mesenchyme` | twist1a, twist1b, prrx1a, prrx1b, pdgfra, nr2f2, nr2f5, msx1b, dlx2a |
   | `chondro` | sox9a, col2a1a, matn1, col9a1a, acana, ucmaa |
   | `ligament_tendon` | scxa, tnmd, mkxa, col12a1a, matn4 |
   | `fibro` | dcn, lum, col6a1, col6a2 |
   | `proliferation` | mki67, pcna, top2a |

   Contaminant modules:

   | Module | Genes |
   |---|---|
   | `muscle` | myod1, myog, mylpfa, acta1a, ttn.1, ttn.2 |
   | `epithelial` | krt4, krt5, epcam, cdh1 |
   | `endothelial` | kdrl, cdh5. Not fli1a, which is also NC mesenchyme. |
   | `blood` | hbaa1, hbba1 |
   | `pigment` | mitfa, pmela, tyrp1b |
   | `glia_neural` | plp1b, mbpa, elavl3 |

   Positional priors:

   | Module | Genes |
   |---|---|
   | `hox2` | hoxa2b, hoxb2a |
   | `hox3` | hoxa3a, hoxb3a, hoxd3a |
   | `dorsal` | jag1b |
   | `intermediate` | dlx3b, dlx5a, dlx6a, msxe |
   | `ventral` | hand2 |

## Phase 1: Pools (R, `analysis/10_build_pools.Rmd`)
**`zscape`**
1. Keep cells with `gene_target != "tbx16-mut"` and timepoint 48, 72 or 96, whose `cell_type_sub` is in this INCLUDE list:
   - osteoblast
   - mineralized tissue, bone
   - jaw chondrocyte
   - chondrocranium
   - pharyngeal arch (contains muscle, early cartilage)
   - pharyngeal arch (early)
   - pharyngeal arch (NC-derived)
   - head mesenchyme/PA cartilage
   - head mesenchyme (maybe ventral, hand2+)
   - head mesenchyme (contains sclerotome)
   - cranial neural crest (ntm, foxd3+)
   - periocular mesenchyme
   - head/eye connective tissue
   - connective tissue-meninges-dermal FB (late)
   - odontoblast

   Exclude fin and pectoral types, cranial muscle types, mesodermal progenitors, RPE, epidermis and vasculature. Expect about 32k cells, with 72 hpf at about 9.5k.
2. Run `build_clean()`, then 2,000 vst HVGs, `embed_features()`, `choose_dims()`, `cluster_sweep()` and `pick_resolution()`, then UMAP.
3. Contaminant filter, on per-cluster mean `AddModuleScore` values z-scored across clusters:
   - drop a cluster when any contaminant z > 1.5 and the mesenchymal (nc_mesenchyme + ob_* + chondro) z < 0;
   - re-embed and check once more;
   - **never drop a cluster that has ≥ 20% annotated osteoblasts**;
   - log every removed cluster with its top 10 markers.

   The main expected removal is the muscle inside "pharyngeal arch (contains muscle…)".

**`fabian3`:** load with `UpdateSeuratObject` and keep `bio.rep`. Same embedding and filter; expect pigment and glia rather than muscle. Report OB-module detection per cluster. Expect it to confirm that no real OB cluster exists.

**`fabian5`:** same as `fabian3`. It is used **only** to give the Fabian trajectory an OB endpoint.
- Its OBs = cluster 14 / `5dpf_RNA_OB.rds` cells.
- Fabian 3 and 5 dpf come from the same lab and protocol, so merging them within Fabian is acceptable for trajectory only.
- If the user objects at the Phase 1 gate, run the Fabian trajectory on 3 dpf alone, with a manually set OB-like terminal state.

**Outputs** (`results/lineage/`):
- `{ds}_pool.rds`
- `pool_composition.csv` (label × timepoint, cluster × module z)
- `removed_clusters.csv`
- html reports

**Done when:**
- all 186 72 hpf OBs and all 92 48 hpf OBs are in the ZSCAPE pool;
- no remaining cluster meets the contaminant rule.

**Gate 1:** show the user the composition tables.

## Phase 2: States and trajectory (Python, `py/20_trajectory.py`; ZSCAPE and Fabian separately, same code)
1. Read the exports (mtx → AnnData):
   - `normalize_total(1e4)`, `log1p`, 2,000 HVGs and PCA(50);
   - neighbors with k = 30 for ZSCAPE (shallow) and k = 15 for Fabian;
   - leiden and diffmap.
2. **Root:** the cell with the highest `nc_mesenchyme − ob_mature` score among the earliest timepoint (ZSCAPE 48 hpf; Fabian 3 dpf). Then compute DPT.
   - Neither dataset has spliced/unspliced layers, so there is no RNA velocity.
3. **CellRank 2:**
   - Primary kernel: `0.8 * PseudotimeKernel(time_key="dpt_pseudotime") + 0.2 * ConnectivityKernel`.
   - Sensitivity kernel: `0.8 * CytoTRACEKernel + 0.2 * ConnectivityKernel`.
   - GPCCA: sweep n_states from 4 to 12, choose by the eigengap, and name terminal states by their dominant module.
   - If no OB terminal state appears, set the top `ob_mature` cells as one with `set_terminal_states`.
   - Compute fate probabilities, then `compute_lineage_drivers(lineages=["OB"])` and save the top 200 with q-values.
4. **`obs["state"]`** (same rubric in both datasets):
   - `OB_mature`: in the OB terminal state.
   - `OB_committed`: OB fate ≥ 0.5 and `ob_early` z > 0.
   - `OB_progenitor`: OB fate is the largest fate and ≥ 0.4.
   - `uncommitted_mesenchyme`: max fate < 0.4.
   - `chondro_lineage`, `ligament_lineage`, `fibro_lineage`, …: max fate ≥ 0.4 to that state.
   - Report how membership changes at thresholds 0.3 and 0.5.
   - **`obs["state_fine"]`** = leiden subclusters within each state, with at least 15 cells (merge smaller ones). spapros uses this.
5. **Positional (descriptive only):** compute kNN-smoothed (15 NN) positional module scores and report, per `state_fine`:
   - detection rates of hox2, hox3 and hand2;
   - any subclusters with a distinct positional signature. Flag them as opercle-candidate.
6. **Cross-dataset concordance (no integration):**
   - Pseudobulk each `state` in each dataset.
   - Correlate state-vs-rest logFCs over shared genes and write `state_concordance.csv`.
   - Matching states should be on the diagonal.

**Done when:**
- an OB terminal state exists in both datasets;
- the order of peak expression along OB-lineage pseudotime is runx2 → sp7 → col10a1a/ifitm5;
- OB-lineage membership Jaccard between the two kernels is reported. Target ≥ 0.6; if it is lower, report it rather than tune it away;
- in ZSCAPE, 48 hpf cells are enriched early in pseudotime (sanity check only).

**Gate 2:** show the user the terminal states and the state table.

## Phase 3: Panel (Python, `py/30_spapros_panel.py`)
1. **Gene universe:**
   - Start from genes shared by ZSCAPE and Fabian 3 dpf.
   - Remove `^mt-`, `^rp[sl]`, hemoglobins and contaminant-module genes.
   - **Trainability filter:** a gene must be detected in ≥ 5% of the ZSCAPE 72 hpf pool cells of at least one `state_fine`. This is the same 5% rule used for Tangram training genes, and it reflects that detection drives add-back gain.
   - Also keep a "Fabian-only detectable" list for reporting.
2. **Cells:**
   - ZSCAPE: the 72 hpf refined pool, all states (≈ 9k cells). Neighbour states such as hyomandibular cartilage and ligament are spatially next to the opercle, so they stay.
   - Fabian: the 3 dpf refined pool.
3. **spapros on each dataset separately:**
   - `ProbesetSelector(adata, celltype_key="state_fine", n=..., preselected_genes=<29 panel genes present>, genes_key=<8k HVG ∪ panel>, save_dir=...)`, then `.select_probeset()`.
   - Check in the installed docs whether `n` counts the preselected genes, and set `n` so that about 70 new genes are chosen.
   - Keep the default selection procedure and pass **no** `prior_genes`, so the selection stays unbiased.
4. **Consensus to 60 new genes** (stay within 50–70):
   - tier 1: selected in both datasets;
   - tier 2: selected in one dataset and detected ≥ 5% in the concordant state of the other, ordered by mean spapros rank;
   - fill to 60 alternating between the datasets, ZSCAPE first because it is the Tangram reference;
   - break ties by ZSCAPE detection.
5. **Evaluation** (`ProbesetEvaluator`, scheme "full", on each dataset). Compare five sets:
   - panel-only (29);
   - panel + consensus;
   - panel + 60 random genes (5 draws, detection-matched);
   - panel + top-60 scanpy DE markers;
   - panel + top-60 HVGs.

   Report forest-classifier F1 overall and **per state**, knn_overlap and cluster_similarity.
6. **Write `results/lineage/panel/hcr_panel_v1.csv`** with these columns:
   - gene, tier, source, states marked;
   - spapros rank in ZSCAPE and in Fabian;
   - % detected per state in each dataset;
   - ZSCAPE 72 hpf OB detection;
   - OB-lineage driver rank;
   - flags for the existing nominee lists and `positional_genes()`;
   - notes. Flag genes that are broadly high in every state, and genes that also mark muscle.

**Done when:**
- the panel has 50–70 genes, none of them already in the panel;
- panel + new beats panel-only and the random sets on overall F1 in both datasets;
- per-state F1 for `OB_progenitor`, `uncommitted_mesenchyme` and the neighbour states improves clearly. That gap puts a number on the original problem.

**Part 1 → Part 2 handoff (`py/31_export_tangram_sets.py`):** write raw-count h5ad files in the `iter0` format (`counts` in X; obs with embryo, gene_target, timepoint, n_umi, state, state_fine, fate probabilities) for:
- `lineage_z72`: ZSCAPE 72 hpf refined pool;
- `lineage_z48_96`: ZSCAPE 48–96 hpf refined pool;
- `lineage_f3`: Fabian 3 dpf refined pool;
- `lineage_z72_sub`: a stratified subsample (below).

## Phase 4: Tangram LOO (`scrna_matching/iter4_lineage_reference/`, using `scrna_matching/env`)
Use the same 28 training genes and hold out ihha, as in Iteration 1. Seed 1, 1000 epochs, lr 0.1.

1. **`tangram_ext.py`:** a small wrapper modelled on `common.tangram_map`:
   - adds `mode` (cells | constrained | clusters), `target_count`, `lambda_f_reg`, `lambda_count` and `cluster_label`;
   - records the effective lambdas for **both** `Mapper` and `MapperConstrained`;
   - returns `F_out` (the filter) in constrained mode.

   Do not edit `common.py`.
2. **Baseline reproduction (gate):**
   - Rerun Iteration 1's LOO (`zscape_ob72`, cells mode) with seed 1.
   - Every gene must have |Δ Pearson| < 0.03 against `iter1_ob_matching/results/tg_loo.csv`, and the median must be ≈ 0.36.
   - If it fails, stop and report.
3. **Cost control:** time one full mapping per arm first and extrapolate to 28 LOO mappings. Run arms in parallel processes, as in `iter2_pool_impute`.
   - `lineage_z72_sub`: all OBs plus up to 150 cells per `state_fine`, a fixed random sample, about 3k cells. It is the default reference for the big arms.
   - Run the full pool only if the timing allows.
4. **Arms** (target = the opercle region, 33,451 voxels, unless stated):

   | Arm | Reference | Mode | Tests |
   |---|---|---|---|
   | A0 | zscape_ob72 | cells | baseline |
   | A1 | lineage_z72_sub | constrained, `target_count` ∈ {186, 500, 1000} | expansion and filter (main arm) |
   | A2 | lineage_z72_sub | cells | adding cells without a filter: are non-opercle cells forced in? |
   | A3 | lineage_z72_sub | clusters (`state_fine`) | cheap, robust comparator |
   | A4 | lineage_z48_96 subsample | constrained | time pooling. The user expects this to hurt; test it. |
   | A5 | lineage_f3 subsample | constrained | Fabian 3 dpf, exploratory (no OBs) |
   | W1 | best of A1/A3 | same | **whole shared mask** (`RUN_DIR/voxel_matrix.csv.gz`), rebinned with `spatial_stats/common.rebin` to about 20–40k voxels; record the factors |
   | W0 | zscape_ob72 | cells | same whole-mask grid; the baseline for W1 |

5. **Metrics** for every arm, using `match_tangram.py` and `run_stage.py` logic:
   - per-gene LOO Pearson;
   - 200-permutation null and genes beating it;
   - markers-top-correct vs null;
   - ihha;
   - rim artifact (`run_check.rim_artifact`: rim call vs UMI and batch), for OB cells only.

   **Peripheral metrics.** The atlas has only about two rim-enriched genes (runx2a, bcl2a), so gene strata alone are too thin. Two primary metrics:
   - **Voxel-band LOO Pearson:** each gene's LOO Pearson computed within E_rim voxels only and within A–D voxels only.
   - **For W arms:** the same, within non-opercle tissue voxels vs opercle voxels.

   Secondary: genes stratified by the rim ratio from `subdomain_profiles_mean.csv` (mean_E / mean_A–D; top vs bottom tertile).
6. **Diagnostics:**
   - state composition of the cells kept by the filter (F_out > 0.5) in A1/A4/A5: which states the atlas region "wants";
   - mapping-mass share per state per subdomain;
   - predicted vs measured mean-over-slices maps for runx2a, bcl2a, sdc2, klf13, mfsd2b and ihha;
   - projection of the 60 panel genes from the best arm as an "expected pattern" sheet for HCR planning, labelled low confidence in the periphery.

**Outputs:**
- `iter4_lineage_reference/results/{gate_summary.csv, loo_by_arm.csv, band_loo.csv, filter_composition.csv, README.md, *.png}`
- mappings on the share

**Done when:**
- the baseline is reproduced;
- every arm has complete LOO;
- the gate table compares each arm with A0 (and W1 with W0) on: median and band Pearson, a paired Wilcoxon test over genes for the E_rim band, central-gene Δ (criterion: does not worsen by more than 0.02), ihha and the rim artifact.

Write the interpretation whichever way the result goes. If expansion doesn't help the rim band, the missing spatial anchors are the bottleneck, and the Phase 3 panel is the next experiment.

## Overall definition of done
1. ZSCAPE and Fabian pools are built and documented, and Gate 1 has been reviewed.
2. States and OB-lineage trajectories are done in both datasets, with the checks above, and Gate 2 has been reviewed.
3. `hcr_panel_v1.csv` holds 50–70 new genes with an evaluation table showing gains over panel-only and random.
4. The Tangram LOO gate table covers all arms, with peripheral (band) metrics and a written interpretation.
5. Write-ups:
   - a "Lineage expansion" section in `identifying_srgs/README.md` (Results style, honest caveats);
   - `iter4_lineage_reference/results/README.md`;
   - a progress row in `scrna_matching/README.md`.
6. `analysis/run_lineage.sh` and the iter4 run commands regenerate everything.
7. **Commit only when the user asks** (in both repos).

## Verification
- Part 1: run `bash analysis/run_lineage.sh` from the `identifying_srgs` root. It renders `10_build_pools.Rmd` per dataset, exports, then runs `uv run --project env python py/20_trajectory.py --dataset {zscape,fabian}`, then `30_spapros_panel.py` and `31_export_tangram_sets.py`.
- Part 2: from the `voxel_measurements_and_analysis` root, `uv run --project scrna_matching/env python scrna_matching/iter4_lineage_reference/run_arms.py <arm>`, then `summarize.py`.
- Built-in asserts:
  - OB counts (186 / 92) are present in the pool;
  - no NA `state`;
  - the panel has 50–70 genes, none already in the panel, and all pass the trainability filter;
  - the baseline matches `tg_loo.csv`;
  - a same-seed rerun of Phase 1 gives clustering ARI = 1;
  - the h5ad exports reproduce R-side detection for the panel genes (as `check_exports.py` does).
- Before the write-up, list every Done criterion and say which ones passed and which failed.

## Risks and caveats to carry into the write-up
- **ZSCAPE depth (~360 UMI):** positional TFs drop out. Positional calls are hypotheses.
- **Other condensations in the pool:** cleithrum, parasphenoid, PA1 and ventral hyoid cells will be in it. Constrained Tangram, the whole-mask arm and HCR are what finally resolve this.
- **Forced mapping:** mapping non-OB cells onto the opercle-only region forces them in. This is why the constrained and whole-mask arms exist; A2 shows the size of the problem.
- **Batch:** 81% of the OBs come from one plate. The rim artifact must not get worse in the expanded arms.
- **Fabian 3 dpf lacks OBs:** its value is in the progenitor and neighbour states, which are exactly the missing ones. It is not a good OB reference.
