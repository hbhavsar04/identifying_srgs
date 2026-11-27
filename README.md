# identifying_srgs

Searching for **spatially restricted genes (SRGs)** in zebrafish osteoblasts from the ZSCAPE
atlas (Saunders et al.). The question is whether osteoblasts split into subclusters that reflect
where they sit in the skeleton, and which genes mark those subclusters.

The analysis compares two ways of choosing the genes that drive subclustering:

- **Seurat (HVG) arm:** the standard workflow. Pick the top 2,000 variable genes, run PCA, build a kNN graph, cluster with Louvain.
- **Festem arm:** [Festem](https://github.com/Edward-Z-Chen/Festem) chooses genes with an EM test for cluster-informative expression. The cells are then reclustered on those genes only. Both arms use the Seurat clusters as Festem's prior, so they start from the same place.

## Cell sets

The pipeline runs on two cell sets so that the results can be compared:

| Cell set | Source | Cells | Genes | Timepoints (hpf) |
|---|---|---|---|---|
| `ob334` | `ob_subset.rds` (`cell_type_sub == "osteoblast"`) | 334 | 5,772 | 24–96, mostly 48 and 72 |
| `obfin` | every cell in `saunders.rds` annotated `osteoblast` or `pectoral fin bud (osteoblast+)` | 6,083 | 18,837 | 18–96, plus 774 with no timepoint |

`obfin` contains all of `ob334` plus 5,749 pectoral fin bud cells.

## Pipeline

```
analysis/00_setup.R            install Festem, clustree, aricode, UpSetR, ...
analysis/00b_build_cellsets.R  build both cell sets from the RDS files -> results/cellset_*.rds
analysis/01_seurat_subcluster.Rmd   HVG arm: PCA, resolution sweep, markers, HVG-count sensitivity
analysis/02_festem_subcluster.Rmd   Festem arm: gene selection, prior sensitivity, recluster, markers
analysis/03_compare.Rmd             side-by-side comparison and SRG candidate table
analysis/run_all.R             render 01-03 for each cell set -> results/<cellset>/*.html
R/utils.R                      shared paths, parameters and helpers
```

Run everything from the repo root (the `/Volumes/kwongroup` share must be mounted; it is only read from):

```sh
Rscript analysis/00_setup.R
Rscript analysis/00b_build_cellsets.R
Rscript analysis/run_all.R            # or: Rscript analysis/run_all.R ob334
```

Shared settings (in `R/utils.R`): seed 321; 2,000 HVGs (vst); a resolution sweep from 0.1 to 1.2;
kNN k = 15 for fewer than 1,000 cells and 20 otherwise. The number of PCs is picked by the elbow rule,
and both cell sets ended up with 8. The script picks the highest resolution at which every cluster
still has at least 15 cells. Festem genes are those at FDR < 0.05 with EM > 0. If fewer than 100 genes
pass, the threshold is relaxed, and failing that, the top genes by Festem rank are used.

`03_compare.Rmd` scores the two arms on:

- overlap between the feature sets;
- correlation of the top PCs with sequencing depth and cell cycle;
- ARI/NMI agreement at matched K;
- silhouette width, in each arm's own PCA and in a shared PCA;
- marker count and specificity;
- stability when reclustering 50 random 80% subsamples;
- Cramér's V of the clusters against timepoint, embryo, gene target and sample;
- a dot plot of osteoblast and positional (hox/dlx/pax/tbx/…) genes.

It ends with the SRG candidate table.

## Results

Outputs are written to `results/<cellset>/` (gitignored): the HTML reports, the marker tables
(`markers_seurat.csv`, `markers_festem.csv`), Festem gene ranks (`festem_genes.csv`),
`partition_agreement.csv`, `cluster_quality.csv` and `srg_candidates.csv`.

### Headline numbers

| | ob334 HVG | ob334 Festem | obfin HVG | obfin Festem |
|---|---|---|---|---|
| Genes used | 2,000 | 298 ¹ | 2,000 | 164 (FDR < 0.05) |
| Clusters (K) | 5 | 7 | 15 | 14 |
| Silhouette, own PCA | 0.093 | 0.110 | 0.107 | 0.157 |
| Silhouette, shared PCA | 0.101 | −0.028 | 0.104 | 0.086 |
| Significant markers | 90 | 95 | 4,008 | 3,777 |
| Median specificity (top 10 per cluster) | 0.31 | 0.36 | 0.36 | 0.42 |
| Subsample stability, median ARI | **0.62** | 0.47 | 0.57 | **0.65** |
| Cramér's V vs timepoint | 0.42 | 0.37 | 0.31 | 0.32 |
| Cramér's V vs embryo | 0.61 | 0.62 | 0.49 | 0.50 |

¹ Only 1 gene reached FDR < 0.05 in `ob334`, so the top 298 genes by Festem rank were used.

### ob334 (334 osteoblasts)

- **Festem has very little power here.** With 334 cells only one gene passed FDR < 0.05, so the Festem gene set is a ranked list and not a set of significant genes. It changed a lot with the prior: the K=5 and K=2 priors gave gene sets with a Jaccard overlap of only 0.25.
- **The two arms disagree.** ARI between them was 0.27–0.43 at K = 2–5. Festem's 298 genes share 124 genes with the HVGs (Jaccard 0.06). None of the positional-panel genes is in the Festem set. The HVG set does contain some, including dlx3b, dlx5a, hoxb3a, hoxc3a, pax3a, pitx2, tbx1 and tbx4.
- **Neither partition is clearly better.** Festem clusters separate slightly better in their own PCA and have somewhat more specific markers. In the shared PCA they do not separate (silhouette −0.03), and they are less stable than the HVG clusters (median ARI 0.47 vs 0.62).
- **There is strong confounding.** Both arms' clusters track embryo (V ≈ 0.6) and timepoint (V ≈ 0.4).
- **SRG candidates:** 34 genes across 5 Festem clusters. Clusters 0 and 1 (14 genes) are flagged as 72 hpf-driven. In those, enam, phex, col1a1b, fgfr3 and abcc6a look like mineralization or maturation markers and not positional ones. The 20 candidates in unflagged clusters are all from 48 hpf-dominated clusters:
  - cluster 3: usta, tanc1b, slc8a4b, smpd3, kremen1, …
  - cluster 5: plod2, bicc1a, samd10a, …
  - cluster 6: col12a1a/b, tnc, col5a1, col5a2a, thbs4b, ror1, matn4. This looks like an ECM or perichondral-type signature.

### obfin (6,083 osteoblast-annotated cells)

- **Festem has enough power.** 164 genes passed FDR < 0.05, though the gene set still depends on the prior (Jaccard 0.22 between the K=15 and K=3 priors).
- **The arms agree at low K** (ARI 0.79 at K=3, 0.62 at K=5) and less so at higher K (0.44–0.49 at K=6–7). Festem clusters separate better in their own PCA (0.157 vs 0.107), have more specific markers, and are more stable (median ARI 0.65 vs 0.57).
- **The main split is cell type, not position.** Both arms put all 334 `osteoblast` cells into one cluster and the fin-bud cells into the rest (Cramér's V vs label = 1.0). In other words, the `obfin` subclusters mostly divide up the pectoral fin bud by stage (clusters are dominated by 18, 24 or 36 hpf). They do not divide the osteoblasts into regions.
- **The positional signal is weak.** The HVGs include many patterning genes, such as hand2, tbx5a, tbx15, lhx9, several fgf genes, hoxb genes, pax3/7 and six1a. Festem kept only one of them (jag1b).
- **SRG candidates:** 288 genes across 14 clusters, none flagged as timepoint- or embryo-driven. 85 of them are not HVGs. Most are fin-bud ECM and signalling genes, for example vcana, frem2a/3, fras1, lama5, lamb1a, notch3, adamts6, mmp9, wnt5b, rspo2 and hmcn1. The osteoblast cluster (cluster 11) gives col1a1b, col1a2, col11a1b, col11a2, col22a1, cdh11, mef2cb, ncam1a and nrp2a.

### Takeaways

1. With these annotations and this depth, neither method finds osteoblast subclusters that are clearly spatial. Most of the structure follows stage, embryo, or the osteoblast/fin-bud split.
2. Festem is useful once there are thousands of cells, where it gives more stable and more specific clusters from about 10× fewer genes. On the 334-cell subset it lacks the power to select genes reliably.
3. The candidate tables in `results/*/srg_candidates.csv` are hypotheses, not validated SRGs. The `ob334` unflagged clusters (3, 5 and 6) are the most direct leads for osteoblast heterogeneity. Validation needs spatial data, such as in situ hybridization or HCR. Next steps could be to restrict `obfin` to the osteoblast cells, to regress out or stratify by timepoint, or to bring in more osteoblasts from later stages.

## Lineage expansion (OB-lineage reference, HCR panel, Tangram test)

Plan: `prompts/ob_lineage_plan_v1.md` (run gate by gate, 2026-11-25/26). Goal: expand the 186-cell ZSCAPE
72 hpf OB reference to the skeletogenic neural-crest cells around it, find OB progenitor / neighbour
states, choose ~60 new HCR genes that resolve them, and test whether the expanded reference improves
the Tangram map of the opercle (Part 2, in `voxel_measurements_and_analysis/scrna_matching/iter4_lineage_reference/`).

```sh
bash analysis/run_lineage.sh     # Part 1, from the repo root (several hours, mostly GPCCA; uv env in env/, config in config/lineage.yml)
```

```
config/lineage.yml            paths, modules, filters, panel settings and every gate decision (with comments)
R/lineage.R                   pool building, module scoring, contaminant filter (sources R/utils.R)
analysis/10_build_pools.Rmd   Phase 1 pools per dataset (rendered by 10_render_pools.R)
analysis/11_export_pools.R    raw-count export for Python
py/20_trajectory.py           Phase 2: CellRank 2 states and OB fate (ZSCAPE and Fabian separately)
py/30_spapros_panel.py        Phase 3: spapros selection, consensus, DE comparators, evaluation, panel table
py/31_export_tangram_sets.py  Part 1 -> Part 2 handoff (h5ad reference sets on the share)
```

Outputs in `results/lineage/` (gitignored): pools, `trajectory/`, `panel/` (`hcr_panel_v1.csv`), `tangram_sets/`.
ZSCAPE and Fabian are never integrated, and nothing is imputed.

### Phase 1: pools

- ZSCAPE: skeletogenic NC pool at 48–96 hpf, 22,547 cells (48 hpf 12,634 / 72 hpf 8,823 / 96 hpf 1,090),
  21 clusters; the 186 (72 hpf) and 92 (48 hpf) annotated OBs are all present and sit in one cluster (96% OB).
  Fabian sox10 lineage: 3 dpf 6,314 cells, 5 dpf 9,742 cells (242 OBs, one 96%-OB cluster).
- The plan's relative contaminant rule removed real mesenchyme and never converged; with an added
  absolute detection floor (user decision) **no cluster was removed** in any dataset.
  "connective tissue-meninges-dermal FB (late)" was dropped from the ZSCAPE label list: it is basal
  epidermis (krt5, krtt1c19e, tp63). Fast-muscle transcripts are ambient across ZSCAPE; one muscle-ambient
  cluster (no myod1/myf5) was kept and flagged.
- **Fabian 3 dpf has no OB cluster** (no cluster has > 10% of cells with any mature-OB gene; runx2a ≤ 2%).

### Phase 2: states and OB fate (CellRank 2)

Kernel 0.8 pseudotime (DPT) + 0.2 connectivity; GPCCA macrostates; Fabian 3 + 5 dpf merged within Fabian.

- OB is a predicted terminal state in both datasets (alongside chondro ×2–3, fibro, mesenchyme, OB_early).
  The OB lineage is 99% annotated OBs (ZSCAPE 266 cells, Fabian 233). CytoTRACE-kernel sensitivity run:
  OB-lineage Jaccard 0.98 (ZSCAPE) and 0.88 (Fabian).
- **No OB progenitors outside the OB cluster** in either dataset: at the 0.4 fate threshold only 8 (ZSCAPE)
  and 19 (Fabian) cells are OB_progenitor; 76% / 68% of cells are uncommitted.
- runx2 → sp7 → col10a1a/ifitm5 ordering along pseudotime passes in Fabian and **fails in ZSCAPE**
  (median DPT 0.080 at 48 hpf vs 0.083 at 72 hpf: pseudotime barely separates the stages).
- State concordance between datasets (state-vs-rest logFC correlation): OB ≈ 0.55, chondro_2 0.68,
  uncommitted 0.60; fibro, OB_early and mesenchyme have no clear match.
- Gate 2 decision (user): head mesenchyme / pharyngeal arch are plausible OB sources, so the panel adds an
  **"OB_proximal" class = top 1% of non-OB cells by OB fate** (86 ZSCAPE 72 hpf, 63 Fabian 3 dpf cells).
  **This is a weakness**: it is a ranking, not a progenitor state, and there is no better option yet.

### Phase 3: HCR panel (`results/lineage/panel/hcr_panel_v1.csv`, 60 genes)

- Cells: ZSCAPE 72 hpf and Fabian 3 dpf pools; classes = fine CellRank states (+ OB_proximal). Universe:
  genes shared by both, minus mt/ribosomal/hemoglobin/contaminant genes, detected in ≥ 5% of one 72 hpf
  class (5,990 genes). The 29 panel genes are preselected.
- spapros (70 new genes per dataset, 10 shared) gave a 60-gene tiered consensus. It also picked
  dissociation-stress, housekeeping and mitochondrial genes (jun, junba, atf3, eef1g, NC-002333.4, …),
  excluded after selection; rhag (erythroid) was excluded on review.
- **A plain DE list beat spapros**, so the final panel is **DE_both** (user decision): wilcoxon DE per class,
  round-robin across classes, alternating the ZSCAPE and Fabian rankings. 25 of its 60 genes are in the
  spapros consensus.

Forest-classifier balanced accuracy (spapros `forest_clfs`: XGBoost, 5-fold × 5 seeds, mean per-class recall;
classes < 40 cells are not scored, so OB_progenitor is never evaluated):

| Gene set (panel + 60 new) | ZSCAPE 72 hpf | Fabian 3 dpf |
|---|---|---|
| panel only (29) | 0.277 | 0.381 |
| best of 5 detection-matched random | 0.451 | 0.512 |
| top-60 HVG | 0.442 | 0.547 |
| spapros consensus | 0.628 | 0.663 |
| DE chosen in the other dataset | 0.489 | 0.580 |
| **DE_both (final panel)** | **0.735** | **0.702** |

- All sets except "DE in the other dataset" were chosen on the cells they are scored on, so their
  numbers are optimistic; the cross-dataset row shows how much is lost in transfer.
- OB_proximal recall stays low (0.30 ZSCAPE, 0.24 Fabian).
- Panel content: matrix/cartilage genes (collagens, matrilins), OB-class genes (enam, col5a2b, plxna3, cntfr,
  ptprfb, kremen1), mesenchyme TFs (foxp4, pbx3b, lmx1ba/bb, foxf2a, her6) and sost/sfrp2/ucmab. Several
  are < 1% detected in ZSCAPE OBs (sost, sfrp2, ucmab, f13a1b); eight are unnamed (si:/zgc:/BX…) genes.
  hmgb2a (proliferation) and tmsb4x (near-ubiquitous) are kept but flagged in `notes`.

### Phase 4: Tangram with the lineage reference (stopped after the sub-sample arms)

Baseline reproduced (seed 1, max |Δ| 0.005). Adding the lineage cells **as cells** makes the opercle map
worse: LOO median 0.35 → 0.24–0.28 in cells and constrained (filter 186/500/1000) modes, E_rim band −0.08 to
−0.10, and ihha → D is lost. The filter keeps the same fraction of every state, so the opercle region shows
no preference for any state, OB_proximal included. Pooling 48–96 hpf hurts; Fabian 3 dpf fails. Clusters
mode (28 state profiles) scores 0.59 with gains in both bands, but sends ihha to the wrong subdomain and
puts ~80% of the mass on non-OB states; it most likely reflects averaging, and the OB-only clusters
control was not run. Details: `scrna_matching/iter4_lineage_reference/results/README.md`.

### Caveats

- ZSCAPE depth (~360 UMI): positional TFs drop out, and pseudotime hardly separates 48 from 72 hpf.
- OB_proximal is a fate ranking, not a validated progenitor state (see Phase 2).
- The pools contain other condensations (cleithrum, parasphenoid, PA1, ventral hyoid); the opercle-only
  target forces them in, which is one reason the expanded arms lose.
- 81% of the 72 hpf OBs come from one plate; the OB rim artifact disappears in the expanded arms only
  because non-OB cells take the rim calls.
- The panel was chosen to separate scRNA states, not for spatial information; whether it helps the map is
  untested until HCR (the add-back result says detection is what predicts gain).

### Done criteria (`prompts/ob_lineage_plan_v1.md`)

| # | Criterion | Status |
|---|---|---|
| 1 | Pools built and documented, Gate 1 reviewed | pass (no cluster removed; filter rule changed by the user) |
| 2 | States and OB trajectories in both datasets, checks, Gate 2 reviewed | pass, with failed checks reported: marker order fails in ZSCAPE; no progenitor state |
| 3 | 50–70 new genes beating panel-only and random | pass (60; DE_both 0.735 / 0.702 vs random 0.451 / 0.512). Per-state gain for OB_progenitor not testable (too few cells) |
| 4 | Tangram gate table for all arms, with band metrics and interpretation | partial: A0–A5 done; W0/W1 not run (stopped by the user) |
| 5 | Write-ups (this section, iter4 README, scrna_matching progress row) | done |
| 6 | `analysis/run_lineage.sh` and the iter4 commands regenerate everything | done (scripts written; not re-run end to end) |
| 7 | Commit only when asked | nothing committed |
