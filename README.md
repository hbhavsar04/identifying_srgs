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
