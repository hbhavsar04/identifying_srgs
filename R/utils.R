# Shared paths, parameters and helpers for the osteoblast subclustering comparison.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---- Paths (from CLAUDE.md; read-only) ---------------------------------------
DATA_DIR  <- "/Volumes/kwongroup/Hayden Bhavsar/R files/RDS_seurat_objects"
OB_PATH   <- file.path(DATA_DIR, "ob_subset.rds")
FULL_PATH <- file.path(DATA_DIR, "saunders.rds")

find_proj_dir <- function(start = getwd()) {
  d <- normalizePath(start)
  while (!file.exists(file.path(d, "CLAUDE.md"))) {
    parent <- dirname(d)
    if (parent == d) stop("Could not find project root (CLAUDE.md) above ", start)
    d <- parent
  }
  d
}
PROJ_DIR    <- find_proj_dir()
RESULTS_DIR <- file.path(PROJ_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Shared analysis parameters ----------------------------------------------
SEED        <- 321
N_PCS       <- 30
K_PARAM     <- 15
RESOLUTIONS <- seq(0.1, 1.2, by = 0.1)
MIN_CLUSTER <- 15

META_KEEP <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent_mito", "timepoint",
               "cell_type_sub", "cell_type_broad", "tissue", "gene_target", "embryo",
               "hash_umis", "log_n_umi", "seurat_clusters")

OB_MARKERS <- c("sp7", "runx2a", "runx2b", "col10a1a", "entpd5a", "spp1", "bglap", "ifitm5",
                "col1a1a", "sparc")

# Cell sets compared side by side:
#   ob334  = ob_subset.rds (cell_type_sub == "osteoblast"), 334 cells
#   obfin  = all osteoblast-annotated cells in saunders.rds, i.e. "osteoblast" +
#            "pectoral fin bud (osteoblast+)", ~6,083 cells
CELLSETS <- c("ob334", "obfin")
OBFIN_TYPES <- c("osteoblast", "pectoral fin bud (osteoblast+)")

cellset_path <- function(cellset) file.path(RESULTS_DIR, paste0("cellset_", cellset, ".rds"))
cellset_dir  <- function(cellset) {
  d <- file.path(RESULTS_DIR, cellset)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# kNN size scaled to dataset size (Seurat default 20; smaller for tiny sets).
k_for <- function(n) if (n < 1000) 15L else 20L

# ---- Preprocessing ------------------------------------------------------------
#' Build a clean, normalized Seurat object from any Seurat object's counts.
#' Drops stale HVGs/scale.data/reductions; keeps atlas clusters as `atlas_cluster`.
build_clean <- function(ob, min_cells = 3) {
  counts <- LayerData(ob, assay = "RNA", layer = "counts")
  meta <- ob@meta.data[, intersect(META_KEEP, colnames(ob@meta.data)), drop = FALSE]
  names(meta)[names(meta) == "seurat_clusters"] <- "atlas_cluster"
  meta$timepoint <- factor(meta$timepoint)
  obj <- CreateSeuratObject(counts, meta.data = meta, min.cells = min_cells)
  obj <- NormalizeData(obj, verbose = FALSE)
  stopifnot(ncol(obj) == ncol(ob))
  obj
}

#' Load a prepared cell set (built by analysis/00b_build_cellsets.R).
load_cellset <- function(cellset) {
  cellset <- match.arg(cellset, CELLSETS)
  p <- cellset_path(cellset)
  if (!file.exists(p)) stop("Missing ", p, "; run analysis/00b_build_cellsets.R first")
  obj <- readRDS(p)
  # Some atlas cells lack timepoint/embryo; keep them as an explicit "unknown" level.
  # Levels are sorted numerically when possible (timepoints), with "unknown" last.
  for (v in c("timepoint", "embryo", "gene_target", "orig.ident", "cellset_label")) {
    x <- as.character(obj[[v]][, 1])
    x[is.na(x)] <- "unknown"
    known <- unique(x[x != "unknown"])
    num <- suppressWarnings(as.numeric(known))
    known <- if (!anyNA(num)) known[order(num)] else sort(known)
    obj[[v]] <- factor(x, levels = c(known, if ("unknown" %in% x) "unknown"))
  }
  obj
}

#' ARI / NMI that tolerate character vectors and NAs (aricode needs factors without NA).
as_lab <- function(x) { x <- as.character(x); x[is.na(x)] <- "unknown"; factor(x) }
safe_ari <- function(a, b) aricode::ARI(as_lab(a), as_lab(b))
safe_nmi <- function(a, b) aricode::NMI(as_lab(a), as_lab(b))

#' Pick Festem genes from feature metadata. Uses FDR-significant genes (EM > 0) if there
#' are at least `min_genes`; otherwise relaxes to `fdr_relaxed`; otherwise takes the top
#' `top_n` by Festem_rank (significant genes first, then EM > 0 ordered by EM).
select_festem_genes <- function(fe, fdr = 0.05, fdr_relaxed = 0.1, min_genes = 100,
                                top_n = 500) {
  fm <- fe[["RNA"]][[]]
  sig <- function(q) rownames(fm)[which(fm$Festem_p_adj < q & fm$Festem_EM > 0)]
  g <- sig(fdr)
  if (length(g) >= min_genes) return(list(genes = g, rule = sprintf("FDR < %.2f", fdr)))
  g <- sig(fdr_relaxed)
  if (length(g) >= min_genes) return(list(genes = g, rule = sprintf("FDR < %.2f", fdr_relaxed)))
  ranked <- rownames(fm)[order(fm$Festem_rank, na.last = NA)]
  ranked <- ranked[which(fm[ranked, "Festem_EM"] > 0)]
  list(genes = head(ranked, top_n),
       rule = sprintf("top %d by Festem_rank (only %d genes at FDR < %.2f)",
                      min(top_n, length(ranked)), length(sig(fdr)), fdr))
}

#' Run Festem with a metadata column as prior (vector priors hit a bug in Festem 1.2.3).
run_festem <- function(obj, prior_col, threads = 8) {
  G <- nlevels(droplevels(factor(obj[[prior_col]][, 1])))
  suppressWarnings(Festem::RunFestem(obj, prior = prior_col, G = G, FDR_level = 0.05,
                                     num.threads = threads, seed = SEED))
}

#' Scale, PCA, neighbors on a given feature set (features already set as VariableFeatures).
embed_features <- function(obj, features = VariableFeatures(obj), npcs = N_PCS) {
  set.seed(SEED)
  obj <- ScaleData(obj, features = features, verbose = FALSE)
  obj <- RunPCA(obj, features = features, npcs = npcs, seed.use = SEED, verbose = FALSE)
  obj
}

#' Run FindClusters at several resolutions; returns object with RNA_snn_res.* columns.
cluster_sweep <- function(obj, dims, resolutions = RESOLUTIONS, k = k_for(ncol(obj)),
                          reduction = "pca", prefix = "res") {
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims, k.param = k, verbose = FALSE)
  for (r in resolutions) {
    obj <- FindClusters(obj, resolution = r, random.seed = SEED, verbose = FALSE)
    obj[[paste0(prefix, "_", r)]] <- Idents(obj)
  }
  obj
}

#' Summarise a resolution sweep: number of clusters and smallest cluster size.
sweep_summary <- function(obj, prefix = "res") {
  cols <- grep(paste0("^", prefix, "_[0-9.]+$"), colnames(obj@meta.data), value = TRUE)
  data.frame(
    column     = cols,
    resolution = as.numeric(sub(paste0("^", prefix, "_"), "", cols)),
    n_clusters = vapply(cols, function(cc) nlevels(droplevels(obj[[cc]][, 1])), integer(1)),
    min_size   = vapply(cols, function(cc) min(table(obj[[cc]][, 1])), integer(1)),
    row.names  = NULL
  )
}

#' Pick the highest resolution whose smallest cluster still has >= min_size cells
#' (and at least 2 clusters). Falls back to the lowest resolution with >= 2 clusters.
pick_resolution <- function(summ, min_size = MIN_CLUSTER) {
  ok <- summ[summ$n_clusters >= 2 & summ$min_size >= min_size, ]
  if (nrow(ok)) return(ok$resolution[which.max(ok$resolution)])
  two <- summ[summ$n_clusters >= 2, ]
  if (nrow(two)) return(min(two$resolution))
  min(summ$resolution)
}

#' Resolution giving exactly K clusters (closest K if none), for matched comparisons.
resolution_for_k <- function(summ, K) {
  hit <- summ[summ$n_clusters == K, ]
  if (nrow(hit)) return(min(hit$resolution))
  summ$resolution[which.min(abs(summ$n_clusters - K))]
}

#' Choose PCs: smallest n where cumulative variance gain per extra PC drops below 0.1%.
choose_dims <- function(obj, reduction = "pca", min_dims = 8, max_dims = 20) {
  sdev <- Stdev(obj, reduction = reduction)
  pct <- sdev^2 / sum(sdev^2) * 100
  n <- which(diff(pct) > -0.1)[1]
  if (is.na(n)) n <- max_dims
  min(max(n, min_dims), max_dims)
}

# ---- Comparison helpers -------------------------------------------------------
#' Subsampling stability: recluster 80% subsamples and compute ARI to full labels.
stability_ari <- function(obj, features, dims, resolution, labels,
                          n_iter = 50, frac = 0.8, k = k_for(ncol(obj))) {
  set.seed(SEED)
  vapply(seq_len(n_iter), function(i) {
    cells <- sample(colnames(obj), floor(frac * ncol(obj)))
    sub <- suppressWarnings(subset(obj, cells = cells))
    sub <- ScaleData(sub, features = features, verbose = FALSE)
    sub <- RunPCA(sub, features = features, npcs = max(dims), seed.use = SEED + i, verbose = FALSE)
    sub <- FindNeighbors(sub, dims = dims, k.param = k, verbose = FALSE)
    sub <- FindClusters(sub, resolution = resolution, random.seed = SEED + i, verbose = FALSE)
    safe_ari(Idents(sub), labels[colnames(sub)])
  }, numeric(1))
}

#' Cramér's V between two categorical vectors.
cramers_v <- function(x, y) {
  tab <- table(x, y)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (min(dim(tab)) < 2) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
}

#' Mean silhouette width of labels in an embedding.
mean_silhouette <- function(labels, emb) {
  lab <- as.integer(factor(labels))
  if (length(unique(lab)) < 2) return(NA_real_)
  s <- cluster::silhouette(lab, dist(emb))
  mean(s[, "sil_width"])
}

#' Positional / regional gene panel present in the object (zebrafish symbols).
positional_genes <- function(obj) {
  g <- rownames(obj)
  pats <- c("^hox[abcd][0-9]+[ab]?a?$", "^dlx[1-6][ab]?$", "^hand2$", "^barx1$",
            "^pax[0-9]+[ab]?$", "^six[0-9]+[ab]?$", "^tbx[0-9]+[ab]?$", "^msx[0-9][a-e]?$",
            "^shha?$", "^en[12][ab]?$", "^pitx[12]$", "^nr2f[0-9]$", "^edn1$", "^jag1b$",
            "^irx[0-9][ab]?$", "^fgf[0-9]+[ab]?$", "^lhx[0-9][ab]?$")
  sort(unique(unlist(lapply(pats, function(p) grep(p, g, value = TRUE)))))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
