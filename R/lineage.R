# Helpers for the OB-lineage expansion (prompts/ob_lineage_plan_v1.md):
# config loader, pool rules, module scoring, contaminant-cluster filter, mtx export.
# Builds on R/utils.R (build_clean, embed_features, choose_dims, cluster_sweep, ...).

source(file.path(rprojroot::find_root("CLAUDE.md"), "R/utils.R"))
suppressPackageStartupMessages(library(Matrix))

# Extra metadata kept by build_clean(): Fabian replicate + original Fabian clusters, pool fields.
META_KEEP <- unique(c(META_KEEP, "bio.rep", "fabian_cluster", "pool_label", "ob_annot", "dataset"))

DATASETS <- c("zscape", "fabian3", "fabian5")

load_lineage_config <- function() {
  cfg <- yaml::read_yaml(file.path(PROJ_DIR, "config", "lineage.yml"))
  cfg$paths$results <- file.path(PROJ_DIR, cfg$paths$results)
  cfg
}

lineage_dir <- function(cfg, ...) {
  d <- file.path(cfg$paths$results, ...)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

#' Module gene lists restricted to genes in `genes`; logs every missing gene (never silent).
#' Returns list(modules = named list of present genes, missing = data.frame).
module_genes <- function(cfg, genes, dataset) {
  mods <- unlist(unname(cfg$modules), recursive = FALSE)
  miss <- do.call(rbind, lapply(names(mods), function(m) {
    g <- setdiff(mods[[m]], genes)
    if (length(g)) data.frame(dataset = dataset, module = m, gene = g) else NULL
  }))
  if (!is.null(miss) && nrow(miss))
    message(sprintf("[%s] module genes missing: %s", dataset,
                    paste(sprintf("%s:%s", miss$module, miss$gene), collapse = ", ")))
  list(modules = lapply(mods, intersect, genes),
       missing = if (is.null(miss)) data.frame(dataset = character(), module = character(),
                                               gene = character()) else miss)
}

# ---- Raw pools -------------------------------------------------------------------
#' Raw (unprocessed) cell pool for a dataset, with pool_label / ob_annot / timepoint set.
#' The ZSCAPE subset is cached in results/lineage/cache/ (the full atlas is 2 GB).
load_raw_pool <- function(ds, cfg = load_lineage_config()) {
  ds <- match.arg(ds, DATASETS)
  if (ds == "zscape") {
    pr <- cfg$pools$zscape
    cache <- file.path(lineage_dir(cfg, "cache"), "zscape_raw_pool.rds")
    if (file.exists(cache)) {
      obj <- readRDS(cache)
      # Cache is reusable when the current include list is a subset of the one it was built with.
      if (all(pr$include_cell_type_sub %in% obj@misc$include) &&
          identical(sort(obj@misc$timepoints), sort(pr$timepoints))) {
        keep <- obj$pool_label %in% pr$include_cell_type_sub
        message(sprintf("[zscape] cached raw pool: keeping %d of %d cells (current include list)",
                        sum(keep), ncol(obj)))
        return(if (all(keep)) obj else subset(obj, cells = colnames(obj)[keep]))
      }
      message("[zscape] cache built with another include list; rebuilding")
    }
    full <- readRDS(cfg$paths$zscape)
    m <- full@meta.data
    gt <- as.character(m$gene_target)
    keep <- !(gt %in% pr$drop_gene_target) & !is.na(gt) &
      as.character(m$timepoint) %in% pr$timepoints &
      as.character(m$cell_type_sub) %in% pr$include_cell_type_sub
    absent <- setdiff(pr$include_cell_type_sub, unique(as.character(m$cell_type_sub)))
    if (length(absent)) message("[zscape] include labels absent from the atlas: ",
                                paste(absent, collapse = "; "))
    obj <- subset(full, cells = colnames(full)[keep])
    rm(full); gc()
    obj$pool_label <- as.character(obj$cell_type_sub)
    obj$ob_annot <- obj$pool_label == "osteoblast"
    obj$dataset <- "zscape"
    obj@misc$include <- pr$include_cell_type_sub
    obj@misc$timepoints <- pr$timepoints
    saveRDS(obj, cache)
    return(obj)
  }
  obj <- suppressWarnings(suppressMessages(UpdateSeuratObject(readRDS(cfg$paths[[ds]]))))
  obj$fabian_cluster <- as.character(obj$RNA_snn_res.0.8)
  obj$timepoint <- if (ds == "fabian3") "3dpf" else "5dpf"
  obj$gene_target <- "sox10_lineage"
  obj$embryo <- paste0("rep", obj$bio.rep)          # replicate stands in for batch
  obj$ob_annot <- FALSE
  if (ds == "fabian5") {
    ob <- colnames(readRDS(cfg$paths$fabian5_ob))
    stopifnot(all(ob %in% colnames(obj)))
    obj$ob_annot <- colnames(obj) %in% ob
  }
  obj$pool_label <- ifelse(obj$ob_annot, "osteoblast (cluster 14)",
                           paste0("fabian_cl", obj$fabian_cluster))
  obj$dataset <- ds
  obj
}

# ---- Embedding ---------------------------------------------------------------------
#' HVG -> PCA -> sweep -> pick resolution -> UMAP. Returns the object with `pool_cluster`.
embed_pool <- function(obj, nfeatures = 2000) {
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  obj <- embed_features(obj)
  dims <- 1:choose_dims(obj)
  obj <- cluster_sweep(obj, dims, prefix = "pool")
  summ <- sweep_summary(obj, "pool")
  res <- pick_resolution(summ)
  obj$pool_cluster <- droplevels(obj[[paste0("pool_", res)]][, 1])
  Idents(obj) <- "pool_cluster"
  obj <- RunUMAP(obj, dims = dims, seed.use = SEED, verbose = FALSE)
  obj@misc$embed <- list(dims = dims, resolution = res, sweep = summ,
                         k = k_for(ncol(obj)))
  obj
}

#' Same-seed rerun of neighbours + clustering at the chosen settings; returns ARI to pool_cluster.
rerun_ari <- function(obj) {
  e <- obj@misc$embed
  o2 <- FindNeighbors(obj, dims = e$dims, k.param = e$k, verbose = FALSE)
  o2 <- FindClusters(o2, resolution = e$resolution, random.seed = SEED, verbose = FALSE)
  safe_ari(Idents(o2), obj$pool_cluster)
}

# ---- Module scores and contaminant filter -------------------------------------------
#' AddModuleScore for every module with >= 1 gene present; columns are named mod_<module>.
score_modules <- function(obj, modules) {
  modules <- modules[lengths(modules) > 0]
  for (m in names(modules)) {
    obj <- AddModuleScore(obj, features = list(modules[[m]]), name = "tmpmod",
                          seed = SEED, ctrl = 100)
    obj[[paste0("mod_", m)]] <- obj$tmpmod1
    obj$tmpmod1 <- NULL
  }
  obj
}

zscore <- function(x) if (sd(x) > 0) (x - mean(x)) / sd(x) else x * 0

#' Per-cluster mean module scores z-scored across clusters, plus the drop rule.
#' mesenchymal = z of the per-cluster sum of the mesenchymal modules' mean scores.
#' A contaminant module counts for a cluster only when its z > contaminant_z AND >= detection_min_frac
#' of the cluster's cells express >= detection_min_genes of its genes AND that fraction is
#' >= detection_fold x the pool median over clusters (the detection floor; see config).
cluster_module_z <- function(obj, cfg, modules, cluster_col = "pool_cluster") {
  cl <- as.character(obj[[cluster_col]][, 1])
  sc <- grep("^mod_", colnames(obj@meta.data), value = TRUE)
  means <- sapply(sc, function(s) tapply(obj[[s]][, 1], cl, mean))
  means <- matrix(means, ncol = length(sc), dimnames = list(sort(unique(cl)), sub("^mod_", "", sc)))
  z <- apply(means, 2, zscore)
  z <- matrix(z, nrow = nrow(means), dimnames = dimnames(means))
  contam <- intersect(names(cfg$modules$contaminant), colnames(z))
  mes_mods <- intersect(cfg$mesenchymal_modules, colnames(means))
  mes_z <- zscore(rowSums(means[, mes_mods, drop = FALSE]))
  ob_frac <- tapply(obj$ob_annot, cl, mean)[rownames(z)]
  n <- as.integer(table(cl)[rownames(z)])
  f <- cfg$pools$filter
  cnt <- LayerData(obj, layer = "counts")
  det <- sapply(contam, function(m) {
    g <- modules[[m]]
    tapply(Matrix::colSums(cnt[g, , drop = FALSE] > 0) >= min(f$detection_min_genes, length(g)), cl, mean)
  })
  det <- matrix(det, ncol = length(contam), dimnames = list(sort(unique(cl)), contam))[rownames(z), , drop = FALSE]
  pool_med <- apply(det, 2, median)
  hit <- z[, contam, drop = FALSE] > f$contaminant_z & det >= f$detection_min_frac &
    sweep(det, 2, f$detection_fold * pool_med, ">=")
  max_c <- apply(z[, contam, drop = FALSE], 1, max)
  top_c <- contam[apply(z[, contam, drop = FALSE], 1, which.max)]
  hit_mods <- apply(hit, 1, function(h) paste(contam[h], collapse = "+"))
  flagged <- rowSums(hit) > 0 & mes_z < f$mesenchymal_z_max
  protected <- flagged & ob_frac >= f$protect_ob_frac
  colnames(det) <- paste0("det_", contam)
  data.frame(cluster = rownames(z), n_cells = n, ob_frac = round(ob_frac, 3),
             round(z, 2), round(det, 3), mesenchymal_z = round(mes_z, 2),
             max_contaminant_z = round(max_c, 2), top_contaminant = top_c,
             contaminant_hits = hit_mods, flagged = flagged, protected = protected,
             drop = flagged & !protected, check.names = FALSE, row.names = NULL)
}

#' Top n positive markers (cluster vs rest) for each cluster in `clusters`.
top_markers <- function(obj, clusters, n = 10, cluster_col = "pool_cluster") {
  Idents(obj) <- cluster_col
  vapply(clusters, function(cl) {
    mk <- FindMarkers(obj, ident.1 = cl, only.pos = TRUE, logfc.threshold = 0.25,
                      min.pct = 0.1, verbose = FALSE)
    paste(head(rownames(mk)[order(-mk$avg_log2FC)], n), collapse = ", ")
  }, character(1))
}

#' Iterated contaminant filter: embed, score, drop flagged clusters, re-embed; at most
#' `max_rounds` removal rounds, and the returned object is always embedded and re-checked
#' after the last removal. Returns list(obj, removed, zs (per round), final_clean).
refine_pool <- function(obj, cfg, modules, max_rounds = 3) {
  removed <- list(); zs <- list()
  for (r in seq_len(max_rounds + 1)) {
    obj <- embed_pool(obj)
    obj <- score_modules(obj, modules)
    z <- cluster_module_z(obj, cfg, modules)
    z$round <- r
    zs[[r]] <- z
    drop <- z$cluster[z$drop]
    message(sprintf("round %d: %d cells, %d clusters, dropping %d (%s)", r, ncol(obj),
                    nrow(z), length(drop), paste(drop, collapse = ",")))
    if (!length(drop) || r > max_rounds) break
    rm_df <- z[z$drop, c("round", "cluster", "n_cells", "ob_frac", "contaminant_hits",
                         "max_contaminant_z", "mesenchymal_z")]
    rm_df$top_markers <- top_markers(obj, drop, cfg$pools$filter$n_markers_log)
    cl <- as.character(obj$pool_cluster)
    rm_df$labels <- vapply(drop, function(d) {
      t <- sort(table(obj$pool_label[cl == d]), decreasing = TRUE)
      paste(sprintf("%s (%d)", names(t)[seq_len(min(3, length(t)))], t[seq_len(min(3, length(t)))]),
            collapse = "; ")
    }, character(1))
    removed[[r]] <- rm_df
    obj <- subset(obj, cells = colnames(obj)[!cl %in% drop])
    obj <- DietSeurat(obj, layers = c("counts", "data"))
  }
  list(obj = obj, removed = do.call(rbind, removed), zs = do.call(rbind, zs),
       final_clean = !any(zs[[length(zs)]]$drop))
}

# ---- Export --------------------------------------------------------------------------
#' Write raw counts (genes x cells, MatrixMarket) + genes.csv + cells.csv, as export_reference.R.
write_mtx_set <- function(counts, meta, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  counts <- as(counts, "CsparseMatrix")
  frac <- max(abs(counts@x - round(counts@x)))
  if (frac > 1e-6) stop(dir, ": counts are not integers (max deviation ", frac, ")")
  writeMM(counts, file.path(dir, "counts.mtx"))
  write.csv(data.frame(gene = rownames(counts)), file.path(dir, "genes.csv"), row.names = FALSE)
  write.csv(cbind(cell = colnames(counts), meta), file.path(dir, "cells.csv"), row.names = FALSE)
  message(sprintf("  wrote %s: %d genes x %d cells", dir, nrow(counts), ncol(counts)))
}
