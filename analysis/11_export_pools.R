# Export refined pools (results/lineage/{ds}_pool.rds) as raw counts for Python:
#   results/lineage/export/<ds>/{counts.mtx (genes x cells), genes.csv, cells.csv}
# Same format as scrna_matching/iter0_data_prep/export_reference.R.
# Also writes panel-gene detection per dataset for check_exports-style verification in Python.
# Run from the repo root: Rscript analysis/11_export_pools.R [zscape fabian3 fabian5]

source("R/lineage.R")
cfg <- load_lineage_config()
args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else DATASETS
out <- lineage_dir(cfg, "export")

det_rows <- list()
for (ds in datasets) {
  p <- file.path(cfg$paths$results, paste0(ds, "_pool.rds"))
  if (!file.exists(p)) stop("missing ", p, "; render analysis/10_build_pools.Rmd first")
  obj <- readRDS(p)
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  m <- obj@meta.data
  um <- Embeddings(obj, "umap")
  meta <- data.frame(
    dataset = ds, embryo = as.character(m$embryo), gene_target = as.character(m$gene_target),
    timepoint = as.character(m$timepoint), n_umi = Matrix::colSums(counts),
    n_genes = Matrix::colSums(counts > 0), pool_label = m$pool_label, ob_annot = m$ob_annot,
    pool_cluster = as.character(m$pool_cluster),
    bio_rep = if ("bio.rep" %in% names(m)) m$bio.rep else NA,
    fabian_cluster = if ("fabian_cluster" %in% names(m)) m$fabian_cluster else NA,
    umap_1 = um[, 1], umap_2 = um[, 2],
    m[, grep("^mod_", names(m)), drop = FALSE], check.names = FALSE)
  write_mtx_set(counts, meta, file.path(out, ds))
  g <- intersect(cfg$panel_genes, rownames(counts))
  det_rows[[ds]] <- data.frame(dataset = ds, gene = g,
                               pct_detected = 100 * Matrix::rowMeans(counts[g, ] > 0))
}
det <- do.call(rbind, det_rows)
f <- file.path(out, "panel_detection_R.csv")
if (file.exists(f)) {
  old <- read.csv(f)
  det <- rbind(old[!old$dataset %in% datasets, ], det)
}
write.csv(det, f, row.names = FALSE)
