# Render analysis/10_build_pools.Rmd for each dataset.
# Run from the repo root: Rscript analysis/10_render_pools.R [zscape fabian3 fabian5]

source("R/lineage.R")

# Outside RStudio, fall back to RStudio's bundled pandoc if none is on PATH.
if (!rmarkdown::pandoc_available("1.12.3")) {
  arch <- if (R.version$arch == "aarch64") "aarch64" else "x86_64"
  rs <- file.path("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools", arch)
  if (file.exists(file.path(rs, "pandoc"))) Sys.setenv(RSTUDIO_PANDOC = rs)
  stopifnot(rmarkdown::pandoc_available("1.12.3"))
}

args <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(args)) args else DATASETS
out <- lineage_dir(load_lineage_config())

for (ds in datasets) {
  message("== 10_build_pools :: ", ds)
  rmarkdown::render("analysis/10_build_pools.Rmd", params = list(dataset = ds),
                    output_file = paste0("10_build_pools_", ds, ".html"), output_dir = out,
                    intermediates_dir = tempdir(), knit_root_dir = PROJ_DIR,
                    envir = new.env(), quiet = TRUE)
}
