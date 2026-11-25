# Render the full comparison for both cell sets.
# Run from the repo root: Rscript analysis/run_all.R [cellset ...]
# Prerequisites: Rscript analysis/00_setup.R && Rscript analysis/00b_build_cellsets.R

source("R/utils.R")

# Outside RStudio, fall back to RStudio's bundled pandoc if none is on PATH.
if (!rmarkdown::pandoc_available("1.12.3")) {
  arch <- if (R.version$arch == "aarch64") "aarch64" else "x86_64"
  rs <- file.path("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools", arch)
  if (file.exists(file.path(rs, "pandoc"))) Sys.setenv(RSTUDIO_PANDOC = rs)
  stopifnot(rmarkdown::pandoc_available("1.12.3"))
}

args <- commandArgs(trailingOnly = TRUE)
cellsets <- if (length(args)) args else CELLSETS

for (cs in cellsets) {
  for (rmd in c("01_seurat_subcluster", "02_festem_subcluster", "03_compare")) {
    message("== ", cs, " :: ", rmd)
    rmarkdown::render(file.path("analysis", paste0(rmd, ".Rmd")),
                      params = list(cellset = cs),
                      output_file = paste0(rmd, ".html"),
                      output_dir = cellset_dir(cs),
                      intermediates_dir = tempdir(),
                      envir = new.env(), quiet = TRUE)
  }
}
