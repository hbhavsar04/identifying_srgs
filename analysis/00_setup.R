# One-time environment setup for the Seurat vs Festem comparison.
# Run with: Rscript analysis/00_setup.R

repos <- "https://cloud.r-project.org"

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = repos)

cran_pkgs <- c("clustree", "aricode", "UpSetR")
missing <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = repos)

if (!requireNamespace("Festem", quietly = TRUE)) pak::pkg_install("XiDsLab/Festem")

needed <- c("Seurat", "Festem", "clustree", "aricode", "UpSetR", "presto",
            "cluster", "pheatmap", "patchwork", "dplyr", "ggplot2", "rmarkdown")
ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("Missing packages: ", paste(needed[!ok], collapse = ", "))
message("Setup complete.")
