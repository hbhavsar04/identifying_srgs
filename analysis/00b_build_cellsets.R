# Build the two cell sets compared side by side and cache them in results/.
# Source data on /Volumes/kwongroup is only read, never written.
# Run from the repo root: Rscript analysis/00b_build_cellsets.R

source("R/utils.R")

# ob334: the existing osteoblast subset
ob <- readRDS(OB_PATH)
ob334 <- build_clean(ob)
ob334$cellset_label <- ob334$cell_type_sub
saveRDS(ob334, cellset_path("ob334"))
message("ob334: ", ncol(ob334), " cells x ", nrow(ob334), " genes")
rm(ob); gc()

# obfin: every osteoblast-annotated cell in the full atlas
full <- readRDS(FULL_PATH)
keep <- colnames(full)[full$cell_type_sub %in% OBFIN_TYPES]
sub <- subset(full, cells = keep)
rm(full); gc()
obfin <- build_clean(sub)
obfin$cellset_label <- obfin$cell_type_sub
stopifnot(all(colnames(ob334) %in% colnames(obfin)))
saveRDS(obfin, cellset_path("obfin"))
message("obfin: ", ncol(obfin), " cells x ", nrow(obfin), " genes")
print(table(obfin$cell_type_sub, obfin$timepoint))
