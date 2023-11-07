suppressPackageStartupMessages(library(SeuratObject))
library(stringr)
library(Matrix)

out_dir <- "pbmc_small"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarning = FALSE)
}

mat_files <- c(
  "matrix.mtx",
  "features.tsv",
  "barcodes.tsv"
)

check_files <- function(suffix) {
  return(all(as.logical(lapply(
    file.path(out_dir, str_c(mat_files, suffix)),
    file.exists
  ))))
}

if (!check_files(".gz")) {
  if (!check_files("")) {
    mat <- pbmc_small@assays$RNA@counts
    writeMM(mat, file.path(out_dir, "matrix.mtx"))
    write.table(
      data.frame("gene_name" = rownames(mat), "feature_type" = "Gene Expression"),
      quote = FALSE,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      file = file.path(out_dir, "features.tsv")
    )
    write.table(
      colnames(mat),
      quote = FALSE,
      sep = "\t",
      col.names = FALSE,
      row.names = FALSE,
      file = file.path(out_dir, "barcodes.tsv")
    )
  }
}
