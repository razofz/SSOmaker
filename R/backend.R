detect_features_format <- function(features) {
  print(head(features))
  # check type of input: {dataframe,character vector,etc}
  # check number of columns (if df)
  # regex to make best-effort at detecting type of columns
  #   e.g. if line starts with "ENS", that is probably a geneID column
}

mat_files <- c(
  "matrix.mtx",
  "features.tsv",
  "barcodes.tsv"
)

check_files <- function(directory, suffix) {
  return(all(as.logical(lapply(
    file.path(directory, stringr::str_c(mat_files, suffix)),
    file.exists
  ))))
}


validate_directory <- function(directory) {
  if (!dir.exists(directory)) {
    return(FALSE)
  }
  # TODO: add check for cellranger v2 `genes.tsv` as well
  if (!check_files(directory, ".gz")) {
    if (!check_files(directory, "")) {
      return(FALSE)
    }
  }
  return(TRUE)
}


read_data <- function(
    chosen_folder,
    project_name = "shiny") {
  print(chosen_folder)
  # read in features file to see what structure it has, and what column is the gene name column
  # and set below accordingly
  so_data <- Seurat::Read10X(chosen_folder, cell.column = NULL, gene.column = 1)
  # handle the case of e.g. CITE-seq, where there are multiple assays.
  # pick out only the RNA assay for now
  if (class(so_data) == "list") {
    print(names(so_data))
    so_data <- so_data[[
      grep(
        pattern = paste(c("RNA", "Gene Expression"), collapse = "|"),
        x = names(so_data),
        ignore.case = TRUE,
        value = TRUE
      )
    ]]
  }
  # so <- CreateSeuratObject(so_data, project = project_name)
  return(so_data)
}

make_seurat_object <- function(
    data_matrix) {
  # print("Dimensions of data uploaded:")
  # print(dim(data_matrix))
  sobj <- Seurat::CreateSeuratObject(counts = data_matrix, project = "SOM") # , min.cells = 3, min.features = 200)
  # print(sobj)
  sobj[["percent_mt"]] <- Seurat::PercentageFeatureSet(sobj, pattern = "^(MT)|(mt)")
  return(sobj)
}

process_seurat_object <- function(seurat_object) {
  # shinyjs::logjs(seurat_object)

  print("Normalising data..")
  shinyjs::logjs("Normalising data..")
  seurat_object <- Seurat::NormalizeData(object = seurat_object, verbose = FALSE)
  n_hvgs <- 2e3
  if (length(Seurat::Cells(seurat_object)) < 2e3) n_hvgs <- length(Seurat::Cells(seurat_object)) * 0.4
  print("Finding HVGs..")
  shinyjs::logjs("Finding HVGs..")
  seurat_object <- Seurat::FindVariableFeatures(object = seurat_object, nfeatures = n_hvgs, verbose = FALSE)
  print("Scaling data..")
  shinyjs::logjs("Scaling data..")
  seurat_object <- Seurat::ScaleData(object = seurat_object, verbose = FALSE)
  print("Calculating PCA..")
  shinyjs::logjs("Calculating PCA..")
  seurat_object <- Seurat::RunPCA(object = seurat_object, verbose = FALSE)
  print("Finding neighbours..")
  shinyjs::logjs("Finding neighbours..")
  seurat_object <- Seurat::FindNeighbors(object = seurat_object, verbose = FALSE)
  print("Clustering..")
  shinyjs::logjs("Clustering..")
  seurat_object <- Seurat::FindClusters(object = seurat_object, verbose = FALSE)
  print("Running UMAP..")
  shinyjs::logjs("Running UMAP..")
  seurat_object <- Seurat::RunUMAP(object = seurat_object, dims = 1:10, verbose = FALSE)

  return(seurat_object)
}

find_degs <- function(seurat_object) {
  print("Identifying DEGs..")
  shinyjs::logjs("Identifying DEGs..")
  markers <- Seurat::FindAllMarkers(seurat_object, only.pos = TRUE)
  markers <- dplyr::arrange(.data = markers, cluster, dplyr::desc(avg_log2FC))
  return(markers)
}


hline <- function(y = 0, color = "black", alpha = .8) {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color),
    alpha = alpha
  )
}

new_somaker_dataobject <- function(x = data.frame()) {
  stopifnot(is.data.frame(x))
  foo <- list(
    "selected_directory" = "",
    "nCount_RNA" = list(
      "low" = 0,
      "high" = Inf
    ),
    "nFeatures_RNA" = list(
      "low" = 0,
      "high" = Inf
    ),
    "percent_mt" = list(
      "low" = 0,
      "high" = 100
    )
  )
  return(
    structure(foo, class = "somaker_dataobject")
  )
}
