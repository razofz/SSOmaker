detect_features_format <- function(features) {
  print(head(features))
  # check type of input: {dataframe,character vector,etc}
  # check number of columns (if df)
  # regex to make best-effort at detecting type of columns
  #   e.g. if line starts with "ENS", that is probably a geneID column
}


check_files <- function(directory, suffix) {
  feature_files <- c(
    "features.tsv",
    "genes.tsv"
  )
  mat_files <- c(
    "matrix.mtx",
    "barcodes.tsv"
  )

  if (
    !any(as.logical(
      lapply(
        file.path(
          directory,
          stringr::str_c(
            feature_files, suffix
          )
        ),
        file.exists
      )
    ))
  ) {
    return(FALSE)
  } else {
    return(
      all(as.logical(
        lapply(
          file.path(directory, stringr::str_c(mat_files, suffix)),
          file.exists
        )
      ))
    )
  }
}


validate_directory <- function(directory) {
  if (!dir.exists(directory)) {
    return(FALSE)
  }
  if (!check_files(directory, ".gz")) {
    if (!check_files(directory, "")) {
      return(FALSE)
    }
  }
  return(TRUE)
}


read_data <- function(
    chosen_folder,
    cell_column = NULL,
    gene_column = 2) {
  so_data <- Seurat::Read10X(
    chosen_folder,
    cell.column = cell_column,
    gene.column = gene_column
  )
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
  return(so_data)
}

make_seurat_object <- function(
    data_matrix,
    project_name = "SOM") {
  # print("Dimensions of data uploaded:")
  # print(dim(data_matrix))
  sobj <- Seurat::CreateSeuratObject(
    counts = data_matrix,
    project = project_name
  )
  # , min.cells = 3, min.features = 200)
  sobj[["percent_mt"]] <- Seurat::PercentageFeatureSet(
    sobj,
    pattern = "^(MT-)|^(mt-)"
  )
  return(sobj)
}

process_seurat_object <- function(
    seurat_object,
    filtering_thresholds = list()) {
  params <- list(
    scale_factor = 10e3,
    normalization_method = "LogNormalize",
    n_hvgs = 2e3,
    clustering_resolution = 0.8,
    n_pcs_to_use = 15,
    features_for_scaling = "hvgs" # c("all", "hvgs")
  )

  print("Normalising data..")
  shinyjs::logjs("Normalising data..")
  seurat_object <- Seurat::NormalizeData(
    object = seurat_object,
    verbose = FALSE
  )
  n_hvgs <- params$n_hvgs

  print("Finding HVGs..")
  shinyjs::logjs("Finding HVGs..")
  if (length(Seurat::Cells(seurat_object)) < params$n_hvgs) {
    n_hvgs <- length(Seurat::Cells(seurat_object)) * 0.4
  }
  seurat_object <- Seurat::FindVariableFeatures(
    object = seurat_object,
    nfeatures = n_hvgs,
    verbose = FALSE
  )

  print("Scaling data..")
  shinyjs::logjs("Scaling data..")
  if (params$features_for_scaling == "all") {
    features <- rownames(seurat_object)
  } else if (params$features_for_scaling == "hvgs") {
    features <- Seurat::VariableFeatures(seurat_object)
  }
  seurat_object <- Seurat::ScaleData(
    object = seurat_object,
    features = features,
    verbose = FALSE
  )

  print("Calculating PCA..")
  shinyjs::logjs("Calculating PCA..")
  seurat_object <- Seurat::RunPCA(
    object = seurat_object,
    features = features,
    verbose = FALSE
  )

  print("Finding neighbours..")
  shinyjs::logjs("Finding neighbours..")
  seurat_object <- Seurat::FindNeighbors(
    object = seurat_object,
    dims = 1:params$n_pcs_to_use,
    verbose = FALSE
  )

  print("Clustering..")
  shinyjs::logjs("Clustering..")
  seurat_object <- Seurat::FindClusters(
    object = seurat_object,
    resolution = params$clustering_resolution,
    verbose = FALSE
  )

  print("Running UMAP..")
  shinyjs::logjs("Running UMAP..")
  seurat_object <- Seurat::RunUMAP(
    object = seurat_object,
    dims = 1:params$n_pcs_to_use,
    verbose = FALSE
  )

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
