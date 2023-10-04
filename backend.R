library(Seurat)
library(ggplot2)

detect_features_format <- function(features) {
  print(head(features))
  # check type of input: {dataframe,character vector,etc}
  # check number of columns (if df)
  # regex to make best-effort at detecting type of columns
  #   e.g. if line starts with "ENS", that is probably a geneID column
}

read_data <- function(
    chosen_folder,
    project_name = "shiny") {
  print(chosen_folder)
  # check that folder exists yada yada
  # read in features file to see what structure it has, and what column is the gene name column
  # and set below accordingly
  so_data <- Read10X(chosen_folder, cell.column = NULL, gene.column = 1)
  # so <- CreateSeuratObject(so_data, project = project_name)
  return(so_data)
}
