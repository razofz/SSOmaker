library(Seurat)
library(ggplot2)
library(stringr)

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

make_qc_slider <- function(x, col) {
  if (is.numeric(x)) {
    col_range <- range(x, na.rm = TRUE)
    sliderInput(
      inputId = str_c("qc_slider_", col),
      label = col,
      min = col_range[1],
      max = col_range[2],
      value = col_range,
      width = "75%"
    )
  } else if (is.factor(x)) {
    levs <- levels(x)
    selectInput(col, col, choices = levs, selected = levs, multiple = TRUE)
  } else {
    # Not supported
    NULL
  }
}