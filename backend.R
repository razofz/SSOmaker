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


hline <- function(y = 0, color = "black") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color)
  )
}


make_qc_plots <- function(sobj, col) {
  if (col %in% colnames(sobj[[]])) {
    fig1 <- sobj[[]] %>%
      plot_ly(
        y = as.formula(str_c(" ~ ", col)),
        type = "violin",
        box = list(visible = T),
        meanline = list(visible = T),
        name = "nCount_RNA",
        x0 = "nCount_RNA"
      ) %>%
      layout(yaxis = list(zeroline = F), shapes = list(hline(6), hline(600)))
    fig2 <- sobj[[]] %>%
      plot_ly(
        y = as.formula(str_c(" ~ log2(", col, ")")),
        type = "violin",
        box = list(visible = T),
        meanline = list(visible = T),
        name = "log2(nCount_RNA)",
        x0 = "nCount_RNA"
      ) %>%
      layout(yaxis = list(zeroline = F), shapes = list(hline(6), hline(9.5)))
    fig <- subplot(fig1, fig2)
    return(fig)
  } else {
    return(NULL)
  }
}
