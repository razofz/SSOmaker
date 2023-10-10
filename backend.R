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

mat_files <- c(
  "matrix.mtx",
  "features.tsv",
  "barcodes.tsv"
)

check_files <- function(directory, suffix) {
  return(all(as.logical(lapply(
    file.path(directory, str_c(mat_files, suffix)),
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
  # check that folder exists yada yada
  # read in features file to see what structure it has, and what column is the gene name column
  # and set below accordingly
  so_data <- Read10X(chosen_folder, cell.column = NULL, gene.column = 1)
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
  print("Dimensions of data uploaded:")
  print(dim(data_matrix))
  sobj <- CreateSeuratObject(counts = data_matrix, project = "SOM") # , min.cells = 3, min.features = 200)
  # print(sobj)
  # print(sobj[[]] %>% head)
  sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  return(sobj)
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


make_qc_plots <- function(
    sobj,
    col,
    input, output, session
    # min_cutoff,
    # max_cutoff
    ) {
  fig <- reactive({
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
        layout(
          yaxis = list(zeroline = F),
          shapes = list(
            # hline(min_cutoff),
            hline(input$qc_slider_nCount_RNA[1]),
            hline(input$qc_slider_nCount_RNA[2]),
            # hline(max_cutoff)
          )
        )
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
  })
  return(fig)
}

qc_plot_ui <- function(id, label = "QC") {
  make_qc_slider(x = pbmc_small$nCount_RNA, col = "nCount_RNA")
}

qc_plot_server <- function(
    id,
    sobj,
    col) {
  moduleServer(
    id,
    function(input, output, session) {
      min_cutoff <- reactiveVal(0)
      observeEvent(input$qc_slider_nCount_RNA, {
        min_cutoff(input$qc_slider_nCount_RNA[1])
      })
    }
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
