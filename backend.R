library(Seurat)
library(ggplot2)
library(stringr)
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
library(plotly)
library(bslib)
library(bsicons)
library(Seurat)
library(stringr)
library(dplyr)
library(DT)
library(htmltools)



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
  # print("Dimensions of data uploaded:")
  # print(dim(data_matrix))
  sobj <- CreateSeuratObject(counts = data_matrix, project = "SOM") # , min.cells = 3, min.features = 200)
  # print(sobj)
  # print(sobj[[]] %>% head)
  sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^(MT)|(mt)")
  return(sobj)
}

process_seurat_object <- function(seurat_object) {
  # shinyjs::logjs(seurat_object)

  print("Normalising data..")
  shinyjs::logjs("Normalising data..")
  seurat_object <- NormalizeData(object = seurat_object, verbose = FALSE)
  n_hvgs <- 2e3
  if (length(Cells(seurat_object)) < 2e3) n_hvgs <- length(Cells(seurat_object)) * 0.4
  print("Finding HVGs..")
  shinyjs::logjs("Finding HVGs..")
  seurat_object <- FindVariableFeatures(object = seurat_object, nfeatures = n_hvgs, verbose = FALSE)
  print("Scaling data..")
  shinyjs::logjs("Scaling data..")
  seurat_object <- ScaleData(object = seurat_object, verbose = FALSE)
  print("Calculating PCA..")
  shinyjs::logjs("Calculating PCA..")
  seurat_object <- RunPCA(object = seurat_object, verbose = FALSE)
  print("Finding neighbours..")
  shinyjs::logjs("Finding neighbours..")
  seurat_object <- FindNeighbors(object = seurat_object, verbose = FALSE)
  print("Clustering..")
  shinyjs::logjs("Clustering..")
  seurat_object <- FindClusters(object = seurat_object, verbose = FALSE)
  print("Running UMAP..")
  shinyjs::logjs("Running UMAP..")
  seurat_object <- RunUMAP(object = seurat_object, dims = 1:10, verbose = FALSE)

  return(seurat_object)
}

find_degs <- function(seurat_object) {
  print("Identifying DEGs..")
  shinyjs::logjs("Identifying DEGs..")
  markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
  markers <- markers %>% arrange(cluster, desc(avg_log2FC))
  return(markers)
}

################################################################################
# shiny modules
################################################################################

qc_slider_ui <- function(id) {
  tagList(
    sliderInput(
      inputId = NS(id, "slider"),
      label = "slider",
      min = 0,
      max = 5000,
      value = c(0, 5000),
      width = "75%"
    ),
  )
}
qc_slider_server <- function(id, col, metadata, start_values) {
  stopifnot(is.reactivevalues(metadata))
  stopifnot(!is.reactive(col))
  stopifnot(!is.reactive(start_values))

  moduleServer(id, function(input, output, session) {
    update_it <- reactiveVal(0)
    observe({
      if (isolate(update_it()) < 3) {
        isolate(update_it(update_it() + 1))
      }
      if (update_it() < 2) {
        # print(str_c("update_it(): ", update_it()))
        invalidateLater(300, session)
      }
      col_range <- range(isolate(metadata$data[[col]]), na.rm = TRUE)
      updateSliderInput(
        session,
        inputId = "slider",
        label = col,
        min = col_range[1],
        max = col_range[2],
        value = start_values
      )
    })
    # shinyjs::logjs(start_values)

    return(reactive(input$slider))
  })
}

qc_plot_ui <- function(id) {
  uiOutput(NS(id, "output"))
}
qc_plot_server <- function(id, col, metadata, ranges) {
  stopifnot(is.reactivevalues(metadata))
  stopifnot(!is.reactive(col))
  stopifnot(is.reactive(ranges))

  moduleServer(id, function(input, output, session) {
    # observe({
    #   print("ranges():")
    #   print(ranges())
    # })
    output$output <- renderUI({
      tagList(
        # output$foo <- renderText(ranges()),
        # textOutput(ranges()),
        renderPlotly({
          metadata$data %>%
            plotly::plot_ly(
              y = as.formula(str_c(" ~ ", col)),
              type = "violin",
              box = list(visible = T),
              meanline = list(visible = T),
              name = col,
              x0 = col
            ) %>%
            layout(
              yaxis = list(zeroline = F),
              shapes = list(
                hline(ranges()[1]),
                hline(ranges()[2])
              )
            )
        })
      )
    })
  })
}

qc_module_UI <- function(id) {
  layout_column_wrap(
    width = "600px",
    div(
      style = "display: flex; justify-content: center;",
      div(
        div(
          style = "display: flex; justify-content: center;",
          tooltip(
            span(
              helpText("Need help? "),
              bs_icon("info-circle")
            ),
            str_c(
              "Hover over the plot to see extra statistics about the metadata column. ",
              "You can also drag to zoom in on a specific area of the plot."
            )
          )
        ),
        qc_plot_ui(NS(id, "qc_plot")),
        # verbatimTextOutput(NS(id, "foo")),
        div(
          style = "display: flex; justify-content: center;",
          qc_slider_ui(NS(id, "qc_slide"))
        ),
        div(
          style = "display: flex; justify-content: center;",
          span(
            helpText("Adjust the sliders to set the filtering cutoffs."),
            tooltip(
              bs_icon("info-circle"),
              "Drag the sliders to set the cutoffs for this metadata column. ",
              "Notice the change in the black lines in the plot above, which correspond to the selected cutoffs. ",
              "Also note that the number of cells displayed up top is updated according to how you change the sliders.",
              placement = "right"
            )
          )
        ),
        style = " width: 60%; max-width: 600px;",
      )
    )
  )
}
qc_module_server <- function(id, col, metadata, start_values) {
  stopifnot(!is.reactive(start_values))
  stopifnot(is.reactivevalues(metadata))
  stopifnot(!is.reactive(col))

  moduleServer(id, function(input, output, session) {
    slider_input_vals <- qc_slider_server(
      "qc_slide",
      col = col,
      metadata = metadata,
      start_values = start_values
    )
    qc_plot_server(
      "qc_plot",
      col = col,
      metadata = metadata,
      ranges = slider_input_vals
    )
    output$foo <- renderText({
      slider_input_vals()
    })

    return(slider_input_vals)
  })
}

violin_plot_ui <- function(id) {
  uiOutput(NS(id, "output"))
}
violin_plot_server <- function(id, col, metadata) {
  stopifnot(is.reactivevalues(metadata))
  stopifnot(!is.reactive(col))

  moduleServer(id, function(input, output, session) {
    output$output <- renderUI({
      tagList(
        div(
          style = "display: flex; justify-content: center;",
          tooltip(
            span(
              helpText("Need help? "),
              bs_icon("info-circle")
            ),
            str_c(
              "Hover over the plot to see extra statistics about the metadata column. ",
              "You can also drag to zoom in on a specific area of the plot."
            )
          )
        ),
        renderPlotly({
          metadata$data %>%
            plotly::plot_ly(
              y = as.formula(str_c(" ~ ", col)),
              type = "violin",
              box = list(visible = T),
              meanline = list(visible = T),
              name = col,
              x0 = col
            ) %>%
            layout(
              yaxis = list(zeroline = F)
            )
        })
      )
    })
  })
}


test_ui <- function(id) {
  uiOutput(NS(id, "output"))
}
test_server <- function(id, foo) {
  moduleServer(id, function(input, output, session) {
    output$text <- DT::renderDataTable({
      DT::datatable(
        foo,
        options = list(
          pageLength = 5
        )
      )
    })

    observeEvent(
      {
        input$n
      },
      {
        output$text <- DT::renderDataTable({
          DT::datatable(
            foo,
            options = list(
              pageLength = input$n
            )
          )
        })
      }
    )

    output$output <- renderUI({
      tagList(
        dataTableOutput(NS(id, "text")),
        numericInput(NS(id, "n"), "n", value = 10)
      )
    })
  })
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
