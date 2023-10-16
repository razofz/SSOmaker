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
  sobj[["percent_mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  return(sobj)
}

qc_slider_ui <- function(id) {
  div(
    style = "display: flex; justify-content: center;",
    sliderInput(
      inputId = NS(id, "slider"),
      label = "slider",
      min = 0,
      max = 5000,
      value = range(0, 5000),
      width = "75%"
    )
  )
  # uiOutput(NS(id, "output"))
}
qc_slider_server <- function(id, col, metadata) {
  # stopifnot(is.reactive(metadata))
  stopifnot(!is.reactive(col))

  moduleServer(id, function(input, output, session) {
    observe({
      col_range <- range(metadata$data[[col]], na.rm = TRUE)
      updateSliderInput(
        inputId = "slider",
        # inputId = NS(id, "slider"),
        label = col,
        min = col_range[1],
        max = col_range[2],
        value = col_range
      )
    })

    return(
      slider = shiny::reactive(input$slider)
    )
  })
}

qc_plot_ui <- function(id) {
  uiOutput(NS(id, "output"))
}
qc_plot_server <- function(id, col, metadata, ranges) {
  # stopifnot(is.reactive(metadata))
  stopifnot(!is.reactive(col))

  moduleServer(id, function(input, output, session) {
    output$output <- renderUI({
      tagList(
        # output$foo <- renderText(ranges()),
        textOutput(ranges()),
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
                # hline(300),
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
        qc_slider_ui(NS(id, "qc_slide")),
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
qc_module_server <- function(id, col, metadata) {
  moduleServer(id, function(input, output, session) {
    slider_input_vals <- qc_slider_server(
      "qc_slide",
      col = col,
      metadata = metadata
    )
    qc_plot_server(
      "qc_plot",
      col = col,
      metadata = metadata,
      ranges = slider_input_vals
    )
    output$foo <- renderText({slider_input_vals()})

    return(slider_input_vals)
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
