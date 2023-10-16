library(Seurat)
library(ggplot2)
library(stringr)
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinycssloaders)
library(plotly)
# library(gridlayout)
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

# make_qc_slider <- function(x, col) {
#   if (is.numeric(x)) {
#     col_range <- range(x, na.rm = TRUE)
#     sliderInput(
#       inputId = str_c("qc_slider_", col),
#       label = col,
#       min = col_range[1],
#       max = col_range[2],
#       value = col_range,
#       width = "75%"
#     )
#   } else if (is.factor(x)) {
#     levs <- levels(x)
#     selectInput(col, col, choices = levs, selected = levs, multiple = TRUE)
#   } else {
#     # Not supported
#     NULL
#   }
# }

make_qc_slider <- function(x, id, col) {
  if (is.numeric(x)) {
    col_range <- range(x, na.rm = TRUE)
    sliderInput(
      inputId = id,
      # label = col,
      min = col_range[1],
      max = col_range[2],
      value = col_range,
      width = "75%"
    )
    # } else if (is.factor(x)) {
    # levs <- levels(x)
    # selectInput(col, col, choices = levs, selected = levs, multiple = TRUE)
  } else {
    # Not supported
    NULL
  }
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
    # output$output <- renderUI({
    #   col_range <- range(metadata$data[[col]], na.rm = TRUE)
    #   div(
    #     style = "display: flex; justify-content: center;",
    #     sliderInput(
    #       # inputId = NS(id, "slider"),
    #       inputId = "slider",
    #       label = col,
    #       min = col_range[1],
    #       max = col_range[2],
    #       value = col_range,
    #       width = "75%"
    #     )
    #   )
    # })

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
    # print(input)
    # print(session)
    # print(reactive({slider_input_vals()}))
    # output$foo <- renderText("aaah")
    # output$foo <- renderText(input$qc_slide)
    output$foo <- renderText({slider_input_vals()})
    # output$foo <- renderText({input[["qc_slide-slider"]]})

    return(slider_input_vals)
    # return(reactive({input[["qc_slide-slider"]]}))
    # return(input[["qc_slide-slider"]])
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

# qc_plot_ui <- function(id, label = "QC") {
# qc_plot_ui <- function(id) {
#   uiOutput(NS(id, "output"))
#   # uiOutput(NS(id, str_c("qc_slider", col)))
#   # make_qc_slider(x = pbmc_small$nCount_RNA, col = "nCount_RNA")
# }
# qc_plot_server <- function(id,
#     col,
#     metadata
#   ) {
#   stopifnot(is.reactive(metadata))
#   stopifnot(!is.reactive(col))

#   moduleServer(id, function(input, output, session, col = col) {
#       col_range <- range(metadata[[col]], na.rm = TRUE)
#       print(col)
#       NS(id, "output") <- renderUI({
#         sliderInput(
#           # inputId = str_c("qc_slider_", col),
#           inputId = NS(id, str_c("qc_slider_", col)),
#           label = col,
#           # min = range(metadata[[col]], na.rm = TRUE)[1],
#           # max = range(metadata[[col]], na.rm = TRUE)[2],
#           # value = range(metadata[[col]], na.rm = TRUE),
#           min = col_range[1],
#           max = col_range[2],
#           value = col_range,
#           width = "75%"
#         )
#       })
#       observeEvent(metadata(), {
#         updateSliderInput(
#           session,
#           str_c("qc_slider_", col),
#           min = col_range[1],
#           max = col_range[2],
#           value = col_range,
#         )
#       })

#       reactive(col_range)
#     }
#   )
# }


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

# histogramUI <- function(id) {
#   tagList(
#     # selectInput(NS(id, "var"), "Variable", choices = names(mtcars)),
#     numericInput(NS(id, "bins"), "bins", value = 10, min = 1)
#     # plotOutput(NS(id, "hist"))
#   )
# }
# histogramServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#     data <- reactive(mtcars[[input$var]])
#     # output$hist <- renderPlot({
#     #   hist(data(), breaks = input$bins, main = input$var)
#     # }, res = 96)
#   })
# }
