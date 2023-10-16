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


source("backend.R")

options(shiny.maxRequestSize = 900 * 1024^2) # 300 MB limit

ui <- tagList(
  useShinyjs(),
  page_navbar(
    id = "nav",
    selected = "load_data",
    # selected = "filtering",
    theme = bs_theme(
      bootswatch = "pulse",
      version = 5
    ),
    title = "Seurat Object Maker",
    fillable = F,
    sidebar = sidebar(
      conditionalPanel(
        "input.nav === 'load_data'",
        markdown(
          mds = c(
            "## Flow of app",
            "1. **Select a directory**",
            "2. Filter the data (QC)",
            "3. Results"
          )
        )
      ),
      conditionalPanel(
        "input.nav === 'filtering'",
        markdown(
          mds = c(
            "## Flow of app",
            "1. Select a directory",
            "2. **Filter the data (QC)**",
            "3. Results"
          )
        ) # ,
        # value_box(
        #   title = "Original number of cells",
        #   value = textOutput(outputId = "original_n_cells_side"),
        #   showcase = bsicons::bs_icon("clipboard-data"),
        #   theme_color = "success"
        # ),
        # value_box(
        #   title = "Number of cells left after filtering",
        #   value = textOutput(outputId = "filtered_n_cells_side"),
        #   showcase = bsicons::bs_icon("receipt-cutoff"),
        #   theme_color = "warning"
        # )
      ),
      conditionalPanel(
        "input.nav === 'results'",
        markdown(
          mds = c(
            "## Flow of app",
            "1. Select a directory",
            "2. Filter the data (QC)",
            "3. **Results**"
          )
        )
      ),
      uiOutput(outputId = "dataobject"),
      fillable = T,
      position = "right"
    ),
    nav_panel(
      id = "load_data",
      title = card_title("Load data"),
      value = "load_data",
      card_body(
        markdown(
          mds = c(
            "## Select a directory",
            "Select the **directory** (**folder**) where the _feature-barcode_ count matrices are. The files should be in the [Matrix Market Exchange format](https://math.nist.gov/MatrixMarket/formats.html) that e.g. the [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) pipeline outputs.",
            "",
            "That format consists of these files:",
            "",
            "```bash",
            "    matrix.mtx",
            "    features.tsv",
            "    barcodes.tsv",
            "```",
            "",
            "The files can also be compressed (in the `gzip` format) and then have the file extension `.gz`, e.g. `matrix.mtx.gz`.",
            "",
            "If the matrices were produced with e.g. a version of Cellranger below v3, `features.tsv` is instead named `genes.tsv`."
          )
        ),
        div(
          style = "display: flex; justify-content: center;",
          shinyDirButton(
            "directory",
            "Folder select",
            "Please select a folder",
            icon = bsicons::bs_icon("folder"),
            buttonType = "primary",
            style = "width: 40%;"
          )
        ),
        shinyjs::hidden(uiOutput("dir_valid_UI")),
        shinyjs::hidden(uiOutput("dir_invalid_UI")),
        # shinyjs::hidden(
        #   div(
        #     id = "data_preview_container",
        #     # shinycssloaders::withSpinner(
        #     #   DT::dataTableOutput(outputId = "data_preview")
        #     # ),
        #     # verbatimTextOutput("sobj_out")
        #   )
        # ),
        shinyjs::hidden(uiOutput("start_processing_UI"))
      ),
    ),
    ##############################################################################
    # Filtering tab
    ##############################################################################
    nav_panel(
      id = "filtering",
      value = "filtering",
      title = card_title("Filtering"),
      card_body(
        markdown(
          mds = c(
            "# Filtering (QC)"
          )
        ),
        # markdown(
        #   mds = c(
        #     "**Selected directory:**"
        #   )
        # ),
        # verbatimTextOutput(outputId = "selected_directory_filtering"),
        # histogramUI("hist1"),
        # qc_plot_ui("qc_test"),
        # test_ui("test"),
        # verbatimTextOutput("qc"),
        # uiOutput("qc"),
        markdown("---"),
        # value_box(
        #   # title = "Selected directory",
        #   value = "",
        #   title = textOutput(outputId = "selected_directory"),
        #   # value = textOutput(outputId = "selected_directory"),
        #   showcase = bsicons::bs_icon("file-earmark-spreadsheet"), # , size = NULL),
        #   theme_color = "success"
        # ),
        layout_column_wrap(
          width = 2,
          value_box(
            title = "Original number of cells",
            value = textOutput(outputId = "original_n_cells"),
            showcase = bsicons::bs_icon("clipboard-data"),
            theme_color = "success"
          ),
          value_box(
            title = "Number of cells left after filtering",
            value = textOutput(outputId = "filtered_n_cells"),
            showcase = bsicons::bs_icon("receipt-cutoff"),
            theme_color = "warning"
          )
        ),
        # ) |>
        #   popover(
        #     "Tooltip text",
        #     title = "Help"
        #   ),
        markdown(
          mds = c(
            "### Violin plots of basic characteristics"
          )
        ),
        # qc_plot_ui("qc_plot"),
        # qc_slider_ui("qc_sliderr"),
        # verbatimTextOutput("tmp"),
        qc_module_UI("qc_nCount_RNA"),
        # qc_slider_ui("qc_test"),
        # layout_column_wrap(
        #   width = "600px",
        #   div(
        #     style = "display: flex; justify-content: center;",
        #     div(
        #       div(
        #         style = "display: flex; justify-content: center;",
        #         tooltip(
        #           span(
        #             helpText("Need help? "),
        #             bs_icon("info-circle")
        #           ),
        #           str_c(
        #             "Hover over the plot to see extra statistics about the metadata column. ",
        #             "You can also drag to zoom in on a specific area of the plot."
        #           )
        #         )
        #       ),
        #       plotlyOutput(outputId = "violin_plot"),
        #       uiOutput("qc_slider"),
        #       div(
        #         style = "display: flex; justify-content: center;",
        #         span(
        #           helpText("Adjust the sliders to set the filtering cutoffs."),
        #           tooltip(
        #             bs_icon("info-circle"),
        #             "Drag the sliders to set the cutoffs for this metadata column. ",
        #             "Notice the change in the black lines in the plot above, which correspond to the selected cutoffs. ",
        #             "Also note that the number of cells displayed up top is updated according to how you change the sliders.",
        #             placement = "right"
        #           )
        #         )
        #       ),
        #       style = " width: 60%; max-width: 600px;",
        #     ),
        #   )
        # ),
        textOutput(outputId = "qc_value"),
        # splitLayout(
        #   make_qc_slider(x = pbmc_small$nCount_RNA, col = "nCount_RNA"),
        #   make_qc_slider(x = pbmc_small$nFeature_RNA, col = "nFeature_RNA")
        # ),
        markdown(
          mds = c(
            "### Choose filtering parameters"
          )
        ),
        textOutput(outputId = "filtered_dimensions"),
        markdown(
          mds = c(
            "Showing the metadata for the dataset, in order to help choose which columns to filter. Some suggestions have been selected in the checkboxes below."
          )
        ),
        DT::dataTableOutput("metadata")
      ) # ,
    ),
    nav_panel(
      id = "results",
      value = "results",
      title = card_title("Results")
    ),
    nav_spacer(),
    nav_item(
      markdown(
        '<img src = "https://www.staff.lu.se/sites/staff.lu.se/files/styles/lu_wysiwyg_full_tablet/public/2021-04/Lunduniversity-horisontal.png.webp?itok=_rp_OxRe" width="200px" />'
      )
    ) # ,
    # tags$head(tags$script(src = "./message_handler.js"))
  )
)

# -------------------------------------------------------------------------------

server <- function(input, output, session) {
  ################################################################################
  # Variable setup
  ################################################################################

  somaker_dataobject <- reactiveValues(
    selected_directory = character(0),
    is_valid_directory = FALSE,
    is_valid_data = FALSE,
    nCount_RNA_low = 6,
    nCount_RNA_high = 600,
    nFeature_RNA_high = 600,
    nFeature_RNA_low = 200,
    percent_mito_high = 20
  )

  # reactive_metadata <- reactiveValues(data = data.frame())
  reactive_metadata <- reactiveValues(data = pbmc_small[[]])
  reactive_sobj <- reactiveValues(data = pbmc_small)

  nav_hide(id = "nav", target = "filtering")
  nav_hide(id = "nav", target = "results")

  # filtering debugging
  nav_show(id = "nav", target = "filtering")
  nav_select(id = "nav", selected = "filtering")

  ################################################################################
  # Sidebar
  ################################################################################

  output$nav_now <- renderText(input$nav)

  output$dataobject <- renderUI({
    bar <- reactiveValuesToList(somaker_dataobject)
    result <- tagList(tags$h3("Data object"))
    for (i in seq_along(bar)) {
      result <- tagList(
        result,
        tags$h6(names(bar)[i]),
        tags$p(bar[[i]])
      )
    }
    return(result)
  })

  ################################################################################
  # Load data tab
  ################################################################################

  volumes <- c(Home = getwd(), "R Installation" = R.home(), getVolumes()()) # TODO: change back to fs::path_home() when done testing

  shinyDirChoose(
    input,
    "directory",
    roots = volumes,
    session = session,
    restrictions = system.file(package = "base"),
    allowDirCreate = FALSE
  )

  # Validate and store the selected directory
  observeEvent(input$directory, {
    if (!is.integer(input$directory)) {
      selected_directory <- parseDirPath(volumes, input$directory)
      somaker_dataobject$selected_directory <- selected_directory
      if (length(somaker_dataobject$selected_directory) > 0) {
        is_valid_dir <- validate_directory(somaker_dataobject$selected_directory)
        somaker_dataobject$is_valid_directory <- is_valid_dir
        if (is_valid_dir) {
          shinyjs::hide("dir_invalid_UI")
          shinyjs::show("dir_valid_UI")
        } else {
          shinyjs::hide("dir_valid_UI")
          shinyjs::show("dir_invalid_UI")
        }
      } else {
        shinyjs::hide("dir_valid_UI")
        shinyjs::hide("dir_invalid_UI")
      }
      # somaker_dataobject$selected_directory <- parseDirPath(roots = c(wd = getwd()), input$directory)
    }
  })

  shinyjs::onclick("directory",
    expr = {
      shinyjs::hide("dir_valid_UI")
      shinyjs::hide("dir_invalid_UI")
      shinyjs::enable("load_into_seurat")
      shinyjs::hide("start_processing_UI")
      somaker_dataobject$is_valid_data <- FALSE
      somaker_dataobject$is_valid_directory <- FALSE
    }
  )

  output$dir_valid_UI <- renderUI({
    div(
      markdown(str_c(
        "`", somaker_dataobject$selected_directory, "`",
        " is a valid count matrix directory.\n\n",
        "Press button below to load the data into `Seurat`."
      )),
      div(
        style = "display: flex; justify-content: center;",
        actionButton(
          "load_into_seurat",
          label = HTML(
            bsicons::bs_icon("box-arrow-in-up-right"),
            " Load data into Seurat"
          ),
          class = "btn-warning",
          style = "width: 40%;"
        )
      )
    )
  })

  output$dir_invalid_UI <- renderUI({
    markdown(str_c(
      "`", somaker_dataobject$selected_directory, "`",
      " is _not_ a valid count matrix directory. ",
      "Please see the description of the format above."
    ))
  })

  shinyjs::onclick("load_into_seurat",
    expr = {
      shinyjs::disable("load_into_seurat")
      start_time <- Sys.time()
      load_time <- system.time(
        showPageSpinner(
          type = 6,
          expr = {
            sobj_data <- read_data(somaker_dataobject$selected_directory)
            sobj <- make_seurat_object(sobj_data)
          },
          caption = HTML(
            "Loading data into Seurat..",
            bsicons::bs_icon("tools")
          )
        )
      )
      end_time <- Sys.time()
      print(end_time - start_time)
      print(load_time["elapsed"])
      somaker_dataobject$is_valid_data <- TRUE
    }
  )

  output$sobj_out <- renderPrint({
    reactive_sobj$data
  })

  observeEvent(somaker_dataobject$is_valid_data, {
    req(somaker_dataobject$is_valid_data)
    # if (somaker_dataobject$is_valid_data == TRUE) {
    reactive_metadata$data <- sobj[[]]
    reactive_sobj$data <- sobj
    shinyjs::show("start_processing_UI")
    # }
  })

  output$start_processing_UI <- renderUI({
    div(
      markdown("Data loaded correctly:"),
      verbatimTextOutput("sobj_out"),
      markdown("Press button below to start processing!"),
      div(
        style = "display: flex; justify-content: center;",
        actionButton(
          "start_processing",
          label = HTML(
            bsicons::bs_icon("pc-display-horizontal"),
            " Start processing"
          ),
          class = "btn-success",
          style = "width: 40%;" # ,
        )
      )
    )
  })

  shinyjs::onclick("start_processing",
    expr = {
      # shinyjs::alert("will move to Filtering tab")
      # hideTab(inputId = "nav", target = "load_data")
      # showTab(inputId = "nav", target = "filtering")
      # updateTabsetPanel(inputId = "nav", selected = "filtering")
      nav_hide(id = "nav", target = "load_data")
      nav_show(id = "nav", target = "filtering")
      nav_select(id = "nav", selected = "filtering")
    }
  )

  ################################################################################
  # Filtering tab
  ################################################################################

  # histogramServer("hist1")
  # qc <- qc_plot_server("qc_test", metadata = reactive_metadata, col = "nCount_RNA")
  # qc <- qc_plot_server("qc_test", metadata = reactive_metadata$data, col = "nCount_RNA")
  # reactive({
  #   qc <- qc_plot_server("qc_test", metadata = reactive_metadata$data, col = "nCount_RNA")
  #   output$qc <- renderPrint(qc)
  # })

  # output$qc <- renderUI({
  #   qc_slider_server("qc", "nCount_RNA", reactive_metadata$data)
  # })

  # output$qc_value <- renderText(input$qc_slider_nCount_RNA)
  # output$qc_value <- renderText(input$qc_test_qc_slider)

  # foo <- reactiveVal("oooh")
  # output$test <- test_server("test", foo = reactive_metadata$data)

  # slider_input_valss <- qc_slider_server(
  #   "qc_sliderr",
  #   # "qc_test",
  #   col = "nCount_RNA",
  #   metadata = reactive_metadata # $data
  # )
  # output$tmp <- renderPrint({input[["qc_sliderr-output"]]})
  # output$tmp <- renderPrint({slider_input_valss()})
  # observe(
  #   # print(slider_input_valss())#,
  #   print(head(reactive_metadata$data))
  #   # input$qc_sliderr
  #   # input[["qc_sliderr-output"]]
  # )
  # observe(
  #   print("aah")
  # )
  # qc_plot_server(
  #   "qc_plot",
  #   # "qc_test",
  #   col = "nCount_RNA",
  #   metadata = reactive_metadata, # $data
  #   ranges = slider_input_vals
  # )
  slider_input_vals <- qc_module_server(
    "qc_nCount_RNA",
    col = "nCount_RNA",
    metadata = reactive_metadata # $data
  )
  observe(
    print(slider_input_vals())
  )
  # output$tmp <- renderPrint({slider_input_vals})


  output$original_n_cells <- renderText({
    nrow(reactive_metadata$data[, ])
  })
  output$filtered_n_cells <- renderText({
    nrow(reactive_metadata$data[
      reactive_metadata$data$nCount_RNA >= slider_input_vals()[1] &
        reactive_metadata$data$nCount_RNA <= slider_input_vals()[2],
      # reactive_metadata$data$nCount_RNA >= input$qc_slider_nCount_RNA[1] &
      #   reactive_metadata$data$nCount_RNA <= input$qc_slider_nCount_RNA[2],
    ])
  })

  # output$filtered_dimensions <- renderText({
  #   str_c(
  #     "Number of cells before filtering: ",
  #     nrow(reactive_metadata$data[reactive_metadata$data$nCount_RNA >= input$qc_slider_nCount_RNA[1], ]),
  #     " cells. ",
  #     "Cells left after filtering: ",
  #     nrow(reactive_metadata$data[reactive_metadata$data$nCount_RNA <= input$qc_slider_nCount_RNA[2], ]),
  #     " cells."
  #   )
  # })

  # output$qc_slider <- renderUI({
  #   div(
  #     style = "display: flex; justify-content: center;",
  #     sliderInput(
  #       inputId = str_c("qc_slider_", "nCount_RNA"),
  #       label = "nCount_RNA",
  #       min = range(reactive_metadata$data$nCount_RNA, na.rm = TRUE)[1],
  #       max = range(reactive_metadata$data$nCount_RNA, na.rm = TRUE)[2],
  #       value = range(reactive_metadata$data$nCount_RNA, na.rm = TRUE),
  #       width = "80%"
  #     )
  #   )
  # })

  # output$violin_plot <- renderPlotly({
  #   reactive_metadata$data %>%
  #     plot_ly(
  #       y = as.formula(str_c(" ~ ", "nCount_RNA")),
  #       type = "violin",
  #       box = list(visible = T),
  #       meanline = list(visible = T),
  #       name = "nCount_RNA",
  #       x0 = "nCount_RNA"
  #     ) %>%
  #     layout(
  #       yaxis = list(zeroline = F),
  #       shapes = list(
  #         # hline(min_cutoff),
  #         hline(input$qc_slider_nCount_RNA[1]),
  #         hline(input$qc_slider_nCount_RNA[2])
  #         # hline(max_cutoff)
  #       )
  #     ) # %>%
  #   # config(
  #   #   selectmode = "lasso",
  #   # )
  # })

  # observeEvent(input$show, {
  #     showNotification(
  #       ui = markdown(c(
  #         "## This is a notification. ",
  #         "_Take care._"
  #       )),
  #       type = "warning"
  #     )
  #   })


  # fig <- make_qc_plots(
  #   sobj = pbmc_small,
  #   col = "nCount_RNA",
  #   min_cutoff = input$qc_slider_nCount_RNA[1],
  #   max_cutoff = input$qc_slider_nCount_RNA[2],
  # )
  # fig <- callModule(make_qc_plots, "backend")
  # output$violin_plot <- renderPlotly(
  #   fig
  # make_qc_plots(
  #   sobj = pbmc_small,
  #   col = "nCount_RNA",
  #   input = input,
  #   output = output,
  #   session = session,
  # )
  # )
  # output$violin_plot <- renderPlotly({
  #   reactive(make_qc_plots(
  #     sobj = pbmc_small,
  #     col = "nCount_RNA",
  #     min_cutoff = input$qc_slider_nCount_RNA[1],
  #     max_cutoff = input$qc_slider_nCount_RNA[2],
  #   ))
  # })
  # input$qc_slider_nCount_RNA
  # output$colnames_output <- renderText(colnames(pbmc_small[[]]))
  output$metadata <- DT::renderDataTable({
    DT::datatable(
      reactive_metadata$data,
      caption = "Metadata",
      options = list(
        pageLength = 5
      ),
      fillContainer = T
      # style = 'bootstrap'
    )
  })
  # print(session$clientData)

  ################################################################################
  # Results tab
  ################################################################################

  # TODO: this tab :-)
}

shinyApp(ui, server)
