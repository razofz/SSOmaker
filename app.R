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


source("backend.R")

options(shiny.maxRequestSize = 900 * 1024^2) # 900 MB limit

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
        )
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
        markdown("---"),
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
        uiOutput(outputId = "qc_UI"),
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
  print("############################################## new run")
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
    percent_mt_high = 20,
    columns_to_filter = list(
      "nCount_RNA",
      "nFeature_RNA"
    )
  )

  # reactive_metadata <- reactiveValues(data = data.frame())
  reactive_metadata <- reactiveValues(data = pbmc_small[[]])
  reactive_sobj <- reactiveValues(data = pbmc_small)

  nav_hide(id = "nav", target = "filtering")
  nav_hide(id = "nav", target = "results")

  # # filtering debugging
  # nav_show(id = "nav", target = "filtering")
  # nav_select(id = "nav", selected = "filtering")

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

  slider_input_vals <- reactiveValues()

  shinyjs::onclick("start_processing",
    expr = {
      # shinyjs::alert("will move to Filtering tab")
      # hideTab(inputId = "nav", target = "load_data")
      # showTab(inputId = "nav", target = "filtering")
      # updateTabsetPanel(inputId = "nav", selected = "filtering")
      nav_hide(id = "nav", target = "load_data")
      nav_show(id = "nav", target = "filtering")
      nav_select(id = "nav", selected = "filtering")
      for (column in colnames(reactive_metadata$data)) {
        if (is.numeric(reactive_metadata$data[[column]])) {
          # print(column)
          values <- stats::quantile(
            reactive_metadata$data[[column]],
            probs = c(0.05, 0.95)
          ) %>%
            round() %>%
            as.vector()
          # slider_input_vals[[column]] <- values
          if (values[1] != values[2]) {
            # print(column)
            if (!column %in% somaker_dataobject$columns_to_filter) {
              somaker_dataobject$columns_to_filter <- c(
                somaker_dataobject$columns_to_filter,
                column
              )
            }
            print(column)
            print(values)
            somaker_dataobject[[str_c(column, "_low")]] <- values[1]
            somaker_dataobject[[str_c(column, "_high")]] <- values[2]
            slider_input_vals[[column]] <- qc_module_server(
              str_c("qc_", column),
              col = column,
              metadata = reactive_metadata,
              start_values = values
            )
          }
        }
        # print(column)
        output$qc_UI <- renderUI({
          tagList(
            lapply(
              somaker_dataobject$columns_to_filter,
              FUN = \(column) {
                # print(column)
                qc_module_UI(str_c("qc_", column))
              }
            )
          )
        })
      }
      # print(slider_input_vals)
    }
  )

  ################################################################################
  # Filtering tab
  ################################################################################


  # slider_input_vals[[column]] <- qc_module_server(
  #   str_c("qc_", column),
  #   col = column,
  #   metadata = reactive_metadata
  # )

  # slider_input_vals[["nCount_RNA"]] <- qc_module_server(
  #   "qc_nCount_RNA",
  #   col = "nCount_RNA",
  #   metadata = reactive_metadata
  # )
  # observe(
  #   print(slider_input_vals[["nCount_RNA"]]())
  # )

  output$original_n_cells <- renderText({
    nrow(reactive_metadata$data[, ])
  })
  output$filtered_n_cells <- renderText({
    if (somaker_dataobject$is_valid_data) {
      nrow(reactive_metadata$data[
        reactive_metadata$data$nCount_RNA >= slider_input_vals[["nCount_RNA"]]()[1] &
          reactive_metadata$data$nCount_RNA <= slider_input_vals[["nCount_RNA"]]()[2],
      ])
    }
    # nrow(reactive_metadata$data[
    #   reactive_metadata$data$nCount_RNA >= slider_input_vals[["nCount_RNA"]]()[1] &
    #     reactive_metadata$data$nCount_RNA <= slider_input_vals[["nCount_RNA"]]()[2],
    # ])
  })


  # observeEvent(input$show, {
  #     showNotification(
  #       ui = markdown(c(
  #         "## This is a notification. ",
  #         "_Take care._"
  #       )),
  #       type = "warning"
  #     )
  #   })


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
