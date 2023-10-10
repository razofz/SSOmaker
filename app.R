library(shiny)
library(shinyFiles)
library(shinyjs)
library(plotly)
# library(gridlayout)
library(bslib)
library(bsicons)
library(Seurat)
library(stringr)
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
        # layout_column_wrap(
        #   width = 1,
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
        # div(
        #   style = "display: flex; justify-content: center;",
        # verbatimTextOutput(outputId = "selected_directory"),
        # actionButton("do", "Click Me")
        # )
        # )
        uiOutput("dir_validation_UI")
      ),
    ),
    nav_panel(
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
              plotlyOutput(outputId = "violin_plot"),
              div(
                style = "display: flex; justify-content: center;",
                sliderInput(
                  inputId = str_c("qc_slider_", "nCount_RNA"),
                  label = "nCount_RNA",
                  min = range(pbmc_small$nCount_RNA, na.rm = TRUE)[1],
                  max = range(pbmc_small$nCount_RNA, na.rm = TRUE)[2],
                  value = range(pbmc_small$nCount_RNA, na.rm = TRUE),
                  width = "80%"
                ),
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
            ),
          )
        ),
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
  somaker_dataobject <- reactiveValues(
    selected_directory = character(0),
    is_valid_directory = FALSE,
    nCount_RNA_low = 6,
    nCount_RNA_high = 600,
    nFeature_RNA_high = 600,
    nFeature_RNA_low = 200,
    percent_mito_high = 20
  )

  # volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
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
      # somaker_dataobject$selected_directory <- parseDirPath(roots = c(wd = getwd()), input$directory)
    }
  })

  # Display new UI components only if a valid directory is selected
  output$dir_validation_UI <- renderUI({
    if (length(somaker_dataobject$selected_directory) > 0) {
      is_valid_dir <- validate_directory(somaker_dataobject$selected_directory)
      somaker_dataobject$is_valid_directory <- is_valid_dir
      if (is_valid_dir) {
        somaker_dataobject$start_seurat_processing <- TRUE
        return(
          div(
            markdown(str_c(
              "`", somaker_dataobject$selected_directory, "`",
              " is a valid count matrix directory.\n\n",
              "Loading data into `Seurat`.."
            ))
          )
        )
      } else if (!is_valid_dir) {
        return(
          div(
            markdown(str_c(
              "`", somaker_dataobject$selected_directory, "`",
              " is _not_ a valid count matrix directory. ",
              "Please see the description of the format above."
            ))
          )
        )
      } else {
        return(NULL)
      }
    } else {
      return(NULL) # Don't show new UI if no directory is selected
    }
  })

  observeEvent(somaker_dataobject$start_seurat_processing, {
    if (somaker_dataobject$start_seurat_processing == TRUE) {
      somaker_dataobject$start_seurat_processing <- FALSE
      print("Starting Seurat processing")
      sobj_data <- read_data(somaker_dataobject$selected_directory)
      print(sobj_data[1:5, 1:5])
      sobj <- make_seurat_object(sobj_data)
      print(sobj)
      print(sobj[[]] %>% head())
    }
  })

  # observeEvent(input$do, {
  #   session$sendCustomMessage(type = 'testmessage',
  #     message = 'Thank you for clicking')
  # })

  # output$selected_directory <- renderPrint({
  #   if (is.integer(input$directory)) {
  #     cat("No directory has yet been selected.")
  #   } else {
  #     selected_directory <- parseDirPath(volumes, input$directory)
  #     somaker_dataobject$selected_directory <- selected_directory
  #     return(selected_directory)
  #   }
  # })

  # output$selected_directory_filtering <- renderPrint({
  #   if (length(somaker_dataobject$selected_directory) > 0) {
  #     somaker_dataobject$selected_directory
  #   } else {
  #     cat("No directory has yet been selected.")
  #   }
  # })
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


  output$nav_now <- renderText(input$nav)
  output$qc_value <- renderText(input$qc_slider_nCount_RNA)

  reactive_metadata <- reactive({
    pbmc_small[[]]
  })

  output$original_n_cells <- renderText({
    nrow(reactive_metadata()[, ])
  })
  output$filtered_n_cells <- renderText({
    nrow(reactive_metadata()[
      pbmc_small[[]]$nCount_RNA >= input$qc_slider_nCount_RNA[1] &
        pbmc_small[[]]$nCount_RNA <= input$qc_slider_nCount_RNA[2],
    ])
  })

  # output$original_n_cells_side <- renderText({
  #   nrow(reactive_metadata()[, ])
  # })
  # output$filtered_n_cells_side <- renderText({
  #   nrow(reactive_metadata()[
  #     pbmc_small[[]]$nCount_RNA >= input$qc_slider_nCount_RNA[1] &
  #       pbmc_small[[]]$nCount_RNA <= input$qc_slider_nCount_RNA[2],
  #   ])
  # })

  output$filtered_dimensions <- renderText({
    str_c(
      "Number of cells before filtering: ",
      nrow(reactive_metadata()[pbmc_small[[]]$nCount_RNA >= input$qc_slider_nCount_RNA[1], ]),
      " cells. ",
      "Cells left after filtering: ",
      nrow(reactive_metadata()[pbmc_small[[]]$nCount_RNA <= input$qc_slider_nCount_RNA[2], ]),
      " cells."
    )
  })


  output$violin_plot <- renderPlotly({
    pbmc_small[[]] %>%
      plot_ly(
        y = as.formula(str_c(" ~ ", "nCount_RNA")),
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
          hline(input$qc_slider_nCount_RNA[2])
          # hline(max_cutoff)
        )
      ) # %>%
    # config(
    #   selectmode = "lasso",
    # )
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
  output$colnames_output <- renderText(colnames(pbmc_small[[]]))
  output$metadata <- DT::renderDataTable({
    DT::datatable(
      pbmc_small[[]],
      caption = "Metadata",
      options = list(
        pageLength = 5
      ),
      fillContainer = T
      # style = 'bootstrap'
    )
  })
  # print(session$clientData)
}

shinyApp(ui, server)
