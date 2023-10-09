library(shiny)
library(shinyFiles)
library(plotly)
# library(gridlayout)
library(bslib)
library(bsicons)
library(Seurat)
library(stringr)
library(DT)
library(htmltools)


source("backend.R")

options(shiny.maxRequestSize = 300 * 1024^2) # 300 MB limit

ui <- page_navbar(
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
      # "input.nav === '<h5>Filtering</h5>'",
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
    # verbatimTextOutput(outputId = "dataobject"),
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
      shinyDirButton("directory", "Folder select", "Please select a folder"),
      # shinyFilesButton(
      #   "files",
      #   label = "File select",
      #   title = "Please select a file",
      #   multiple = FALSE
      # ),
      verbatimTextOutput(outputId = "selected_directory")
      # actionButton(
      #   inputId = "select_files_button",
      #   label = "Select a directory",
      #   width = "40%"
      # )
    ),
  ),
  nav_panel(
    value = "filtering",
    title = card_title("Filtering"),
    card_body(
      markdown(
        mds = c(
          "**Selected directory:**"
        )
      ),
      verbatimTextOutput(outputId = "selected_directory_filtering"),
      # value_box(
      #   # title = "Selected directory",
      #   value = "",
      #   title = textOutput(outputId = "selected_directory"),
      #   # value = textOutput(outputId = "selected_directory"),
      #   showcase = bsicons::bs_icon("file-earmark-spreadsheet"), # , size = NULL),
      #   theme_color = "success"
      # ),
      splitLayout(
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
      markdown(
        mds = c(
          "### Violin plots of basic characteristics"
        )
      ),
      plotlyOutput(outputId = "violin_plot"),
      # make_qc_slider(x = pbmc_small$nCount_RNA, col = "nCount_RNA"),
      sliderInput(
        inputId = str_c("qc_slider_", "nCount_RNA"),
        label = "nCount_RNA filtering",
        min = range(pbmc_small$nCount_RNA, na.rm = TRUE)[1],
        max = range(pbmc_small$nCount_RNA, na.rm = TRUE)[2],
        value = range(pbmc_small$nCount_RNA, na.rm = TRUE),
        width = "75%"
      ),
      helpText("Adjust the sliders to set the filtering cutoffs."),
      textOutput(outputId = "qc_value"),
      # actionButton("show", "Show"),
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
  )
)


server <- function(input, output, session) {
  # somaker_dataobject <- new_somaker_dataobject()
  somaker_dataobject <- reactiveValues(
      selected_directory = "",
      nCount_RNA = 6
    )
  # somaker_dataobject <- isolate(reactiveValuesToList(somaker_dataobject))

  # volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  volumes <- c(Home = getwd(), "R Installation" = R.home(), getVolumes()()) # TODO: change back to fs::path_home() when done testing

  shinyDirChoose(input, "directory",
    roots = volumes, session = session,
    restrictions = system.file(package = "base"),
    allowDirCreate = FALSE
  )

  # placeholder_dir <- file.path("/Users/johndoe/Downloads/pbmc_small")
  output$selected_directory <- renderPrint({
    if (is.integer(input$directory)) {
      cat("No directory has yet been selected.")
    } else {
      selected_directory <- parseDirPath(volumes, input$directory)
      somaker_dataobject$selected_directory <- selected_directory
      # somaker_dataobject(selected_directory)
      return(selected_directory)
    }
  })
  output$selected_directory_filtering <- renderPrint({
    if (is.integer(input$directory)) {
      cat("No directory has yet been selected.")
    } else {
      selected_directory <- parseDirPath(volumes, input$directory)
      somaker_dataobject$selected_directory <- selected_directory
      return(selected_directory)
    }
  })
  # output$dataobject <- renderUI({somaker_dataobject()})
  # output$dataobject <- renderUI({somaker_dataobject$dir})
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
  # output$dataobject <- renderText({
  #     bar <- reactiveValuesToList(somaker_dataobject)
  #     result <- ""
  #     for (i in seq_along(bar)) {
  #       result <- paste(result, names(bar)[i], ": ", bar[[i]], "\n")
  #     }
  #     return(result)
  # })
  # output$dataobject <- renderUI({
  #   values <- somaker_dataobject
  #   if (length(values) > 0) {
  #     # tagList(
  #     #   p(reactiveValuesToList(somaker_dataobject))
  #     #   # lapply(values, function(value) {
  #     #   #   p(typeof(value))
  #     #   # })
  #     # )
  #   } else {
  #     p("List is empty.")
  #   }
  # })
  # print(somaker_dataobject)
  

  output$nav_now <- renderText(input$nav)
  output$qc_value <- renderText(input$qc_slider_nCount_RNA)

  reactive_metadata <- reactive({
    pbmc_small[[]] # [
    # pbmc_small[[]]$nCount_RNA > input$qc_slider_nCount_RNA[1] |
    # pbmc_small[[]]$nCount_RNA < input$qc_slider_nCount_RNA[2]
    # ,]
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


  # input$qc_slider_nCount_RNA
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
