#' Initialise the SSOMaker in a browser
#' @name SeuratObjectMaker
#' @description Start the SSOmaker in a browser so you can import data and process it.
#' @export SeuratObjectMaker
SeuratObjectMaker <- function(
    appDir = getwd(),
    port = getOption("shiny.port"),
    launch.browser = getOption("shiny.launch.browser", interactive()),
    host = getOption("shiny.host", "127.0.0.1"),
    workerId = "",
    quiet = FALSE,
    display.mode = c("auto", "normal", "showcase"),
    test.mode = getOption("shiny.testmode", FALSE)) {
  ui <- htmltools::tagList(
    shinyjs::useShinyjs(),
    bslib::page_navbar(
      id = "nav",
      selected = "load_data",
      # selected = "filtering",
      theme = bslib::bs_theme(
        bootswatch = "pulse",
        version = 5
      ),
      title = "Seurat Object Maker",
      fillable = FALSE,
      sidebar = bslib::sidebar(
        shiny::conditionalPanel(
          "input.nav === 'load_data'",
          shiny::markdown(
            mds = c(
              "## Flow of app",
              "1. **Select a directory**",
              "2. Filter the data (QC)",
              "3. Results"
            )
          )
        ),
        shiny::conditionalPanel(
          "input.nav === 'filtering'",
          shiny::markdown(
            mds = c(
              "## Flow of app",
              "1. Select a directory",
              "2. **Filter the data (QC)**",
              "3. Results"
            )
          )
        ),
        shiny::conditionalPanel(
          "input.nav === 'results'",
          shiny::markdown(
            mds = c(
              "## Flow of app",
              "1. Select a directory",
              "2. Filter the data (QC)",
              "3. **Results**"
            )
          )
        ),
        shiny::uiOutput(outputId = "dataobject"),
        fillable = TRUE,
        open = "closed",
        position = "right"
      ),
      bslib::nav_panel(
        id = "load_data",
        title = bslib::card_title("Load data"),
        value = "load_data",
        bslib::card_body(
          shiny::markdown(
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
          htmltools::div(
            style = "display: flex; justify-content: center;",
            shinyFiles::shinyDirButton(
              "directory",
              "Folder select",
              "Please select a folder",
              icon = bsicons::bs_icon("folder"),
              buttonType = "primary",
              style = "width: 40%;"
            )
          ),
          shinyjs::hidden(shiny::uiOutput("dir_valid_UI")),
          shinyjs::hidden(shiny::uiOutput("dir_invalid_UI")),
          shinyjs::hidden(shiny::uiOutput("start_processing_UI"))
        ),
      ),
      ##############################################################################
      # Filtering tab
      ##############################################################################
      bslib::nav_panel(
        id = "filtering",
        value = "filtering",
        title = bslib::card_title("Filtering"),
        bslib::card_body(
          shiny::markdown(
            mds = c(
              "# Filtering (QC)"
            )
          ),
          shiny::markdown("---"),
          bslib::layout_column_wrap(
            width = 2,
            bslib::value_box(
              title = "Original number of cells",
              value = shiny::textOutput(outputId = "original_n_cells"),
              showcase = bsicons::bs_icon("clipboard-data"),
              theme_color = "success"
            ),
            bslib::value_box(
              title = "Number of cells left after filtering",
              value = shiny::textOutput(outputId = "filtered_n_cells"),
              showcase = bsicons::bs_icon("receipt-cutoff"),
              theme_color = "warning"
            )
          ),
          shiny::markdown(
            mds = c(
              "### Violin plots of basic characteristics"
            )
          ),
          shiny::uiOutput(outputId = "qc_UI"),
          shiny::markdown(
            mds = c(
              "### Percent mitochondrial genes vs. number of transcripts"
            )
          ),
          htmltools::div(
            style = "display: flex; justify-content: center;",
            bslib::tooltip(
              htmltools::span(
                shiny::helpText("Current limitations of this plot "),
                bsicons::bs_icon("info-circle")
              ),
              stringr::str_c(
                "For the moment the selection done in this plot does ",
                "not update the sliders for the violin plots above. ",
                "Nevertheless, the selection will affect the final filtering. ",
                "The narrowest filtering thresholds will be used."
              )
            )
          ),
          htmltools::div(
            style = "display: flex; justify-content: center;",
            htmltools::div(
              style = "width: 50%",
              plotly::plotlyOutput("perc_mt_scatter") # ,
            )
          ),
          # bslib::layout_column_wrap(
          #   width = 4,
          shiny::uiOutput(outputId = "filtering_thresholds"),
          # ),
          shiny::uiOutput(outputId = "confirm_filtering_UI"),
          # shiny::markdown(
          #   mds = c(
          #     "### Choose filtering parameters"
          #   )
          # ),
          # shiny::textOutput(outputId = "filtered_dimensions"),
          # shiny::markdown(
          #   mds = c(
          #     "Showing the metadata for the dataset, in order to help choose which columns to filter. Some suggestions have been selected in the checkboxes below."
          #   )
          # ),
          # DT::dataTableOutput("metadata")
        ) # ,
      ),
      ##############################################################################
      # Results tab
      ##############################################################################
      bslib::nav_panel(
        id = "results",
        value = "results",
        title = bslib::card_title("Results"),
        shiny::markdown(
          mds = c(
            "# Results",
            "---"
          )
        ),
        shiny::markdown(stringr::str_c("## QC results")),
        htmltools::br(),
        shiny::uiOutput("final_n_cells"),
        htmltools::br(),
        htmltools::h3("Distribution of metadata columns post filtering"),
        bslib::layout_column_wrap(
          width = "600px",
          shiny::uiOutput("results_violin_plots"),
          # htmltools::div(
          #   style = "display: flex; justify-content: center;",
          #   htmltools::div(
          #     violin_plot_ui("vln_nCount_RNA"),
          #     style = "width: 60%; max-width: 600px;",
          #   )
          # )
        ),
        htmltools::br(),
        shiny::markdown(stringr::str_c(
          "### Elbow plot\n\n",
          "**PCA evaluation**"
        )),
        bslib::tooltip(
          htmltools::span(
            shiny::helpText("What does an elbow plot show? "),
            bsicons::bs_icon("info-circle")
          ),
          stringr::str_c(
            "An elbow plot in this context shows the PCs' (principal components') ",
            "standard deviation. We do PCA (principal component analysis) to ",
            "reduce the number of dimensions in the dataset. PCA captures the directions in ",
            "the dataset that have the most variance. The elbow plot shows this variation, ",
            "and we specifically want to look for an \"elbow\" in the plot, which is ",
            "the point where the \"line\" tapers off. This is the point where increasing the ",
            "number of PCs doesn't add much information, and we can stop there."
          )
        ),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          htmltools::div(
            shiny::plotOutput("elbow_plot"),
            style = "width: 100%; max-width: 800px;"
          )
        ),
        shiny::markdown(stringr::str_c(
          "### Volcano plots of highly variable genes (HVGs)\n\n"
        )),
        # plotOutput("hvg_plot"),
        # htmltools::div(
        #   style = "display: flex; justify-content: center;",
        #   htmltools::div(
        #     plotOutput("hvg_plot"),
        #     # plotOutput("pca_plot"),
        #     style = "width: 60%; max-width: 1200px;",
        #   )
        # ),
        bslib::layout_column_wrap(
          width = 2,
          # width = "400px",
          # style = "display: flex; justify-content: center;",
          htmltools::div(
            style = "display: flex; justify-content: center;",
            htmltools::div(
              shiny::plotOutput("hvg_plot_unlabelled"),
              style = "width: 100%; max-width: 600px;"
            )
          ),
          htmltools::div(
            style = "display: flex; justify-content: center;",
            htmltools::div(
              shiny::plotOutput("hvg_plot_labelled"),
              style = "width: 100%; max-width: 600px;"
            )
          )
        ),
        htmltools::br(),
        shiny::markdown(stringr::str_c(
          "### Scatter plot of the first two PCs\n\n"
        )),
        # shiny::plotOutput("pca_plot"),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          htmltools::div(
            shiny::plotOutput("pca_plot"),
            style = "width: 60%; max-width: 600px;"
          )
        ),
        shiny::markdown(stringr::str_c(
          "---\n\n",
          "## More informative plots",
          "\n\n",
          "### UMAP"
        )),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          htmltools::div(
            shiny::plotOutput("umap_plot"),
            style = "width: 60%; max-width: 600px;"
          )
        ),
        shiny::markdown(stringr::str_c(
          "### Differentially expressed genes (DEGs)"
        )),
        bslib::tooltip(
          htmltools::span(
            shiny::helpText("Note on the following tables "),
            bsicons::bs_icon("info-circle")
          ),
          stringr::str_c(
            "The values in the following DEG tables are rounded to 8 decimal points. ",
            "This is for clarity reasons in this web interface. ",
            "The exported tables will have the full precision."
          )
        ),
        htmltools::br(),
        htmltools::br(),
        shiny::uiOutput("deg"),
        htmltools::br(),
        shiny::markdown(stringr::str_c(
          "---\n\n",
          "### Filtered metadata"
        )),
        DT::dataTableOutput("metadata_filtered"),
        htmltools::br(),
        shiny::markdown(stringr::str_c(
          "---\n\n",
          "# Download data"
        )),
        htmltools::br(),
        bslib::layout_column_wrap(
          width = 2,
          htmltools::div(
            style = "display: flex; justify-content: center;",
            shiny::downloadButton(
              "download_seurat_object",
              "Download Seurat Object",
              icon = shiny::icon("cloud-download"),
              class = "btn btn-primary",
              style = "width: 50%"
            )
          ),
          htmltools::div(
            style = "display: flex; justify-content: center;",
            shiny::downloadButton(
              "download_degs",
              "Download DEGs",
              icon = shiny::icon("cloud-download"),
              class = "btn btn-primary",
              style = "width: 50%"
            )
          )
        )
      ),
      ##############################################################################
      # Rest of the top bar
      ##############################################################################
      bslib::nav_spacer(),
      bslib::nav_item(
        shiny::markdown(
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

    somaker_dataobject <- shiny::reactiveValues(
      selected_directory = character(0),
      is_valid_directory = FALSE,
      is_valid_data = FALSE,
      nCount_RNA_low = 6,
      nCount_RNA_high = 600,
      nFeature_RNA_high = 600,
      nFeature_RNA_low = 200,
      percent_mt_low = 0,
      percent_mt_high = 20,
      columns_to_filter = list(
        "nCount_RNA",
        "nFeature_RNA"
      )
    )

    # reactive_metadata <- reactiveValues(data = data.frame())
    reactive_metadata <- shiny::reactiveValues(data = SeuratObject::pbmc_small[[]])
    reactive_sobj <- shiny::reactiveValues(data = SeuratObject::pbmc_small)

    bslib::nav_hide(id = "nav", target = "filtering")
    bslib::nav_hide(id = "nav", target = "results")

    # # debugging
    # bslib::nav_show(id = "nav", target = "results")
    # bslib::nav_select(id = "nav", selected = "results")

    ################################################################################
    # Sidebar
    ################################################################################

    output$nav_now <- shiny::renderText(input$nav)

    output$dataobject <- shiny::renderUI({
      bar <- shiny::reactiveValuesToList(somaker_dataobject)
      result <- htmltools::tagList(tags$h3("Data object"))
      for (i in seq_along(bar)) {
        result <- htmltools::tagList(
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

    volumes <- c(
      Home = getwd(),
      "R Installation" = R.home(),
      shinyFiles::getVolumes()()
    ) # TODO: change back to fs::path_home() when done testing

    withr::with_options(
      new = list(shiny.maxRequestSize = 900 * 1024^2), # 900 MB limit
      code = shinyFiles::shinyDirChoose(
        input,
        "directory",
        roots = volumes,
        session = session,
        restrictions = system.file(package = "base"),
        allowDirCreate = FALSE
      )
    )



    reactive_features <- shiny::reactiveValues(
      data = SeuratObject::pbmc_small[[]]
    )

    # Validate and store the selected directory
    shiny::observeEvent(input$directory, {
      if (!is.integer(input$directory)) {
        selected_directory <- shinyFiles::parseDirPath(volumes, input$directory)
        somaker_dataobject$selected_directory <- selected_directory
        if (length(somaker_dataobject$selected_directory) > 0) {
          is_valid_dir <- validate_directory(
            somaker_dataobject$selected_directory
          )
          if (is_valid_dir) {
            filenames <- c("features.tsv", "genes.tsv")
            for (suffix in c("", ".gz")) {
              for (fn in stringr::str_c(filenames, suffix)) {
                if (file.exists(
                  file.path(somaker_dataobject$selected_directory, fn)
                )) {
                  reactive_features$data <- data.table::fread(
                    file = file.path(somaker_dataobject$selected_directory, fn),
                    header = FALSE,
                    nrows = 5
                  )
                }
              }
            }
            shinyjs::hide("dir_invalid_UI")
            shinyjs::show("dir_valid_UI")
          } else {
            shinyjs::hide("dir_valid_UI")
            shinyjs::show("dir_invalid_UI")
          }
          somaker_dataobject$is_valid_directory <- is_valid_dir
        } else {
          shinyjs::hide("dir_valid_UI")
          shinyjs::hide("dir_invalid_UI")
        }
        # somaker_dataobject$selected_directory <- shinyFiles::parseDirPath(roots = c(wd = getwd()), input$directory)
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

    output$features_tsv <- DT::renderDataTable({
      if (somaker_dataobject$is_valid_directory) {
        return(DT::datatable(
          reactive_features$data,
          options = list(
            dom = "t"
          )
        ))
      }
    })

    output$feature_column_selection <- shiny::renderUI({
      if (somaker_dataobject$is_valid_directory) {
        if (shiny::is.reactivevalues(reactive_features) &
          !is.null(reactive_features$data)
        ) {
          return(
            htmltools::tagList(
              shiny::radioButtons(
                inputId = "feat_radiobuttons",
                label = shiny::markdown(stringr::str_c(
                  "Select which column contains the gene names:"
                )),
                choices = colnames(reactive_features$data),
                selected = colnames(reactive_features$data)[2],
                inline = TRUE
              )
            )
          )
        }
      }
    })

    output$dir_valid_UI <- shiny::renderUI({
      htmltools::div(
        htmltools::br(),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          bslib::card(
            bslib::card_header(htmltools::h3("Success!")),
            shiny::markdown(
              stringr::str_c(
                "`", somaker_dataobject$selected_directory, "`",
                " is a valid count matrix directory."
              )
            ),
            class = "card text-bg-success",
            style = "width: 30%"
          ),
        ),
        htmltools::br(),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          bslib::card(
            bslib::card_header("Feature file has these columns:"),
            htmltools::div(
              DT::dataTableOutput("features_tsv"),
              # style = "width: 60%"
            ),
            # ),
            # htmltools::div(
            #   style = "display: flex; justify-content: center;",
            htmltools::div(
              shiny::uiOutput("feature_column_selection"),
              # style = "width: 60%"
            ),
            style = "width: 60%"
          ),
        ),
        htmltools::br(),
        shiny::markdown(stringr::str_c(
          "Press button below to load the data into `Seurat`."
        )),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          shiny::actionButton(
            "load_into_seurat",
            label = htmltools::HTML(
              bsicons::bs_icon("box-arrow-in-up-right"),
              " Load data into Seurat"
            ),
            class = "btn-warning",
            style = "width: 40%;"
          )
        )
      )
    })

    output$dir_invalid_UI <- shiny::renderUI({
      htmltools::br()
      htmltools::div(
        style = "display: flex; justify-content: center;",
        bslib::card(
          bslib::card_header(htmltools::h4("Did not succeed.")),
          shiny::markdown(
            stringr::str_c(
              "`", somaker_dataobject$selected_directory, "`",
              " is _not_ a valid count matrix directory. \n\n",
              "Please see the description of the format above."
            )
          ),
          class = "card text-bg-danger",
          style = "width: 30%"
        ),
      )
    })

    shinyjs::onclick("load_into_seurat",
      expr = {
        shinyjs::disable("load_into_seurat")
        load_time <- system.time(
          shinycssloaders::showPageSpinner(
            type = 6,
            expr = {
              gene_column <- which(
                colnames(reactive_features$data) == input$feat_radiobuttons
              )
              shinyjs::logjs(gene_column)
              sobj_data <- read_data(
                somaker_dataobject$selected_directory,
                gene_column = gene_column
              )
              sobj <- make_seurat_object(sobj_data)
            },
            caption = htmltools::HTML(
              "Loading data into Seurat..",
              bsicons::bs_icon("tools")
            )
          )
        )
        print(load_time["elapsed"])
        if (!is.null(sobj)) {
          somaker_dataobject$is_valid_data <- TRUE
        }
      }
    )

    output$sobj_out <- shiny::renderPrint({
      reactive_sobj$data
    })

    shiny::observeEvent(somaker_dataobject$is_valid_data, {
      shiny::req(somaker_dataobject$is_valid_data)
      # if (somaker_dataobject$is_valid_data == TRUE) {
      reactive_metadata$data <- sobj[[]]
      reactive_sobj$data <- sobj
      shinyjs::show("start_processing_UI")
      # }
    })

    output$start_processing_UI <- shiny::renderUI({
      htmltools::div(
        shiny::markdown("Data loaded correctly:"),
        shiny::verbatimTextOutput("sobj_out"),
        shiny::markdown("Press button below to start processing!"),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          shiny::actionButton(
            "start_processing",
            label = htmltools::HTML(
              bsicons::bs_icon("pc-display-horizontal"),
              " Start processing"
            ),
            class = "btn-success",
            style = "width: 40%;" # ,
          )
        )
      )
    })

    slider_input_vals <- shiny::reactiveValues()

    shinyjs::onclick("start_processing",
      expr = {
        bslib::nav_hide(id = "nav", target = "load_data")
        bslib::nav_show(id = "nav", target = "filtering")
        bslib::nav_select(id = "nav", selected = "filtering")
        for (column in colnames(reactive_metadata$data)) {
          if (is.numeric(reactive_metadata$data[[column]])) {
            values <- as.vector(
              round(
                stats::quantile(
                  reactive_metadata$data[[column]],
                  probs = c(0.05, 0.95)
                )
              )
            )
            # slider_input_vals[[column]] <- values
            if (values[1] != values[2]) {
              if (!column %in% somaker_dataobject$columns_to_filter) {
                somaker_dataobject$columns_to_filter <- c(
                  somaker_dataobject$columns_to_filter,
                  column
                )
              }
              print(column)
              print(values)
              somaker_dataobject[[stringr::str_c(column, "_low")]] <- values[1]
              somaker_dataobject[[stringr::str_c(column, "_high")]] <- values[2]
              slider_input_vals[[column]] <- qc_module_server(
                stringr::str_c("qc_", column),
                col = column,
                metadata = reactive_metadata,
                start_values = values
              )
            }
          }
          somaker_dataobject$original_n_cells <- nrow(reactive_metadata$data)
          output$qc_UI <- shiny::renderUI({
            htmltools::tagList(
              lapply(
                somaker_dataobject$columns_to_filter,
                FUN = \(column) {
                  qc_module_UI(stringr::str_c("qc_", column))
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

    output$original_n_cells <- shiny::renderText({
      nrow(reactive_metadata$data[, ])
    })
    output$filtered_n_cells <- shiny::renderText({
      if (somaker_dataobject$is_valid_data) {
        cell_filters <- list()
        for (column in somaker_dataobject$columns_to_filter) {
          # print(column)
          # reactive_metadata$data[[column]] %>% head %>% print
          # slider_input_vals[[column]]()[1] %>% print
          ## TODO: more reasonable check for type of column.
          ## Should probably have a better data model for filterable columns.
          ## E.g. "Is it supposed to only have an upper cutoff? Or a lower? Or both?"
          if (length(slider_input_vals[[column]]()) > 1) {
            cell_filters <- c(
              cell_filters,
              list(reactive_metadata$data[[column]] >= slider_input_vals[[column]]()[1])
            )
          }
          cell_filters <- c(
            cell_filters,
            list(reactive_metadata$data[[column]] <= slider_input_vals[[column]]()[2])
          )
        }
        # print(cell_filters)
        master_filter <- Reduce(`&`, cell_filters)
        nrow(reactive_metadata$data[
          master_filter,
        ])
      }
    })


    output$confirm_filtering_UI <- shiny::renderUI({
      htmltools::div(
        # shiny::markdown("Data loaded correctly:"),
        # shiny::verbatimTextOutput("sobj_out"),
        # shiny::markdown("Press button below to start processing!"),
        # shiny::verbatimTextOutput("filtering_thresholds"),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          shiny::actionButton(
            "confirm_filtering",
            label = htmltools::HTML(
              bsicons::bs_icon("bar-chart-fill"),
              "<strong>",
              " Confirm filtering thresholds",
              "</strong>"
            ),
            class = "btn-success",
            style = "width: 40%;" # ,
          )
        )
      )
    })

    # event_data <- shiny::reactiveVal(0)

    output$perc_mt_scatter <- plotly::renderPlotly({
      # shiny::req(reactive_metadata$data)
      # colours <- shiny::reactive({
      #   ifelse(
      #     reactive_metadata$data$percent_mt >=
      #       input$`qc_percent_mt-qc_slide-slider`[1] &
      #       reactive_metadata$data$percent_mt <=
      #         input$`qc_percent_mt-qc_slide-slider`[2],
      #     "red", "blue"
      #   )
      # })
      p <- plotly::plot_ly(
        data = reactive_metadata$data,
        x = ~nCount_RNA,
        y = ~percent_mt,
        type = "scatter",
        customdata = ~ rownames(reactive_metadata$data),
        source = "perc_mt_plot",
        mode = "markers",
        # color = ~ colours(),
        # colors = "identity",
        name = "perc_mt_vs_ncount_rna" # ,
        # x0 = col
      )
      # plotly::event_register(p, "plotly_selected")
      # removing lasso selection for now, not compatible with the mental model
      # from the violin plots and already set up info boxes for the min/max
      # values per metadata
      plotly::config(p, modeBarButtonsToRemove = list("lasso2d"))

      # this chunk is probably better to have in another place,
      # having it here resulted in extremely long response times in the UI
      # results <- plotly::event_data(
      #     "plotly_selected",
      #     source = "perc_mt_plot"
      #   )
      # event_data(results)

      # if (!is.null(shiny::isolate(event_data())$y)) {
      # # if (max(event_data()$y) != Inf & min(event_data()$y) != -Inf) {
      #   y_min <- min(event_data()$y)
      #   y_max <- max(event_data()$y)
      #   shiny::updateSliderInput(
      #     session,
      #     inputId = "qc_percent_mt-qc_slide-slider",
      #     value = range(y_min, y_max)
      #   )
      # }

      plotly::layout(
        p = p,
        yaxis = list(zeroline = FALSE),
        dragmode = "select"
      )
    })

    # output$foo <- shiny::renderPrint({
    #   # shiny::req(wait_for_scatter_plot)
    #   # if (wait_for_scatter_plot() == TRUE) {
    #   if (!is.null(event_data())) {
    #     df <- reactive_metadata$data[event_data()$customdata, ]
    #     print(dim(df))
    #     # slider_input_vals$percent_mt <<- range(y_min, y_max)
    #     print(head(df))
    #     print(head(event_data()))
    #   }
    #   # }
    # })

    output$filtering_thresholds <- shiny::renderUI({
      if (somaker_dataobject$is_valid_data) {
        cell_filters <- list()
        for (column in somaker_dataobject$columns_to_filter) {
          if (length(slider_input_vals[[column]]()) > 1) {
            cell_filters <- c(
              cell_filters,
              list(
                reactive_metadata$data[[column]] >=
                  slider_input_vals[[column]]()[1]
              )
            )
          }
          cell_filters <- c(
            cell_filters,
            list(
              reactive_metadata$data[[column]] <=
                slider_input_vals[[column]]()[2]
            )
          )
        }
        results <- plotly::event_data(
          "plotly_selected",
          source = "perc_mt_plot"
        )
        if (!is.null(results$y)) {
          cell_filters <- c(
            cell_filters,
            list(
              reactive_metadata$data[["percent_mt"]] <=
                max(results$y)
            ),
            list(
              reactive_metadata$data[["percent_mt"]] >=
                min(results$y)
            ),
            list(
              reactive_metadata$data[["nCount_RNA"]] <=
                max(results$x)
            ),
            list(
              reactive_metadata$data[["nCount_RNA"]] >=
                min(results$x)
            )
          )
        }

        master_filter <- Reduce(`&`, cell_filters)

        to_return <- htmltools::tagList()
        j <- 1
        for (i in seq_len(length(somaker_dataobject$columns_to_filter))) {
          high_value <- slider_input_vals[[
            somaker_dataobject$columns_to_filter[[i]]
          ]]()[2]
          low_value <- slider_input_vals[[
            somaker_dataobject$columns_to_filter[[i]]
          ]]()[1]
          if (somaker_dataobject$columns_to_filter[[i]] == "percent_mt") {
            if (!is.null(results$y)) {
              high_value <- min(
                max(results$y),
                slider_input_vals[[
                  somaker_dataobject$columns_to_filter[[i]]
                ]]()[2]
              )
              low_value <- max(
                min(results$y),
                slider_input_vals[[
                  somaker_dataobject$columns_to_filter[[i]]
                ]]()[1]
              )
            }
          } else if (
            somaker_dataobject$columns_to_filter[[i]] == "nCount_RNA"
          ) {
            if (!is.null(results$x)) {
              high_value <- min(
                max(results$x),
                slider_input_vals[[
                  somaker_dataobject$columns_to_filter[[i]]
                ]]()[2]
              )
              low_value <- max(
                min(results$x),
                slider_input_vals[[
                  somaker_dataobject$columns_to_filter[[i]]
                ]]()[1]
              )
            }
          } else {
            low_value <- slider_input_vals[[
              somaker_dataobject$columns_to_filter[[i]]
            ]]()[1]
            high_value <- slider_input_vals[[
              somaker_dataobject$columns_to_filter[[i]]
            ]]()[2]
          }
          to_return[[j]] <- bslib::value_box(
            title = shiny::markdown(stringr::str_c(
              "**", somaker_dataobject$columns_to_filter[[i]], "_low", "**"
            )),
            value = low_value,
            showcase = bsicons::bs_icon("filter-left"),
            theme_color = "primary"
            # theme_color = "secondary"
          )
          j <- j + 1

          to_return[[j]] <- bslib::value_box(
            title = shiny::markdown(stringr::str_c(
              "**", somaker_dataobject$columns_to_filter[[i]], "_high", "**"
            )),
            value = high_value,
            showcase = bsicons::bs_icon("filter-right"),
            theme_color = "primary"
          )
          j <- j + 1
        }

        return(
          htmltools::tagList(
            bslib::layout_column_wrap(
              width = "240px",
              !!!to_return,
              heights_equal = "all",
              fill = TRUE
            ),
            bslib::layout_column_wrap(
              width = 2,
              bslib::value_box(
                title = "Original number of cells",
                value = nrow(reactive_metadata$data),
                showcase = bsicons::bs_icon("clipboard-data"),
                theme_color = "success"
              ),
              bslib::value_box(
                title = "Number of cells left after filtering",
                value = nrow(reactive_metadata$data[master_filter, ]),
                showcase = bsicons::bs_icon("receipt-cutoff"),
                theme_color = "warning"
              ),
              heights_equal = "all",
              fill = TRUE
            )
          )
        )
      }
    })

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

    ############################################################################
    # Results tab
    ############################################################################

    shinyjs::onclick("confirm_filtering",
      expr = {
        shinyjs::disable("confirm_filtering")

        bslib::nav_hide(id = "nav", target = "load_data")
        bslib::nav_show(id = "nav", target = "results")
        bslib::nav_select(id = "nav", selected = "results")
        bslib::nav_hide(id = "nav", target = "filtering")


        ## NB: duplicated code, here be dragons
        cell_filters <- list()
        for (column in somaker_dataobject$columns_to_filter) {
          ## TODO: more reasonable check for type of column.
          ## Should probably have a better data model for filterable columns.
          ## E.g. "Is it supposed to only have an upper cutoff? Or a lower? Or both?"
          ## TODO: maybe just have a lower and a higher threshold for non-complexity's sake,
          ## and if only one is applicable set the other to 0 or Inf, respectively.
          if (length(slider_input_vals[[column]]()) > 1) {
            cell_filters <- c(
              cell_filters,
              list(
                reactive_metadata$data[[column]] >=
                  slider_input_vals[[column]]()[1]
              )
            )
          }
          cell_filters <- c(
            cell_filters,
            list(
              reactive_metadata$data[[column]] <=
                slider_input_vals[[column]]()[2]
            )
          )
        }
        results <- plotly::event_data(
          "plotly_selected",
          source = "perc_mt_plot"
        )
        if (!is.null(results$y)) {
          cell_filters <- c(
            cell_filters,
            list(
              reactive_metadata$data[["percent_mt"]] <=
                max(results$y)
            ),
            list(
              reactive_metadata$data[["percent_mt"]] >=
                min(results$y)
            ),
            list(
              reactive_metadata$data[["nCount_RNA"]] <=
                max(results$x)
            ),
            list(
              reactive_metadata$data[["nCount_RNA"]] >=
                min(results$x)
            )
          )
        }
        master_filter <- Reduce(`&`, cell_filters)
        somaker_dataobject$filtered_n_cells <- nrow(reactive_metadata$data[
          master_filter,
        ])
        reactive_metadata$data <- reactive_metadata$data[
          master_filter,
        ]
        shinycssloaders::showPageSpinner(
          type = 6,
          expr = {
            reactive_sobj$data <- subset(
              reactive_sobj$data,
              cells = rownames(reactive_metadata$data)
            )
          },
          caption = htmltools::HTML(
            "Subsetting data..",
            bsicons::bs_icon("scissors")
          )
        )

        for (md_column in somaker_dataobject$columns_to_filter) {
          somaker_dataobject[[
            stringr::str_c(md_column, "_low")
          ]] <- slider_input_vals[[md_column]]()[1]
          somaker_dataobject[[
            stringr::str_c(md_column, "_high")
          ]] <- slider_input_vals[[md_column]]()[2]
        }

        shinycssloaders::showPageSpinner(
          type = 6,
          expr = {
            reactive_sobj$data <- process_seurat_object(reactive_sobj$data)
            degs <- shiny::reactiveValues(data = find_degs(reactive_sobj$data))
          },
          caption = htmltools::HTML(
            "Processing data, hold on..",
            bsicons::bs_icon("clipboard2-pulse"),
            bsicons::bs_icon("hourglass-split")
          )
        )
        reactive_metadata$data <- reactive_sobj$data@meta.data
      }
    )

    output$original_n_cells_results <- shiny::renderText({
      somaker_dataobject$original_n_cells
    })
    output$filtered_n_cells_results <- shiny::renderText({
      somaker_dataobject$filtered_n_cells
    })
    output$final_n_cells <- shiny::renderUI({
      htmltools::tagList(
        bslib::layout_column_wrap(
          width = 2,
          bslib::value_box(
            title = "Original number of cells",
            value = somaker_dataobject$original_n_cells,
            showcase = bsicons::bs_icon("clipboard-data"),
            theme_color = "success"
          ),
          bslib::value_box(
            title = "Number of cells left after filtering",
            value = somaker_dataobject$filtered_n_cells,
            showcase = bsicons::bs_icon("receipt-cutoff"),
            theme_color = "warning"
          ),
          heights_equal = "all",
          fill = TRUE
        )
      )
    })

    output$results_violin_plots <- shiny::renderUI({
      to_return <- htmltools::tagList()
      for (i in seq_len(length(somaker_dataobject$columns_to_filter))) {
        violin_plot_server(
          stringr::str_c("vln_", somaker_dataobject$columns_to_filter[[i]]),
          somaker_dataobject$columns_to_filter[[i]],
          reactive_metadata
        )
        to_return[[i]] <- htmltools::div(
          style = "display: flex; justify-content: center;",
          htmltools::div(
            shiny::markdown(stringr::str_c(
              "**", somaker_dataobject$columns_to_filter[[i]], "**"
            )),
            violin_plot_ui(
              stringr::str_c("vln_", somaker_dataobject$columns_to_filter[[i]])
            ),
            style = "width: 60%; max-width: 600px;",
          )
        )
      }
      return(to_return)
    })

    output$elbow_plot <- shiny::renderPlot({
      Seurat::ElbowPlot(reactive_sobj$data)
    })

    output$hvg_plot_unlabelled <- shiny::renderPlot({
      top10 <- head(Seurat::VariableFeatures(reactive_sobj$data), 10)
      return(Seurat::VariableFeaturePlot(reactive_sobj$data))
    })
    output$hvg_plot_labelled <- shiny::renderPlot({
      top10 <- head(Seurat::VariableFeatures(reactive_sobj$data), 10)
      plot1 <- Seurat::VariableFeaturePlot(reactive_sobj$data)
      plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
      return(plot2)
    })

    output$pca_plot <- shiny::renderPlot({
      Seurat::DimPlot(reactive_sobj$data, reduction = "pca", group.by = "orig.ident")
    })

    output$umap_plot <- shiny::renderPlot({
      Seurat::DimPlot(reactive_sobj$data, reduction = "umap", group.by = "seurat_clusters")
    })

    output$deg <- shiny::renderUI({
      stopifnot(!is.null(degs$data))
      print(head(degs$data))
      # shiny::verbatimTextOutput(degs %>% head())
      to_return <- list()
      # print(reactive_metadata$data$seurat_clusters %>% levels() %>% as.numeric())
      panels <- lapply(
        levels(
          reactive_metadata$data$seurat_clusters
        ),
        FUN = \(i) {
          print(i)
          out_id <- stringr::str_c("deg_", i)
          deg_cluster <- degs$data[degs$data$cluster == as.character(i), ]
          deg_cluster[, sapply(deg_cluster, is.numeric)] <- round(deg_cluster[, sapply(deg_cluster, is.numeric)], 8)
          # output[[out_id]] <- renderPrint(head(degs$data[degs$data$cluster == as.character(i),]))
          return(bslib::accordion_panel(
            stringr::str_c("Genes identifying cluster ", i, ":"),
            # shiny::verbatimTextOutput(out_id),
            DT::renderDataTable(
              DT::datatable(
                deg_cluster,
                caption = stringr::str_c("Cluster ", i),
                fillContainer = FALSE,
                width = "100%"
              )
            ),
            icon = bsicons::bs_icon("table")
            # icon = bsicons::bs_icon("sign-intersection-fill")
          ))
        }
      )
      bslib::accordion(
        !!!panels,
        id = "deg_accordion",
        title = "The top DEGs for the different clusters",
        multiple = FALSE
      )
    })

    # output$deg <- shiny::renderUI({
    #   stopifnot(!is.null(degs))
    #   print(head(degs))
    #   # shiny::verbatimTextOutput(degs %>% head())
    #   to_return <- tagList()
    #   j <- 1
    #   print(reactive_metadata$data$seurat_clusters %>% levels %>% as.numeric)
    #   for (i in reactive_metadata$data$seurat_clusters %>% levels %>% as.numeric) {
    #     print(i)
    #     out_id <- stringr::str_c("deg_", i)
    #     print(out_id)
    #     print(head(degs[degs$cluster == as.character(i),]))
    #     output[[out_id]] <- renderPrint(head(degs[degs$cluster == as.character(i),]))
    #     to_return[[j]] <- shiny::markdown(stringr::str_c(
    #       "DEGs for cluster ", i, ":"
    #     ))
    #     j <- j + 1
    #     to_return[[j]] <- shiny::verbatimTextOutput(out_id)
    #     j <- j + 1
    #   }
    #   return(to_return)
    # })

    output$metadata_filtered <- DT::renderDataTable({
      shinyjs::logjs(
        dim(
          reactive_sobj$data@meta.data
        )
      )
      shinyjs::logjs(
        dim(
          reactive_metadata$data
        )
      )
      DT::datatable(
        reactive_sobj$data@meta.data,
        caption = "Metadata",
        # options = list(
        #   pageLength = 5
        # ),
        fillContainer = F
        # style = 'bootstrap'
      )
    })


    output$download_seurat_object <- shiny::downloadHandler(
      filename = function() {
        stringr::str_c("SOM_output_", Sys.Date(), ".rds")
      },
      content = function(filename) {
        saveRDS(reactive_sobj$data, filename)
      }
    )

    output$download_degs <- shiny::downloadHandler(
      filename = function() {
        stringr::str_c("SOM_output_DEGs_", Sys.Date(), ".csv")
      },
      content = function(filename) {
        write.table(x = degs$data, file = filename)
      }
    )
  }

  return(
    shiny::shinyApp(
      ui,
      server,
      options = list(
        appDir = appDir,
        port = port,
        launch.browser = launch.browser,
        host = host,
        workerId = workerId,
        quiet = quiet,
        display.mode = display.mode,
        test.mode = test.mode
      )
    )
  )
}
