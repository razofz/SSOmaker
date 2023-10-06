library(shiny)
library(shinyFiles)
library(plotly)
library(gridlayout)
library(bslib)
library(Seurat)
library(stringr)
library(DT)
library(htmltools)


source("../backend.R")

ui <- page_navbar(
  id = "nav",
  selected = "filtering",
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
      actionButton(
        inputId = "select_files_button",
        label = "Select a directory",
        width = "40%"
      )
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
      textOutput(outputId = "selected_directory"),
      # value_box(
      #   # title = "Selected directory",
      #   value = "",
      #   title = textOutput(outputId = "selected_directory"),
      #   # value = textOutput(outputId = "selected_directory"),
      #   showcase = bsicons::bs_icon("file-earmark-spreadsheet"), # , size = NULL),
      #   theme_color = "success"
      # ),
      markdown(
        mds = c(
          "### Violin plots of basic characteristics"
        )
      ),
      plotlyOutput(outputId = "violin_plot"),
      make_qc_slider(x = pbmc_small$nCount_RNA, col = "nCount_RNA"),
      splitLayout(
        make_qc_slider(x = pbmc_small$nCount_RNA, col = "nCount_RNA"),
        make_qc_slider(x = pbmc_small$nFeature_RNA, col = "nFeature_RNA")
      ),
      markdown(
        mds = c(
          "### Choose filtering parameters"
        )
      ),
      markdown(
        mds = c(
          "Showing the metadata for the dataset, in order to help choose which columns to filter. Some suggestions have been selected in the checkboxes below."
        )
      ),
      DTOutput(outputId = "metadata", width = "100%")
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


server <- function(input, output) {
  # shinyFilesButton("files", label = "File select", title = "Please select a file", multiple = FALSE)

  # pbmc_small

  placeholder_dir <- file.path("/Users/johndoe/Downloads/pbmc_small")
  output$selected_directory <- renderText(placeholder_dir)
  output$nav_now <- renderText(input$nav)

  fig <- make_qc_plots(pbmc_small, "nCount_RNA")
  output$violin_plot <- renderPlotly(fig)
  # output$colnames_output <- renderText(colnames(pbmc_small[[]]))
  output$metadata <- renderDataTable(pbmc_small[[]] %>% head(10))
}

shinyApp(ui, server)
