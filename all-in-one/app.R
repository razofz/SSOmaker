library(shiny)
library(shinyFiles)
library(plotly)
library(gridlayout)
library(bslib)


ui <- grid_page(
  layout = c(
    "logo    header         ",
    "meatandpotatoes meatandpotatoes"
  ),
  row_sizes = c(
    "100px",
    "1fr"
  ),
  col_sizes = c(
    "250px",
    "1fr"
  ),
  gap_size = "1rem",
  grid_card(area = "sidebar", card_header("Overview")),
  grid_card_text(
    area = "header",
    content = "Seurat Object Maker (all-in-one)",
    alignment = "center",
    is_title = TRUE
  ),
  grid_card(
    area = "meatandpotatoes",
    card_header(strong("Content section")),
    card_body(
      tabPanel(
        title = "My Shiny App",
        tabsetPanel(
          tabPanel(
            title = "Load data",
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
            )
          ),
          tabPanel(
            title = "Filtering",
            card_body(
              markdown(
                mds = c(
                  "**Selected directory:**"
                )
              ),
              textOutput(outputId = "selected_directory"),
              markdown(
                mds = c(
                  "### Violin plots of basic characteristics"
                )
              ),
              plotlyOutput(outputId = "violin_plot"),
              markdown(
                mds = c(
                  "### Choose filtering parameters"
                )
              ),
              sliderInput(
                inputId = "rna_count",
                label = "nCount_RNA",
                min = 0,
                max = 1000,
                value = 5,
                width = "75%"
              ),
              numericInput(
                inputId = "num_input_rna_count",
                label = "nCount_RNA",
                value = 500,
                min = 0,
                max = 1000
              )
            ),
            card_body()
          ),
          tabPanel(title = "Results")
        )
      )
    )
  ),
  grid_card(
    area = "logo",
    card_body(
      markdown(
        mds = c(
          "<img src='https://www.staff.lu.se/sites/staff.lu.se/files/styles/lu_wysiwyg_full_tablet/public/2021-04/Lunduniversity-horisontal.png.webp?itok=_rp_OxRe' width='100%' />"
        )
      )
    )
  )
)


server <- function(input, output) {
  # shinyFilesButton("files", label = "File select", title = "Please select a file", multiple = FALSE)
}

shinyApp(ui, server)
