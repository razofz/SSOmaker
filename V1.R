# Load the required libraries
library(shiny)
library(shinyFiles)
library(Seurat)

# Define the UI
ui <- fluidPage(
  titlePanel("Folder Selection App"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("folder", "Choose a folder", "Please select a folder", title = "Select a folder")
    ),
    mainPanel(
      h4("Selected Folder Path:"),
      verbatimTextOutput("folder_path")
    )
  )
)

# Define the server
server <- function(input, output, session) {
  shinyDirChoose(input, "folder", roots = c(home = "~"))

  observe({
    if (!is.null(input$folder)) {
      chosen_folder_path <- parseDirPath(roots = c(home = "~"), input$folder)
      output$folder_path <- renderText({
        chosen_folder_path
        print(chosen_folder_path)
        so.data <- Read10X(data.dir = chosen_folder_path)
      })
    }
  })
}

# Run the Shiny app
shinyApp(ui, server)
