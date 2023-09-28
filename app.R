# Load the required libraries
library(shiny)
library(shinyFiles)
library(Seurat)
library(ggplot2)

# Define the UI
ui <- fluidPage(
  titlePanel("Shamit's Seurat Object Maker"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("folder", "Choose a folder", "Please select a folder", title = "Select a folder"),
      numericInput("mtpc_max", "% mito max:", min = 0, max = 100, value = 5),
      numericInput("mtpc_min", "% mito min:", min = 0, max = 100, value = 0),
      numericInput("nfrna_max", "nFeature_RNA max:", min = 0, max = 1000000, value = 10000),
      numericInput("nfrna_min", "nFeature_RNA min:", min = 0, max = 1000000, value = 200),
      numericInput("ncrna_max", "nCount_RNA max:", min = 0, max = 1000000, value = 10000),
      numericInput("ncrna_min", "nCount_RNA min:", min = 0, max = 1000000, value = 100),
      numericInput("scale.fac", "Scale factor:", min = 0, max = 1000000, value = 10000),
      numericInput("hvgs", "Number of HVGs:", min = 0, max = 10000, value = 2000),
      numericInput("ncomp", "Number of PCAs:", min = 3, max = 100, value = 10),
      numericInput("res", "Cluster resolution:", min = 0.1, max = 3, value = 0.5),
      downloadButton("download", "Download Seurat Object")
    ),
    mainPanel(
      h4("Selected Folder Path:"),
      verbatimTextOutput("folder_path"),
      h4("Violin plots of basic characteristics:"),
      plotOutput("vlnplot"),
      h4("Basic QC:"),
      plotOutput("qc"),
      h4("Basic QC post filter:"),
      plotOutput("qc.pf"),
      h4("HVGs:"),
      plotOutput("hvg.plot"),
      h4("PCA plot:"),
      plotOutput("pca.plot"),
      h4("UMAP:"),
      plotOutput("umap")
    )
  )
)

# Define the server
server <- function(input, output, session) {
  shinyDirChoose(input, "folder", roots = c(home = "~"))

  observe({
    if (!is.null(input$folder)) {
      chosen_folder_path <- parseDirPath(roots = c(home = "~"), input$folder)
      output$folder_path <- renderPrint({
        chosen_folder_path
        print(chosen_folder_path)
        so.data <- Read10X(data.dir = chosen_folder_path)
        print("Dimensions of data uploaded:")
        print(dim(so.data))
        so <- CreateSeuratObject(counts = so.data, project = "shiny", min.cells = 3, min.features = 200)
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

        qc1 <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

        output$vlnplot <- renderPlot({
          qc1
        })

        plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        output$qc <- renderPlot({
          plot1 + plot2
        })

        so <- subset(so, subset = nFeature_RNA > input$nfrna_min & nFeature_RNA < input$nfrna_max & percent.mt < input$mtpc_max)
        print(dim(so))

        plot1qc <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2qc <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        output$qc.pf <- renderPlot({
          plot1qc + plot2qc
        })

        so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = input$scale.fac)

        so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = input$hvgs)

        top10 <- head(VariableFeatures(so), 10)

        # plot variable features with and without labels
        plot1h <- VariableFeaturePlot(so)
        plot2h <- LabelPoints(plot = plot1h, points = top10, repel = TRUE)
        output$hvg.plot <- renderPlot({
          plot1h + plot2h
        })

        all.genes <- rownames(so)
        so <- ScaleData(so, features = all.genes)

        so <- RunPCA(so, features = VariableFeatures(object = so))

        dimplot <- DimPlot(so, reduction = "pca")

        output$pca.plot <- renderPlot({
          dimplot
        })

        so <- FindNeighbors(so, dims = 1:input$ncomp)
        so <- FindClusters(so, resolution = input$res)
        so <- RunUMAP(so, dims = 1:input$ncomp)
        umap <- DimPlot(so, reduction = "umap")

        output$umap <- renderPlot(umap)

        output$download <- downloadHandler(
          filename = function() {
            paste("Seurat_Object_SSOMoutput_", Sys.Date(), ".rds", sep = "")
          },
          content = function(filename) {
            saveRDS(so, filename)
          }
        )
      })
    }
  })
}

# Run the Shiny app
shinyApp(ui, server)
