qc_slider_ui <- function(id) {
  shiny::tagList(
    shiny::sliderInput(
      inputId = shiny::NS(id, "slider"),
      label = "slider",
      min = 0,
      max = 5000,
      value = c(0, 5000),
      width = "75%"
    ),
  )
}
qc_slider_server <- function(id, col, metadata, start_values) {
  stopifnot(shiny::is.reactivevalues(metadata))
  stopifnot(!shiny::is.reactive(col))
  stopifnot(!shiny::is.reactive(start_values))

  shiny::moduleServer(id, function(input, output, session) {
    update_it <- shiny::reactiveVal(0)
    shiny::observe({
      if (shiny::isolate(update_it()) < 3) {
        shiny::isolate(update_it(update_it() + 1))
      }
      if (update_it() < 2) {
        # print(stringr::str_c("update_it(): ", update_it()))
        shiny::invalidateLater(700, session)
      }
      col_range <- range(shiny::isolate(metadata$data[[col]]), na.rm = TRUE)
      shiny::updateSliderInput(
        session,
        inputId = "slider",
        label = col,
        min = col_range[1],
        max = col_range[2],
        value = start_values
      )
    })
    return(shiny::reactive(input$slider))
  })
}

qc_plot_ui <- function(id) {
  shiny::uiOutput(shiny::NS(id, "output"))
}
qc_plot_server <- function(id, col, metadata, ranges) {
  stopifnot(shiny::is.reactivevalues(metadata))
  stopifnot(!shiny::is.reactive(col))
  stopifnot(shiny::is.reactive(ranges))

  shiny::moduleServer(id, function(input, output, session) {
    output$output <- shiny::renderUI({
      htmltools::tagList(
        plotly::renderPlotly({
          p <- plotly::plot_ly(
            data = metadata$data,
            y = as.formula(stringr::str_c(" ~ ", col)),
            type = "violin",
            box = list(visible = T),
            meanline = list(visible = T),
            name = col,
            x0 = col
          )
          return(
            plotly::layout(
              p = p,
              yaxis = list(zeroline = F),
              shapes = list(
                hline(ranges()[1]),
                hline(ranges()[2])
              )
            )
          )
        })
      )
    })
  })
}

qc_module_UI <- function(id) {
  bslib::layout_column_wrap(
    width = "600px",
    htmltools::div(
      style = "display: flex; justify-content: center;",
      htmltools::div(
        htmltools::div(
          style = "display: flex; justify-content: center;",
          bslib::tooltip(
            htmltools::span(
              shiny::helpText("Need help? "),
              bsicons::bs_icon("info-circle")
            ),
            stringr::str_c(
              "Hover over the plot to see extra statistics about the metadata column. ",
              "You can also drag to zoom in on a specific area of the plot."
            )
          )
        ),
        qc_plot_ui(shiny::NS(id, "qc_plot")),
        # verbatimTextOutput(NS(id, "foo")),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          qc_slider_ui(shiny::NS(id, "qc_slide"))
        ),
        htmltools::div(
          style = "display: flex; justify-content: center;",
          htmltools::span(
            shiny::helpText("Adjust the sliders to set the filtering cutoffs."),
            bslib::tooltip(
              bsicons::bs_icon("info-circle"),
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
qc_module_server <- function(id, col, metadata, start_values) {
  stopifnot(!shiny::is.reactive(start_values))
  stopifnot(shiny::is.reactivevalues(metadata))
  stopifnot(!shiny::is.reactive(col))

  shiny::moduleServer(id, function(input, output, session) {
    slider_input_vals <- qc_slider_server(
      "qc_slide",
      col = col,
      metadata = metadata,
      start_values = start_values
    )
    qc_plot_server(
      "qc_plot",
      col = col,
      metadata = metadata,
      ranges = slider_input_vals
    )
    output$foo <- shiny::renderText({
      slider_input_vals()
    })

    return(slider_input_vals)
  })
}

violin_plot_ui <- function(id) {
  shiny::uiOutput(shiny::NS(id, "output"))
}
violin_plot_server <- function(id, col, metadata) {
  stopifnot(shiny::is.reactivevalues(metadata))
  stopifnot(!shiny::is.reactive(col))

  shiny::moduleServer(id, function(input, output, session) {
    output$output <- shiny::renderUI({
      htmltools::tagList(
        htmltools::div(
          style = "display: flex; justify-content: center;",
          bslib::tooltip(
            htmltools::span(
              shiny::helpText("Need help? "),
              bsicons::bs_icon("info-circle")
            ),
            stringr::str_c(
              "Hover over the plot to see extra statistics about the metadata column. ",
              "You can also drag to zoom in on a specific area of the plot."
            )
          )
        ),
        plotly::renderPlotly({
          p <- plotly::plot_ly(
            data = metadata$data,
            y = as.formula(stringr::str_c(" ~ ", col)),
            type = "violin",
            box = list(visible = T),
            meanline = list(visible = T),
            name = col,
            x0 = col
          )
          return(
            plotly::layout(
              p = p,
              yaxis = list(zeroline = F)
            )
          )
        })
      )
    })
  })
}
