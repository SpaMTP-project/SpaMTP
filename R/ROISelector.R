#' Launch an Interactive ROI Annotation App for SpatialExperiment
#'
#' This function opens a Shiny app to allow users to manually annotate Regions of Interest (ROIs)
#' on spatial transcriptomics or metabolomics data stored in a SpatialExperiment. Users can interactively
#' select regions using lasso selection, assign custom names to each ROI, and save the results as
#' binary (0/1) columns in colData.
#'
#' @param object A SpatialExperiment with spatial coordinates and metadata.
#' @param sampleId The sample to select regions from; required for multiple samples.
#' @param launch Launch the app; FALSE returns a Shiny app for embedding or testing.
#'
#' @return A SpatialExperiment with ROI columns (1 = selected, 0 = unselected
#'   within sample, NA = other samples), or a Shiny app if launch is FALSE.
#'
#' @details
#' The app includes options to:
#' \itemize{
#'   \item Choose a metadata column to display.
#'   \item Adjust the spot size for display.
#'   \item Use lasso to select spots and name each ROI.
#'   \item Save each ROI to colData without overwriting existing columns.
#'   \item Return the final SpatialExperiment upon clicking "Finish".
#' }
#'
#' Selection uses the exact points returned by plotly's lasso tool.
#'
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny selectInput radioButtons sliderInput textInput actionButton verbatimTextOutput
#' @importFrom shiny plotOutput renderPlot reactiveVal runApp showNotification stopApp
#' @importFrom plotly plotlyOutput renderPlotly plot_ly add_trace layout event_data
#' @importFrom sf st_polygon st_sfc st_make_valid st_as_sf st_within
#' @importFrom dplyr %>%
#' @examples
#' utils::str(formals(selectROIs))
#' @export
selectROIs <- function(object, sampleId = NULL, launch = TRUE) {
  object <- .nativeSpatialObject(object)
  coords <- .nativeCoordinates(object)
  if (is.null(sampleId)) {
    if (length(unique(coords$sample_id)) != 1L)
      stop("Choose sampleId for ROI selection on multiple samples.", call. = FALSE)
    sampleId <- unique(coords$sample_id)
  }
  if (length(sampleId) != 1L || !sampleId %in% coords$sample_id)
    stop("sampleId must identify one sample in colData.", call. = FALSE)
  selected <- which(coords$sample_id == sampleId)
  coords <- coords[selected, , drop = FALSE]
  meta_cols <- colnames(.cellMetadata(object))

  app <- shinyApp(
    ui = fluidPage(
      titlePanel("Select Regions of Interest"),
      sidebarLayout(
        sidebarPanel(
          selectInput("meta_col", "Select Metadata Column to Plot", choices = meta_cols),
          radioButtons("plot_type", "Plot Type",
                       choices = c("Feature (Continuous)" = "feature", "Categorical" = "categorical")),
          sliderInput("pt_size", "Spot Size", min = 1, max = 20, value = 5),
          textInput("roi_name", "Name for ROI", value = "ROI_1"),
          actionButton("save_roi", "Save ROI"),
          actionButton("reset_roi", "Reset ROI Selection"),
          actionButton("done", "Finish & Return Object"),
          verbatimTextOutput("status")
        ),
        mainPanel(
          plotlyOutput("spatial_plot", height = "600px")
        )
      )
    ),
    server = function(input, output, session) {
      rv <- reactiveVal(object)
      roi_mask <- reactiveVal(rep(0, nrow(coords)))

      output$spatial_plot <- renderPlotly({
        meta_col <- input$meta_col
        plot_type <- input$plot_type
        meta_data <- .cellMetadata(rv())[[meta_col]][selected]
        spot_size <- input$pt_size

        plot_ly() %>%
          add_trace(
            x = coords$x,
            y = coords$y,
            type = "scattergl",
            mode = "markers",
            marker = if (plot_type == "feature") {
              list(color = meta_data, colorscale = "Viridis", showscale = TRUE, size = spot_size)
            } else {
              list(color = as.factor(meta_data), showscale = FALSE, size = spot_size)
            },
            text = rownames(coords),
            hoverinfo = "text"
          ) %>%
          layout(
            title = paste("Spatial Plot:", meta_col),
            xaxis = list(title = "X Coordinate"),
            yaxis = list(title = "Y Coordinate"),
            dragmode = "lasso"
          )
      })

      observeEvent(event_data("plotly_selected"), {
        sel_data <- event_data("plotly_selected")
        if (!is.null(sel_data)) {
          sel_points <- sel_data$pointNumber + 1
          sel_points <- sel_points[sel_points >= 1L & sel_points <= nrow(coords)]
          if (length(sel_points)) {
            new_mask <- roi_mask()
            new_mask[sel_points] <- 1
            roi_mask(new_mask)
            output$status <- renderText("Points selected for ROI.")
          } else {
            showNotification("Select at least one point", type = "warning")
          }
        }
      })

      observeEvent(input$save_roi, {
        name <- input$roi_name
        if (name == "") {
          showNotification("Please enter a name for the ROI.", type = "error")
          return()
        }

        updated <- rv()
        if (name %in% colnames(SummarizedExperiment::colData(updated))) {
          showNotification("Choose a new column name; existing metadata is preserved.",
                           type = "error")
          return()
        }
        value <- rep(NA_integer_, ncol(updated))
        value[selected] <- roi_mask()
        SummarizedExperiment::colData(updated)[[name]] <- value
        rv(updated)
        roi_mask(rep(0, nrow(coords)))
        output$status <- renderText(paste0("Saved ROI to metadata column: ", name))
      })

      observeEvent(input$reset_roi, {
        roi_mask(rep(0, nrow(coords)))
        output$status <- renderText("ROI selection reset.")
      })

      observeEvent(input$done, {
        stopApp(rv())
      })
    }
  )
  if (launch) runApp(app) else app
}
