.alignmentMatrix <- function(angle, scaleX, scaleY, shiftX, shiftY, center) {
  values <- c(angle, scaleX, scaleY, shiftX, shiftY, center)
  if (length(center) != 2L || any(!is.finite(values)) || scaleX == 0 || scaleY == 0)
    stop("Alignment parameters must be finite with non-zero scales.", call. = FALSE)
  angle <- angle * pi / 180
  linear <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2L) %*%
    diag(c(scaleX, scaleY))
  output <- diag(3L)
  output[1:2, 1:2] <- linear
  output[1:2, 3] <- center + c(shiftX, shiftY) - as.numeric(linear %*% center)
  output
}

.singleAlignmentSample <- function(object) {
  samples <- unique(.nativeCoordinates(object)$sample_id)
  if (length(samples) != 1L)
    stop("Align one tissue sample at a time by subsetting the object.", call. = FALSE)
  samples[[1L]]
}

.alignmentApp <- function(source, target, imageId = NULL, multiplier = 1,
                           showTarget = TRUE) {
  sourceFrame <- .nativeCoordinates(source)
  targetFrame <- .nativeCoordinates(target)
  sampleId <- .singleAlignmentSample(target)
  .singleAlignmentSample(source)
  image <- if (is.null(imageId)) NULL else .nativeImage(target, imageId, sampleId)
  initial <- diag(c(multiplier, multiplier, 1))
  center <- colMeans(sourceFrame[, c("x", "y")]) * multiplier
  extent <- max(diff(range(c(sourceFrame$x * multiplier, targetFrame$x))),
                diff(range(c(sourceFrame$y * multiplier, targetFrame$y))), 1)
  ui <- shiny::fluidPage(
    shiny::titlePanel("Align spatial measurements"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::sliderInput("angle", "Rotation (degrees)", -180, 180, 0),
        shiny::sliderInput("shiftX", "Horizontal shift", -extent, extent, 0),
        shiny::sliderInput("shiftY", "Vertical shift", -extent, extent, 0),
        shiny::sliderInput("scaleX", "Horizontal scale", 0.1, 4, 1, step = 0.01),
        shiny::sliderInput("scaleY", "Vertical scale", 0.1, 4, 1, step = 0.01),
        shiny::checkboxInput("flipX", "Flip horizontally", FALSE),
        shiny::checkboxInput("flipY", "Flip vertically", FALSE),
        shiny::actionButton("apply", "Return aligned data"),
        shiny::actionButton("cancel", "Cancel")),
      shiny::mainPanel(shiny::plotOutput("preview", height = "700px"))))
  server <- function(input, output, session) {
    transform <- shiny::reactive({
      .alignmentMatrix(input$angle,
        input$scaleX * if (isTRUE(input$flipX)) -1 else 1,
        input$scaleY * if (isTRUE(input$flipY)) -1 else 1,
        input$shiftX, input$shiftY, center) %*% initial
    })
    output$preview <- shiny::renderPlot({
      xy <- cbind(as.matrix(sourceFrame[, c("x", "y")]), 1) %*% t(transform())
      frame <- data.frame(x = xy[, 1], y = xy[, 2])
      plot <- ggplot2::ggplot()
      if (!is.null(image)) plot <- plot + .imageLayer(image)
      if (showTarget) plot <- plot + ggplot2::geom_point(
        data = targetFrame, ggplot2::aes(x = x, y = y), colour = "steelblue", alpha = 0.5)
      plot + ggplot2::geom_point(data = frame, ggplot2::aes(x = x, y = y),
        colour = "tomato", alpha = 0.7) + ggplot2::coord_fixed() +
        ggplot2::scale_y_reverse() + ggplot2::theme_void()
    })
    shiny::observeEvent(input$apply, shiny::stopApp(transform()))
    shiny::observeEvent(input$cancel, shiny::stopApp(NULL))
  }
  shiny::shinyApp(ui, server)
}

#' Align two spatial experiments using an affine transformation
#'
#' A supplied homogeneous matrix is applied reproducibly without a GUI. With
#' transformation = NULL, an interactive preview allows translation, rotation,
#' scaling and reflection. Coordinates are always stored in their full spatial
#' units; an imgData scaleFactor affects only display of the reference image.
#' @param sm.data,st.data Source and reference SpatialExperiment objects, each
#'   containing one tissue sample.
#' @param transformation A 3 by 3 homogeneous matrix mapping source coordinates
#'   to reference coordinates, or NULL for interactive alignment.
#' @param image.slice Optional reference image ID in imgData(st.data).
#' @param msi.pixel.multiplier Explicit initial scaling of source coordinates.
#'   Defaults to one. Applied before transformation, including GUI transforms.
#' @param shiny.host,shiny.port Host and optional port for the interactive app.
#' @param verbose Show alignment messages.
#' @return The source SpatialExperiment with updated spatialCoords and
#'   alignment provenance recorded by applySpatialAlignment.
#' @export
alignSpatialOmics <- function(
    sm.data, st.data, transformation = NULL, image.slice = NULL,
    msi.pixel.multiplier = 1, shiny.host = "127.0.0.1", shiny.port = NULL,
    verbose = FALSE
) {
  sm.data <- .nativeSpatialObject(sm.data)
  st.data <- .nativeSpatialObject(st.data)
  .singleAlignmentSample(sm.data)
  .singleAlignmentSample(st.data)
  if (length(msi.pixel.multiplier) != 1L || !is.finite(msi.pixel.multiplier) ||
      msi.pixel.multiplier <= 0)
    stop("msi.pixel.multiplier must be one positive number.", call. = FALSE)
  if (is.null(transformation)) {
    transformation <- shiny::runApp(
      .alignmentApp(sm.data, st.data, image.slice, msi.pixel.multiplier),
      host = shiny.host, port = shiny.port)
    if (is.null(transformation)) stop("Alignment cancelled.", call. = FALSE)
  } else {
    if (!is.matrix(transformation) || !identical(dim(transformation), c(3L, 3L)) ||
        any(!is.finite(transformation)))
      stop("transformation must be a finite 3 by 3 matrix.", call. = FALSE)
    transformation <- transformation %*% diag(c(msi.pixel.multiplier,
                                                msi.pixel.multiplier, 1))
  }
  applySpatialAlignment(sm.data, ST.data = st.data,
    alignment = list(transformation = list(matrix = transformation)),
    verbose = verbose)
}

#' Inspect the alignment of two spatial experiments
#' @param SM.data,ST.data SpatialExperiment objects in the same coordinate units.
#' @param names Two labels for the source and target measurements.
#' @param cols Two point colours.
#' @param image.slice Optional image ID in ST.data. With multiple samples,
#'   subset to one sample before drawing an image overlay.
#' @param size Point size.
#' @return A ggplot; multiple samples are shown in separate panels.
#' @export
checkAlignment <- function(SM.data, ST.data, names = c("SM", "ST"), cols = NULL,
                            image.slice = NULL, size = 0.5) {
  source <- .nativeCoordinates(SM.data)
  target <- .nativeCoordinates(ST.data)
  if (length(names) != 2L || (!is.null(cols) && length(cols) != 2L))
    stop("Provide two dataset labels and two colours.", call. = FALSE)
  source$dataset <- names[[1L]]
  target$dataset <- names[[2L]]
  frame <- rbind(source, target)
  plot <- ggplot2::ggplot(frame, ggplot2::aes(x = x, y = y, colour = dataset))
  if (!is.null(image.slice)) {
    .singleAlignmentSample(SM.data)
    sampleId <- .singleAlignmentSample(ST.data)
    plot <- plot + .imageLayer(.nativeImage(ST.data, image.slice, sampleId))
  }
  plot <- plot + ggplot2::geom_point(size = size) +
    ggplot2::scale_colour_manual(values = cols %||% c("#F8766D", "#00BFC4")) +
    ggplot2::scale_y_reverse() + ggplot2::coord_fixed() + ggplot2::theme_void()
  if (is.null(image.slice) && length(unique(frame$sample_id)) > 1L)
    plot <- plot + ggplot2::facet_wrap(~sample_id)
  plot
}

#' Attach and optionally align an optical image to MSI
#' @param image_path Local image path.
#' @param SpaMTP A SpatialExperiment.
#' @param imageId Identifier for the attached image.
#' @param sampleId Sample identifier; required for multi-sample objects.
#' @param scaleFactor Image pixels per spatial coordinate unit.
#' @param transformation Optional 3 by 3 matrix mapping MSI to image coordinates
#'   in full spatial units.
#' @param interactive Open the manual alignment app when no matrix is supplied.
#' @param ... Additional arguments to shiny::runApp when interactive is TRUE.
#' @return A SpatialExperiment with an imgData entry and optionally transformed
#'   spatialCoords. Coordinate alignment requires one tissue sample.
#' @export
addSMImage <- function(image_path, SpaMTP, imageId = "optical", sampleId = NULL,
                       scaleFactor = 1, transformation = NULL,
                       interactive = TRUE, ...) {
  object <- .nativeSpatialObject(SpaMTP)
  sampleId <- sampleId %||% .singleAlignmentSample(object)
  reference <- addSpatialImage(object, imageSource = image_path,
    sampleId = sampleId, imageId = imageId, scaleFactor = scaleFactor)
  if (interactive && is.null(transformation)) {
    transformation <- shiny::runApp(
      .alignmentApp(object, reference, imageId, showTarget = FALSE), ...)
    if (is.null(transformation)) stop("Alignment cancelled.", call. = FALSE)
  }
  if (!is.null(transformation)) {
    object <- alignSpatialOmics(object, reference, transformation = transformation)
    SpatialExperiment::imgData(object) <- SpatialExperiment::imgData(reference)
    return(object)
  }
  reference
}
