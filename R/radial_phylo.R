#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
# allow users to set viewer.suppress to FALSE to see the thing in RStudio
radial_phylo <- function(dataObject, radius='fill', autoResize=FALSE, width=NULL, height=NULL) {

  # Determine what the parameter 'radius' is
  err = "The argument 'radius' is invalid"
  radius_options <- c('fill')
  tryCatch({
    if (!is.element(radius, radius_options) && (radius != floor(radius))) {
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(err, call.=FALSE)
  }
  )
  
  # Warn user if setting 'radius' and autoResize=TRUE at same time
  warn = "Using autoResize=TRUE will override the 'radius' parameter"
  tryCatch({
    if (is.numeric(radius) && autoResize == TRUE) {
      warning(warn, call.=FALSE)
    }
  })
  
  # forward options using x
  x = list(
    dataObject = dataObject,
    radius = radius,
    autoResize = autoResize
  )

  # create widget
  htmlwidgets::createWidget(
    name = 'radial_phylo',
    x,
    width = width,
    height = height,
    htmlwidgets::sizingPolicy(
      viewer.padding = 0,
      #viewer.suppress = TRUE,
      viewer.paneHeight = 500,
      browser.fill = TRUE
    ),
    package = 'receptormarker'
  )
}

#' Widget output function for use in Shiny
#'
#' @export
radial_phyloOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'radial_phylo', width, height, package = 'receptormarker')
}

#' Widget render function for use in Shiny
#'
#' @export
renderRadial_phylo <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, radial_phyloOutput, env, quoted = TRUE)
}
