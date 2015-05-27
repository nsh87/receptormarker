#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
# allow users to set viewer.suppress to FALSE to see the thing in RStudio
radial_phylo <- function(dataObject, width = NULL, height = NULL) {

  # forward options using x
  x = list(
    dataObject = dataObject
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
