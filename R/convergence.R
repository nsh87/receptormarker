#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
convergence <- function(isLabel=TRUE,browser = FALSE) {

  # forward options using x
  x = list(
    isLabel = isLabel
  )
  
#  convergencexml <- htmltools::htmlDependency(
#    name = 'convergencexml',
#    version = 1.0,
#    src = "htmlwidgets/lib/convergencexml",
#    attachment = list(xml="basename(xml_file)")
#  )

  # create widget
  htmlwidgets::createWidget(
    name = 'convergence',
    x,
    htmlwidgets::sizingPolicy(
      padding = 22,
      viewer.suppress = browser,
      browser.fill = TRUE
    ),
    width = NULL,
    height = NULL,
    package = 'receptormarker'
  )
}

#' Widget output function for use in Shiny
#'
#' @export
convergenceOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'convergence', width, height, package = 'receptormarker')
}

#' Widget render function for use in Shiny
#'
#' @export
renderConvergence <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, convergenceOutput, env, quoted = TRUE)
}
