
#' SinCE
#'
#' A function to lauch shinyApp SinCE
#' @param
#' @keyword mean
#' @export
#' @examples
#' SinCe()

SinCE <-function() {
  library(Seurat)
  shiny::runApp(file.path(system.file("SinCEapp", package = "SinCE")))
}
