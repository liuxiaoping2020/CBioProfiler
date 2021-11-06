#' CBioExplorer: a web and standalone application for screening, validation, and annotation of cancer survival related biomarkers from molecular level to clinical settings
#'
#'
#'
#'
#' @export
#'
CBioExplorer <- function(run = TRUE) {
  if(run) suppressMessages(shiny::runApp(system.file("app", package = "CBioExplorer"),launch.browser=TRUE))
}
