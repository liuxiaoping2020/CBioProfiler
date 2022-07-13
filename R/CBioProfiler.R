#' CBioProfiler: a web and standalone pipeine for cancer biomarker and subtype characteration from molecular level to clinical settings
#'
#'
#'
#'
#' @param run
#'
#' @export
#'
CBioProfiler <- function(run = TRUE) {
  if(run) suppressMessages(shinyParallel::runApp(system.file("app", package = "CBioProfiler"),launch.browser=TRUE))
}
