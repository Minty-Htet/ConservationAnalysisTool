#' Run the Shiny App
#' @export
run_app <- function(browser = FALSE) {
  shiny::runApp(
    list(ui = app_ui(), server = app_server),
    launch.browser = browser
  )
}