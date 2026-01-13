#' Run the Shiny App
#' @export
app <- function(browser = TRUE) {
  shiny::runApp(
    list(ui = app_ui(), server = app_server),
    launch.browser = launch.browser
  )
}