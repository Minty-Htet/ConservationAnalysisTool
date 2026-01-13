#' Run the Shiny App
#' @export
app <- function(browser = TRUE) {
  shiny::shinyApp(
    ui = app_ui(),
    server = app_server,
    launch.browser = browser
  )
}