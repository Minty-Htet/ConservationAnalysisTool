#' Run the Shiny App
#' @importFrom shiny shinyApp
#' @importFrom shiny runApp
#' @importFrom pkgload load_all
#' @param browser TRUE or FALSE. Wherether you want to open browser or not.
#' @param version Default is 'One'. The two way is 'Two'
#' @return The shiny app
#' @export
run_app <- function(browser = FALSE, version = 'One') {
  if (version == 'Two') {
    shiny::runApp(
    list(ui = app_ui(), server = app_server),
    launch.browser = browser
    )
  } else {
    shiny::runApp(
    list(ui = app_ui_oneway(), server = app_server_oneway),
    launch.browser = browser
    )
  }
  
}