#' Run the Shiny App
#' @export
run_app <- function(browser = FALSE, version = 'Two') {
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