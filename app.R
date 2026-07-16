# app.R
library(shiny)
library(pkgload)
pkgload::load_all(".")
shiny::shinyApp(ui = app_ui_oneway(), server = app_server_oneway)