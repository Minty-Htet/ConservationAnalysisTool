# app.R
# Manually load libraries so the server doesn't have to "think" as hard
library(shiny)
pkgload::load_all(".")

shinyApp(ui = app_ui(), server = app_server)