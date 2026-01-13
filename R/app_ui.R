#' @keywords internal

app_ui <- function() {
  fluidPage(
    # 1. Add busy spinner and global styles
    add_busy_spinner(spin = "fading-circle"),
    
    tags$head(
      tags$style(HTML("
        .shiny-notification {
          font-size: 24px !important;
          position: fixed;
          top: 50%;
          left: 50%;
          transform: translate(-50%, -50%);
        }
        .well { 
          width: 200px !important; 
          min-width: 200px !important;
          max-width: 220px !important; 
        }
        #go_button {
          font-size: 14px !important;
          padding: 4px 12px !important;
          height: 32px !important;
          min-width: 70px !important;
          width: auto !important;
          display: inline-block !important;
        }
      "))
    ),
    
    titlePanel("Ortholog Dashboard"),
    
    sidebarLayout(
      sidebarPanel(
        width = 2,
        shinyjs::useShinyjs(),
        
        # Inputs
        selectInput('species', 'Species', c('Mmusculus', 'Hsapiens'), width = "100%"),
        
        textInput('list_in', 'Input list of genes (format: <gene name>) separated by commas or spaces.', width = "100%"),
        textInput('list_in2', 'Input list of genes (format: <gene name>) separated by commas or spaces.', width = "100%"),
        
        conditionalPanel(
          condition = "input.species === 'Mmusculus'",
          selectInput('mmver', 'Mouse Genome Version for Peaks (to mm39)',
                      c('mm10', 'mm39'), width = "100%"),
          # Unique ID for Mouse Chromosomes
          selectInput('chromosome_mm', 'Chromosome (Mouse)', 
                      choices = c("All", chr_order), multiple = TRUE, selected = "All", width = "100%")
        ),
        
        conditionalPanel(
          condition = "input.species === 'Hsapiens'",
          selectInput('hsver', 'Human Genome Version for Peaks (to mm39)', 
                      c('hg38'), width = "100%"),
          # Unique ID for Human Chromosomes
          selectInput('chromosome_hs', 'Chromosome (Human)', 
                      choices = c("All", hs_chr_order), multiple = TRUE, selected = "All", width = "100%")
        ),
        
        textInput('list_in_peak', 'Input list of peaks (chr.start.end)', width = "100%"),
        textInput('list_in_peak2', 'Input list of peaks (chr.start.end)', width = "100%"),
        
        div(style = "margin-top: 10px;", actionButton("go_button", "GO"))
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Risky Regions",
            #plotOutput('risk_plot'),
            DT::dataTableOutput('risk_render1'),
            DT::dataTableOutput('risk_render2')
          )
        ) # end tabsetPanel
      ) # end mainPanel
    ) # end sidebarLayout
  ) # end fluidPage
}