app_ui_oneway <- function() {
  fluidPage(
    # 1. Add busy spinner and global styles
    #add_busy_spinner(spin = "fading-circle"),
    
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
      ## Side bar Panel for on the fly Chromosome Selection
      sidebarPanel(
        width = 2,
        shinyjs::useShinyjs(),
        conditionalPanel(
          condition = "input.species === 'Mmusculus'",
          # Unique ID for Mouse Chromosomes
          selectInput(
            'chromosome', 'Chromosome (Mouse)', 
            choices = c("All", chr_order), multiple = TRUE,
            selected = "All", width = "100%"
          )
        ),
        conditionalPanel(
          condition = "input.species === 'Hsapiens'",
          # Unique ID for Human Chromosomes
          selectInput(
            'chromosome', 'Chromosome (Human)', 
            choices = c("All", hs_chr_order), multiple = TRUE,
            selected = "All", width = "100%"
          )
        )
      ),
      
      mainPanel(
        tabsetPanel(
          ## Page for Inputs and Settings
          tabPanel(
            "Inputs and Settings",
            fluidRow(
              column(
                width = 6,
                selectInput(
                  'species', 'Species',
                  c('Mmusculus', 'Hsapiens'),
                  width = "100%"
                )
              ),
              column(
                width = 6,
                conditionalPanel(
                  condition = "input.species === 'Mmusculus'",
                  selectInput(
                    'mmver', 'Mouse Genome Version for Peaks',
                    c('mm10', 'mm39'), width = "100%"
                  )
                ),
                conditionalPanel(
                  condition = "input.species === 'Hsapiens'",
                  selectInput(
                    'hsver', 'Human Genome Version for Peaks', 
                    c('hg38'), width = "100%"
                  )
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                textInput(
                  'list_in',
                  'List 1. Input list of genes (format: <gene name>) separated by commas or spaces.'
                )
              ),
              column(
                width = 6,
                textInput(
                  'list_in_peak',
                  'List 1. Input list of peaks (chr.start.end)'
                )
              )
            ),
            ## Adjustable settings
            fluidRow(
              column(
                width = 6,
                sliderInput(
                  inputId = "exp_rocr",
                  label = "Expected size of Regions of Coordinated Regulation:",
                  min = 0, max = 1000000, value = 500000, step = 250000
                ),
                selectInput(
                  'split_choice', 
                  'Maximum distance between the orthologous regions of a group you tolerate before considereing it to be split. Default is same value is Expected Size of the Regions of Coordinated Regulation with the floor being 500kb:', 
                  choices = c("Default", 'Manual'), multiple = FALSE,
                  selected = "Default",
                ),
                conditionalPanel(
                  condition = "input.split_choice === 'Manual'",
                  textInput('split_dis', 'Manually Input Split Distance')
                )
              ),
              column(
                width = 6,
                sliderInput(
                  inputId = "threshold",
                  label = "Peak Conversion. Min ratio of bases that must remap:",
                  min = 0, max = 1, value = 0.5
                ),
                sliderInput(
                  inputId = "dist_to_group",
                  label = "Maximum distance from peak to edge of group:",
                  min = 0, max = 100000, value = 50000, step = 10000
                ),
              )
            ),
            fluidRow(
              actionButton("go_button", "GO")
            )
          ),
          
          tabPanel(
            "Risky Regions Plot",
            fluidRow(
              column(
                width = 12,
                height = 10,
                plotOutput('risk_plot', height = "800px")
              )
            ),
            fluidRow(
              column(
                width = 4,
                sliderInput(
                  inputId = "disttobreak",
                  label = "Distance at which regulatory groups are considered AT RISK of crossing over EVOLUTIONARY BREAKPOINT.",
                  min = 0, max = 1000000, value = 100000, step = 10000
                )
              ),
              column(
                width = 4,
                selectInput(
                  'show_choice', 
                  'Showing:', 
                  choices = c("All", 'Paralogs'), multiple = FALSE,
                  selected = "All",
                )
              ),
              column(
                width = 4,
                conditionalPanel(
                  condition = "input.show_choice === 'Paralogs'",
                  textInput('dup_group', 'Input Group to show Paralogous Genes')
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                tags$div(
                  h4("List 1"),
                  DT::DTOutput('risk_render1')
                )
              )
            )
          ),
          
          tabPanel(
            "Grouping Heatmaps",
            fluidRow(
              column(
                width = 4,
                tags$div(
                  h4("Distance Heatmap"),
                  plotOutput('distance_plot', hover = hoverOpts(id = 'heatmap_hover', delay = 100))
                ),
                tags$div(
                  h4('Grouping Status'),
                  plotOutput('grouping_plot', hover = hoverOpts(id = 'heatmap_hover', delay = 100))
                )
              ),
              column(
                width = 4,
                tags$div(
                  h4('Sequence Conservation Status'),
                  plotOutput('sequence_plot', hover = hoverOpts(id = 'heatmap_hover', delay = 100))
                ),
                tags$div(
                  h4('Location Conservation Status'),
                  plotOutput('location_plot', hover = hoverOpts(id = 'heatmap_hover', delay = 100))
                )
              ),
              column(
                width = 4,
                DT::DTOutput('position_groups_table'),
                sliderInput(
                  inputId = "pos_range",
                  label   = "Position range:",
                  min     = 1,
                  max     = 100,
                  value   = c(1, 20),   # default range
                  step    = 1
                ),
                tags$div(
                  h4("Hover over the heatmaps to see group IDs"),
                  DT::DTOutput('hovered_groups_table')
                )
              )
            )
          ),
          tabPanel(
            "Grouping Statistics",
            fluidRow(
                  selectInput(
                  'group_gene_mode', 
                  'Gene or Group statistic:', 
                  choices = c('Group', 'Gene'), multiple = FALSE,
                  selected = "Group",
                )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('group_by_break1'),
                plotOutput('group_split1')
              ),
              column(
                width = 6,
                plotOutput('group_split1')
              )
            )
          )
          )
        )
      )
    )
}