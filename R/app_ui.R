#' @importFrom shinyjs useShinyjs
#' @importFrom DT DTOutput

app_ui <- function() {
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
    
    titlePanel("Conservation Analysis Tool-2 way"),
    
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
                ),
                textInput(
                  'list_in_peak',
                  'List 1. Input list of peaks (chr.start.end)'
                )
              ),
              column(
                width = 6,
                textInput(
                  'list_in2',
                  'List 2. Input list of genes (format: <gene name>) separated by commas or spaces.'
                ),
                textInput(
                  'list_in_peak2',
                  'List 2. Input list of peaks (chr.start.end)'
                )
              )
            ),
            ## Adjustable settings
            fluidRow(
              column(
                width = 6,
                sliderInput(
                  inputId = "dist_to_group",
                  label = "Maximum distance from peak to edge of group:",
                  min = 0, max = 100000, value = 50000, step = 10000
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
                  inputId = "exp_rocr",
                  label = "Expected size of Regions of Coordinated Regulation:",
                  min = 0, max = 1000000, value = 500000, step = 250000
                )
              )
            ),
            fluidRow(
              actionButton("go_button", "GO")
            )
          ),
          tabPanel(
            "Instructions",
            fluidRow(
              column(
                width = 12,
                markdown(
                  "
                  ### Please feel free to read our paper and check the [GitHub](https://github.com/Minty-Htet/ConservationAnalysisTool_Shiny).
                  "
                )
              )
            ),
            fluidRow( 
              column(
                width = 6,
                plotOutput('risk_example1')
              ),
              column( 
                width = 6,
                markdown(
                  "
                  # Risky Regions Plot
                  On the X-axis coordinates: <br>
                  0.0-0.5: The lines extending from the x-axis represent well annotated one-to-one orthologous gene pairs from mice and humans. The position of the lines on the x-axis represents the genomic coordinates of the gene in the reference genome and the height of the lines represent the genomic coordinates of the ortholog in the comparator genome. The color of the line represents the chromosome the ortholog is located on in the comparator genome, with colors for each chromosome indicated in the legend. <br>
                  1.0, 1.5, 2.0: The arrows represent genes or peaks from the input datasets that have a one-to-one ortholog in the comparator species. 
                  (1.0, 1.5) The smaller arrows represent peaks, and (2.0) the larger arrows represent genes. 
                  Peaks unattached to a gene group will be represented as a small gray arrow since these cannot have a group status. <br>
                  2.5, 3.0: The open triangles represent genes, and the open diamonds represent peaks that do not have an identifiable one-to-one ortholog in the comparator species. <br>
                  The colors for the arrows, triangles, and diamonds represent the status of the groups as indicated in the legend.
                  Not applicable indicates a single gene that does not meet the definition of a group of two or more features. <br>
                  To view individual chromosomes, use the drop-down menu on the left sidebar.
                  "
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('grouping_example1')
              ),
              column(
                width = 6,
                markdown(
                  "
                  # Grouping Heatmaps
                  The triangles on the diagonal axis with bold colors represent genes or peaks with an ortholog in the comparator genome, while those with faded colors indicate there is no one-to-one ortholog. <br>
                  The faint color larger size triangle encompassing more than one triangle on the diagonal axis denote the genes and peaks considered to be grouped. <br>
                  The example shows a group in blue made up of 2 genes and 2 peaks, with one gene and one peak having orthologs in the comparator genome while the second gene and peak does not have an ortholog. Also shown is a second group made up of a single gene by itself, with this gene having an ortholog. <br>
                  The meaning for the colors of the triangles is denoted in the legends for each plot. <br>
                  Clicking on the features in the heatmaps will bring up the gene and peak IDs for all members in the group. <br>
                  The position range bar can be used to zoom in on specific areas of the heatmaps. 
                  "
                )
              )
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
                  choices = c("All", 'Shared', 'Paralogs'), multiple = FALSE,
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
                width = 6,
                tags$div(
                  h4("List 1"),
                  downloadButton("downloadTable1", "Download CSV"),
                  DT::DTOutput('risk_render1')
                )
              ),
              column(
                width = 6,
                tags$div(
                  h4("List 2"),
                  downloadButton("downloadTable2", "Download CSV"),
                  DT::DTOutput('risk_render2')
                )
              )
            )
          ),
          
          tabPanel(
            "Grouping Heatmaps",
            fluidRow(
              column(
                width = 6,
                sliderInput(
                  inputId = "pos_range",
                  label   = "Position range:",
                  min     = 1,
                  max     = 100,
                  value   = c(1, 20),   # default range
                  step    = 1
                )
              ),
              column(
                width = 6,
                h4("Click on the heatmaps to see group IDs in the table below:"),
                DT::DTOutput('hovered_groups_table')
              )
            ),
            fluidRow(
              column(
                width = 6,
                tags$div(
                  h4("Distance Heatmap"),
                  plotOutput('distance_plot', click = clickOpts(id = 'heatmap_click'))
                ),
                tags$div(
                  h4('Grouping Status'),
                  plotOutput('grouping_plot', click = clickOpts(id = 'heatmap_click'))
                )
              ),
              column(
                width = 6,
                tags$div(
                  h4('Sequence Conservation Status'),
                  plotOutput('sequence_plot', click = clickOpts(id = 'heatmap_click'))
                ),
                tags$div(
                  h4('Location Conservation Status'),
                  plotOutput('location_plot', click = clickOpts(id = 'heatmap_click'))
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                DT::DTOutput('position_groups_table')
              )
            )
          ),
          tabPanel(
            "Grouping Summaries",
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
                plotOutput('group_split1'),
                plotOutput('perc_attached_peaks1'),
                downloadButton("Downloadable_peaktable1", "Download CSV"),
                DT::DTOutput("attached_peaks1")
              ),
              column(
                width = 6,
                plotOutput('group_by_break2'),
                plotOutput('group_split2'),
                plotOutput('perc_attached_peaks2'),
                downloadButton("Downloadable_peaktable2", "Download CSV"),
                DT::DTOutput("attached_peaks2")
              )
            )
          ),
          tabPanel(
            "Orthologous Summaries",
            fluidRow(
              width = 6,
              selectInput(
                'chosen_ortho_stat',
                'Choose to see summary for genes or peaks',
                choices = c('Genes', 'Peaks'), multiple = FALSE,
                selected = "Genes"
              )
            ),
            fluidRow(
              column(
                width = 6,
                h4("List 1"),
                downloadButton("Downloadable_orthotable1", "Download CSV"),
                DT::DTOutput('orthologous_genepeak_table1'),
                h4("List 2"),
                downloadButton("Downloadable_orthotable2", "Download CSV"),
                DT::DTOutput('orthologous_genepeak_table2')
              ),
              column(
                width = 6,
                h4("List 1"),
                plotOutput('orthologous_genepeak_plot1'),
                h4("List 2"),
                plotOutput('orthologous_genepeak_plot2')
              )
            )
          ),
          tabPanel(
            "Region Sizes Summaries",
            fluidRow(
              column(
                width = 6
              ),
              column(
                width = 6,
                selectInput(
                'chosen_region_stat',
                'Choose to see summary for groups or genes',
                choices = c('Groups', 'Genes'), multiple = FALSE,
                selected = "Groups"
              )
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('singlet_summaryplot1')
              ),
              column(
                width = 6,
                plotOutput('groupsize_plot1')
                
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('singlet_summaryplot2')
              ),
              column(
                width = 6,
                plotOutput('groupsize_plot2')
                
              )
            )
          )

          )
        )
      )
    )
}