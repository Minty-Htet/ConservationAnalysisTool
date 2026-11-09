#' @importFrom shinyjs useShinyjs
#' @importFrom DT DTOutput

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
    
    titlePanel("Conservation Analysis Tool"),
    
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
                  'Maximum distance between the orthologous regions of a group you tolerate before considering it to be split. Default is same value is Expected Size of the Regions of Coordinated Regulation with the floor being 500kb:', 
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
                  inputId = "dist_to_group",
                  label = "Maximum distance from peak to edge of group:",
                  min = 0, max = 100000, value = 50000, step = 10000
                ),
                sliderInput(
                  inputId = "threshold",
                  label = "Peak Conversion. Min ratio of bases that must remap:",
                  min = 0, max = 1, value = 0.5
                )
              )
            ),
            fluidRow(
              actionButton("go_button", "GO")
            )
          ),
          tabPanel(
            "Background Concepts",
            fluidRow(
              column(
                width = 12,
                markdown(
                  "
                  ### Please feel free to read our paper and check the [GitHub](https://github.com/Minty-Htet/ConservationAnalysisTool).
                  "
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                markdown(
                  "
                  # Settings Explanations
                  "
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('edgesandvertices')
              ),
              column(
                width = 6,
                markdown(
                  "
                  # Expected size of Regions of Coordinated Regulation
                  Groupings are done based on a undirected graph. 
                  Each gene is ordered based on their genomic location and they will be represented as vertices (Vi) in the graph.
                  The distance from each gene to the next gene is represented as an edge (Ei,i+1).
                  All edges longer than the Expected size of Regions of Coordinated Regulation are cut.
                  Vertices connected by edges are grouped together.
                  "
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('peaktogroup')
              ),
              column(
                width = 6,
                markdown(
                  "
                  # Maximum distance from peak to edge of group
                  This is the cut-off distance from the boundaries of a group before a peak is no longer assigned to a group.
                  Here is how it is used: <br>
                  Suppose the red dots represent the genes within a group and the lines represent peaks.
                  The bounds (Bi,1 , Bi,2) of a group (Gi) is determined by the most 5' start (Si) and the most 3' end (Ei)
                  and the cut-off distance (d). <br>
                  5' bound Bi,1 is the most 5' start minus cut-off distance (Bi,1 = Si - d) <br>
                  3' bound Bi,2 is the most 3' end minus cut-off distance (Bi,2 = Ei - d) <br>
                  Peaks (Pi) are determined to be within the bounds by checking if there is any overlap between the coordinates of
                  bounds and the coordinates of the peaks.
                  "
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                plotOutput('remapbases')
              ),
              column(
                width = 6,
                markdown(
                  "
                  # Minimum ratio of bases that must remap
                  This is the ratio (R) of the orthologous region length (Oi,j) to a given reference peak length (Li). 
                  The length of the region and the peak is measured in bases. 
                  Passing this threshold will designate the orthologous region as the orthologous peak of the reference peak. <br>
                  Note that the representative orthologous region is chosen by combining orthologous fragments 
                  (Fi,j,k) which are: <br> 
                  1) of the same orthologous chromosome and within the reference peak length (Li) of each other. <br>
                  2) taking the 5' start and the 3' end of the largest combined fragment <br>
                  The resulting orthologous region is considered an orthologous peak if the ratio of the ortholgous region length 
                  to the reference peak length is at least the chosen ratio. (Oi,j * R >= Li) <br>
                  The example to the left has Reference peak (Pi) queried at 0.5 R. Let's say the query gives us 3 fragments.
                  2 fragments map to orthoglous chromsome a (CHRa) and one to orthologous chromosome b (CHRb). 
                  The 2 fragments on CHRa (Fi,1,a ; Fi,2,a) are combined into orthologous region 1 (Oi,1) 
                  which covers 50% of the reference peak length. The 1 fragment on orthologous chromosome b (Fi,3,b) 
                  covers 40% of the reference peak length and is designated orthologous region 2 (Oi,2). 
                  The longest of the orthoglous region is chosen to be the representative region which is ortholgous 
                  region 1 (Oi,1 > Oi,2) and since orthologous region 1 meets the 0.5 threshold, it is designated as the 
                  orthologous peak.
                  "
                )
              )
            ),
            fluidRow(
              column(
                width = 6,
                markdown("# Plot Explanations")
              ),
              column(
                width = 6
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
                  The squares on the diagonal axis with bold colors represent genes or peaks with an ortholog in the comparator genome, while those with faded colors indicate there is no one-to-one ortholog. <br>
                  The faint color larger size squares encompassing more than one square on the diagonal axis denote the genes and peaks considered to be grouped. <br>
                  The example shows a group in blue made up of 2 genes and 2 peaks, with one gene and one peak having orthologs in the comparator genome while the second gene and peak does not have an ortholog. Also shown is a second group made up of a single gene by itself, with this gene having an ortholog. <br>
                  The meaning for the colors of the squares is denoted in the legends for each plot. <br>
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
                ## Choices are filled in by the server from the Group Status
                ## values currently present in the table and plot.
                selectInput(
                  'status_filter',
                  'Group Status to show in the plot and table: (You can unselect by clicking on the choice and hitting delete. Reselect by clicking on the menu below)',
                  choices = NULL, multiple = TRUE, width = "100%"
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
                  downloadButton("downloadTable1", "Download CSV"),
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
                  plotOutput('distance_plot', click = clickOpts(id = 'heatmap_click'))
                ),
                tags$div(
                  h4('Grouping Status'),
                  plotOutput('grouping_plot', click = clickOpts(id = 'heatmap_click'))
                )
              ),
              column(
                width = 4,
                tags$div(
                  h4('Sequence Conservation Status'),
                  plotOutput('sequence_plot', click = clickOpts(id = 'heatmap_click'))
                ),
                tags$div(
                  h4('Location Conservation Status'),
                  plotOutput('location_plot', click = clickOpts(id = 'heatmap_click'))
                )
              ),
              column(
                width = 4,
                tags$div(
                  h4("Click on the heatmaps to see group IDs. See Table at bottom of page for more info."),
                  DT::DTOutput('hovered_groups_table')
                ),
                sliderInput(
                  inputId = "pos_range",
                  label   = "Position range:",
                  min     = 1,
                  max     = 100,
                  value   = c(1, 20),   # default range
                  step    = 1
                )
              )
            ),
            fluidRow(
              column(
                width = 12,
                h4("Genomic Locations Table"),
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
                plotOutput('perc_attached_peaks1')
              ),
              column(
                width = 6,
                plotOutput('group_split1'),
                downloadButton("Downloadable_peaktable1", "Download CSV"), 
                DT::DTOutput("attached_peaks1")
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
                downloadButton("Downloadable_orthotable1", "Download CSV"),
                DT::DTOutput('orthologous_genepeak_table1')
              ),
              column(
                width = 6,
                plotOutput('orthologous_genepeak_plot1')
              )
            )
          ),
          tabPanel(
            "Region Size Summaries",
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
                plotOutput('singlet_summaryplot')
              ),
              column(
                width = 6,
                plotOutput('groupsize_plot')
                
              )
            )
          )
          )
        )
      )
    )
}