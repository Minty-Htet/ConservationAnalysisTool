#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table fifelse
#' @importFrom data.table copy
#' @importFrom stringr str_split
#' @keywords internal

app_server <- function(input, output, session) {

  # Message shown when starting the thing  
  observeEvent(input$go_button, {
    showNotification("Running. Will take a minute.", type = "message", duration = NULL)
  })

  # Choose chromosomes of interest
  chosen_chr <- reactive({
    if ("All" %in% input$chromosome) {
      if (input$species == "Mmusculus") chr_order else hs_chr_order
    } else {
      input$chromosome
    }
  })

  # Choose genes of interest.
  chosen_inputs1 <- eventReactive(input$go_button, {
  # Split your input based on commas or spaces
    genes <- str_split(input$list_in, pattern = "[, ]+")[[1]]
    return(genes)
  })

  # Choose genes of interest.
  chosen_inputs2 <- eventReactive(input$go_button, {
  # Split your input based on commas or spaces
    genes <- str_split(input$list_in2, pattern = "[, ]+")[[1]]
    return(genes)
  })

  # Choose genome version. Manily for mouse mm10 vs mm39
  chosen_genome <- reactive({
    if (input$species == "Mmusculus") {
      input$mmver
    } else {
      input$hsver
    }
  })

  ## Spit peak inputs into chromosome start end
  ## Should give back in a table format.
  chosen_peaks <- eventReactive(input$go_button, {
    req(input$list_in_peak)
    tryCatch({
      check1 <- data.table(
        peak = unlist(strsplit(input$list_in_peak, split = "[, ]+"))
      )
      check1[, c("chromosome", "start", "end") := tstrsplit(
        peak,
        split = "[.,]+",
        fixed = FALSE
      )]
      check1[
          , start := as.integer(start)
        ][
          , end := as.integer(end)
        ][
          , symbol := paste(chromosome,start,end,sep = "_")
        ]
      check1 <- unique(na.omit(check1[, .(chromosome, start, end, symbol)]))
      # reorder start/end if necessary
      check2 <- copy(check1)
      check2[
        , start2 := start 
      ][
        , end2 := end
      ]
      check2[
        , start := fifelse(start2 > end2, end, start)
      ][
        , end := fifelse(start2 > end2, start2, end2)
      ]
      check2 <- check2[, .(chromosome, start, end, symbol)]
      return(unique(check2))
    }, error = function(e) {
      showNotification("Invalid peak format! Use: chr.start.end", type = "error")
      NULL
    })
  })

  ## Spit peak inputs into chromosome start end
  ## Should give back in a table format.
  chosen_peaks2 <- eventReactive(input$go_button, {
    req(input$list_in_peak2)
    tryCatch({
      check1 <- data.table(
        peak = unlist(strsplit(input$list_in_peak2, split = "[, ]+"))
      )
      check1[, c("chromosome", "start", "end") := tstrsplit(
        peak,
        split = "[.,]+",
        fixed = FALSE
      )]
      check1[
          , start := as.integer(start)
        ][
          , end := as.integer(end)
        ][
          , symbol := paste(chromosome,start,end,sep = "_")
        ]
      check1 <- unique(na.omit(check1[, .(chromosome, start, end, symbol)]))
      # reorder start/end if necessary
      check2 <- copy(check1)
      check2[
        , start2 := start 
      ][
        , end2 := end
      ]
      check2[
        , start := fifelse(start2 > end2, end, start)
      ][
        , end := fifelse(start2 > end2, start2, end2)
      ]
      check2 <- check2[, .(chromosome, start, end, symbol)]
      return(unique(check2))
    }, error = function(e) {
      showNotification("Invalid peak format! Use: chr.start.end", type = "error")
      NULL
    })
  })

  ## Start Computations
  ## Get the tables ready
  input1_tables <- reactive({
    req(chosen_inputs1())
    cat1 <- run_CAT(
    species = input$species,
    genome = chosen_genome(),
    genes = chosen_inputs1()
    )
    return(cat1)
  })

  ## Get the tables ready
  input2_tables <- reactive({
    req(chosen_inputs2())
    cat2 <- run_CAT(
    species = input$species,
    genome = chosen_genome(),
    genes = chosen_inputs2()
    )
    return(cat2)
  })
  ## End Computations

  ## Rendering Tables
  output$risk_render1 <- DT::renderDataTable({
    req(input1_tables())
    return(plot_tab_CAT(chosen_chr(), input1_tables()[1] ))
  },
  options = list(
  dom = 't',
  scrollX = TRUE,    # horizontal scrolling
  scrollY = "400px", # vertical scrolling
  paging = FALSE   # disables pagination if you want scroll only,
  ),
  filter = "top"
  )

  output$risk_render2 <- DT::renderDataTable({
    req(input2_tables())
    return(plot_tab_CAT(chosen_chr(), input2_tables()[[1]]))
  },
  options = list(
  dom = 't',
  scrollX = TRUE,    # horizontal scrolling
  scrollY = "400px", # vertical scrolling
  paging = FALSE   # disables pagination if you want scroll only,
  ),
  filter = "top"
  )

  ## Rendering Plots
  output$risk_plot <- renderPlot({
    req(input1_tables())
    req(input2_tables())
    return(risky_plot(input$species, chosen_chr(), input1_tables(), input2_tables()))
  })

}