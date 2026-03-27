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
  chosen_peaks1 <- eventReactive(input$go_button, {
    req(input$list_in_peak)
    tryCatch({
      check1 <- data.table(
        peak = unlist(strsplit(input$list_in_peak, split = "[, ]+"))
      )
      ## Parse in order to check if start > end
      ## Need to make sure everything is on the forward strand
      check1[, c("chromosome", "start", "end") := tstrsplit(
        peak,
        split = "[,\\._]+",
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
      check2 <- check2[, .(chromosome, start, end, symbol)][symbol != ""]
      return(unique(check2$symbol))
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
      ## Parse in order to check if start > end
      ## Need to make sure everything is on the forward strand
      check1[, c("chromosome", "start", "end") := tstrsplit(
        peak,
        split = "[,\\._]+",
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
      check2 <- check2[, .(chromosome, start, end, symbol)][symbol != ""]
      return(unique(check2$symbol))
    }, error = function(e) {
      showNotification("Invalid peak format! Use: chr.start.end", type = "error")
      NULL
    })
  })

  ## Start Computations
  ## Get the tables ready
  input1_genepeak <- reactive({
    req(chosen_inputs1())
    if (!is.null(input$list_in_peak) && nzchar(input$list_in_peak)) {
      distance_tables <- comb_genepeak(
        input$species, 
        chosen_genome(), 
        chosen_inputs1(), 
        chosen_peaks1(), 
        dist_to_group = input$dist_to_group, 
        exp_rocr = input$exp_rocr, 
        threshold = input$threshold)
      return(distance_tables)
    } else {
      distance_tables <- comb_genepeak(
        input$species, 
        chosen_genome(), 
        chosen_inputs1(), 
        dist_to_group = input$dist_to_group, 
        exp_rocr = input$exp_rocr, 
        threshold = input$threshold)
      return(distance_tables)
    }
  })
  input1_tables <- reactive({
    req(input1_genepeak())
    if (input$split_choice == 'Manual') {
      split_distance = input$split_dis
    } else {
      if (input$exp_rocr < 500000) {
        split_distance = 500000
      } else {
        split_distance = input$exp_rocr
      }
    }
    risk_tables <- conservation_table(
      input1_genepeak()[[2]], 
      input1_genepeak()[[1]],
      split_distance
      )
    final_tables <- break_pts(
      input$species, 
      risk_tables, 
      disttobreak = input$disttobreak)
    return(final_tables)
  })

  ## Get the tables ready
  input2_genepeak <- reactive({
    req(chosen_inputs2())
    if (!is.null(input$list_in_peak2) && nzchar(input$list_in_peak2)) {
      distance_tables <- comb_genepeak(
        input$species, 
        chosen_genome(), 
        chosen_inputs2(), 
        chosen_peaks2(), 
        dist_to_group = input$dist_to_group, 
        exp_rocr = input$exp_rocr, 
        threshold = input$threshold)
      return(distance_tables)
    } else {
      distance_tables <- comb_genepeak(
        input$species, 
        chosen_genome(), 
        chosen_inputs2(), 
        dist_to_group = input$dist_to_group, 
        exp_rocr = input$exp_rocr, 
        threshold = input$threshold)
      return(distance_tables)
    }
  })
  input2_tables <- reactive({
    req(input2_genepeak())
    if (input$split_choice == 'Manual') {
      split_distance = input$split_dis
    } else {
      if (input$exp_rocr < 500000) {
        split_distance = 500000
      } else {
        split_distance = input$exp_rocr
      }
    }
    risk_tables <- conservation_table(
      input2_genepeak()[[2]], 
      input2_genepeak()[[1]],
      split_distance
      )
    final_tables <- break_pts(
      input$species, 
      risk_tables, 
      disttobreak = input$disttobreak)
    return(final_tables)
  })
  ## End Computations

  ## Heatmap matrices
  heatmap_matrices <- reactive({
    req(input1_tables())
    req(input2_tables())
    heatmap_for_plot <- matrices_for_plot(
      input$species,
      chosen_chr(), 
      input1_tables(), 
      input2_tables(), 
      mode = 'grouped')
    return(heatmap_for_plot)
  })

  ## Reactive Table for Positions
  position_groups <- reactive({
    req(input1_tables())
    req(input2_tables())
    return(search_up(input1_tables()[[1]], input2_tables()[[1]], chosen_chr()))
  })

  ## Reactive slider for positions
  observe({
  req(position_groups())
  vals <- position_groups()
  
  # If you want a numeric slider & your reactive range is numeric:
  min_val <- min(vals$position)
  max_val <- max(vals$position)
  
  updateSliderInput(
    session,
    inputId = "pos_range",
    label   = "Position range:",
    min     = min_val,
    max     = max_val,
    value   = c(min_val, max_val)  # or any default within [min_val, max_val]
  )
  })

  ## Reactive positions for the groups within the chosen range of positions
  ## I hate to do this on the server.R, but putting this as a standalone function
  ## Might be overkill. It's pretty much only for the shiny app.
  zoom_groups <- reactive({
    req(input$pos_range)
    req(position_groups)
    ## Get the min and max positions
    min_val <- input$pos_range[1]
    max_val <- input$pos_range[2]
    new_positions <- position_groups()[
      position >= min_val & position <= max_val
    ]
    ## Get the boundary group id1 and group id2 within positions
    new_group1_min <- min(na.omit(new_positions$List1_group))
    new_group1_max <- max(na.omit(new_positions$List1_group))
    new_group2_min <- min(na.omit(new_positions$List2_group))
    new_group2_max <- max(na.omit(new_positions$List2_group))
    ## Return the boundaries as a vector
    return(c(new_group1_min, new_group1_max, new_group2_min, new_group2_max))
  })

  ## Rendering Tables
  ## Risk Tables
  output$risk_render1 <- DT::renderDT({
    req(input1_tables())
    return( 
      DT::datatable(plot_tab_CAT(chosen_chr(), input1_tables()[1] ), rownames = FALSE,
  options = list(
  dom = 't',
  scrollX = TRUE,    # horizontal scrolling
  scrollY = "400px", # vertical scrolling
  paging = FALSE   # disables pagination if you want scroll only,
  ),
  filter = "top"
  )
    )
  }
  )
  output$risk_render2 <- DT::renderDT({
    req(input2_tables())
    return(
      DT::datatable(plot_tab_CAT(chosen_chr(), input2_tables()[[1]]), rownames = FALSE,
  options = list(
  dom = 't',
  scrollX = TRUE,    # horizontal scrolling
  scrollY = "400px", # vertical scrolling
  paging = FALSE   # disables pagination if you want scroll only,
  ),
        filter = "top")
    )
  }
  )
  
  ## Position and Group ID
  output$position_groups_table <- DT::renderDT({
    return(
      DT::datatable(position_groups(), rownames = FALSE, options = list(
        dom = 't',
        scrollX = TRUE,    # horizontal scrolling
        scrollY = "400px", # vertical scrolling
        paging = FALSE   # disables pagination if you want scroll only,
        ),
        filter = "top")
    )
    })


  ## Rendering Plots
  ## Risky Plot
  output$risk_plot <- renderPlot({
    return(risky_plot(input$species, chosen_chr(), input1_tables()[[1]], input2_tables()[[1]], way = 'Two' ) )
  })

  ## Heatmap Plots
  output$distance_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'distance'))
    return(heatmap_plot(heatmap_matrices(), type = 'distance', mode = 'shiny', bounds = zoom_groups()))
  })
  output$grouping_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'grouped'))
    return(heatmap_plot(heatmap_matrices(), type = 'grouped', mode = 'shiny', bounds = zoom_groups()))
  })
  output$sequence_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'sequence'))
    return(heatmap_plot(heatmap_matrices(), type = 'sequence', mode = 'shiny', bounds = zoom_groups()))
  })
  output$location_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'location'))
    return(heatmap_plot(heatmap_matrices(), type = 'location', mode = 'shiny', bounds = zoom_groups()))
  })
  
  ## Group Statistics Plots
  ## By Breakpoint
  output$group_by_break1 <- renderPlot({
    req(input1_tables()[[1]])
    req(chosen_chr())
    return(plot_group_stats(input1_tables()[[1]], input$species, 'Risk', "List1", chosen_chr(), input$group_gene_mode))
  })
  output$group_by_break2 <- renderPlot({
    req(input2_tables()[[1]])
    req(chosen_chr())
    return(plot_group_stats(input2_tables()[[1]], input$species, 'Risk',"List2", chosen_chr(), input$group_gene_mode))
  })
  ## By Duplications/Insertions
  output$group_split1 <- renderPlot({
    req(input1_tables()[[1]])
    req(chosen_chr())
    return(plot_group_stats(input1_tables()[[1]], input$species, 'Split',"List1", chosen_chr(), input$group_gene_mode))
  })
  output$group_split2 <- renderPlot({
    req(input2_tables()[[1]])
    req(chosen_chr())
    return(plot_group_stats(input2_tables()[[1]], input$species, 'Split',"List2", chosen_chr(), input$group_gene_mode))
  })

}