#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table fifelse
#' @importFrom data.table copy
#' @importFrom stringr str_split
#' @keywords internal

app_server_oneway <- function(input, output, session) {

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
  ## Heatmap matrices
  heatmap_matrices <- reactive({
    req(input1_tables())
    heatmap_for_plot <- matrices_for_plot(
      input$species,
      chosen_chr(), 
      input1_tables(), 
      input1_tables(), 
      mode = 'grouped')
    return(heatmap_for_plot)
  })
  ## End Computations

  ## Reactive Table for Positions
  position_groups <- reactive({
    req(input1_tables())
    return(search_up(input1_tables()[[1]], input1_tables()[[1]], chosen_chr(), version = 'One'))
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
    ## Return the boundaries as a vector
    return(c(new_group1_min, new_group1_max, new_group1_min, new_group1_max))
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
  
  ## Position and Group ID
  output$position_groups_table <- DT::renderDT({
    pos_groups <- copy(position_groups())[
      , List1_group := as.character(List1_group)
    ]
    return(
      DT::datatable(pos_groups, 
      rownames = FALSE, options = list(
        dom = 't',
        #scrollX = TRUE,      
        scrollY = "200px", 
        #scrollCollapse = TRUE,
        paging = FALSE       
        #pageLength = 5       
        ),
        filter = "top")
    )
    })


  ## Rendering Plots
  ## Risky Plot
  output$risk_plot <- renderPlot({
    return(risky_plot(input$species, chosen_chr(), input1_tables()[[1]], input1_tables()[[1]]))
  })

  ## Heatmap Plots
  output$distance_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'distance'))
    return(heatmap_plot(heatmap_matrices(), type = 'distance', mode = 'shiny', bounds = zoom_groups(), version = 'One' ))
  })
  output$grouping_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'grouped'))
    return(heatmap_plot(heatmap_matrices(), type = 'grouped', mode = 'shiny', bounds = zoom_groups(), version = 'One'))
  })
  output$sequence_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'sequence'))
    return(heatmap_plot(heatmap_matrices(), type = 'sequence', mode = 'shiny', bounds = zoom_groups(), version = 'One'))
  })
  output$location_plot <- renderPlot({
    req(zoom_groups())
    req(heatmap_matrices())
    #return(heatmap_plot(heatmap_matrices(), type = 'location'))
    return(heatmap_plot(heatmap_matrices(), type = 'location', mode = 'shiny', bounds = zoom_groups(), version = 'One'))
  })
  
  ## Group Statistics Plots
  ## By Breakpoint
  output$group_by_break1 <- renderPlot({
    req(input1_tables()[[1]])
    req(chosen_chr())
    return(plot_group_stats(input1_tables()[[1]], input$species, 'Risk', "List1", chosen_chr(), input$group_gene_mode))
  })
  ## By Duplications/Insertions
  output$group_split1 <- renderPlot({
    req(input1_tables()[[1]])
    req(chosen_chr())
    return(plot_group_stats(input1_tables()[[1]], input$species, 'Split',"List1", chosen_chr(), input$group_gene_mode))
  })

  ## Click on heatmaps to show group IDs
  hovered_groups <- reactive({
    req(input$heatmap_click)
    req(heatmap_matrices())
    
    click <- input$heatmap_click
    
    # Get upper_group (list_matrices[[5]]) and lower_group (list_matrices[[6]])
    upper_groups <- heatmap_matrices()[[5]][, source := "Upper"]
    lower_groups <- heatmap_matrices()[[6]][, source := "Lower"]
    
    # Combine both data frames
    all_groups <- rbind(upper_groups, lower_groups)
    all_groups <- all_groups[
        , min_x := min(x), by = .(id, source)
      ][
        , max_x := max(x), by = .(id, source)
      ][
        , min_y := min(y), by = .(id, source)
      ][
        , max_y := max(y), by = .(id, source)
      ][ 
        , .( id, source, min_x, max_x, min_y, max_y)
      ]
    all_groups <- unique(all_groups)
    
    # Find groups where click coordinates fall within the polygon
    # The x and y values in the data frame define polygons
    click_x <- click$x
    click_y <- click$y

    all_groups <- all_groups[
      (click_x >= min_x & click_x <= max_x) &
      (click_y >= min_y & click_y <= max_y)
    ]
    display_table <- position_groups()[List1_group %in% all_groups[source == "Upper"]$id]
    return(display_table)
    
  })
  
  ## Render hovered groups table
  output$hovered_groups_table <- DT::renderDT({
    hovered_data <- hovered_groups()
    return(DT::datatable(hovered_data, 
                  rownames = FALSE,
                  options = list(
                    dom = 't',
                    scrollX = TRUE,
                    scrollY = "400px",
                    paging = FALSE
                  ))
    )
  })

}