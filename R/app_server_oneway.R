#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table fifelse
#' @importFrom data.table copy
#' @importFrom stringr str_split
#' @importFrom ggnewscale new_scale_color
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @keywords internal

app_server_oneway <- function(input, output, session) {
  output$risk_example1 <- renderPlot({
    ex_back <- data.table(
      x = seq(0, 0.5, length.out = 10), 
      y = seq(0,1, length.out = 10),
      Ortholog_Chromosome = c(rep("chr1",10))
    )
    ex_risk <- data.table(
      x = c(1, 1.5, 2),
      y = c(1.5, 1.5, 2),
      y2 = c(1.1, 1.1, 1.1),
      Status = c("Unattached Peak", rep("By Breakpoint", 2))
    )
    ex_lost <- data.table(
      x = c(2.5, 3),
      y = c(3, 4),
      Status = rep("By Breakpoint", 2),
      ty = c("Peak", "Gene")
    )
    
    ex_plot <- ggplot() +
    geom_segment(data = ex_back, aes(x = x, xend = x, y = 0, yend = y, color = Ortholog_Chromosome)) +
    new_scale_color() +
    geom_segment(
      data = ex_risk, 
      aes(x = x, xend = x, y = y, yend = y2, color = Status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
    geom_point(
      data = ex_lost,
      aes(x = x, y = y, shape = ty, color = Status)
    ) +
    scale_color_manual(values=c(
      "By Breakpoint and Split" = "#d73027",
      "By Breakpoint" = "orange", 
      "Group Split" = "darkorchid2",
      'Group Together' = 'blue',
      'Not Applicable' = 'black',
      'Unattached Peak' = 'grey'
      ),
      breaks = c("By Breakpoint and Split", "By Breakpoint", "Group Split", 'Group Together', 'Not Applicable', 'Unattached Peak')
      ) +
    scale_shape_manual(
      name = 'Missing Ortholog',
      values = c("Gene" =  5, "Peak" = 6)
    )
    return(ex_plot)
  })

  output$grouping_example1 <- renderPlot({
    ex_features <- data.table(
      x = c(0, 0, 1, 1, 1, 1, 3, 3, 3, 3, 4, 4, 4, 4, 6, 6, 6, 6, 8, 8),
      y = c(0, 1, 1, 0, 1, 3, 3, 1, 3, 4, 4, 3, 4, 6, 6, 4, 6, 8, 8, 6),
      id = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5),
      Has_ortholog = c(rep('YES', 8), rep('NO', 8), rep('YES', 4)),
      grouped = c(rep('Grouped', 16), rep('Singlet group', 4))
    )
    ex_group <- data.table(
      x = c(0, 0, 6, 6),
      y = c(0,6, 6, 0),
      id = c(1, 1, 1, 1),
      grouped = c('Grouped', 'Grouped', 'Grouped', 'Grouped')
    )
    ex_plot <- ggplot() +
      geom_polygon(data = ex_features, aes(x = x ,y = y, group = id, alpha = Has_ortholog, fill = grouped), color = 'black', linewidth = 0.1) +
      geom_polygon(data = ex_group, aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5) +
      scale_fill_manual(values = c(
      "Grouped" = "#00BFC4",
      "Singlet group" = "#F8766D"
      ))
    return(ex_plot)
  })

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
  
  ## Must transpose the distance comparison table
  input2_tables <- reactive({
    req(input1_tables())
    list_of_tables <- copy(input1_tables())
    transpose_table <- list_of_tables[[2]][
      !is.na(distance)
    ][
      , .(symbol, symbol2, distance)
    ]
    colnames(transpose_table) <- c('symbol2', 'symbol', 'distance')
    list_of_tables[[2]] <- list_of_tables[[2]][ 
      transpose_table, on = .(symbol, symbol2)
    ][ 
      , distance := fifelse(is.na(distance), i.distance, distance)
    ][ 
      , i.distance := NULL
    ]
    return(list_of_tables)
  })
  
  ## Heatmap matrices
  heatmap_matrices <- reactive({
    req(input1_tables())
    heatmap_for_plot <- matrices_for_plot(
      input$species,
      chosen_chr(), 
      input2_tables(), 
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

  ## Orthologous tables
  ortho_table1 <- reactive({
    req(input1_tables())
    return(ortho_summary(chosen_inputs1(), chosen_peaks1(), input1_tables()[[1]], input$chosen_ortho_stat))
  })

  ## Attached Unattached Table
  peak_attach_table1 <- reactive({
    req(input1_tables()[[1]])
    req(chosen_chr())
    return(table_attached_peaks(input1_tables()[[1]], chosen_chr()))
  })

  ## Rendering Tables
  ## Risk Tables
  render1 <- reactive({
    req(input1_tables())
    render1_table <- plot_tab_CAT(chosen_chr(), input1_tables()[1] )
    return(render1_table)
  })
  output$risk_render1 <- DT::renderDT({
    req(render1())
    return( 
      DT::datatable(
        render1(), 
        rownames = FALSE,
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
  ## Downloadable Table
  output$downloadTable1 <- downloadHandler(
    filename = function() {
      paste("List1table", Sys.Date(), ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(render1(), file, row.names = FALSE)
    }
  )

  ## Position and Group ID
  output$position_groups_table <- DT::renderDT({
    pos_groups <- input1_tables()[[1]][ 
      , List1_group := as.character(group_id)
    ][ 
      chromosome %in% chosen_chr()
    ][ 
      , .(symbol, symbol_o, List1_group, chromosome, start, end)
    ]
    return(
      DT::datatable(pos_groups, 
      rownames = FALSE, options = list(
        dom = 't',
        #scrollX = TRUE,      
        scrollY = "300px", 
        #scrollCollapse = TRUE,
        paging = FALSE       
        #pageLength = 5       
        ),
        filter = "top")
    )
    })

  ## Render hovered groups table
  output$hovered_groups_table <- DT::renderDT({
    hovered_data <- hovered_groups()
    hovered_data <- input1_tables()[[1]][, .(symbol)][
      hovered_data, on = .(symbol)
    ]
    return(DT::datatable(hovered_data,     
                  rownames = FALSE,
                  options = list(
                    dom = 't',
                    scrollX = TRUE,
                    scrollY = "300px",
                    paging = FALSE
                  ))
    )
  })

  ## Render Attach Unattach Table
  output$attached_peaks1 <- DT::renderDT({
    return(peak_attach_table1())
  })
  ## Downloadable Table
  output$Downloadable_peaktable1 <- downloadHandler(
    filename = function() {
      paste("Peak_attach", Sys.Date(), ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(peak_attach_table1(), file, row.names = FALSE)
    }
  )

  ## Render orthologous tables
  output$orthologous_genepeak_table1 <- DT::renderDT({
    req(ortho_table1())
    return(
      DT::datatable(ortho_table1(),
      rownames = FALSE, options = list(
        dom = 't',
        scrollX = TRUE,
        scrollY = "300px",
        paging = FALSE
        ),
        filter = "top")
    )
  })
  ## Downloadable Table
  output$Downloadable_orthotable1 <- downloadHandler(
    filename = function() {
      paste("Ortho_status", Sys.Date(), ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(ortho_table1(), file, row.names = FALSE)
    }
  )

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
  ## Attach Unattached Peaks
  output$perc_attached_peaks1 <- renderPlot({
    req(peak_attach_table1())
    return(plot_attached_peaks(input1_tables()[[1]], chosen_chr()))
  })

# ortho_summary_prop <- function(ortho_table, type) {
  output$orthologous_genepeak_plot1 <- renderPlot({
    req(ortho_table1())
    return(ortho_summary_prop(ortho_table1(), input$chosen_ortho_stat))
  })

  ## Group Size Summaries
  output$singlet_summaryplot <- renderPlot({
    req(input1_tables()[[1]])
    return(singlet_summary(input1_tables()[[1]], chosen_chr() ) )
  })

  ## Group Size Summaries
  output$groupsize_plot <- renderPlot({
    req(input1_tables()[[1]])
    return(group_size_plot(input1_tables()[[1]], chosen_chr(), input$chosen_region_stat))
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

}