#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom grid arrow unit
#' @importFrom ggnewscale new_scale_color
#' @importFrom data.table as.data.table
#' @keywords internal

## Plot for Risk region
risky_plot <- function(species, chromosome_vec, front_df1, front_df2, background = mouse_to_human_genes) {
  species = 'Mmusculus'
  back_df <- copy(mouse_to_human_genes)
  back_df <- as.data.table(back_df)[
    , .(mmusculus_gene_name, hsapiens_gene_name)
  ][ 
    , symbol := if(species == "Mmusculus") mmusculus_gene_name else hsapiens_gene_name
  ][ 
    , symbol_o := if(species == "Mmusculus") hsapiens_gene_name else mmusculus_gene_name
  ][ 
    , .(symbol, symbol_o)
  ]
  if (species == "Mmusculus") {
    back_df <- back_df[
        mouse_all_startend, on = 'symbol'
    ][ 
        human_all_startend, on = c('symbol_o' = 'symbol')
    ]
    cumulative <- copy(mm_cumulative_df)
    setnames(cumulative, c('chromosome', 'position'))
    cumulative[
        , position := as.numeric(position)
    ][
      , position := data.table::shift(position, n = 1, type = 'lag')
    ][ 
      , position := fifelse(is.na(position), 0, position)
    ]
  } 
  if (species == 'Hsapiens') {
    back_df <- back_df[
        human_all_startend, on = 'symbol'
    ][ 
        mouse_all_startend, on = c('symbol_o' = 'symbol')
    ]
    cumulative <- copy(hs_cumulative_df)
    setnames(cumulative, c('chromosome', 'position'))
    cumulative[
        , position := as.numeric(position)
    ][
      , position := data.table::shift(position, n = 1, type = 'lag')
    ][ 
      , position := fifelse(is.na(position), 0, position)
    ]
  }
  old_names <- names(back_df)           
  has_prefix <- grepl("i\\.", old_names)
  names_to_fix <- old_names[has_prefix]
  new_names <- gsub("i\\.", "", names_to_fix)  
  new_names <- paste0(new_names, '_o')
  old_names[has_prefix] <- new_names
  setnames(back_df, names(back_df), old_names)
  back_df <- na.omit(back_df)
  back_df <- back_df[
   cumulative, on = 'chromosome'
  ][
    , mid := (start + end)/2 + position
  ][ 
    , mid_o := (start_o + end_o)/2
  ]

  up_df1 <- front_df1[
    , .(group_id, symbol, chromosome, start, end, chromosome_o, start_o, end_o, In_Risk, Risk_status, Group_status)
  ][ 
    , status := fifelse(In_Risk == TRUE, 'In Risk',
        fifelse(In_Risk == FALSE & Risk_status == 'By Breakpoint', 'Group in Risk', 
        fifelse(Risk_status == 'Group Split', 'Group Split', 'Not in Risk')))
  ][
    cumulative, on = 'chromosome'
  ][
    , mid := (start + end)/2 + position
  ][ 
    , mid_o := (start_o + end_o)/2
  ][ 
    , type := fifelse(grepl("chr", symbol), 0.75, 1)
  ][ 
    , show := fifelse(is.na(status), 0.75, 1)
  ][
    , type_spacing := fifelse(grepl("chr", symbol), 0.9, 1)
  ]
  down_df1 <- front_df2[
    , .(group_id, symbol, chromosome, start, end, chromosome_o, start_o, end_o, In_Risk, Risk_status, Group_status)
  ][ 
    , status := fifelse(In_Risk == TRUE, 'In Risk',
        fifelse(In_Risk == FALSE & Risk_status == 'By Breakpoint', 'Group in Risk', 
        fifelse(Risk_status == 'Group Split', 'Group Split', 'Not in Risk')))
  ][
    cumulative, on = 'chromosome'
  ][
    , mid := (start + end)/2 + position
  ][ 
    , mid_o := (start_o + end_o)/2
  ][ 
    , type := fifelse(grepl("chr", symbol), 0.75, 1)
  ][ 
    , show := fifelse(is.na(status), 0.75, 1)
  ][
    , type_spacing := fifelse(grepl("chr", symbol), 0.9, 1)
  ]
  
  fil_back <- back_df[chromosome %in% chromosome_vec & !grepl("chr", symbol)]
  up_df1 <- up_df1[chromosome %in% chromosome_vec]
  down_df1 <- down_df1[chromosome %in% chromosome_vec]
  noto_space <- max(fil_back$mid_o) * 1.3
  yoto_space <- max(fil_back$mid_o) * 0.1
  
  risk_plot <- ggplot(data = fil_back) +
    geom_segment(
      aes(x = mid, xend = mid, y=0, yend = mid_o, color = chromosome_o), 
      alpha = 1
      ) +
    geom_segment(
      aes(x = mid, xend = mid, y=0, yend = -mid_o, color = chromosome_o), 
      alpha = 1
      ) +
    labs(
      x = paste(species, "Gene Midpoint Position"),
      y = "Ortholog Midpoint Position",
      color = "Ortholog Chromosome"
    ) +
    scale_color_manual(values=matching_colors) +
    new_scale_color() +
    geom_segment(
      data = up_df1[!is.na(chromosome_o)],
      aes(x = mid, xend = mid, yend = mid_o + yoto_space, y = mid_o + 2 * yoto_space * type, color = status, alpha = show),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) + # for one to one orthos
    geom_point(data = up_df1[is.na(chromosome_o)], aes(x = mid, y = noto_space/type_spacing , color = status, alpha = show), shape = 6) + # for no orthos
    geom_segment(
      data = down_df1[!is.na(chromosome_o)], 
      aes(x = mid, xend = mid, yend = -mid_o - yoto_space, y = -mid_o - 2 * yoto_space * type, color = status, alpha = show),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) +
    geom_point(data = down_df1[is.na(chromosome_o)], aes(x = mid, y = -noto_space/type_spacing, color = status, alpha = show), shape = 2) +
    scale_color_manual(values=c("In Risk" = "#FF5733", "Not in Risk" = "blue", "Group in Risk" = "#E68613", "Group Split" = "darkorchid2")) +
    guides(alpha = "none")
  return(risk_plot)
}

## Helpers for heatmaps
get_order_pluspeaks <- function(cat_obj1, cat_obj2, mode = 'all') {
    if (mode == 'grouped') {
    gene_order <- unique(rbind(
        cat_obj1[[1]][!is.na(group_id)][ , .(symbol, chromosome, start)],
        cat_obj2[[1]][!is.na(group_id)][, .(symbol, chromosome, start)]
    ))
    }
    if (mode == 'all') {
    gene_order <- unique(rbind(
        cat_obj1[[1]][ , .(symbol, chromosome, start)],
        cat_obj2[[1]][, .(symbol, chromosome, start)]
    ))
    }
    
    gene_order <- gene_order[ 
        order(start)
    ]
    gene_order[
        , coord := 1
    ][ 
        , coord := cumsum(ifelse(grepl('chr', symbol), 1L, 2L))
    ][ 
        , diff := coord - shift(coord)
    ][ 
        , diff := fifelse(is.na(diff) & grepl('chr', symbol), 1, fifelse(is.na(diff) & !grepl('chr', symbol), 2, diff))
    ][
        , n_coord := shift(coord) + 0.5*diff
    ][
        , n_coord := fifelse(is.na(n_coord), 0.5*coord, n_coord)
    ][ 
        , coord := n_coord
    ][ 
        , n_coord := NULL
    ]
  return(gene_order)
}

## Plot for Groupings
## Prepare the heatmap matrices
matrices_for_plot <- function(species, chromosome_vec, cat_obj1, cat_obj2, mode = 'all') {
  order_dt <- get_order_pluspeaks(cat_obj1, cat_obj2, mode)[
    , .(symbol,diff,coord)
  ]
  if (mode == 'all') {
    up_filtr <- cat_obj1[[2]]
    down_filtr <- cat_obj2[[2]]
  }
  if (mode == 'grouped') {
    up_filtr <- cat_obj1[[2]][!is.na(group_id) & !is.na(group_id2)]
    down_filtr <- cat_obj2[[2]][!is.na(group_id) & !is.na(group_id2)]
  }

  ## Top half of the distances
  in_group_df <- up_filtr[
    chromosome == chromosome2
  ][
    (chromosome %in% chromosome_vec) & (chromosome2 %in% chromosome_vec)
  ][
    , .(symbol, symbol2, symbol_o, distance, grouped, Group_status, one_to_one, group_id, group_id2)
  ]
  in_group_df <- order_dt[ 
   in_group_df, on = .(symbol)
  ]

  in_group_df <- order_dt[ 
    in_group_df, on = list(symbol = symbol2)
  ]
  in_group_df <- in_group_df[
    , symbol2 := i.symbol
  ][
    , i.symbol := NULL
  ][ 
    , distance := abs(distance)
  ][
    , .(symbol,symbol2,coord,diff,i.coord,i.diff,distance, group_id, group_id2, grouped, Group_status, one_to_one)
  ][
    order(coord, i.coord)
  ][
    i.coord >= coord
  ]
  # Diagonals
  in_group_diag <- in_group_df[
    symbol == symbol2
  ][, .(
    x = c(coord - 0.5*diff, coord - 0.5*diff, coord + 0.5*diff),
    y = c(i.coord - 0.5*i.diff, i.coord + 0.5*i.diff, i.coord + 0.5*i.diff),
    distance = distance,
    grouped = grouped,
    Group_status = Group_status,
    one_to_one = one_to_one,
    id = paste(coord, i.coord, "A", sep = "_")
    ), by = .(coord, i.coord)
  ]
  # The Blocks
  upper_blocks <- in_group_df[
    symbol != symbol2
  ][
    , .(
    x = c(coord - 0.5*diff, coord + 0.5*diff, coord + 0.5*diff, coord - 0.5*diff),
    y = c(i.coord - 0.5*i.diff, i.coord - 0.5*i.diff, i.coord + 0.5*i.diff, i.coord + 0.5*i.diff),
    distance = distance,
    id = paste(coord, i.coord, "A", sep = "_")
    ), by = .(coord, i.coord)
  ]
  # Groups
  upper_group <- in_group_df[
        group_id == group_id2
    ][
        , min_x := min(coord), by = group_id
    ][ 
        , max_x := max(coord), by = group_id
    ][ 
        , min_y := min(i.coord), by = group_id
    ][ 
        , max_y := max(i.coord), by = group_id
    ][
    , .SD[c(1L, .N)], by = group_id
  ][, left_diff := {
    # Check the first row's symbol in this group
    first_symbol <- symbol[1L]

    # Choose min or max of diff based on whether it starts with 'chr'
    if (startsWith(first_symbol, "chr")) {
      min(diff, na.rm = TRUE)
    } else {
    max(diff, na.rm = TRUE)
  }
  }, by = group_id
  ][, right_diff := {
  # Check the first row's symbol in this group
  first_symbol <- symbol[2L]

  # Choose min or max of diff based on whether it starts with 'chr'
  if (startsWith(first_symbol, "chr")) {
    min(diff, na.rm = TRUE)
  } else {
    max(diff, na.rm = TRUE)
  }
  }, by = group_id
  ][, bottom_i.diff := {
  # Check the first row's symbol in this group
  first_symbol <- symbol[1L]

  # Choose min or max of diff based on whether it starts with 'chr'
  if (startsWith(first_symbol, "chr")) {
    min(i.diff, na.rm = TRUE)
  } else {
    max(i.diff, na.rm = TRUE)
  }
  }, by = group_id
  ][, top_i.diff := {
  # Check the first row's symbol in this group
  first_symbol <- symbol[2L]

  # Choose min or max of diff based on whether it starts with 'chr'
  if (startsWith(first_symbol, "chr")) {
    min(i.diff, na.rm = TRUE)
  } else {
    max(i.diff, na.rm = TRUE)
  }
  }, by = group_id
  ]
  upper_group <- unique(upper_group[
    , .(group_id,grouped,Group_status,one_to_one,min_x,max_x,min_y,max_y,left_diff,right_diff,bottom_i.diff,top_i.diff)
  ])[, .(
    x = c(min_x - 0.5*left_diff, min_x - 0.5*left_diff, max_x + 0.5*right_diff),
    y = c(min_y - 0.5*bottom_i.diff, max_y + 0.5*top_i.diff, max_y + 0.5*top_i.diff),
    grouped = grouped,
    Group_status = Group_status,
    one_to_one = one_to_one,
    id = group_id
    ), by = .(group_id)
  ]

  ## Bottom half of the distances
  in_group_df2 <- down_filtr[
    (chromosome %in% chromosome_vec) & (chromosome2 %in% chromosome_vec)
  ][
    chromosome == chromosome2
  ][
    , .(symbol, symbol2, symbol_o, distance, grouped, Group_status, one_to_one, group_id, group_id2)
  ]
  in_group_df2 <- order_dt[ 
   in_group_df2, on = .(symbol)
  ]
  in_group_df2 <- order_dt[ 
    in_group_df2, on = list(symbol = symbol2)
  ]
  in_group_df2 <- in_group_df2[
    , symbol2 := i.symbol
  ][
    , i.symbol := NULL
  ][ 
    , distance := abs(distance)
  ][
    , .(symbol,symbol2,coord,diff,i.coord,i.diff,distance, group_id, group_id2, grouped, Group_status, one_to_one)
  ][
    order(coord, i.coord)
  ][
     i.coord <= coord
  ]
  
  in_group_diag2 <- in_group_df2[
    symbol == symbol2
  ][, .(
    x = c(coord - 0.5*diff, coord + 0.5*diff, coord + 0.5*diff),
    y = c(i.coord - 0.5*i.diff, i.coord - 0.5*i.diff, i.coord + 0.5*i.diff),
    distance = distance,
    grouped = grouped,
    Group_status = Group_status,
    one_to_one = one_to_one,
    id = paste(coord, i.coord, "A", sep = "_")
    ), by = .(coord, i.coord)
  ]

  lower_blocks <- in_group_df2[
    symbol != symbol2
  ][, .(
    x = c(coord - 0.5*diff, coord + 0.5*diff, coord + 0.5*diff, coord - 0.5*diff),
    y = c(i.coord - 0.5*i.diff, i.coord - 0.5*i.diff, i.coord + 0.5*i.diff, i.coord + 0.5*i.diff),
    distance = distance,
    id = paste(coord, i.coord, "A", sep = "_")
    ), by = .(coord, i.coord)
  ]

  # Groups
  lower_group <- in_group_df2[
        group_id == group_id2
  ][
    , min_x := min(coord), by = group_id
  ][ 
    , max_x := max(coord), by = group_id
  ][ 
    , min_y := min(i.coord), by = group_id
  ][ 
    , max_y := max(i.coord), by = group_id
  ][
    , .SD[c(1L, .N)], by = group_id
  ][, left_diff := {
    # Check the first row's symbol in this group
    first_symbol <- symbol[1L]

    # Choose min or max of diff based on whether it starts with 'chr'
    if (startsWith(first_symbol, "chr")) {
      min(diff, na.rm = TRUE)
    } else {
    max(diff, na.rm = TRUE)
  }
  }, by = group_id
  ][, right_diff := {
  # Check the first row's symbol in this group
  first_symbol <- symbol[2L]

  # Choose min or max of diff based on whether it starts with 'chr'
  if (startsWith(first_symbol, "chr")) {
    min(diff, na.rm = TRUE)
  } else {
    max(diff, na.rm = TRUE)
  }
  }, by = group_id
  ][, bottom_i.diff := {
  # Check the first row's symbol in this group
  first_symbol <- symbol[1L]

  # Choose min or max of diff based on whether it starts with 'chr'
  if (startsWith(first_symbol, "chr")) {
    min(i.diff, na.rm = TRUE)
  } else {
    max(i.diff, na.rm = TRUE)
  }
  }, by = group_id
  ][, top_i.diff := {
  # Check the first row's symbol in this group
  first_symbol <- symbol[2L]

  # Choose min or max of diff based on whether it starts with 'chr'
  if (startsWith(first_symbol, "chr")) {
    min(i.diff, na.rm = TRUE)
  } else {
    max(i.diff, na.rm = TRUE)
  }
  }, by = group_id
  ]
  lower_group <- unique(lower_group[
    , .(group_id,grouped,Group_status,one_to_one,min_x,max_x,min_y,max_y,left_diff,right_diff,bottom_i.diff,top_i.diff)
  ])[, .(
    x = c(min_x - 0.5*left_diff, max_x + 0.5*right_diff, max_x + 0.5*right_diff),
    y = c(min_y - 0.5*bottom_i.diff, min_y - 0.5*bottom_i.diff, max_y + 0.5*top_i.diff),
    grouped = grouped,
    Group_status = Group_status,
    one_to_one = one_to_one,
    id = group_id
    ), by = .(group_id)
  ]

  return(list(in_group_diag, in_group_diag2, upper_blocks, lower_blocks, upper_group, lower_group, in_group_df, in_group_df2))
}

## Plot for Distance
heatmap_plot <- function(list_matrices, type = 'all') {
  #list_matrices <- matrices_for_plot(species, chromosome_vec, cat_obj1, cat_obj2, mode)
  if (type == 'distance') {
    plot <- ggplot() +
    geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = log10(abs(distance) + 1))) +
    geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
    geom_polygon(data = list_matrices[[3]], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
    geom_polygon(data = list_matrices[[4]], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
    geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, color = grouped), fill = NA, linewidth = 0.5) +
    geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, color = grouped), fill = NA, linewidth = 0.5) +
    scale_fill_gradient(low = "navyblue", high = "lightgoldenrodyellow") +
    coord_fixed() +
    theme(legend.position = 'bottom') +
    labs(
        x = "List2",   
        y = "List1" 
    )
  }
  if (type == 'grouped') {
    plot <- ggplot() +
    geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = grouped)) +
    geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = grouped)) +
    geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5, linewidth = 0.5) +
    geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5, linewidth = 0.5) +
    coord_fixed() +
    theme(legend.position = 'bottom') +
    labs(
        x = "List2",   
        y = "List1" 
    )
  }
  if (type == 'sequence') {
    plot <- ggplot() +
    geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = one_to_one), color = 'black', linewidth = 0.1) +
    geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = one_to_one), color = 'black', linewidth = 0.1) +
    geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, fill = one_to_one), alpha = 0.5) +
    geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, fill = one_to_one), alpha = 0.5) +
    coord_fixed() +
    theme(legend.position = 'bottom') +
    labs(
        x = "List2",   
        y = "List1" 
    )
  }
  if (type == 'location') {
    plot <- ggplot() +
    geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = Group_status), color = "black", linewidth = 0.1) +
    geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = Group_status), color = "black", linewidth = 0.1) +
    geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, fill = Group_status), alpha = 0.5) +
    geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, fill = Group_status), alpha = 0.5) +
    coord_fixed() +
    theme(legend.position = 'bottom') +
    labs(
        x = "List2",   
        y = "List1" 
    )
  }
  return(plot)
}

## Plot for Group Statistics
plot_group_stats <- function(input_table, species, type, label) {
    #Group_breaks <- c(
    #'Became Singlet',
    #'Lost',
    #'No Loss',
    #'Partial Loss',
    #'Singlet')
    summary_table <- unique(input_table[
        !grepl("chr", symbol)
    ][ 
    , .(group_id, one_to_one, Risk_status, Group_status)
    ])
    summary_table <- unique(summary_table[
        , Group_count := .N, by = .(one_to_one)
    ][ 
        , Risk_count := .N, by = .(one_to_one, Risk_status)
    ][
        , Split_count := .N, by = .(one_to_one, Group_status)
    ][
        , group_id := NULL
    ][
        order(one_to_one)
    ])
    if (type == 'Risk') {
        summary_table <- unique(summary_table[
            , .(one_to_one, Risk_status, Group_count, Risk_count)
        ])
        summary_table[ 
            , one_to_one_x := as.numeric(factor(summary_table$one_to_one))
        ][
            , y_cumulative := cumsum(Risk_count), by = one_to_one_x
        ][
            , y_base := shift(y_cumulative), by = one_to_one_x
        ][
            , y_base := fifelse(is.na(y_base), 0, y_base)
        ][
            , how_many := .N, by = one_to_one_x
        ][
            , Status := Risk_status
        ][ 
            , Risk_status := NULL
        ][ 
            , Count := Risk_count
        ][ 
            , Risk_count := NULL
        ]
        phrase <- "\nGroups by Breakpoint"
    } 
    if (type == 'Split') {
        summary_table <- unique(summary_table <- summary_table[
            , .(one_to_one, Group_status, Group_count, Split_count)
        ])
        summary_table[ 
            , one_to_one_x := as.numeric(factor(summary_table$one_to_one))
        ][
            , y_cumulative := cumsum(Split_count), by = one_to_one_x
        ][
            , y_base := shift(y_cumulative), by = one_to_one_x
        ][
            , y_base := fifelse(is.na(y_base), 0, y_base)
        ][
            , how_many := .N, by = one_to_one_x
        ][
            , Status := Group_status
        ][ 
            , Group_status := NULL
        ][ 
            , Count := Split_count
        ][ 
            , Split_count := NULL
        ]
        phrase <- "\nGroups by Split"
    }
    Group_breaks <- unique(summary_table$one_to_one)
    plot <- ggplot()  +
    geom_rect(
        data = summary_table,
        aes(xmin = one_to_one_x - 0.5, xmax = one_to_one_x + 0.5, ymin = y_base, ymax = y_cumulative, fill = Status)
    ) + 
    geom_text(
        data = unique(summary_table[,.(one_to_one_x, Group_count)]),
        aes(
            x = one_to_one_x, 
            y = Group_count * 1.1, 
            label = paste("Group total:", Group_count, sep = '\n')
        )
    ) +
    geom_text(
        data = summary_table,
        aes(
            x = one_to_one_x - 0.25,
            y = y_cumulative,
           label = ifelse(how_many == 1, NA, Count)
        )
    ) +
    scale_x_continuous(
        breaks = rep(1:max(summary_table$one_to_one_x) ),
        labels = Group_breaks
    ) +
    labs(
        x = "Group Status",   
        y = "Group composition" 
    ) +
    ggtitle(paste(species, label, phrase, sep = ' '))
    return(plot)
}