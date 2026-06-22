#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 geom_violin
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
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 expand_limits
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 position_dodge
#' @importFrom forcats fct_rev
#' @importFrom grid arrow unit
#' @importFrom ggnewscale new_scale_color
#' @importFrom ggnewscale new_scale
#' @importFrom data.table as.data.table
#' @keywords internal

## Plot for Risk region
risky_plot <- function(species, chromosome_vec, front_df1, front_df2, background = mouse_to_human_genes, mode = 'All', way = 'One') {
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
    , status := fifelse(Risk_status == "By Breakpoint" & Group_status == "Group Split", "By Breakpoint and Split",
    fifelse(Risk_status == "By Breakpoint" & Group_status == "Group Together", "By Breakpoint",
    fifelse(Risk_status == 'By Breakpoint' & Group_status == 'Not Applicable', 'By Breakpoint', 
    fifelse(Risk_status == "Not by Breakpoint" & Group_status == "Group Split", "Group Split",
    Group_status))
    ))
  ][
    , status := fifelse(is.na(status), "Unattached Peak", status)
  ][
    cumulative, on = 'chromosome'
  ][
    , mid := (start + end)/2 + position
  ][ 
    , mid_o := (start_o + end_o)/2
  ][ 
    , type := fifelse(grepl("chr", symbol), 0.75, 1)
  ][
    , type_spacing := fifelse(grepl("chr", symbol), 0.9, 1)
  ][ 
    , dotshape := fifelse(type == 1, "Gene", "PeakFir")
  ]

  down_df1 <- front_df2[
    , .(group_id, symbol, chromosome, start, end, chromosome_o, start_o, end_o, In_Risk, Risk_status, Group_status)
  ][ 
    , status := fifelse(Risk_status == "By Breakpoint" & Group_status == "Group Split", "By Breakpoint and Split",
    fifelse(Risk_status == "By Breakpoint" & Group_status == "Group Together", "By Breakpoint",
    fifelse(Risk_status == 'By Breakpoint' & Group_status == 'Not Applicable', 'By Breakpoint', 
    fifelse(Risk_status == "Not by Breakpoint" & Group_status == "Group Split", "Group Split",
    Group_status))
    ))
  ][
    , status := fifelse(is.na(status), "Unattached Peak", status)
  ][
    cumulative, on = 'chromosome'
  ][
    , mid := (start + end)/2 + position
  ][ 
    , mid_o := (start_o + end_o)/2
  ][ 
    , type := fifelse(grepl("chr", symbol), 0.75, 1)
  ][
    , type_spacing := fifelse(grepl("chr", symbol), 0.9, 1)
  ][ 
    , dotshape := fifelse(type == 1, "Gene", "PeakSec")
  ]
  
  fil_back <- back_df[chromosome %in% chromosome_vec & !grepl("chr", symbol)]
  up_df1 <- up_df1[chromosome %in% chromosome_vec]
  down_df1 <- down_df1[chromosome %in% chromosome_vec]

  noto_space <- max(fil_back$mid_o) * 1.3
  yoto_space <- max(fil_back$mid_o) * 0.1

  ## Change when more than one species.
  if (species == "Mmusculus") {
    ortholog_species = "Hsapiens"
  } else {
    ortholog_species = "Mmusculus"
  }
  
  if (mode == 'All') {
    ## For labels of groups by breakpoints
    up_label <- up_df1[
      !is.na(group_id) & !is.na(chromosome_o) & (Risk_status == 'By Breakpoint' | Group_status == 'Group Split') 
    ][
      , overall := mean(mid), by = .(group_id)
    ][
      , overall_o := mean(mid_o), by = .(group_id)
    ][ 
      , .(group_id, chromosome, overall, overall_o)
    ]
    up_label <- unique(up_label)
    down_label <- down_df1[
      !is.na(group_id) & !is.na(chromosome_o) & (Risk_status == 'By Breakpoint' | Group_status == 'Group Split')
    ][
      , overall := mean(mid), by = .(group_id)
    ][
      , overall_o := mean(mid_o), by = .(group_id)
    ][ 
      , .(group_id, chromosome, overall, overall_o)
    ]
    down_label <- unique(down_label)
    ## Risky plot
    if (way == 'One') {
      risk_plot <- ggplot(data = fil_back) +
      geom_segment(
      aes(x = mid, xend = mid, y=0, yend = mid_o, color = chromosome_o), 
      alpha = 1
      ) +
      labs(
      x = paste(species, "Cumulative Gene Midpoint Position from 5prime of Chromosome 1 (bp)"),
      y = paste(ortholog_species, "Midpoint Position (bp)"),
      color = "Ortholog Chromosome"
      ) +
      scale_color_manual(values=matching_colors) +
      new_scale_color() +
      geom_segment(
      data = up_df1[!is.na(chromosome_o)],
      aes(x = mid, xend = mid, yend = mid_o + yoto_space, y = mid_o + 2 * yoto_space * type, color = status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) + # for one to one orthos
      geom_text(data = up_label, aes(x = overall, y = overall_o + 2 * yoto_space, label = group_id)) +
      geom_point(
      data = data.frame(
      x = NA, y = NA,
      flag = c('Has ortholog', 'No ortholog')
      ),
      aes(x = x, y = y, shape = flag),
      size = 4
      ) +
      scale_shape_manual(
      name   = "One to One ortholog",
      values = c("Has ortholog" = 17, "No ortholog" =  2),
      guide = guide_legend(order = 1)
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
      new_scale('shape') +
      geom_point(data = up_df1[is.na(chromosome_o)], aes(x = mid, y = noto_space/type_spacing , color = status, shape = dotshape)) + # for no orthos
      scale_shape_manual(
      name   = "Missing Ortholog",
      values = c("Gene" =  5, "PeakFir" = 6),
      labels = c("Gene" = "Gene", "PeakFir" = "Peak"),
      guide = guide_legend(order = 2)
      ) +
      labs(color = 'Group Status') +
      guides(alpha = "none") +
      theme(text = element_text(size = 14))

    } else {
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
      x = paste(species, "Cumulative Gene Midpoint Position from 5prime of Chromosome 1 (bp)"),
      y = paste(ortholog_species, "Midpoint Position (bp)"),
      color = "Ortholog Chromosome"
      ) +
      scale_color_manual(values=matching_colors) +
      new_scale_color() +
      geom_segment(
      data = up_df1[!is.na(chromosome_o)],
      aes(x = mid, xend = mid, yend = mid_o + yoto_space, y = mid_o + 2 * yoto_space * type, color = status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) + # for one to one orthos
      geom_text(data = up_label, aes(x = overall, y = overall_o + 2 * yoto_space, label = group_id)) +
      geom_segment(
      data = down_df1[!is.na(chromosome_o)], 
      aes(x = mid, xend = mid, yend = -mid_o - yoto_space, y = -mid_o - 2 * yoto_space * type, color = status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) +
      geom_text(data = down_label, aes(x = overall, y = -overall_o - 2 * yoto_space, label = group_id)) +
      geom_point(
      data = data.frame(
      x = NA, y = NA,
      flag = c('Has ortholog', 'No ortholog')
      ),
      aes(x = x, y = y, shape = flag),
      size = 4
      ) +
      scale_shape_manual(
      name   = "One to One ortholog",
      values = c("Has ortholog" = 17, "No ortholog" =  2),
      guide = guide_legend(order = 1)
      ) +
      new_scale("shape") +
      geom_point(data = up_df1[is.na(chromosome_o)], aes(x = mid, y = noto_space/type_spacing , color = status, shape = dotshape)) + # for no orthos
      geom_point(data = down_df1[is.na(chromosome_o)], aes(x = mid, y = -noto_space/type_spacing, color = status, shape = dotshape)) +
      scale_shape_manual(
      name   = "Missing Ortholog",
      values = c("Gene" =  5, "PeakFir" = 6, "PeakSec" = 2),
      labels = c("Gene" = "Gene", "PeakFir" = "Peak List1", "PeakSec" = "Peak List2"),
      guide = guide_legend(order = 2)
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
      labs(color = 'Group Status') +
      guides(alpha = "none") +
      theme(text = element_text(size = 14))
    }
    
  }
  if (mode == 'Shared') {
    # For labels of shared 
    # Do I do this for peaks??
    shared_genes <- unique(c(up_df1[grepl('chr', symbol)]$symbol, down_df1[grepl('chr', symbol)]$symbol))
    up_label <- up_df1[
      !is.na(group_id) & (symbol %in% shared_genes$symbol) 
    ]
    down_label <- down_df1[
      !is.na(group_id) & (symbol %in% shared_genes$symbol)
    ]

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
      x = paste(species, "Gene Midpoint Position (bp)"),
      y = paste(ortholog_species, "Midpoint Position (bp)"),
      color = "Ortholog Chromosome"
    ) +
    scale_color_manual(values=matching_colors) +
    new_scale_color() +
    geom_segment(
      data = up_label[!is.na(chromosome_o)],
      aes(x = mid, xend = mid, yend = mid_o + yoto_space, y = mid_o + 2 * yoto_space * type, color = status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) + 
    geom_text(
      data = up_label[!is.na(chromosome_o)],
      aes(x = mid, xend = mid, yend = mid_o + yoto_space, y = mid_o + 2 * yoto_space * type, color = status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) + # for one to one orthos
    #geom_point(data = up_df1[is.na(chromosome_o)], aes(x = mid, y = noto_space/type_spacing , color = status), shape = 6) + # for no orthos
    #geom_text(data = up_label, aes(x = overall, y = overall_o + 2 * yoto_space, label = group_id)) +
    geom_segment(
      data = down_label[!is.na(chromosome_o)], 
      aes(x = mid, xend = mid, yend = -mid_o - yoto_space, y = -mid_o - 2 * yoto_space * type, color = status),
      arrow = arrow(length = unit(0.1, "cm"), type = "closed"), 
      linewidth = 0.2) +
    #geom_point(data = down_df1[is.na(chromosome_o)], aes(x = mid, y = -noto_space/type_spacing, color = status), shape = 2) +
    #geom_text(data = down_label, aes(x = mid, y = -mid_o - 2 * yoto_space * type, label = group_id)) +
    geom_point(
    data = data.frame(
      x = NA, y = NA,
      flag = c('Has ortholog', 'No ortholog')
    ),
    aes(x = x, y = y, shape = flag),
    size = 4
    ) +
    scale_color_manual(values=c(
      "By Breakpoint and Split" = "#FF5733",
      "By Breakpoint" = "orange", 
      "Group Split" = "darkorchid2",
      'Group Together' = 'blue',
      'Not Applicable' = 'grey'
      ),
      breaks = c("By Breakpoint and Split", "By Breakpoint", "Group Split", 'Group Together', 'Not Applicable')
      ) +
    scale_shape_manual(
      name   = "One to One ortholog",
      values = c("Has ortholog" = 17, "No ortholog" =  2)
    ) +
    labs(color = 'Group Status') +
    guides(alpha = "none")
  }
  return(risk_plot)
}

## Helpers for heatmaps
get_order_pluspeaks <- function(in_obj1, in_obj2, mode = 'grouped') {
    if (mode == 'grouped') {
    gene_order <- unique(rbind(
        in_obj1[!is.na(group_id)][ , .(symbol, chromosome, start)],
        in_obj2[!is.na(group_id)][, .(symbol, chromosome, start)]
    ))
    }
    if (mode == 'all') {
    gene_order <- unique(rbind(
        in_obj1[ , .(symbol, chromosome, start)],
        in_obj2[, .(symbol, chromosome, start)]
    ))
    }
    
    gene_order <- gene_order[ 
        order(chromosome,start)
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
matrices_for_plot <- function(species, chromosome_vec, in_obj1, in_obj2, mode = 'grouped') {
  order_dt <- get_order_pluspeaks(in_obj1[[1]], in_obj2[[1]], mode)[
    , .(symbol,diff,coord)
  ]
  if (mode == 'all') {
    up_filtr <- in_obj1[[2]]
    down_filtr <- in_obj2[[2]]
  }
  if (mode == 'grouped') {
    up_filtr <- in_obj1[[2]][!is.na(group_id) & !is.na(group_id2)]
    down_filtr <- in_obj2[[2]][!is.na(group_id) & !is.na(group_id2)]
  }

  ## Top half of the distances
  in_group_df <- up_filtr[
    chromosome == chromosome2
  ][
    (chromosome %in% chromosome_vec) & (chromosome2 %in% chromosome_vec)
  ][
    , .(symbol, symbol2, symbol_o, distance, grouped, Group_status, one_to_one, group_id, group_id2, chromosome)
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
    , has_ortho := fifelse(is.na(symbol_o) | symbol_o == "NA.NA.NA", FALSE, TRUE)
  ][
    , .(symbol,symbol2,coord,diff,i.coord,i.diff,distance, group_id, group_id2, grouped, Group_status, one_to_one, chromosome, has_ortho)
  ][
    order(coord, i.coord)
  ]
  
  fix <- in_group_df[,.(symbol2, symbol, distance)]
  colnames(fix) <- c('symbol', 'symbol2', 'distance2')
  in_group_df <- fix[
    in_group_df, on = .(symbol,symbol2)
  ][ 
    , distance := fifelse(is.na(distance), distance2, distance)
  ][
    i.coord >= coord
  ][ 
    , distance2 := NULL
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
    group_id = group_id,
    chromosome = chromosome,
    has_ortholog = has_ortho,
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
    group_id = group_id,
    group_id2 = group_id2,
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
    , .(symbol, symbol2, symbol_o, distance, grouped, Group_status, one_to_one, group_id, group_id2, chromosome)
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
    , has_ortho := fifelse(is.na(symbol_o) | symbol_o == "NA.NA.NA", FALSE, TRUE)
  ][
    , .(symbol,symbol2,coord,diff,i.coord,i.diff,distance, group_id, group_id2, grouped, Group_status, one_to_one, chromosome, has_ortho)
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
    group_id = group_id, 
    chromosome = chromosome,
    has_ortholog = has_ortho,
    id = paste(coord, i.coord, "A", sep = "_")
    ), by = .(coord, i.coord)
  ]

  lower_blocks <- in_group_df2[
    symbol != symbol2
  ][, .(
    x = c(coord - 0.5*diff, coord + 0.5*diff, coord + 0.5*diff, coord - 0.5*diff),
    y = c(i.coord - 0.5*i.diff, i.coord - 0.5*i.diff, i.coord + 0.5*i.diff, i.coord + 0.5*i.diff),
    distance = distance,
    group_id = group_id2, ## Note that the group_id here should be group_id2. We want the lower squares.
    group_id2 = group_id,
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
  return(list(in_group_diag, in_group_diag2, upper_blocks, 
  lower_blocks, upper_group, lower_group))
}

## Plot for Distance
heatmap_plot <- function(list_matrices, type, mode = 'interactive', bounds = c(), version = 'Two') {
  x_label <- ifelse(version == 'Two', 'List2', 'Features ordered by position')
  y_label <- ifelse(version == 'Two', 'List1', 'Features ordered by position')
  ortho_label <- ifelse(version == 'Two', '(Diagonal Triangle) Has Ortholog', '(Diagonal Square) Has Ortholog')
  if (mode == 'interactive') {
    if (type == 'distance') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[3]], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[4]], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, color = grouped), fill = NA, linewidth = 0.5) +
      geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, color = grouped), fill = NA, linewidth = 0.5) +
      scale_fill_gradient(name = 'log10 of Distance', low = "navyblue", high = "lightgoldenrodyellow") +
      scale_fill_manual(values = c(
      "Grouped" = "#00BFC4",
      "Singlet group" = "#F8766D"
      )) +
      coord_fixed() +
      theme(legend.position = 'bottom') +
      labs(
        x = x_label,   
        y = y_label
      )
    }
    if (type == 'grouped') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = grouped, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = grouped, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5, linewidth = 0.5) +
      geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5, linewidth = 0.5) +
      coord_fixed() +
      scale_fill_manual(values = c(
      "Grouped" = "#00BFC4",
      "Singlet group" = "#F8766D"
      )) +
      theme(legend.position = 'bottom') +
      labs(
        x = x_label,   
        y = y_label 
      )
    }
    if (type == 'sequence') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = one_to_one, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = one_to_one, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, fill = one_to_one), alpha = 0.5) +
      geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, fill = one_to_one), alpha = 0.5) +
      coord_fixed() +
      scale_fill_manual(values = c(
      "No Loss" = "#00BFC4",
      "Partial Loss" = "#1b9e77",
      "Lost" = "#7570b3",
      "Became Singlet" = "brown",
      "Singlet group" = "#F8766D"
      )) +
      theme(legend.position = 'bottom') +
      labs(
        x = x_label,   
        y = y_label 
      )
    }
    if (type == 'location') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]], aes(x = x ,y = y, group = id, fill = Group_status, alpha = has_ortholog), color = "black", linewidth = 0.1) +
      geom_polygon(data = list_matrices[[2]], aes(x = x, y = y, group = id, fill = Group_status, alpha = has_ortholog), color = "black", linewidth = 0.1) +
      geom_polygon(data = list_matrices[[5]], aes(x = x, y = y, group = id, fill = Group_status), alpha = 0.5) +
      geom_polygon(data = list_matrices[[6]], aes(x = x, y = y, group = id, fill = Group_status), alpha = 0.5) +
      coord_fixed() +
      scale_fill_manual(values = c(
      "Group Together" = "#00BFC4",
      "Group Split" = "#1b9e77",
      "Not Applicable" = "#7570b3",
      "Became Singlet" = "brown",
      "Singlet group" = "#F8766D"
      )) +
      theme(legend.position = 'bottom') +
      labs(
        x = x_label,   
        y = y_label,
        alpha = "Has Ortholog" 
      )
    }
    return(plot)
  }
  if (mode == 'shiny') {
    # return(list(in_group_diag, in_group_diag2, upper_blocks, lower_blocks, upper_group, lower_group, in_group_df, in_group_df2))
    # return(c(new_group1_min, new_group1_max, new_group2_min, new_group2_max))
    side1 <- list_matrices[[1]][
      group_id >= bounds[1] & group_id <= bounds[2]
    ]
    side2 <- list_matrices[[2]][
      group_id >= bounds[3] & group_id <= bounds[4]
    ]
    both_sides <- rbind(side1,side2)[, .(x, y, chromosome)][
      , lower_x := min(x), by = .(chromosome)
    ][
      , upper_x := max(x), by = .(chromosome)
    ]
    chrom_outline <- unique(both_sides)[
    , .(
      x = c(lower_x, lower_x, upper_x, upper_x),
      y = c(lower_x, upper_x, upper_x, lower_x),
      chromosome = chromosome,
      by = .(chromosome)
      )
    ]
    if (type == 'distance') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]][
          group_id >= bounds[1] & group_id <= bounds[2]
        ], aes(x = x ,y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[2]][
          group_id >= bounds[3] & group_id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[3]][
          group_id >= bounds[1] & group_id <= bounds[2] & group_id2 >= bounds[1] & group_id2 <= bounds[2]
        ], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[4]][
          group_id >= bounds[3] & group_id <= bounds[4] & group_id2 >= bounds[3] & group_id2 <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1))) +
      geom_polygon(data = list_matrices[[5]][
          id >= bounds[1] & id <= bounds[2]
        ], aes(x = x, y = y, group = id, color = grouped), fill = NA, linewidth = 0.5) +
      geom_polygon(data = list_matrices[[6]][
          id >= bounds[3] & id <= bounds[4]
        ], aes(x = x, y = y, group = id, color = grouped), fill = NA, linewidth = 0.5) +
      geom_polygon(data = chrom_outline, aes(x = x, y = y, group = chromosome), color = 'black', fill = NA, linewidth = 0.5) +
      scale_fill_gradient(name = 'log10 of Distance', low = "navyblue", high = "lightgoldenrodyellow") +
      scale_color_manual(name = 'Group Status', 
        values = c('Singlet group' = '#F8766D', 'Grouped' = '#00BFC4'),
        labels = c('Grouped' = 'Multi-gene group', 'Singlet group' = 'Singlet group')
        ) +
      coord_fixed() +
      theme(legend.position = 'bottom', legend.box="vertical", text = element_text(size = 12)) +
      labs(
        x = x_label,   
        y = y_label,
        fill = 'Log10 of distance',
        color = 'Group Status'
      ) +
      guides(color=guide_legend(nrow=2,byrow=TRUE)) 
    }
    if (type == 'grouped') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]][
          group_id >= bounds[1] & group_id <= bounds[2]
        ], aes(x = x ,y = y, group = id, fill = grouped, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[2]][
          group_id >= bounds[3] & group_id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = grouped, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[5]][
          id >= bounds[1] & id <= bounds[2]
        ], aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5, linewidth = 0.1) +
      geom_polygon(data = list_matrices[[6]][
          id >= bounds[3] & id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = grouped), alpha = 0.5, linewidth = 0.1) +
      geom_polygon(data = chrom_outline, aes(x = x, y = y, group = chromosome), color = 'black', fill = NA, linewidth = 0.5) +
      scale_fill_manual(name = 'Group Status', 
      values = c('Singlet group' = '#F8766D', 'Grouped' = '#00BFC4'),
      labels = c('Grouped' = 'Multi-gene group', 'Singlet group' = 'Singlet group')) +
      scale_alpha_manual(name = ortho_label, 
        values = c('FALSE' = 0.1, 'TRUE' = 1),
        labels = c('FALSE' = 'NO', 'TRUE' = 'YES')
      ) +
      coord_fixed() +
      theme(legend.position = 'bottom', legend.box="vertical", text = element_text(size = 12)) +
      labs(
        x = x_label,   
        y = y_label
      ) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
    }
    if (type == 'sequence') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]][
          group_id >= bounds[1] & group_id <= bounds[2]
        ], aes(x = x ,y = y, group = id, fill = one_to_one, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[2]][
          group_id >= bounds[3] & group_id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = one_to_one, alpha = has_ortholog), color = 'black', linewidth = 0.1) +
      geom_polygon(data = list_matrices[[5]][
          id >= bounds[1] & id <= bounds[2]
        ], aes(x = x, y = y, group = id, fill = one_to_one), alpha = 0.5) +
      geom_polygon(data = list_matrices[[6]][
          id >= bounds[3] & id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = one_to_one), alpha = 0.5) +
      geom_polygon(data = chrom_outline, aes(x = x, y = y, group = chromosome), color = 'black', fill = NA, linewidth = 0.5) +
      scale_alpha_manual(name = ortho_label, 
        values = c('FALSE' = 0.1, 'TRUE' = 1),
        labels = c('FALSE' = 'NO', 'TRUE' = 'YES')
      ) +
      coord_fixed() +
      scale_fill_manual(values = c(
      "No Loss" = "#00BFC4",
      "Partial Loss" = "#1b9e77",
      "Lost" = "#7570b3",
      "Became Singlet" = "brown",
      "Singlet group" = "#F8766D"
      )) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
      theme(legend.position = 'bottom', legend.box="vertical", text = element_text(size = 12)) +
      labs(
        x = x_label,   
        y = y_label,
        fill = 'Sequence Status'
      )
    }
    if (type == 'location') {
      plot <- ggplot() +
      geom_polygon(data = list_matrices[[1]][
          group_id >= bounds[1] & group_id <= bounds[2]
        ], aes(x = x ,y = y, group = id, fill = Group_status, alpha = has_ortholog), color = "black", linewidth = 0.1) +
      geom_polygon(data = list_matrices[[2]][
          group_id >= bounds[3] & group_id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = Group_status, alpha = has_ortholog), color = "black", linewidth = 0.1) +
      geom_polygon(data = list_matrices[[5]][
          id >= bounds[1] & id <= bounds[2]
        ], aes(x = x, y = y, group = id, fill = Group_status), alpha = 0.5) +
      geom_polygon(data = list_matrices[[6]][
          id >= bounds[3] & id <= bounds[4]
        ], aes(x = x, y = y, group = id, fill = Group_status), alpha = 0.5) +
      geom_polygon(data = chrom_outline, aes(x = x, y = y, group = chromosome), color = 'black', fill = NA, linewidth = 0.5) +
      scale_alpha_manual(name = ortho_label, 
        values = c('FALSE' = 0.1, 'TRUE' = 1),
        labels = c('FALSE' = 'NO', 'TRUE' = 'YES')
      ) +
      coord_fixed() +
      scale_fill_manual(values = c(
      "Group Together" = "#00BFC4",
      "Group Split" = "#1b9e77",
      "Not Applicable" = "#7570b3",
      "Became Singlet" = "brown",
      "Singlet group" = "#F8766D"
      )) +
      theme(legend.position = 'bottom', legend.box="vertical", text = element_text(size = 12)) +
      labs(
        x = x_label,   
        y = y_label,
        fill = 'Location Status'
      ) +
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
    }
    return(plot)
  }
}

plot_group_stats <- function(input_table, species, type, label, in_chr, mode = 'Group') {
    #Group_breaks <- c(
    #'Became Singlet',
    #'Lost',
    #'No Loss',
    #'Partial Loss',
    #'Singlet')
    if (mode == 'Group') {
        summary_table <- unique(input_table[
            !grepl("chr", symbol) & chromosome %in% in_chr
        ][ 
            , .(group_id, one_to_one, Risk_status, Group_status)
        ])
    }
    if (mode == 'Gene') {
        summary_table <- unique(input_table[
            !grepl("chr", symbol) & chromosome %in% in_chr
        ][ 
            , .(symbol, one_to_one, Risk_status, Group_status)
        ])
    }
    summary_table <- unique(summary_table[
        , Group_count := .N, by = .(one_to_one)
    ][ 
        , Risk_count := .N, by = .(one_to_one, Risk_status)
    ][
        , Split_count := .N, by = .(one_to_one, Group_status)
    ][
        order(one_to_one)
    ])
    if (type == 'Risk') {
        summary_table <- unique(summary_table[
            , .(one_to_one, Risk_status, Group_count, Risk_count)
        ])
        summary_table[ 
            , one_to_one_x := as.numeric(factor(summary_table$one_to_one))
        ]
        summary_table <- summary_table[
          order(one_to_one_x, Risk_status)
        ]
        summary_table[
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
        ][ 
          , sep := fifelse(Status == 'By Breakpoint', 0.25, -0.25)
        ]
        phrase <- paste0('\n',mode,' by Breakpoint')
        color_vec <- c(
          "By Breakpoint" = "#F8766D",
          "Not by Breakpoint" = "#00BFC4")
    } 
    if (type == 'Split') {
        summary_table <- unique(summary_table <- summary_table[
            , .(one_to_one, Group_status, Group_count, Split_count)
        ])
        summary_table[ 
            , one_to_one_x := as.numeric(factor(summary_table$one_to_one))
        ]
        summary_table <- summary_table[
          order(one_to_one_x, Group_status)
        ]
        summary_table[
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
        ][ 
          , sep := fifelse(Status == 'Not Applicable' , 0,fifelse(Status == 'Group Split', 0.25, -0.25))
        ]
        phrase <- paste0('\n',mode,' by Split')
        color_vec <- c(
          "Group Together" = "#00BFC4",
          "Group Split" = "#1b9e77",
          "Not Applicable" = "#7570b3")
    }
    buffer <- min(summary_table$Group_count) * 0.33
    summary_table$buffer <- buffer
    Group_breaks <- unique(summary_table$one_to_one)
    plot <- ggplot()  +
    geom_rect(
        data = summary_table,
        aes(xmin = one_to_one_x - 0.375, xmax = one_to_one_x + 0.375, ymin = y_base, ymax = y_cumulative, fill = Status)
    ) + 
    geom_text(
        data = unique(summary_table[,.(one_to_one_x, Group_count, buffer)]),
        aes(
            x = one_to_one_x, 
            y = Group_count * 1.1 + buffer, 
            label = paste("Total:", Group_count, sep = ' ')
        )
    ) +
    geom_text(
        data = summary_table,
        aes(
            x = one_to_one_x + sep,
            y = y_cumulative,
           label = ifelse(how_many == 1, NA, Count)
        )
    ) +
    scale_x_continuous(
        breaks = rep(1:max(summary_table$one_to_one_x) ),
        labels = Group_breaks
    ) +
    labs(
        x = "Status",   
        y = "Composition" 
    ) +
    scale_fill_manual(name = 'Status', 
      values = color_vec
    ) +
    ggtitle(paste(species, label, phrase, sep = ' ')) +
    theme(text = element_text(size = 14))
    return(plot)
}

plot_attached_peaks <- function(input_table, chromosome_vec) {
  plot_dat <- input_table[ 
    chromosome %in% chromosome_vec
  ][
    grepl('chr', symbol)
  ][ 
    , .(symbol, group_id)
  ]
  not_attached <- length(plot_dat[is.na(group_id)]$symbol)
  attached <- length(plot_dat[!is.na(group_id)]$symbol)
  plot_dat2 <- data.table(
    Status = c('Attached to Group', 'Not Attached to Group'),
    Count = c(attached, not_attached)
  )
  plot <- ggplot(plot_dat2, aes(x = Status, fill = Status)) +
    geom_bar(aes(y = Count), stat = 'identity') +
    geom_text(aes(y = Count * 1.1, label = paste("Total:", Count))) +
    scale_fill_manual(values = c('Attached to Group' = '#00BFC4', 'Not Attached to Group' = '#F8766D')) +
    labs(
      x = "Status",   
      y = "Count" 
    ) +
    theme(text = element_text(size = 14)) +
    ggtitle("Peaks Attached to Groups")
  return(plot)

}

ortho_summary_prop <- function(ortho_table, type) {
  if (type == 'Peaks') {
    summary_table <- ortho_table[ 
      , total := .N
    ][ 
      , Has_Ortholog_Num := .N, by = .(Has_Ortholog)
    ][ 
      , Has_Ortholog_Proportion := Has_Ortholog_Num / total
    ]
    summary_table <- unique(summary_table[, .(Has_Ortholog, Has_Ortholog_Proportion)])
    plot <- ggplot(summary_table, aes(x = Has_Ortholog, fill = Has_Ortholog)) +
      geom_bar(aes(y = Has_Ortholog_Proportion), stat = 'identity') +
      geom_text(aes(y = Has_Ortholog_Proportion + 0.1, label = round(Has_Ortholog_Proportion, 2))) +
      scale_fill_manual(values = c('Yes' = '#00BFC4', 'No' = '#F8766D')) +
      labs(
        x = "Has Ortholog",   
        y = "Proportion" 
      ) +
      theme(legend.position = 'none', text = element_text(size = 14)) +
      ggtitle(paste(type, "Ortholog Proportions", sep = ' '))
    return(plot)
  } else if (type == 'Genes') {
    summary_table <- ortho_table[ 
      , total := .N
    ][ 
      , Has_Ortholog_Num := .N, by = .(Has_Ortholog)
    ][ 
      , Has_Ortholog_Proportion := Has_Ortholog_Num / total
    ][ 
      , In_Annotation_Num := .N, by = .(In_Annotation)
    ][ 
      , In_Annotation_Proportion := In_Annotation_Num / total
    ]
    summary_table <- unique(summary_table[, .(Has_Ortholog, Has_Ortholog_Proportion, In_Annotation, In_Annotation_Proportion)])
    summary_ortho <- unique(summary_table[,.(Has_Ortholog, Has_Ortholog_Proportion)][,Type := 'Has_Ortholog'])
    summary_annot <- unique(summary_table[,.(In_Annotation, In_Annotation_Proportion)][,Type := 'In_Annotation'])
    colnames(summary_ortho) <- c('Status', 'Proportion', 'Type')
    colnames(summary_annot) <- c('Status', 'Proportion', 'Type')
    melted_summary <-  rbind(summary_ortho,summary_annot)[
      , Status := factor(Status, levels = c('Yes', 'No'))
    ][ 
      order(Type,Status)
    ][ 
      , cumulative_prop := cumsum(Proportion), by = .(Type)
    ]
    plot <- ggplot(melted_summary, aes(x = Type, fill = forcats::fct_rev(Status))) +
      geom_bar(aes(y = Proportion), stat = 'identity', position = 'dodge') +
      geom_text(aes(y = 1, label = round(Proportion, 2)), position = position_dodge(width = 0.9) ) +
      scale_fill_manual(name = 'Status', values = c('Yes' = '#00BFC4', 'No' = '#F8766D')) +
      labs(
        x = "Status Type",   
        y = "Proportion" 
      ) +
      theme(
        legend.position = 'bottom', 
        text = element_text(size = 14),
        fill = 'Status'
        ) +
      ggtitle(paste(type, "Ortholog Proportions", sep = ' ')) +
      labs(x = 'Status')
  }
  return(plot)
}

singlet_summary <- function(input_table, chromosome_vec) {
  plot_dat <- unique(input_table[ 
    chromosome %in% chromosome_vec
  ][
    !is.na(group_id)
  ][
    , .(group_id, in_group)
  ])
  plot_dat <- plot_dat[
    , total_count := .N
  ][
    , group_count := .N, by = .(in_group)
  ][
    , group_proportion := group_count/total_count
  ]

  plot <- ggplot(plot_dat) +
  geom_bar(aes(x = in_group, y = group_proportion, fill = in_group), stat = 'identity', position = 'dodge') +
  geom_text(data = unique(plot_dat[,.(in_group, group_proportion, group_count)]), aes(x = in_group, y = group_proportion + 0.1, label = paste("Total", group_count, sep = ' ')), position = position_dodge(width = 0.9) ) +
  scale_fill_manual(
    name = 'Grouped Status', 
    values = c('TRUE' = '#00BFC4', 'FALSE' = '#F8766D'),
    labels = c('TRUE' = 'Multi-gene group', 'FALSE' = 'Singlet group')
    ) +
  scale_x_discrete(
    "Grouped Status",
    labels = c(
      "TRUE" = "Multi-gene group",
      "FALSE" = "Singlet group"
    )
  ) +
  labs(x = "Grouped", y = "Proportion") +
  theme(text = element_text(size = 14))
  return(plot)
  
}

group_size_plot <- function(input_table, chromosome_vec, stat_type) {
  plot_dat <- input_table[
    chromosome %in% chromosome_vec
  ][
    !is.na(group_id)
  ][
    , .(symbol, group_id, in_group, start, end)
  ]
  if (stat_type == 'Groups') {
    plot_dat <- plot_dat[
    , group_start := min(start), by = .(group_id)
    ][
    , group_end := max(end), by = .(group_id)
    ][
    , length := log10(group_end - group_start)
    ][
    , log_medianlength := median(length), by = .(in_group)
    ][
    ,.(in_group, log_medianlength, length)
    ]
  }
  if (stat_type == 'Genes') {
    plot_dat <- plot_dat[
    , length := log10(end - start)
    ][
    , log_medianlength := median(length), by = .(in_group)
    ][
    ,.(in_group, log_medianlength, length)
    ]
  }
  plot <- ggplot(plot_dat, aes(x = in_group)) +
  geom_violin(
    aes(y = length, fill = in_group)
    ) +
  geom_point(
    data = unique(plot_dat[,.(in_group, log_medianlength)]), 
    aes(y = log_medianlength), color = 'black'
    ) +
  expand_limits(y = 0) +
  scale_fill_manual(
    name = 'Grouped Status', 
    values = c('TRUE' = '#00BFC4', 'FALSE' = '#F8766D'),
    labels = c('TRUE' = 'Multi-gene group', 'FALSE' = 'Singlet group')
    ) +
  scale_x_discrete(
    "Grouped Status",
    labels = c(
      "TRUE" = "Multi-gene group",
      "FALSE" = "Singlet group"
    )
  ) +
  labs(y = "Log of Length") +
  theme(text = element_text(size = 14))
  return(plot)
}