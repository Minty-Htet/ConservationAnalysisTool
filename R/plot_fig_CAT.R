#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_alpha_manual
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom grid arrow unit
#' @importFrom ggnewscale new_scale_color
#' @importFrom data.table as.data.table
#' @keywords internal

## Helpers
get_order <- function(cat_obj1, cat_obj2) {
  gene_order <- unique(rbind(
    cat_obj1[[1]][ , .(symbol, chromosome, start)],
    cat_obj2[[1]][, .(symbol, chromosome, start)]
  ))[
    !grepl('chr', symbol)
  ][ 
    order(start)
  ]
  gene_order <- unique(
    gene_order[ 
      , symbol := factor(symbol, levels = unique(symbol))
    ][ 
      , coord := as.numeric(symbol)
    ][ 
      , .(symbol, coord)
    ])
  return(gene_order)
}

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

## Plot for Group
## Prepare the heatmap matrices
heatmap_for_plot <- function(species, chromosome_vec, cat_obj1, cat_obj2) {
  order_dt <- get_order(cat_obj1, cat_obj2)
  ## Top half of the distances
  in_group_df <- cat_obj1[[2]][
    !grepl('chr', symbol) & !grepl('chr', symbol2)
  ][
    (chromosome %in% chromosome_vec) & (chromosome2 %in% chromosome_vec)
  ][
    , .(symbol, symbol2, symbol_o, distance, grouped, Group_status, one_to_one, group_id, group_id2)
  ][ 
    order_dt, on = .(symbol)
  ][ 
    order_dt, on = list(symbol2 = symbol)
  ][ 
    , distance := abs(distance)
  ][ 
    !is.na(grouped) | symbol == symbol2 
  ][
     i.coord >= coord
  ]
  ## For the diagnoal triangles
  in_group_diag <- in_group_df[
    symbol == symbol2
  ][
    , coord := coord - 0.5
  ][ 
    , i.coord := i.coord - 0.5
  ][, .(
    x = c(coord + 0.5, coord - 0.5, coord - 0.5),
    y = c(i.coord + 0.5, i.coord + 0.5, i.coord - 0.5),
    grouped = grouped,
    Group_status = Group_status,
    one_to_one = one_to_one,
    half = "lower_left",
    id = paste(coord, i.coord, "A", sep = "_"),
    group_id = group_id,
    group_id2 = group_id2,
    distance = distance,
    symbol_o = symbol_o
    ), by = .(coord, i.coord)
  ]
  ## Adjust the positioning so that the corner starts and 0,0
  in_group_df <- in_group_df[
    symbol != symbol2
  ][
    !is.na(coord) & !is.na(i.coord)
  ][
    , coord := coord - 0.5
  ][ 
    , i.coord := i.coord - 0.5
  ]
  ## Bottom half of the distances
  in_group_df2 <- cat_obj2[[2]][
    !grepl('chr', symbol) & !grepl('chr', symbol2)
  ][
    (chromosome %in% chromosome_vec) & (chromosome2 %in% chromosome_vec)
  ][
    , .(symbol, symbol2, symbol_o, distance, grouped, Group_status, one_to_one, group_id, group_id2)
  ][ 
    order_dt, on = .(symbol)
  ][ 
    order_dt, on = list(symbol2 = symbol)
  ][ 
    , distance := abs(distance)
  ][ 
    !is.na(grouped) | symbol == symbol2 
  ][
     i.coord <= coord
  ]
  ## Bottom Diagonals
  in_group_diag2 <- in_group_df2[
    symbol == symbol2
  ][
    , coord := coord - 0.5
  ][ 
    , i.coord := i.coord - 0.5
  ][, .(
    x = c(coord + 0.5, coord + 0.5, coord - 0.5),
    y = c(i.coord - 0.5, i.coord + 0.5, i.coord - 0.5),
    grouped = grouped,
    Group_status = Group_status,
    one_to_one = one_to_one,
    half = "lower_left",
    id = paste(coord, i.coord, "A", sep = "_"),
    group_id = group_id,
    group_id2 = group_id2,
    distance = distance,
    symbol_o = symbol_o
  ), by = .(coord, i.coord)
  ]
  ## Adjust positioning
  in_group_df2 <- in_group_df2[
    symbol != symbol2
  ][
    !is.na(coord) & !is.na(i.coord)
  ][
    , coord := coord - 0.5
  ][ 
    , i.coord := i.coord - 0.5
  ]
  ## Top part of the Group Boundary
  group_up_df <- rbind(
    in_group_df[ 
        , x := coord
    ][ 
        , y := i.coord
    ][ 
        , .(x, y, grouped, Group_status, one_to_one, group_id, group_id2)
    ],
    in_group_diag[
        , .(x ,y, grouped, Group_status, one_to_one, group_id, group_id2)
    ]
  )[ 
    order(group_id, group_id2, x,y)
  ][ 
    grouped != FALSE
  ][
    , min_x := min(x), by = group_id
  ][ 
    , max_x := max(x), by = group_id
  ][ 
    , min_y := min(y), by = group_id
  ][ 
    , max_y := max(y), by = group_id
  ]
  ## For the larger triangles
  group_diag <- unique(group_up_df[
    , .(grouped, Group_status, one_to_one, group_id, group_id2, min_x, max_x, min_y, max_y)
  ])[
    , .(
        x = c(min_x, min_x, max_x),
        y = c(min_y, max_y, max_y),
        id = paste(group_id,group_id2, sep = "_")
    ), by = .(group_id, group_id2, grouped, Group_status, one_to_one)
  ]
  ## Bottom part of the Group Boundary
  group_down_df <- rbind(
    in_group_df2[ 
        , x := coord
    ][ 
        , y := i.coord
    ][ 
        , .(x, y, grouped, Group_status, one_to_one, group_id, group_id2)
    ],
    in_group_diag2[
        , .(x ,y, grouped, Group_status, one_to_one, group_id, group_id2)
    ]
  )[ 
    order(group_id, group_id2, x,y)
  ][ 
    grouped != FALSE
  ][
    , min_x := min(x), by = group_id
  ][ 
    , max_x := max(x), by = group_id
  ][ 
    , min_y := min(y), by = group_id
  ][ 
    , max_y := max(y), by = group_id
  ]
  group_diag2 <- unique(group_down_df[
    , .(grouped, Group_status, one_to_one, group_id, group_id2, min_x, max_x, min_y, max_y)
  ])[
    , .(
        x = c(min_x, max_x, max_x),
        y = c(min_y, min_y, max_y),
        id = paste(group_id,group_id2, sep = "_")
    ), by = .(group_id, group_id2, grouped, Group_status, one_to_one)
  ]
  return(list(in_group_df, in_group_diag, group_up_df, group_diag, in_group_df2, in_group_diag2, group_down_df, group_diag2))
}

## Plot for Distance
distance_plot <- function(species, chromosome_vec, cat_obj1, cat_obj2) {
  list_matrices <- heatmap_for_plot(species, chromosome_vec, cat_obj1, cat_obj2)
  in_group_df = list_matrices[[1]]
  in_group_diag = list_matrices[[2]]
  group_up_df = list_matrices[[3]]
  group_diag = list_matrices[[4]]
  in_group_df2 = list_matrices[[5]]
  in_group_diag2 = list_matrices[[6]]
  group_up_df = list_matrices[[7]]
  group_diag2 = list_matrices[[8]]
  plot <- ggplot() +
  geom_tile(data = in_group_df, aes(x = coord, y = i.coord, fill = log10(abs(distance) + 1)), width = 1, height = 1 ) +
  geom_tile(data = in_group_df2, aes(x = coord, y = i.coord, fill = log10(abs(distance) + 1)), width = 1, height = 1) +
  geom_polygon(data = in_group_diag, aes(x = x ,y = y, group = id, fill = log10(abs(distance) + 1))) +
  geom_polygon(data = in_group_diag2, aes(x = x, y = y, group = id, fill = log10(abs(distance) + 1) )) +
  geom_polygon(data = group_diag, aes(x = x, y = y, color = grouped, group = id), linewidth = 0.5, fill = NA) +
  geom_polygon(data = group_diag2, aes(x = x, y = y, color = grouped, group = id), linewidth = 0.5, fill = NA) +
  coord_fixed() +
  scale_fill_gradient(low = "navyblue", high = "lightgoldenrodyellow") +
  theme(legend.position = 'bottom') +
  labs(
    x = "Gene List 2",
    y = "Gene List 1",
    color = 'Grouping'
  )
  return(plot)
}

## Plot for Group
group_plot <- function(species, chromosome_vec, cat_obj1, cat_obj2) {
  list_matrices <- heatmap_for_plot(species, chromosome_vec, cat_obj1, cat_obj2)
  in_group_df = list_matrices[[1]]
  in_group_diag = list_matrices[[2]]
  group_up_df = list_matrices[[3]]
  group_diag = list_matrices[[4]]
  in_group_df2 = list_matrices[[5]]
  in_group_diag2 = list_matrices[[6]]
  group_up_df = list_matrices[[7]]
  group_diag2 = list_matrices[[8]]
  plot <- ggplot() +
  geom_polygon(data = group_diag, aes(x = x, y = y, fill = grouped, group = id), color = 'black', linewidth = 0.25) +
  geom_polygon(data = group_diag2, aes(x = x, y = y, fill = grouped, group = id), color = 'black', linewidth = 0.25) +
  coord_fixed() +
  theme(legend.position = 'bottom') +
  labs(
    x = "Gene List 2",
    y = "Gene List 1",
    fill = "Grouping"
  )
  return(plot)
}

## Plot for Ortholog Conservation
ortho_plot <- function(species, chromosome_vec, cat_obj1, cat_obj2) {
  list_matrices <- heatmap_for_plot(species, chromosome_vec, cat_obj1, cat_obj2)
  in_group_df = list_matrices[[1]]
  in_group_diag = list_matrices[[2]][
    , Has_ortho := fifelse(is.na(symbol_o), FALSE, TRUE)
  ]
  group_up_df = list_matrices[[3]]
  group_diag = list_matrices[[4]]
  in_group_df2 = list_matrices[[5]]
  in_group_diag2 = list_matrices[[6]][
    , Has_ortho := fifelse(is.na(symbol_o), FALSE, TRUE)
  ]
  group_up_df = list_matrices[[7]]
  group_diag2 = list_matrices[[8]]
  plot <- ggplot() +
  geom_polygon(data = group_diag, aes(x = x, y = y, fill = one_to_one, group = id), alpha = 0.2, color = 'black', linewidth = 0.25) +
  geom_polygon(data = group_diag2, aes(x = x, y = y, fill = one_to_one, group = id), alpha = 0.2, color = 'black', linewidth = 0.25) +
  geom_polygon(data = in_group_diag, aes(x = x ,y = y, group = id, fill = one_to_one, alpha = Has_ortho), color = 'black', linewidth = 0.1) +
  geom_polygon(data = in_group_diag2, aes(x = x, y = y, group = id, fill = one_to_one, alpha = Has_ortho), color = 'black', linewidth = 0.1) +
  scale_alpha_manual(values = c("FALSE" = 0.2, "TRUE" = 1)) +
  coord_fixed() +
  theme(legend.position = 'bottom') +
  labs(
    x = "Gene List 2",
    y = "Gene List 1",
    fill = "Grouping",
    alpha = 'Has ortho'
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    alpha = guide_legend(nrow = 2, byrow = TRUE)
  )
  return(plot)
}

## Plot for location Conservation
location_plot <- function(species, chromosome_vec, cat_obj1, cat_obj2) {
  list_matrices <- heatmap_for_plot(species, chromosome_vec, cat_obj1, cat_obj2)
  in_group_df = list_matrices[[1]]
  in_group_diag = list_matrices[[2]][
    , Has_ortho := fifelse(is.na(symbol_o), FALSE, TRUE)
  ]
  group_up_df = list_matrices[[3]]
  group_diag = list_matrices[[4]]
  in_group_df2 = list_matrices[[5]]
  in_group_diag2 = list_matrices[[6]][
    , Has_ortho := fifelse(is.na(symbol_o), FALSE, TRUE)
  ]
  group_up_df = list_matrices[[7]]
  group_diag2 = list_matrices[[8]]
  plot <- ggplot() +
  geom_polygon(data = group_diag, aes(x = x, y = y, fill = Group_status, group = id), alpha = 0.2, color = 'black', linewidth = 0.25) +
  geom_polygon(data = group_diag2, aes(x = x, y = y, fill = Group_status, group = id), alpha = 0.2, color = 'black', linewidth = 0.25) +
  geom_polygon(data = in_group_diag, aes(x = x ,y = y, group = id, fill = Group_status, alpha = Has_ortho), color = 'black', linewidth = 0.1) +
  geom_polygon(data = in_group_diag2, aes(x = x, y = y, group = id, fill = Group_status, alpha = Has_ortho), color = 'black', linewidth = 0.1) +
  scale_alpha_manual(values = c("FALSE" = 0.2, "TRUE" = 1)) +
  coord_fixed() +
  theme(legend.position = 'bottom') +
  labs(
    x = "Gene List 2",
    y = "Gene List 1",
    fill = "Grouping",
    alpha = 'Has ortho'
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    alpha = guide_legend(nrow = 2, byrow = TRUE)
  )
  return(plot)
}