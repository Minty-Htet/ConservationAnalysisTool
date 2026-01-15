#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom grid arrow unit
#' @importFrom ggnewscale new_scale_color
#' @keywords internal
## Run all the plots
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
  fil_back <- back_df[chromosome %in% chromosome_vec & !grepl("chr", symbol)]
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
      x = "Gene Midpoint Position",
      y = "Ortholog Midpoint Position",
      color = "Ortholog Chromosome"
    ) +
    scale_color_manual(values=matching_colors)
  return(risk_plot)
}