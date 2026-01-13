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

#' Create output plots after doing run_CAT()
#' @param chromosome_vec chromosomes of interest to view in plots. Default is all chromosomes.
#' @param back_df Dataframe with the start and end positions of known coding genes with one to one ortholog
#' @param front_df Dataframe with the start and end positions of genes/regions of interest
#' @return A manhattan plot comparing the conservation patterns of 2 sets of gene inputs
#' @keywords internal
## Run all the plots
risky_plot <- function(front_df1, front_df2, chromosome_vec = "All", back_df = background) {
  fil_back <- back_df[chromosome %in% chromosome_vec & !grepl("chr", symbol)]
  risk_plot <- ggplot(data = back_df) +
    geom_segment(
      aes(x = mid, xend = mid, y=0, yend = start_o, color = chromosome_o), 
      alpha = 1
      ) +
    geom_segment(
      aes(x = mid, xend = mid, y=0, yend = -start_o, color = chromosome_o), 
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