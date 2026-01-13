#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table fifelse
#' @importFrom data.table setkey
#' @importFrom data.table melt
#' @importFrom data.table copy
#' @importFrom stringr str_replace
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components
#' @keywords internal

assign_groupid <- function(distance_matrix, exp_rocr = 500000) {
  distance_matrix <- distance_matrix[ 
    , in_group := fifelse(abs(distance) < exp_rocr & !symbol == symbol2 & chromosome == chromosome2, TRUE, #Must be within 500kb, and same chromosome
    fifelse(abs(distance) > exp_rocr & !symbol == symbol2 & chromosome == chromosome2, FALSE, 
    fifelse(symbol == symbol2, NA, NA)))
    ][
        , chromosome := factor(chromosome, levels = hs_chr_order)
    ][
        , chromosome2 := factor(chromosome2, levels = hs_chr_order)
    ][ 
    order(chromosome, start, chromosome2, start2)
    ]
  all_genes <- unique(distance_matrix$id)
  matrix_half <- distance_matrix[in_group == TRUE & id != id2]
  # Groups are based on chains of genes within 500kb
  graph <- igraph::graph_from_data_frame(
    d = matrix_half,
    directed = FALSE,
    vertices = all_genes
    )
  components <- igraph::components(graph)
  groups <- data.table(
    id = as.numeric(names(components$membership)),
    group_id = components$membership
    )
  distance_matrix <- distance_matrix[
    groups, on = "id"
    ]
  colnames(groups) <- str_replace(colnames(groups), "$", "2")
  distance_matrix <- distance_matrix[groups, on = "id2"]
  return(distance_matrix)
}