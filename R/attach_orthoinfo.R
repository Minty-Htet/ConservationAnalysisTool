#' @importFrom data.table as.data.table
#' @importFrom data.table fifelse
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom data.table copy
#' @importFrom stringr str_replace
#' @keywords internal

attach_orthoinfo <- function(distance_matrix, peaks, species) {
  ## Get the ortholog info of known genes
  if (species == "Mmusculus") {
    cols <- c("mmusculus_gene_name", "hsapiens_gene_name")
    ortho1 <- copy(mouse_to_human_genes)[, .SD, .SDcols = cols]
    setnames(ortho1, cols, c("symbol", "symbol_o"))
    to_join <- human_all_startend[, .(symbol, chromosome, start, end)]
    colnames(to_join) <- paste(colnames(to_join), "o", sep = "_")
    ortho1 <- to_join[
      ortho1, on = .(symbol_o)
    ]
  } else {
    cols <- c("hsapiens_gene_name", "mmusculus_gene_name")
    ortho1 <- copy(mouse_to_human_genes)[, .SD, .SDcols = cols]
    setnames(ortho1, cols, c("symbol", "symbol_o"))
    to_join <- mouse_all_startend[, .(symbol, chromosome, start, end)]
    colnames(to_join) <- paste(colnames(to_join), "o", sep = "_")
    ortho1 <- to_join[
      ortho1, on = .(symbol_o)
    ]
  }
  ## Get the ortholog info of peaks
  peaks_part <- peaks[, .(symbol, symbol_o, chromosome_o, start_o, end_o)]
  ortho1 <- as.data.table(ortho1)
  ortho1 <- rbind(peaks_part, ortho1)
  distance_matrix <- ortho1[distance_matrix, on = .(symbol)]
  colnames(ortho1) <- str_replace(colnames(ortho1), "$", "2")
  distance_matrix <- ortho1[distance_matrix, on = .(symbol2)]
  return(distance_matrix)
}