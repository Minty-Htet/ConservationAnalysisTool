#' @importFrom data.table as.data.table
#' @importFrom data.table fifelse
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom data.table copy
#' @importFrom stringr str_replace
#' @keywords internal

attach_orthoinfo <- function(distance_matrix, species) {
  if (species == "Mmusculus") {
    cols <- c("mm_id", "hsapiens_gene_name", "hsapiens_chromosome", "hsapiens_start", "hsapiens_end")
    ortho1 <- copy(mouse_to_human_genes)[, cols]
    setnames(ortho1, cols, c("id", "symbol_o", "chromosome_o", "start_o", "end_o"))
  } else {
    cols <- c("hs_id", "mmusculus_gene_name", "mmusculus_chromosome", "mmusculus_start", "mmusculus_end")
    ortho1 <- copy(mouse_to_human_genes)[, cols]
    setnames(ortho1, cols, c("id", "symbol_o", "chromosome_o", "start_o", "end_o"))
  }
  ortho1 <- as.data.table(ortho1)
  ortho2 <- copy(ortho1)
  colnames(ortho2) <- str_replace(colnames(ortho2), "$", "2")
  setkey(ortho1, id)
  setkey(ortho2, id2)
  setkey(distance_matrix, id)
  distance_matrix <- ortho1[distance_matrix]
  setkey(distance_matrix, id2)
  distance_matrix <- ortho2[distance_matrix]
  return(distance_matrix)
}