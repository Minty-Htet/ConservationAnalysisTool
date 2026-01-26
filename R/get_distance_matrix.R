#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table fifelse
#' @importFrom data.table setkey
#' @importFrom data.table melt
#' @importFrom data.table copy
#' @importFrom stringr str_replace
#' @keywords internal

## Helper function to merge processed ensembl id's to gene names
merge_df <- function(species, input_long) {
  input_long <- as.data.table(input_long)
  if (species == "Mmusculus") {
    all_start_end <- as.data.table(mouse_all_startend)
    setkey(all_start_end, id)    
    result <- merge(input_long, all_start_end, by = "id", all.x = TRUE)
    result <- as.data.table(result)
    setkey(result, id)
    result <- merge(result, all_start_end, by.x = "id2", by.y = "id", all.x = TRUE, suffixes = c("", ".id2"))
    
  } else {
    all_start_end <- as.data.table(human_all_start_end)
    setkey(all_start_end, id)    
    result <- merge(input_long, all_start_end, by = "id", all.x = TRUE)
    result <- as.data.table(result)
    setkey(result, id)
    result <- merge(result, all_start_end, by.x = "id2", by.y = "id", all.x = TRUE, suffixes = c("", ".id2"))
  }}

## Function to get distance between genes
get_distance_matrix <- function(species, genes) {
  chosen_genes <- if (species == "Mmusculus") {
    all_start_end <- as.data.table(mouse_all_startend)
    all_start_end[symbol %in% genes]
  } else {
    all_start_end <- as.data.table(human_all_start_end)
    all_start_end[symbol %in% genes]
  }
  all_vec1 <- (chosen_genes$start + chosen_genes$end) * 0.5 # Use midpoints to calculate distance
  names(all_vec1) <- chosen_genes$id
  comp_matrix_all <- outer(all_vec1, all_vec1, "-")
  rownames(comp_matrix_all) <- colnames(comp_matrix_all)
  long_matrix_all <- as.data.table(comp_matrix_all, keep.rownames = "id")
  long_matrix_all <- melt(long_matrix_all, 
                       id.vars = "id", 
                       variable.name = "id2", 
                       value.name = "distance")
  long_matrix_all[, ':='(
   id = as.numeric(id),
   id2 = as.numeric(as.character(id2))  # in case id2 is a factor
  )]

  matrix_long <- merge_df(species, long_matrix_all) 
  colnames(matrix_long) <- c("id2", str_replace(colnames(matrix_long[, .SD, .SDcols = !("id2")]), ".id", ""))
  matrix_long <- matrix_long[chromosome == chromosome2]
  return(matrix_long)
}