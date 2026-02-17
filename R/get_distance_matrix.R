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
    all_start_end <- as.data.table(human_all_startend)
    setkey(all_start_end, id)    
    result <- merge(input_long, all_start_end, by = "id", all.x = TRUE)
    result <- as.data.table(result)
    setkey(result, id)
    result <- merge(result, all_start_end, by.x = "id2", by.y = "id", all.x = TRUE, suffixes = c("", ".id2"))
  }}

## Function to get distance between genes
get_distance_matrix <- function(species, genes, peaks) {
  chosen_genes <- if (species == "Mmusculus") {
    all_start_end <- as.data.table(mouse_all_startend)
    all_start_end[symbol %in% genes][,.(chromosome,start,end,symbol)]
  } else {
    all_start_end <- as.data.table(human_all_startend)
    all_start_end[symbol %in% genes][,.(chromosome,start,end,symbol)]
  }
  chosen_genes <- rbind(chosen_genes, peaks[, .(chromosome,start,end,symbol)])[
    order(chromosome, start)
  ]
  all_vec1 <- (as.numeric(chosen_genes$start) + as.numeric(chosen_genes$end)) * 0.5 # Use midpoints to calculate distance
  names(all_vec1) <- chosen_genes$symbol
  comp_matrix_all <- outer(all_vec1, all_vec1, "-")
  rownames(comp_matrix_all) <- colnames(comp_matrix_all)
  long_matrix_all <- as.data.table(comp_matrix_all, keep.rownames = "id")
  long_matrix_all <- melt(long_matrix_all, 
                       id.vars = "id", 
                       variable.name = "id2", 
                       value.name = "distance")
  colnames(long_matrix_all) <- c('symbol', 'symbol2', 'distance')
  if (species == "Mmusculus") {
    all_start_end <- mouse_all_startend[, .(chromosome,start,end,symbol)]
  } 
  if (species == "Hsapiens") {
    all_start_end <- human_all_startend[, .(chromosome,start,end,symbol)]
  }
  matrix_long <- all_start_end[
    long_matrix_all, on = .(symbol = symbol)
  ]
  matrix_long <- all_start_end[
    matrix_long, on = .(symbol = symbol2)
  ]
  colnames(matrix_long)[grepl("i\\.", colnames(matrix_long))] <- gsub("$", "2", colnames(matrix_long)[grepl("i\\.", colnames(matrix_long))])
  colnames(matrix_long) <- gsub("i\\.", "", colnames(matrix_long))
  return(matrix_long)
}