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

distance_per_chrom <- function(input_table, chromosome_input) {
    chr1_dt <- input_table[chromosome == chromosome_input][order(start,end)]
    start_vec <- chr1_dt$start
    end_vec <- chr1_dt$end
    symbol_vec <- chr1_dt$symbol

    # build the difference matrix using outer
    comparison_matrix_chr1 <- outer(start_vec, end_vec, FUN = function(s, e) s - e)

    colnames(comparison_matrix_chr1) <- symbol_vec
    rownames(comparison_matrix_chr1) <- symbol_vec
    lower_mat <- comparison_matrix_chr1
    lower_mat[!lower.tri(lower_mat, diag = TRUE)] <- NA
    diag(lower_mat) <- 0
    lower_mat[lower_mat <= 0] <- 0
    lower_mat <- as.data.table(lower_mat, keep.rownames = 'id')
    lower_mat <- melt(
        lower_mat,
        id.vars = 'id',
        variable.name = 'id2',
        value.name = 'distance'
    )
}

## Function to get distance between all genes and in long format
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
  chrom_vec <- unique(chosen_genes$chromosome)
  result_list <- lapply(chrom_vec, function(s) distance_per_chrom (chosen_genes, s))
  long_matrix_all <- do.call(rbind, result_list)
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