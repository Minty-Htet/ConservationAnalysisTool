#' @importFrom data.table data.table
#' @importFrom data.table setnames
#' @keywords internal

plot_tab_CAT <- function(chromosome_vec, cat_list) {
    show_table <- cat_list[
      cat_list$chromosome %in% chromosome_vec & cat_list$symbol != ""
    ][
      , group_id := as.character(group_id)
    ][
      , .(symbol, group_id, flagged, symbol_o, chromosome, show_start, end, chromosome_o, start_o, end_o)
    ]
    setnames(show_table, old = c("show_start"), new = c("start"))
    return(as.data.frame(show_table))
}