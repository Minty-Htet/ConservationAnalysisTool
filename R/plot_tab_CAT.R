#' @importFrom data.table as.data.table
#' @keywords internal

plot_tab_CAT <- function(chromosome_vec, cat_list) {
    show_table <- as.data.table(cat_list)[
      chromosome %in% chromosome_vec & symbol != ""
    ][
      , group_id := as.character(group_id)
    ][
      , .(group_id, symbol, symbol_o, chromosome, start, end, chromosome_o, start_o, end_o, Risk_status)
    ]
    return(as.data.frame(show_table))
}