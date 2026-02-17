#' @importFrom data.table as.data.table
#' @keywords internal

plot_tab_CAT <- function(chromosome_vec, cat_list) {
    show_table <- as.data.table(cat_list)[
      chromosome %in% chromosome_vec & symbol != ""
    ][
      , group_id := as.character(group_id)
    ][ 
      , status := fifelse(Risk_status == "By Breakpoint" & Group_status == "Group Split", "By Breakpoint and Split",
        fifelse(Risk_status == "By Breakpoint" & Group_status == "Group Together", "By Breakpoint",
        fifelse(Risk_status == "Not by Breakpoint" & Group_status == "Group Split", "Group Split",
        "Group Together")))
    ][
      , .(symbol, symbol_o, group_id, chromosome_o, start_o, end_o, chromosome, start, end, status, Risk_status, Group_status)
    ]
    colnames(show_table) <- c('gene_or_peak_symbol', 'ortholog_symbol', 'group_id', 
    'ortholog_chromosome', 'ortholog_start', 'ortholog_end', 'chromosome', 'start', 'end', 
    'Group_Status', 'Breakpoint_status', 'Split_status')
    return(as.data.frame(show_table))
}

search_up <- function(cat_obj1, cat_obj2, chromosome_vec, type = 'grouped') {
    order <- get_order_pluspeaks(cat_obj1, cat_obj2, type)[
      chromosome %in% chromosome_vec
    ]
    group <- cat_obj1[
        !is.na(group_id)
    ][,.(symbol,group_id)][ 
        order, on = .(symbol)
    ]
    group <- cat_obj2[
        !is.na(group_id)
    ][,.(symbol,group_id)][ 
        group, on = .(symbol)
    ][
        ,.(symbol,group_id,i.group_id,chromosome,start,coord)
    ]
    group <- group[
      chromosome %in% chromosome_vec
    ][ 
        , Which_List := fifelse(is.na(group_id), 'List2', fifelse(is.na(i.group_id), 'List1', 'Both'))
    ][ 
      , .(symbol, i.group_id, group_id, Which_List, coord)
    ]
    colnames(group) <- c('symbol', 'List1_group', 'List2_group', 'Which_List', 'position')
    return(group)
}