#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table uniqueN
#' @importFrom data.table fifelse
#' @importFrom data.table setkey
#' @importFrom data.table copy
#' @importFrom data.table setcolorder
#' @importFrom data.table shift
#' @keywords internal

conservation_table <- function(matrix, table, split_dis) {

  ## Identify singlet groups
  no_singlet <- unique(copy(table)[
    !grepl("^chr", symbol)
    ][
    , in_group := .N > 1, by = group_id
    ][ 
      , .(group_id, in_group)
    ])

  ## Sequence conservation check
  ortho_dt <- unique(table[
    !grepl("^chr", symbol)
    ][ 
    , .(symbol, symbol_o, chromosome_o, start_o, group_id)
  ][  
    , n_unique_symbols := data.table::uniqueN(symbol), by = group_id
  ][ 
    ,  n_oto := data.table::uniqueN(symbol_o[!is.na(symbol_o)]), by = group_id
  ][, one_to_one := fifelse(n_unique_symbols == 1 & n_oto == 1, "Singlet group",
        fifelse(n_oto == 0, "Lost",
          fifelse(n_unique_symbols != 1 & n_oto == 1, "Became Singlet",
            fifelse(n_unique_symbols != 1 & n_oto != 1 & n_oto != n_unique_symbols, "Partial Loss", "No Loss"))))
  ][
    , .(group_id, one_to_one)
  ])

  ## Location conservation check
  location_dt <- table[
    ortho_dt, on = "group_id"
  ][ 
    , .(symbol, symbol_o, chromosome_o, start_o, group_id)
  ][
    order(group_id, start_o)
  ][ 
    , lag_start := shift(start_o), by = group_id
  ][ 
    , lag_start := fifelse(is.na(lag_start), start_o, lag_start)
  ][ 
    , lag_chr := shift(chromosome_o), by = group_id
  ][ 
    , lag_chr := fifelse(is.na(lag_chr), chromosome_o, lag_chr)
  ][
    , Disrupt := fifelse(abs(lag_start - start_o) <= split_dis & lag_chr == chromosome_o, FALSE, TRUE) , by = group_id
  ][
    , location_disrupted := ifelse(any(Disrupt, na.rm = TRUE), "TRUE", "FALSE"), by = group_id
  ][
    no_singlet, on = .NATURAL
  ][
    , location_disrupted := fifelse(is.na(group_id), NA, location_disrupted)
  ]

  ## only one ortholog left
  only_one_left <- table[
    !is.na(group_id)
  ][
    , .(group_id, symbol, chromosome_o)
  ][ 
    , left := sum(!is.na(chromosome_o)) , by = .(group_id)
  ]

  ## Not Applicable list
  non_appl_vec <- only_one_left[left < 2]$group_id

  ## Summary of group status
  status_dt <- unique(location_dt[
    ortho_dt, on = .NATURAL
  ][ 
    , Group_status := fifelse(location_disrupted == TRUE, 'Group Split',
      fifelse(group_id %in% non_appl_vec, 'Not Applicable', 'Group Together')
    )
  ][
    , Group_status := fifelse(group_id %in% ortho_dt[one_to_one == 'Lost']$group_id, 'Not Applicable', Group_status)
  ][ 
    , .(group_id, one_to_one, location_disrupted, Group_status, in_group)
  ])[
    !is.na(group_id)
  ]

  ## Build the final table
  table_dt <- status_dt[
    table, on = "group_id"
  ][
    , i.in_group := NULL ## quick fix IDK where this comes from
  ]
  
  ## Build the heatmap matrix
  heatmap_dt <- status_dt[
    matrix, on = "group_id"
    ]
  colnames(status_dt) <- str_replace(colnames(status_dt), "$", "2")
  heatmap_dt <- status_dt[heatmap_dt, on = "group_id2"]
  heatmap_dt <- heatmap_dt[
        , one_to_one := ifelse(group_id == group_id2, one_to_one, "FALSE")
    ][ 
        , Group_status := ifelse(group_id == group_id2, Group_status, "FALSE")
    ][, 
    grouped := fifelse(group_id != group_id2, "FALSE", 
    fifelse(in_group == TRUE & in_group2 == TRUE, "Grouped", "Singlet group")
  )]
  heatmap_dt[, symbol := factor(symbol, levels = unique(symbol))]
  heatmap_dt[, symbol2 := factor(symbol2, levels = unique(symbol2))]

  ## table and matrix within a list
  final_tables <- list(
    table_dt[
         , chromosome := factor(chromosome, levels = hs_chr_order)
    ][ 
      order(chromosome, start)
    ],
    heatmap_dt[
         , chromosome := factor(chromosome, levels = hs_chr_order)
    ][ 
      order(chromosome, chromosome2, start, start2)
    ],
    status_dt,
    only_one_left
  )
  return(final_tables)
}