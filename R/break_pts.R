#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table fifelse
#' @importFrom data.table setnames
#' @keywords internal

## Ortholog Percentage Section
## Ortholog percentage tables
perc_table <- function(final_table, section = "genes") {
  if (section == "genes") {
    df_perc <- unique(final_table[
        !grepl("^chr", symbol)
    ][
      , ortho_found := fifelse(is.na(symbol_o), FALSE, TRUE)
    ][ 
      , Count := .N , by = ortho_found
    ][
      , Total := .N
    ][ 
      , .(Count, Total, ortho_found)
    ])[
      , Percentage := signif((Count/Total) * 100, 3)
    ]
  } else {
    df_perc <- unique(final_table[
        grepl("^chr", symbol)
    ][
      , ortho_found := fifelse(is.na(symbol_o), FALSE, TRUE)
    ][ 
      , Count := .N , by = ortho_found
    ][
      , Total := .N
    ][
      , .(Count, Total, ortho_found)
    ])[
      , Percentage := signif((Count/Total) * 100, 3)
    ]
  }
  return(df_perc)
  }

## Risk Region Table
risky_table <- function(long_table, flag_table) {
    flagged_table <- long_table[
        !grepl("chr", symbol)
    ]
    flagged_table <- unique(flagged_table[ 
        , lower_start_o := min(start_o), by = group_id

    ][
      , lower_start := min(start), by = group_id
    ][ 
        , upper_end_o := max(end_o), by = group_id
    ][
      , upper_end := max(end), by = group_id
    ][
        , .(group_id, lower_start, upper_end, chromosome)
    ])
    hits <- as.data.table(flag_table)[ 
        flagged_table, on = .(chromosome), allow.cartesian = TRUE
    ][
        (start <= lower_start & lower_start <= end) | 
        (start <= upper_end & upper_end <= end) 
    ]
    hits_vec <- unique(hits$group_id)
    long_table[
        , In_Risk := fifelse(group_id %in% hits_vec, TRUE, FALSE)
    ][ 
        , Risk_status := fifelse(In_Risk == TRUE, "By Breakpoint", Group_status)
    ]
    return(long_table)
}

## Run all the plots
break_pts <- function(species, final_table, disttobreak) {

  ## Filter long matrix to chosen chromosome
  ## Set the order
  final_table[[2]] <- final_table[[2]][
    order(chromosome, chromosome2, start, start2)
  ]
  final_table[[2]][
    , symbol := factor(symbol, levels=unique(final_table[[2]]$symbol))
  ][ 
    , symbol2 := factor(symbol2, levels=unique(final_table[[2]]$symbol2))
  ]

  ## Choose the right cumulative positions
  flag_dt <- data.table::copy(flag_df)
  ## For mice relabel which is reference and ortholog species
  if (species == "Mmusculus") {
    cumulative_df <- mm_cumulative_df
    colnames(cumulative_df) <- c("chromosome", "position")
    back_dt <- as.data.table(mouse_to_human_genes)[
      mmusculus_chromosome %in% hs_chr_order
    ]
  for (loop in list(back_dt, flag_dt)) {
     setnames(loop, 
    c(
      "mmusculus_chromosome", 
      "hsapiens_chromosome",
      "mmusculus_start",
      "hsapiens_start",
      "mmusculus_end",
      "hsapiens_end"
      ), 
    c(
      "chromosome", 
      "chromosome_o",
      "start",
      "start_o",
      "end",
      "end_o"
      )
      )
  }
  }
  ## For humans relabel which is reference and ortholog species
  if (species == "Hsapiens") {
    cumulative_df <- hs_cumulative_df
    colnames(cumulative_df) <- c("chromosome", "position")
    back_dt <- as.data.table(mouse_to_human_genes)[
      hsapiens_chromosome %in% hs_chr_order
    ]
  for (loop in list(back_dt, flag_dt)) {
     setnames(loop, 
    c(
      "mmusculus_chromosome", 
      "hsapiens_chromosome",
      "mmusculus_start",
      "hsapiens_start",
      "mmusculus_end",
      "hsapiens_end"
      ), 
    c(
      "chromosome_o", 
      "chromosome",
      "start_o",
      "start",
      "end_o",
      "end"
      )
      )
  }
  }
  ## Add in disttobreak <distance to breakpoint until the group is in danger zone>
  flag_dt <- as.data.table(flag_dt)[
    , start := start + disttobreak
  ][
    , end := end + disttobreak
  ]

  if (!species %in% c("Hsapiens", "Mmusculus")) {
    print("Type Hsapiens or Mmusculus for species")
  }

  ## Combine the cumulative positions to the start and end positions
  back_dt <- cumulative_df[
    back_dt, on = "chromosome"
  ][
    , mid := (start + end)/2 + as.numeric(position)
  ][
    , mid_o := (start_o + end_o)/2
  ]
    adj_dt <- cumulative_df[
        risky_table(final_table[[1]], flag_dt), on = "chromosome"
    ][
        , mid := (start + end)/2 + as.numeric(position)
    ][
        , mid_o := (start_o + end_o)/2
    ]

    list_of_tables <- list()
    ## Long table
    list_of_tables[[1]] <- risky_table(final_table[[1]], flag_dt)
    ## Wide table in Long Format
    list_of_tables[[2]] <- final_table[[2]]
    #list_of_tables[[2]] <- perc_table(final_table[[1]], "genes")
    #list_of_tables[[3]] <- perc_table(final_table[[1]], "peaks")
    # Combined the cumulative positions to the start and end positions
    #list_of_tables[[4]] <- back_dt
    #list_of_tables[[5]] <- adj_dt
    return(list_of_tables)
}