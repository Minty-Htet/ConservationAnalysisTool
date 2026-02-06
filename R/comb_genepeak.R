#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom data.table tstrsplit
#' @importFrom data.table fifelse
#' @importFrom data.table setkey
#' @importFrom data.table melt
#' @importFrom data.table copy
#' @importFrom data.table setcolorder
#' @keywords internal

empty_peak <- matrix(nrow = 0, ncol = 0)

comb_genepeak <- function(species, genome, gene_vec, peak_vec = empty_peak, dist_to_group = 50000, exp_rocr = 500000, threshold = 0.5) {
    dim_peak <- c(NROW(peak_vec), NCOL(peak_vec))
    ## If no peaks provided, just combine with empty peak
    if (dim_peak[1] == 0 | dim_peak[2] == 0) {
        chr_rows <- data.table(
            chromosome = character(),
            start = numeric(),
            end = numeric(),
            symbol = character(),
            pmid = numeric(),
            chromosome_o = character(),
            start_o = numeric(),
            end_o = numeric(),
            symbol_o = character()
        )
    }
    ## If the input is a character vector, convert to data frame
    if (typeof(peak_vec) == "character") {
        peak_vec <- as.data.table(peak_vec)
        peak_vec[
                , c("chromosome", "start", "end") := data.table::tstrsplit(peak_vec, split = "[,\\.[:space:]_]+", fixed = FALSE)
            ][
                , start := as.numeric(start)
            ][
                , end := as.numeric(end)
            ]
        colnames(peak_vec) <- c("symbol", "chromosome", "start", "end")
        peak_vec <- peak_vec[, .(chromosome,start,end,symbol)]
        peak_vec <- as.data.frame(peak_vec)
        ## Convert peaks using liftover
        peaks_part <- convert_peaks(peak_vec, genome, threshold)
        ## Prepare for assigning group ids to peaks
        chr_rows <- peaks_part[[3]][
            , pmid := (start + end) / 2
        ]
    }
    ## Long form of pairwise distance matrix
    both_vec <- get_distance_matrix(species, gene_vec, chr_rows)
    ## Need to join adjusted peak locations
    both_vec[
        chr_rows[, .(chromosome, start, end, symbol)], on = .(symbol = symbol), chromosome := i.chromosome
    ][
        chr_rows[, .(chromosome, start, end, symbol)], on = .(symbol = symbol), start := i.start
    ][
        chr_rows[, .(chromosome, start, end, symbol)], on = .(symbol = symbol), end := i.end
    ][
        chr_rows[, .(chromosome, start, end, symbol)], on = .(symbol2 = symbol), chromosome2 := i.chromosome
    ][
        chr_rows[, .(chromosome, start, end, symbol)], on = .(symbol2 = symbol), start2 := i.start
    ][
        chr_rows[, .(chromosome, start, end, symbol)], on = .(symbol2 = symbol), end2 := i.end
    ]
    ## Only take within chromosome distances
    both_vec <- both_vec[chromosome == chromosome2]
    ## Join ortholog info
    both_vec <- attach_orthoinfo(both_vec, chr_rows, species)
    ## Form groups using just genes 
    chr_part <- both_vec[grepl('chr', symbol) & grepl('chr', symbol2)]
    genes_part <- both_vec[!grepl('chr', symbol) & !grepl('chr', symbol2)]
    genes_part <- assign_groupid(genes_part,exp_rocr)

    ## Get the bounds of the groups based on the boundaries
    ## Start of the 5' most gene and end of 3' most gene
    ref_cols <- grep("2$", names(genes_part), value = TRUE, invert = TRUE)
    nonchr_rows <- unique(genes_part[
         , ..ref_cols][
            , distance := NULL
         ][ 
             , in_group := any(in_group == TRUE), by = symbol
         ]
         )[ 
             , in_group := fifelse(is.na(in_group), FALSE, TRUE)
         ]
    
    bound_rows <- nonchr_rows[
            , .(chromosome, start, end, group_id)
        ]
    bound_rows <- unique(bound_rows[
        , start := min(start) , by = group_id
        ][
        , end := max(end), by = group_id
        ])

    ## Join matches within same chromosome
    matches <- bound_rows[chr_part,
        on = .(chromosome),
        allow.cartesian = TRUE,
        nomatch = 0L
    ]

    ## matches now has gene and peak rows
    ## distance is minimum distance of peak midpoint to group start or peak start to group end
    ## Keep only distances less than dist_to_group
    ## And keep peaks between the group start and group end
    matches[
            , start_dist := abs(i.end - start)
        ][
            , end_dist := abs(i.start - end)
        ][ 
            , dist := pmin(end_dist, start_dist)
        ]
    matches <- matches[dist < dist_to_group | (i.end >= start & i.start <= end)]
    ## Keep the nearest
    nearest <- matches[
        matches[, .I[which.min(dist)], by = .(symbol)]$V1
    ][
        , .(symbol, group_id)
    ][
        chr_rows, on = .NATURAL
    ][ 
        , start_o := as.integer(start_o) 
    ][ 
        , end_o := as.integer(end_o)
    ]

    ## Combine back with the peaks part
    ## Remove the boundary rows
    comb_dt <- merge(nonchr_rows, nearest, all = TRUE)[
         , chromosome := factor(chromosome, levels = hs_chr_order)
    ][ 
        order(chromosome, start)
    ][
        !is.na(symbol)
    ][
        , pmid := NULL
    ]

    ## Combine the peak and group_id stuff to the full peak and gene matrix
    groups_df <- comb_dt[,.(symbol,group_id)]
    both_vec <- both_vec[
        groups_df, on = .(symbol)
    ]
    colnames(groups_df) <- c("symbol2", "group_id2")
    both_vec <- both_vec[
        groups_df, on = .(symbol2)
    ]
    
    list_of_dt <- list(
        comb_dt,
        both_vec
    )
    return(list_of_dt)
}