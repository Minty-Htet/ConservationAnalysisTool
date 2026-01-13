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

comb_genepeak <- function(species, genome, gene_vec, peak_vec = empty_peak, dist_to_group = 50000, exp_rocr = 500000) {

    ## Get the gene part
    genes_part <-  attach_orthoinfo(
        assign_groupid(
            get_distance_matrix(species, genome, gene_vec),
            exp_rocr
        ), 
        species
        )
    ref_cols <- grep("2$", names(genes_part), value = TRUE, invert = TRUE)
    dim_peak <- c(NROW(peak_vec), NCOL(peak_vec))
    ## If no peaks provided, just combine with empty peak
    if (dim_peak[1] == 0 & dim_peak[2] == 0) {
        nonchr_rows <- unique(genes_part[
         , ..ref_cols][
            , distance := NULL
         ][ 
             , in_group := any(in_group == TRUE), by = id
         ]
         )[ 
             , in_group := fifelse(is.na(in_group), FALSE, TRUE)
         ]
        list_of_dt <- list(
            nonchr_rows,
            genes_part
        )
        return(list_of_dt)
    }
    ## If the input is a character vector, convert to data frame
    if (typeof(peak_vec) == "character") {
            peak_vec <- as.data.table(peak_vec)
            peak_vec[
                , c("chromosome", "start", "end") := data.table::tstrsplit(peak_vec, split = "[,\\.[:space:]]+", fixed = FALSE)
            ][
                , start := as.numeric(start)
            ][
                , end := as.numeric(end)
            ]
            colnames(peak_vec) <- c("symbol", "chromosome", "start", "end")
            setcolorder(peak_vec, c("chromosome", "start", "end", "symbol"))
            peak_vec <- as.data.frame(peak_vec)
    }
    ## Convert peaks using liftover
    peaks_part <- convert_peaks(peak_vec, genome)

    ## Prepare for assigning group ids to peaks
    chr_rows <- peaks_part[[3]][
            , pmid := (start + end) / 2
        ]
    ## Get the bounds of the groups based on the boundaries
    ## Start of the 5' most gene and end of 3' most gene
    nonchr_rows <- unique(genes_part[
         , ..ref_cols][
            , distance := NULL
         ][ 
             , in_group := any(in_group == TRUE), by = id
         ]
         )[ 
             , in_group := fifelse(is.na(in_group), FALSE, TRUE)
         ]
    bound_rows <- unique(nonchr_rows[
        , start := min(start) , by = group_id
        ][
        , end := max(end), by = group_id
        ][
            , .(chromosome, start, end, group_id)
        ])

    ## Join matches within same chromosome
    matches <- bound_rows[chr_rows,
        on = .(chromosome),
        allow.cartesian = TRUE,
        nomatch = 0L
    ]
    ## matches now has gene and peak rows
    ## distance is minimum distance of peak midpoint to group start or peak start to group end
    ## Keep only distances less than dist_to_group
    ## And keep peaks between the group start and group end
    matches[
            , start_dist := abs(pmid - start)
        ][
            , end_dist := abs(pmid - end)
        ][ 
            , dist := pmin(end_dist, start_dist)
        ]
    matches <- matches[dist < dist_to_group | (pmid >= start & pmid <= end)]
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
        , pmid := NULL
    ][
        !is.na(symbol)
    ]
    
    list_of_dt <- list(
        comb_dt,
        genes_part
    )
     return(list_of_dt)
}