#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table setcolorder
#' @importFrom data.table setnames
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer liftOver
#' @importFrom rtracklayer import.chain
#' @keywords internal
# Choose chains
# Import chain files
choose_chains <-  function(genome, chain_path) {
    chain <- list()
    if(genome == "mm39") {
        chain[[1]] <- import.chain(get_chain_file("mm39", "mm10", chain_path))
        chain[[2]] <- import.chain(get_chain_file("mm10", "hg38", chain_path))
    }
    if(genome == "mm10") {
        chain[[1]] <- import.chain(get_chain_file("mm10", "mm39", chain_path))
        chain[[2]] <- import.chain(get_chain_file("mm10", "hg38", chain_path))
    }
    if(genome == "hg38") {
        chain[[1]] <- import.chain(get_chain_file("hg38", "mm39", chain_path))
    }
    return(chain)
}

# liftover with rtracklayer
# Must do manual smoothing
# liftover with rtracklayer
# Must do manual smoothing
#' @keywords internal
run_liftover <- function(peaks, chain, threshold = 0.5) {
    liftover_dt <- data.table(start = "", end = "", chromosome = "", symbol = "")
    peaks <- as.data.frame(peaks)
    peaks$total_len <- peaks$end - peaks$start
    lift_result <- as.data.table(
        as.data.frame(liftOver(
        GRanges(seqnames = peaks[,1], ranges = IRanges(start = peaks[,2], end = peaks[,3])),
        chain
    )))[ 
        , .(seqnames, start, end)
    ][
        order(seqnames, start)
    ]
    if(nrow(lift_result) > 0) { # If there are ortholog hits
        lift_result <- lift_result[ 
            , adj_start := min(start), by = .(seqnames) # Smooth the start sites
            ][ 
            , adj_end := max(end), by = .(seqnames) # Smooth the end sites
            ][  
            , adj_chromosome := seqnames
            ][ 
            , adj_symbol := peaks$symbol # Label the peak with ID
            ][ 
                as.data.table(peaks)[, .(symbol, total_len)], on = .(adj_symbol = symbol)
            ][ 
                ,len := adj_end - adj_start
            ][ 
                len > threshold*total_len
            ][
                , .(adj_chromosome, adj_start, adj_end, adj_symbol)
            ]
    } else { # If no ortholog hits
        lift_result <- lift_result[ 
            , adj_start := start # Rename
            ][ 
            , adj_end := end 
            ][  
            , adj_chromosome := seqnames
            ][ 
            , adj_symbol := peaks$symbol 
            ]
    }
    adj_cols <- grep("^adj_", names(lift_result), value = TRUE)
    lift_result <- unique(lift_result[, ..adj_cols])
    data.table::setnames(lift_result, sub("^adj_", "", names(lift_result)))
    lift_result <- lift_result[
        , .(chromosome, start, end, symbol)
    ]
    ## Take only the first entry if there ever was a multi-hit.
    if (length(lift_result$symbol > 1)) {
        lift_result <- lift_result[1] 
        return(lift_result)
    } else {
        return(lift_result)
    } 
}

#' Convert genomic positions formatted in
#' bed file like data frame between genomes
#' @param peaks A bed file like dataframe. Requires a "chromosome", "start", "end", "symbol". The symbol is just a unique id for the peak
#' @param genome A genome version. Mmusculus has mm39 and mm10 options. Hsapiens has hg38 option.
#' @param chain_path A path to where you store chain files if you have them
#' @return A list containing tables of 1) the mouse to human coordinates, 2) the mm10 to mm39 or mm39 to mm10 coordinates, 3) the original table with new columns showing the other species' coodinates, 4) the coordinates which failed to map to the other species
#' @keywords export
#' @export
convert_peaks <- function(peaks, genome, threshold = 0.5, chain_path = NULL) {
    # Choose the chains needed
    chain <- choose_chains(genome, chain_path)
    # Make a list of liftover results
    # The first is the within species conversions
    # The second is the between species conversions
    # Or if its human to mouse its just between species
    liftover_list <- list(
        data.frame(chromosome = "", start = "", end = "", symbol = ""),
        data.frame(chromosome = "", start = "", end = "", symbol = "")
    )
    peaks <- as.data.frame(peaks)
    peaks$start <- as.numeric(peaks$start)
    peaks$end <- as.numeric(peaks$end)
    if (genome == "mm39") {
        for (x in 1:length(peaks[,4])) {
            mm39ToMm10 <- run_liftover(peaks[x,], chain[[1]], threshold) #mm39 to mm10
            mm10ToHg38 <- run_liftover(mm39ToMm10, chain[[2]], threshold) #mm10 to hg38
            liftover_list[[2]] <- unique(rbind(liftover_list[[2]], mm39ToMm10))
            liftover_list[[1]] <- unique(rbind(liftover_list[[1]], mm10ToHg38))
            }
        }
    if (genome == "mm10") {
        for (x in 1:length(peaks[,4])) {
            mm10ToMm39 <- run_liftover(peaks[x,], chain[[1]], threshold) #mm10 to mm39
            mm10ToHg38 <- run_liftover(peaks[x,], chain[[2]], threshold) #mm10 to hg38
            liftover_list[[2]] <- unique(rbind(liftover_list[[2]], mm10ToMm39))
            liftover_list[[1]] <- unique(rbind(liftover_list[[1]], mm10ToHg38))
            }
        }
    if (genome == "hg38") { 
        for (x in 1:length(peaks[,4])) {
            hg38ToMm39 <- run_liftover(peaks[x,], chain[[1]], threshold) #hg38 to mm39
            liftover_list[[1]] <- unique(rbind(liftover_list[[1]], hg38ToMm39))
            }
        }
    colnames(liftover_list[[1]]) <- c("chromosome_o", "start_o", "end_o", "symbol") # label the ortholog results
    if (genome == 'mm10') {
        liftover_list[[3]] <- as.data.table(liftover_list[[1]])[  
            as.data.table(liftover_list[[2]]), on = .NATURAL # Join with the original peaks table
        ][ 
            , symbol_o := paste(chromosome_o,start_o,end_o, sep = ".")
        ][
            , start := as.numeric(start)
        ][ 
            , end := as.numeric(end)
        ][
            , start_o := as.numeric(start_o)
        ][ 
            , end_o := as.numeric(end_o)
        ][
            symbol != ""
        ]
    } else {
        liftover_list[[3]] <- as.data.table(liftover_list[[1]])[  
            as.data.table(peaks), on = .NATURAL # Join with the original peaks table
        ][ 
            , symbol_o := paste(chromosome_o,start_o,end_o, sep = ".")
        ][
            , start := as.numeric(start)
        ][ 
            , end := as.numeric(end)
        ][
            , start_o := as.numeric(start_o)
        ][ 
            , end_o := as.numeric(end_o)
        ][
            symbol != ""
        ]
    }
    return(liftover_list)
}