#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table setcolorder
#' @importFrom data.table setnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer liftOver
#' @importFrom rtracklayer import.chain
#' @importFrom data.table :=
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

#' @keywords internal
choose_chains_rds <- function(genome) {
    load_chain <- function(name) {
        readRDS(system.file("extdata", paste0(name, ".rds"),
                            package = "ConservationAnalysisTool", mustWork = TRUE))
    }
    chain <- list()
    if (genome == "mm39") {
        chain[[1]] <- load_chain("mm39_to_mm10")
        chain[[2]] <- load_chain("mm10_to_hg38")
    }
    if (genome == "mm10") {
        chain[[1]] <- load_chain("mm10_to_mm39")   
        chain[[2]] <- load_chain("mm10_to_hg38")
    }
    if (genome == "hg38") {
        chain[[1]] <- load_chain("hg38_to_mm39")
    }
    return(chain)
}

#' @keywords internal
run_liftover <- function(peaks, chain, threshold = 0.5, max_gap = Inf) {
  empty <- data.table::data.table(
    chromosome = character(), start = integer(),
    end = integer(), symbol = character()
  )
 
  peaks <- as.data.frame(peaks)
  n <- nrow(peaks)
  if (n == 0L) return(empty)
 
  sym_vec <- if ("symbol" %in% names(peaks)) {
    as.character(peaks$symbol)
  } else {
    as.character(seq_len(n))
  }
  chr <- as.character(peaks[[1]])
  st  <- as.integer(peaks[[2]])
  en  <- as.integer(peaks[[3]])
  len_vec <- en - st                       # original width (0-based, half-open)
 
  # One liftOver call for every peak. +1L converts 0-based start to 1-based.
  gr     <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = st + 1L, end = en))
  lifted <- rtracklayer::liftOver(gr, chain)
  out    <- unlist(lifted, use.names = FALSE)
  if (length(out) == 0L) return(empty)
 
  # One row per lifted piece, tied to its input peak by index (not by symbol).
  pieces <- data.table::data.table(
    peak       = rep(seq_along(lifted), lengths(lifted)),
    chromosome = as.character(GenomeInfoDb::seqnames(out)),
    start      = BiocGenerics::start(out) - 1L,   # back to 0-based
    end        = BiocGenerics::end(out)
  )
  data.table::setorder(pieces, peak, chromosome, start)
 
  # Gap to the previous piece, computed within each peak and chromosome only,
  # so nothing leaks across peaks or chromosomes.
  pieces[, gap := start - data.table::shift(end), by = .(peak, chromosome)]
 
  # A new fragment starts at the first piece on each chromosome, or after a
  # gap larger than max_gap.
  pieces[, frag := cumsum(is.na(gap) | gap > max_gap), by = peak]
 
  frags <- pieces[, .(chromosome = chromosome[1L],
                      start = min(start),
                      end   = max(end)),
                  by = .(peak, frag)]
  frags[, frag_len := end - start]
 
  # Largest fragment per peak (ties: first by chromosome/position).
  data.table::setorder(frags, peak, -frag_len)
  best <- unique(frags, by = "peak")
 
  best[, `:=`(symbol = sym_vec[peak], total_len = len_vec[peak])]
  best <- best[frag_len > threshold * total_len]
  data.table::setorder(best, peak)
 
  best[, .(chromosome, start, end, symbol)]
}

#' Convert genomic positions formatted in
#' bed file like data frame between genomes
#' @param peaks A bed file like dataframe. Requires a "chromosome", "start", "end", "symbol". The symbol is just a unique id for the peak
#' @param genome A genome version. Mmusculus has mm39 and mm10 options. Hsapiens has hg38 option.
#' @param threshold Ratio (0 to 1) of resulting chain to peak query length which must be met in order to be considered an orthologous peak
#' @param chain_path A path to where you store chain files if you have them
#' @return A list containing tables of 1) the mouse to human coordinates, 2) the mm10 to mm39 or mm39 to mm10 coordinates, 3) the original table with new columns showing the other species' coodinates, 4) the coordinates which failed to map to the other species
#' @keywords export
#' @export
convert_peaks <- function(peaks, genome = c("mm39", "mm10", "hg38"),
                          threshold = 0.5, chain_path = NULL) {
  genome <- match.arg(genome)
  chain  <- choose_chains_rds(genome)
 
  # copy() so we never modify the caller's table by reference
  peaks <- data.table::copy(data.table::as.data.table(peaks))
  if (!"symbol" %in% names(peaks)) stop("`peaks` needs a `symbol` column.")
  if (anyDuplicated(peaks$symbol)) {
    warning("Duplicate symbols in `peaks`; the join on symbol will multiply rows.")
  }
  peaks[, `:=`(start = as.numeric(start), end = as.numeric(end))]
 
  empty <- data.table::data.table(
    chromosome = character(), start = integer(),
    end = integer(), symbol = character()
  )
  within  <- empty   # within-species conversion (mm10 <-> mm39)
  between <- empty   # between-species conversion (to hg38 / from hg38)
 
  if (genome == "mm39") {
    within  <- run_liftover(peaks,  chain[[1]], threshold)  # mm39 -> mm10
    between <- run_liftover(within, chain[[2]], threshold)  # mm10 -> hg38
  } else if (genome == "mm10") {
    within  <- run_liftover(peaks, chain[[1]], threshold)   # mm10 -> mm39
    between <- run_liftover(peaks, chain[[2]], threshold)   # mm10 -> hg38
  } else {                                                  # hg38
    between <- run_liftover(peaks, chain[[1]], threshold)   # hg38 -> mm39
  }
 
  data.table::setnames(between, c("chromosome_o", "start_o", "end_o", "symbol"))
 
  # For mm10 the "original" side is the mm39 conversion; otherwise the input.
  base <- if (genome == "mm10") within else peaks
 
  # Right join: keeps every row of `base`; peaks with no ortholog get NA in *_o.
  result <- between[base, on = "symbol"]
 
  # Build the ID while start_o/end_o are still integers, so large values
  # are not formatted in scientific notation (e.g. 1e+08).
  result[, symbol_o := paste(chromosome_o, start_o, end_o, sep = ".")]
  result[, `:=`(start   = as.numeric(start),   end   = as.numeric(end),
                start_o = as.numeric(start_o), end_o = as.numeric(end_o))]
 
  return(result)
}