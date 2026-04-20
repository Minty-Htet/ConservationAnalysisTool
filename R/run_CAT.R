#' Convert genomic positions formatted in
#' bed file like data frame between genomes
#' @param species Which species did the genes and genomic positions come from? Hsapiens or Mmusculus
#' @param genome A genome version. Mmusculus has mm39 and mm10 options. Hsapiens has hg38 option.
#' @param genes A vector of gene symbols. Be aware the capitalization matters.
#' @param peaks A bed file like dataframe. Requires a "chromosome", "start", "end", "symbol". The symbol is just a unique id for the peak. Requires internet for the first time you run this to download graph files from UCSC
#' @param disttobreak How close should the group be to an evolutionary breakpoint to be considered at risk?
#' @param dist_to_group How close a peak must be to the outside edges of a group to be considered part of said group
#' @param exp_rocr Expected size of regions of co-regulation (Usually the mean size of the TADs)
#' @param threshold Size of the region in the comparator genome when querying regions in liftover at which you say there is an orthologous region
#' @param chain_path A path to a directory containing chain files if you have them
#' @return Tables summarizing 1) groups and proximity to breakpoints, 2-3) percentage of input genes we found one to one orthologs for, 4) the 1st table with cumulative position adjustment for plotting purposes, 6) long form matrix of distances between genes
#' @keywords export
#' @export

## To run the Comparative Analysis Tool
run_CAT <- function(species, genome, genes, peaks = NULL, disttobreak = 100000, dist_to_group = 50000, exp_rocr = 500000, threshold = 0.5, chain_path = NULL) {
    ex1 <- comb_genepeak(species, genome, genes, peaks, dist_to_group, exp_rocr, threshold, chain_path)
    risk_tables <- conservation_table(ex1[[2]], ex1[[1]], exp_rocr)
    final_table <- break_pts(species, risk_tables, disttobreak)
    return(final_table)
}
