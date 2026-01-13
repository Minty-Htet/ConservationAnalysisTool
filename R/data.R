#' Order of the mice chromosomes which will appear on plots
#' @format ## `chr_order`
#' A vector with 21 entries:
#' \describe{
#'      \item{chr1}{Chromosome 1}
#'   ...
#' }
"chr_order"

#' Which colors we used to represent the chromosomes
#' @format ## `color_values`
#' A vector with 25 entries:
#' \describe{
#'      \item{#FF4040}{red}
#'   ...
#' }
"color_values"

#' Which genes are by evolutionary breakpoints 
#' Automated screening by looking at ortholog gene coordinates
#' Followed by manual check: going through the UCSC genome browser by hand
#' @format ## `flag_df`
#' A data frame with 1,104 rows and 8 columns:
#' \describe{
#'      \item{mmusculus_gene_name}{Mouse gene name}
#'      \item{mmusculus_chromosome}{Mouse chromosome}
#'      \item{mmusculus_start}{Mouse start site}
#'   ...
#'      \item{hsapiens_gene_name}{Human gene name}
#'      \item{hsapiens_chromosome}{Human chromosome}
#'      \item{hsapiens_start}{Human start site}
#' }
"flag_df"

#' Order of the mice chromosomes which will appear on plots
#' @format ## `hs_chr_order`
#' A vector with 21 entries:
#' \describe{
#'      \item{chr1}{Chromosome 1}
#'   ...
#' }
"hs_chr_order"

#' The cumulative length of the genome by chromosome
#' Based on Gencode v36 annotation
#' @format ## `hs_cumulative_df`
#' A data frame with 25 rows and 2 columns:
#' \describe{
#'      \item{hsapiens_chromosome}{Human chromosome}
#'      \item{cumulative_position}{Cumulative position}
#'   ... 
#' }
#' @source https://www.gencodegenes.org/human/release_36.html
"hs_cumulative_df"

#' Positions of human genes
#' Based on Gencode V36 annotation
#' @format ## `human_all_startend`
#' A data frame with 70717 columns by 6 rows:
#' \describe{
#'      \item{chromosome}{Chromosome}
#'      \item{start}{Start position}
#'      \item{end}{End position}
#'      \item{strand}{Strand info}
#'      \item{id}{Ensembl id}
#'      \item{symbol}{Gene symbol}
#' }
#' @source https://www.gencodegenes.org/human/release_36.html
"human_all_startend"

#' Numeric id's matched to gene symbols
#' Based on Gencode V36 annotation for Humans
#' Based on Genocome VM35 annotation for Mice
#' @format ## `keys`
#' A list with:
#' A data frame with 718092 rows by 2 columns:
#' \describe{
#'      \item{id}{Numeric id}
#'      \item{symbol}{Mouse gene symbol}
#' }
#' A data frame with 707172 rows by 2 columns:
#' \describe{
#'      \item{id}{Numeric id}
#'      \item{symbol}{Human gene symbol}
#' }
#' @source https://www.gencodegenes.org/human/release_36.html
#' @source https://www.gencodegenes.org/mouse/release_M35.html
"keys"

#' An indexed vector of colors matched to chromsomes for plots
#' @format ## `matching colors`
#' A list with:
#' A named vector with 25 entries:
#' \describe{
#'      \item{chromosome}{Chromosome}
#'      \item{color}{Matching color}
#' }
"matching_colors"

#' The cumulative length of the genome by chromosome
#' Based on Gencode vM35 annotation
#' @format ## `mm_cumulative_df`
#' A data frame with 21 rows and 2 columns:
#' \describe{
#'      \item{mmusculus_chromosome}{Mouse Chromosome}
#'      \item{cumulative_position}{Cumulative position}
#'   ... 
#' }
#' @source https://www.gencodegenes.org/mouse/release_M35.html
"mm_cumulative_df"

#' Positions of mouse genes
#' Based on Gencode VM35 annotation
#' @format ## `mouse_all_startend`
#' A data frame with 71809 columns by 6 rows:
#' \describe{
#'      \item{chromosome}{Chromsome}
#'      \item{start}{Start}
#'      \item{end}{End}
#'      \item{strand}{Strand}
#'      \item{id}{Ensembl id}
#'      \item{symbol}{Gene Symbol}
#' }
#' @source https://www.gencodegenes.org/mouse/release_M35.html
"mouse_all_startend"

#' The list of one to one ortholog genes between mice and human
#' @format ## `mouse_to_human_genes`
#' A data frame with 1626514 columns by 14 rows:
#' #' \describe{
#'      \item{mmusculus_ensembl_id}{Mouse Ensembl id}
#'      \item{mmusculus_gene_name}{Mouse Gene Symbol}
#'   ...
#'      \item{hsapiens_ensembl_id}{Human Ensembl id}
#'      \item{hsapiens_gene_name}{Human Gene Symbol}
#' }
#' @source http://useast.ensembl.org/info/about/publications.html
"mouse_to_human_genes"