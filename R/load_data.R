.onLoad <- function(libname, pkgname) {
  ds <- c(
    "chr_order",
    "color_values",
    "flag_df",
    "hs_chr_order",
    "hs_cumulative_df",
    "human_all_startend",
    "keys",
    "matching_colors",
    "mm_cumulative_df",
    "mouse_all_startend",
    "mouse_to_human_genes"
  )

  for (d in ds) {
    utils::data(list = d, package = pkgname, envir = parent.env(asNamespace(pkgname)))
  }
}