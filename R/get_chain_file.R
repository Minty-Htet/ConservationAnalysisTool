#' @importFrom rappdirs user_data_dir
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @keywords internal

get_chain_file <- function(from, to, chain_dir) {
  chain_filename <- paste0(from,"To", toupper(substr(to,1,1)), substr(to,2, nchar(to)) , ".over.chain")
  cache_dir <- user_data_dir("ComparativeAnalysisTool")
  chain_path <- if (!is.null(chain_dir)) {
      print(paste("Searching for",paste(chain_dir, chain_filename, sep = '/')))
      file.path(chain_dir, chain_filename)
    } else {
      print(paste("Searching for",paste(cache_dir, chain_filename, sep = '/')))
      file.path(cache_dir, chain_filename)
    }

  # Download if file doesn't already exist
  if (!file.exists(chain_path)) {
    print('Downloading from UCSC')
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    gz_url <- sprintf(
      "https://hgdownload.soe.ucsc.edu/goldenPath/%s/liftOver/%s.gz",
      from, paste0(from,"To", toupper(substr(to,1,1)), substr(to,2, nchar(to)) , ".over.chain")
    )
    gz_path <- paste0(chain_path, ".gz")
    download.file(url = gz_url, destfile = gz_path, quiet = FALSE)
    gunzip(gz_path, destname = chain_path, remove = TRUE) ## Need R.utils
    sprintf("Downloaded from UCSC: %s", gz_path)
  } 
  return(chain_path)
}
