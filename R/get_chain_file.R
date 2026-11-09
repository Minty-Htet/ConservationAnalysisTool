#' @importFrom rappdirs user_data_dir
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @keywords internal

get_chain_file <- function(from, to, chain_dir) {
  chain_filename <- paste0(from,"To", toupper(substr(to,1,1)), substr(to,2, nchar(to)) , ".over.chain")
  cache_dir <- user_data_dir("ConservationAnalysisTool")
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

#' @keywords internal
parse_chain <- function(path) {
  x <- data.table(line = readLines(path))      # readLines reads .gz directly
  x <- x[line != ""]
  x[, is_hdr := startsWith(line, "chain")]
  x[, chain_row := cumsum(is_hdr)]             # 1 for every line in chain 1, etc.

  hdr <- x[is_hdr == TRUE, tstrsplit(line, "[ \t]+", keep = 2:13)]
  setnames(hdr, c("score","tName","tSize","tStrand","tStart","tEnd",
                  "qName","qSize","qStrand","qStart","qEnd","id"))
  num <- c("score","tSize","tStart","tEnd","qSize","qStart","qEnd","id")
  hdr[, (num) := lapply(.SD, as.numeric), .SDcols = num]
  hdr[, chain_row := .I]

  body <- x[is_hdr == FALSE]
  body[, c("size","dt","dq") := tstrsplit(line, "[ \t]+", type.convert = TRUE)]
  body[is.na(dt), `:=`(dt = 0L, dq = 0L)]      # last line of each chain has only size
  body[, line := NULL]

  # join header info onto blocks
  blocks <- hdr[body, on = "chain_row"]

  # running offsets within each chain
  blocks[, `:=`(
    t_start = tStart + cumsum(size + dt) - (size + dt),
    q_start = qStart + cumsum(size + dq) - (size + dq)
  ), by = chain_row]
  blocks[, `:=`(t_end = t_start + size, q_end = q_start + size)]

  blocks[, .(id, score, tName, t_start, t_end, qName, qSize, qStrand, q_start, q_end)]
}