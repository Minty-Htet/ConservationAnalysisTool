#' Find Transcription Factor Binding Sites in Peaks
#' Uses JASPAR 2022 database and motifmatchr at default 0.00005 p-value cut-off against background genome
#' @import data.table
#' @importFrom TFBSTools getMatrixSet ID name matrixClass
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom JASPAR2022 JASPAR2022
#' @param peaks_df A bed file like dataframe. Requires a "chromosome", "start", "end", "symbol". The symbol is just a unique id for the peak. Requires internet for the first time you run this to download graph files from UCSC
#' @param species_id NCBI species id
#' @param genome_id UCSC genome
#' @param tf_names A list of transcription factors whose motifs you want to search for in peaks. Include "All Hs TF JASPAR 2022" (or leave NULL) to use all motifs.
#' @return Table showing the peak symbol, the location of the transcription factor binding site, and the type of transcription factor
#' @keywords export
#' @export

find_tf <- function(peaks_df, species_id, genome_id, tf_names = NULL) {
    pfms <- getMatrixSet(JASPAR2022,
                     list(collection = "CORE", tax_group = "vertebrates", species = species_id))
    # "All Hs TF JASPAR 2022": use every motif in the set
    if ("All Hs TF JASPAR 2022" %in% tf_names) tf_names <- NULL
    if (!is.null(tf_names)) {
        keep_idx <- toupper(name(pfms)) %in% toupper(tf_names)
        if (!any(keep_idx)) {
            stop("None of the requested TFs found in JASPAR for species ", species_id,
                 ": ", paste(tf_names, collapse = ", "))
        }
        pfms <- pfms[keep_idx]
    }

    collapse_tag <- function(x, tag) {
        v <- TFBSTools::tags(x)[[tag]]
        if (is.null(v) || length(v) == 0) NA_character_ else paste(v, collapse = "; ")
    }
    motif_info <- data.frame(
        motif_id = ID(pfms),
        tf       = name(pfms),
        family   = vapply(pfms, collapse_tag, character(1), tag = "family"),
        class    = vapply(pfms, function(x) {
            mc <- matrixClass(x)
               if (is.null(mc) || length(mc) == 0) NA_character_ else paste(mc, collapse = "; ")
            }, character(1)),
        row.names = NULL
    )

    peaks_gr <- makeGRangesFromDataFrame(
        peaks_df,
        seqnames.field = "chromosome",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
    )
    hits <- matchMotifs(pfms, peaks_gr,
                    genome   = genome_id,
                    out      = "positions",
                    p.cutoff = 5e-5)
    hits_gr <- unlist(hits, use.names = TRUE)
    hits_gr$motif_id <- names(hits_gr)

    ov <- findOverlaps(hits_gr, peaks_gr)   
    hits_gr <- hits_gr[queryHits(ov)]
    hits_gr$symbol <- peaks_gr$symbol[subjectHits(ov)]

    hits_dt <- as.data.table(hits_gr)[motif_info, on = "motif_id"]

    keep <- c("seqnames", "start", "end", "symbol", "strand", "tf", "motif_id")
    hits_dt <- hits_dt[, keep, with = FALSE]
    setnames(hits_dt, keep, c("tf_chromosome", "tf_start", "tf_end", "symbol", "tf_strand", "tf", "motif_id"))
    
    return(hits_dt)
}

#' wrapper to help find transcription factor binding motifs within reference peaks and comparator peaks
#' @keywords internal
wrapper_tf <- function(genome, render_dt, species_id_ref=9606, species_id_ortho=9606, tf_names=NULL) {
    if (genome == 'mm39' | genome == 'mm10') {
        #species_id_ref = 10090
        species_genome = BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39
        #species_id_ortho = 9606
        ortho_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    }
    if (genome == 'hg38') {
        #species_id_ref = 9606
        species_genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        #species_id_ortho = 10090
        ortho_genome = BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39
    }
    peaks_dt = render_dt[
        grepl("chr", symbol)
    ]
    ref_cols = c("chromosome", "start", "end", "symbol")
    ortho_cols = c("chromosome_o", 'start_o', 'end_o', 'symbol')
    simplref_cols = c("symbol", "tf", "motif_id")
    hits_ref = find_tf(
        peaks_df = peaks_dt[, ref_cols, with = FALSE],
        species_id_ref,
        species_genome,
        tf_names
    )[
        tf_strand == '+'
    ][
        , simplref_cols, with = FALSE
    ]
    hits_ref = unique(hits_ref)
    ortho_dt = peaks_dt[, ortho_cols, with = FALSE]
    ortho_dt = na.omit(ortho_dt)
    setnames(ortho_dt, ortho_cols, c("chromosome", "start", "end", "symbol"))
    hits_ortho = find_tf(
        peaks_df = ortho_dt,
        species_id_ortho,
        ortho_genome,
        tf_names
    )[
        tf_strand == '+'
    ][
        , simplref_cols, with = FALSE
    ]
    setnames(hits_ortho, simplref_cols, c("symbol", "tf_o", "motif_id_o"))
    hits_ortho <- unique(hits_ortho)
    return(list(hits_ref, hits_ortho))
}