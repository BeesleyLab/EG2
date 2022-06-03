#' Import a BED
#'
#' Import a BED file as a tibble
#'
#' @param bedfile A tab-delimited file in bed format with columns: chrom, start, stop, ...(metadata cols)...
#' @param metadata_cols A character vector of the names of metadata columns (bedfile cols 4+) in the BED file, to be imported as metadata in the tibble BED. These must be in the same order as they appear in the file, starting from column 4.
#'
#' @return A tibble representation of the BED file, optionally with metadata columns. If no metadata column names are provided, only the chrom, start and stop columns will be returned.
import_BED <- function(bedfile, metadata_cols = NULL) {

  # silence "no visible binding" NOTE for data variables
  . <- NULL

  read_tibble(bedfile) %>%
    dplyr::rename_with(., ~ c("chrom", "start", "end"), 1:3) %>%
    { if(!is.null(metadata_cols)) dplyr::rename_with(., ~ metadata_cols, 4:(3 + length(metadata_cols))) else . } %>%
    dplyr::select(c("chrom", "start", "end", metadata_cols))

}
