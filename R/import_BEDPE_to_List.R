#' Import a BEDPE
#'
#' Import a paired-end BEDPE file as a list object c(tibble BED of first end, tibble BED of last end).
#' The direction / grouping of first and last ends is arbitrary.
#'
#' @param bedpefile A tab-delimited file in bedpe format with columns: chromA, startA, endA, chromB, startB, endB, ...(metadata cols)...
#' @param metadata_cols A character vector of the names of metadata columns (bedfile cols 4+) in the BED file
#'
#' @return A paired (first + last) list representation of the BEDPE file as two valr-compatible tibble BEDs with an ID column matching pairs, optionally with metadata columns
#' @export
import_BEDPE_to_List <- function(bedpefile, metadata_cols = NULL, prefix_InteractionID = NULL) {

  # silence "no visible binding" NOTE for data variables
  . <- NULL

  cat("Importing paired-end BED file", basename(bedpefile), "from", bedpefile, "\n")

  bedpe <- read_tibble(bedpefile)

  if(ncol(bedpe) < (3 + length(metadata_cols))){stop("Number of metadata column names given (", length(metadata_cols), " ) is more than the number of metadata columns that exist in the file")}

  PEList <- list()
  PEList[["first"]] <- bedpe[,-c(4:6)]
  PEList[["last"]] <- bedpe[,-c(1:3)]
  lapply(PEList,
         function(x)  x %>%
           dplyr::rename_with(., ~ c("chrom", "start", "end"), 1:3) %>%
           {
             if(!is.null(metadata_cols))
               dplyr::rename_with(., ~ metadata_cols, 4:(3 + length(metadata_cols)))
             else
               .
            } %>%
           dplyr::mutate(InteractionID = dplyr::row_number() %>% paste0("i.",.)) %>%
           {
             if(!is.null(prefix_InteractionID))
               dplyr::mutate(., InteractionID = paste(prefix_InteractionID, InteractionID, sep = "_"))
             else
               .
           } %>%
           dplyr::select(c("chrom", "start", "end", "InteractionID", metadata_cols))
  )

}
