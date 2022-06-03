#' Left join intersection of two tibble BEDs
#'
#' Intersect two tibble BEDs and retain the full intervals from the lefthand tibble BED that have any overlap with the righthand tibble BED.
#' Optionally retain the interval coordinates for the righthand tibble in the metadata.
#' If the two input tibble BEDs share any common metadata column names, or if keepBcoords = TRUE, then unique suffixes must be provided.
#' This function is built as a modification of the valr::bed_intersect() function.
#' Equivalent to bedtools -a A.bed -b B.bed -wa
#'
#' @param bedA A tibble BED with columns: chrom, start, end, ...(metadata cols)...
#' @param bedB A tibble BED with columns: chrom, start, end, ...(metadata cols)...
#' @param suffix A character vector of suffixes to add to metadata column names. (If the two input tibble BEDs share any common metadata column names, suffixes must be provided.)
#' @param keepBcoords If TRUE, coordinates from bedB are retained in the metadata. Must provide unique suffix for bedB to avoid start/end column duplication. Default is TRUE.
#' @param keepBmetadata If FALSE, bedB columns are dropped and only the original bedA columns are returned, filtered to those intervals that intersect with any bedB intervals.
#' @param suffixMetadataCols If TRUE, suffixes are also added to the end of metadata column names, not just the coordinate columns. Default is FALSE.
#'
#' @return A tibble BED of all intervals from bedA that contain intersects with bedB, plus all metadata columns from both
bed_intersect_left <- function(bedA, bedB, suffix = c("",""), keepBcoords = T, keepBmetadata = T, suffixMetadataCols = F){

  # silence "no visible binding" NOTE for data variables
  . <- start <- start.x <- end <- end.x <- .overlap <- NULL

  # Make sure suffix has been provided if there are common metadata columns
  check_A_vs_B_metadata_cols(bedA, bedB, suffix, keepBcoords, suffixMetadataCols)

  # Left join intersect
  bedA %>%
    valr::bed_intersect(bedB) %>%
    # only allow intersections with zero-overlaps if the feature is an indel (start == end)
    dplyr::filter(.overlap > 0 | start.x == end.x | start.y == end.y) %>%
    dplyr::select(-.overlap) %>%
    dplyr::rename(start = start.x, end = end.x) %>%
    { if(!keepBcoords) dplyr::select(., -c("start.y", "end.y")) else . } %>%
    { if(suffixMetadataCols)
        dplyr::rename_with(., ~ gsub(".x$", suffix[1], .)) %>%
        dplyr::rename_with(., ~ gsub(".y$", suffix[2], .))
      else
        dplyr::rename_with(., ~ gsub("(?<=start|end).y$", suffix[2], ., perl = T)) %>%
        dplyr::rename_with(., ~ gsub(".[xy]$", "", .))
    } %>%
    { if(!keepBmetadata & !keepBcoords)
        dplyr::select(., dplyr::any_of(c(names(bedA), paste0(names(bedA), suffix[1]))))
      else if (!keepBmetadata & keepBcoords)
        dplyr::select(., dplyr::any_of(c(names(bedA), paste0(names(bedA), suffix[1]), paste0(c("start", "end"), suffix[2]))))
      else
        .
    }  %>%
    dplyr::distinct()

}

check_A_vs_B_metadata_cols <- function(bedA, bedB, suffix, keepBcoords, suffixMetadataCols){
  metaA <- names(bedA) %>% setdiff(c("chrom","start","end"))
  metaB <- names(bedB) %>% setdiff(c("chrom","start","end"))
  if(length(intersect(metaA, metaB) > 0) && (suffix[1] == suffix[2])){
    stop("bedA and bedB have common metadata column names but unique suffixes have not been provided.",
         "\nPlease provide two unique suffixes to be appended to the metadata columns in the output.",
         "\nShared metadata column name(s): ", intersect(metaA, metaB),
         "\nProvided suffixes: ", suffix)
  }
  if(length(intersect(metaA, metaB) > 0) && !suffixMetadataCols){
    stop("bedA and bedB have common metadata column names, but suffixes are not being appended (suffixMetadataCols = F), so bedA and bedB metadata columns are duplicated.",
         "\nPlease rename shared metadata columns or set suffixMetadataCols = TRUE and provide unique suffixes to be appended.",
         "\nShared metadata column name(s): ", intersect(metaA, metaB))
  }
  if(keepBcoords && suffix[1] == suffix[2]){
    stop("Coordinates for bedB are being retained (keepBcoords = TRUE), but a unique bedB suffix has not been provided, so bedA and bedB coordinate columns (start, end) are duplicated",
         "\nPlease provide a suffix for bedB.",
         "\nProvided suffixes: ", suffix)
  }
}
