#' Intersect BEDPE list
#'
#' Intersect BEDPE list object with SNPs and TSSs
#'
#' @param bedpe A tibble BEDPE list object representing chromosomal interaction looping data
#' @param SNPend A tibble BED of trait SNPs
#' @param TSSend A tibble BED of gene TSSs
#'
#' @return A tibble BED of all interaction loops in which one end intersects with a SNP and the other end intersects with a TSS. The coordinates for the SNP
#' @export
intersect_BEDPE <- function(bedpe, SNPend, TSSend) {

  # silence "no visible binding" NOTE for data variables
  . <- InteractionID <- origin <- SNP_origin <- ensg <- enst <- NULL

  # Unlist BEDPE list object (InteractionID and first/last origin can still trace loops)
  PE <- bedpe %>%
    dplyr::bind_rows(.id = "origin")

  # Find loops in which one end overlaps a SNP
  SNP_end <- SNPend %>%
    bed_intersect_left(PE, keepBcoords = F) %>%
    dplyr::rename(SNP_origin = origin)

  # Intersect the opposite ends of these SNP-intersected loops with TSSs...
  TSS_end <- SNP_end %>%
    # Get IDs and origins of the _opposite_ end of all SNP-intersected loops
    dplyr::transmute(InteractionID,
                     SNP_origin,
                     TSS_origin = dplyr::case_when(SNP_origin == "first" ~ "last",
                                                   SNP_origin == "last" ~ "first")) %>%
    # Get coordinates of the _opposite_ end of all SNP-intersected loops
    dplyr::left_join(PE,
                     by = c("InteractionID" = "InteractionID",
                            "TSS_origin" = "origin")) %>%
    # Intersect with TSSs
    bed_intersect_left(TSSend, suffix = c(".SNP", ".TSS")) %>%
    # Extract intersected gene metadata
    dplyr::select(dplyr::ends_with(".TSS"),
                  ensg:enst,
                  InteractionID,
                  SNP_origin)

  # Combine SNP and TSS intersection data
  SNP_end %>%
    dplyr::as_tibble() %>%
    dplyr::right_join(TSS_end,
                      by = c("InteractionID", "SNP_origin")) %>%
    dplyr::select(-SNP_origin) %>%
    dplyr::distinct()
}
