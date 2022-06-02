get_vxt_level_annotations <- function(variants,
                                      DHSs,
                                      vxt_master,
                                      annotations) {

  vxt <- list()

  # TSS distance scoring method (these pairings are the only ones to consider)
  distance <- vxt_master %>%
    dplyr::group_by(variant) %>%
    # Calculate inverse of the absolute bp distance for each variant-transcript pair (avoid infinite values)
    dplyr::mutate(inv_distance = dplyr::case_when(distance == 0 ~ 1,
                                                  TRUE ~ 1 / distance),
                  # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                  inv_distance_rank = 1 / rank(distance, ties.method = "min"))

  # intersect loop ends, by cell type, with enhancer variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  vxt$HiChIP_scores <- annotations$HiChIP %>% names %>%
    sapply(function(celltype){
        # Intersect with the HiChIP data
        annotations$HiChIP[[celltype]] %>%
          intersect_BEDPE(
            # ! For mutually exclusive intersection with ranges, make variant intervals 1bp long, equal to the end position
            SNPend = variants %>% dplyr::mutate(start = end),
            TSSend = TSSs,
            bedpe = .) %>%
          dplyr::transmute(cs, variant, enst,
                           InteractionID,
                           value = score) %>%
          dplyr::inner_join(distance %>% dplyr::select(cs, variant, enst), by = c("cs", "variant", "enst"))
      }, simplify = F, USE.NAMES = T) %>%
        dplyr::bind_rows(.id = "celltype") %>%
        tidyr::pivot_wider(id_cols = c(cs, variant, enst),
                           names_from = celltype,
                           values_from = value)

  # get variants at promoters
  vxt$promoter <- variants %>%
    dplyr::select(chrom:end, variant, cs) %>%
    # Get variants within promoter regions
    bed_intersect_left(
      promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(cs, variant, enst)

  # get variants at promoters - score by sum of signal and specificity per celltype at the DHS (~promoter activity)
  vxt$promoter_H3K27ac_bins_sum <- variants %>%
    dplyr::select(chrom:end, variant, cs) %>%
    # Get variants within promoter regions
    bed_intersect_left(
      promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    # Get DHS bins at those promoter variants
    intersect_H3K27ac(list(),
                   query = .,
                   DHSs,
                   H3K27ac = annotations$H3K27ac,
                   cs, variant, enst) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    # Sum specificity + signal bin per promoter variant per cell type
    dplyr::group_by(cs, variant, enst) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
    dplyr::ungroup()

  # variant-TSS inverse distance score method
  inv_distance <- distance %>%
    dplyr::transmute(cs, variant, enst,
                     value = inv_distance)
  
  # exonic (coding) variants
  exon <- variants %>%
    # Get variants within exons
    bed_intersect_left(
      exons, .,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(cs, variant, enst)

  # exonic or inv distance
  vxt$exon_or_inv_distance <- dplyr::full_join(inv_distance, exon %>% dplyr::mutate(exon = T),
                                               by = c("cs", "variant", "enst")) %>%
    dplyr::transmute(cs, variant, enst,
                     value = dplyr::case_when(exon ~ 1,
                                           TRUE ~ value))
  
  # TADs
  TADs_w_ID <- annotations$TADs %>%
    dplyr::bind_rows(.id = "celltype") %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(TAD = paste0(celltype, "_", dplyr::row_number())) %>%
    dplyr::ungroup()
  vxt$TADs <- dplyr::inner_join(
    bed_intersect_left(
      variants, TADs_w_ID, keepBcoords = F) %>%
      dplyr::select(cs, variant, TAD, celltype),
    bed_intersect_left(
      TSSs, TADs_w_ID, keepBcoords = F) %>%
      dplyr::select(enst, TAD, celltype),
    by = c("TAD", "celltype")
  ) %>%
    dplyr::transmute(cs, variant, enst, celltype, value = 1) %>%
    tidyr::pivot_wider(id_cols = c("cs", "variant", "enst"),
                       names_from = celltype,
                       values_from = value)

  # return
  names(vxt) <- paste0("vxt_", names(vxt))
  return(vxt)
}

# culled annotations: #
# # variant-TSS inverse distance score method
# vxt$inv_distance <- distance %>%
#   dplyr::transmute(cs, variant, enst,
#                    value = inv_distance)
# 
# # exonic (coding) variants
# vxt$exon <- variants %>%
#   # Get variants within exons
#   bed_intersect_left(
#     exons, .,
#     keepBcoords = F, keepBmetadata = T) %>%
#   dplyr::transmute(cs, variant, enst)
# # intronic variants
# vxt$intron <- variants %>%
#   # Get variants within introns
#   bed_intersect_left(
#     introns, .,
#     keepBcoords = F, keepBmetadata = T) %>%
#   dplyr::transmute(cs, variant, enst)
# # variant-TSS distance rank method
# vxt$inv_distance_rank <- distance %>%
#   dplyr::transmute(cs, variant, enst,
#                    value = inv_distance_rank)
# # closest variant-TSS method
# vxt$closest <- distance %>%
#   dplyr::filter(inv_distance_rank == 1) %>%
#   dplyr::transmute(cs, variant, enst)
# # HiChIP binary
# vxt$HiChIP_binary <- vxt$HiChIP_scores %>%
#   dplyr::mutate(., dplyr::across(where(is.numeric), ~ dplyr::case_when(is.na(.) ~ 0, TRUE ~ 1)))
# # exonic or inv distance rank
# vxt$exon_or_inv_distance_rank <- vxt$exon_or_inv_distance %>%
#   dplyr::group_by(cs, variant) %>%
#   dplyr::mutate(value = 1 / rank(-value, ties.method = "min"))
# # sum of best vxt features
# vxt$exon_or_inv_distance_plus_HiChIP_plus_TADs