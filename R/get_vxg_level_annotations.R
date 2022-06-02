get_vxg_level_annotations <- function(variants,
                                      vxt_master,
                                      enriched){

  vxg <- list()

  # coding_mutations #
  intersect_coding_mutations <- function(variants, coding_mutations_df, ...){
    variants %>%
      dplyr::inner_join(coding_mutations_df, by = c("chrom" = "chrom", "end" = "position")) %>%
      tidyr::separate_rows(ensgs) %>%
      dplyr::filter(...) %>%
      dplyr::transmute(cs, variant, ensg = ensgs)
  }

  # missense variants
  vxg$missense <- variants %>% intersect_coding_mutations(missense, score > 0.5)

  # nonsense variants
  vxg$nonsense <- variants %>% intersect_coding_mutations(nonsense)

  # splice site variants
  vxg$splicesite <- variants %>% intersect_coding_mutations(splicesite)

  # prefix names
  names(vxg) <- paste0("vxg_", names(vxg))
  return(vxg)
}

# culled annotations: #
# # variant-gene distance
# # (between the variant and all of the gene's transcripts' TSSs)
# distance <- vxt_master %>%
#   dplyr::group_by(cs, variant, symbol) %>%
#   dplyr::summarise(distance = min(distance)) %>%
#   dplyr::group_by(variant) %>%
#   dplyr::transmute(cs, variant, symbol,
#                    inv_distance = dplyr::case_when(distance == 0 ~ 1,
#                                                    TRUE ~ 1 / distance),
#                    inv_distance_rank = 1 / rank(distance, ties.method = "min")) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct()
# vxg$closest <- distance %>%
#   dplyr::filter(inv_distance_rank == 1) %>%
#   dplyr::transmute(cs, variant, symbol)
# vxg$inv_distance <- distance %>%
#   dplyr::transmute(cs, variant, symbol,
#                    value = inv_distance)
# vxg$inv_distance_rank <- distance %>%
#   dplyr::transmute(cs, variant, symbol,
#                    value = inv_distance_rank)