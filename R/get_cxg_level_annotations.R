get_cxg_level_annotations <- function(vxt_master){

  cxg <- list()

  # CS-gene distance
  # (among all of the CS's variants and all of the gene's transcripts' TSSs)
  distance <- vxt_master %>%
    dplyr::group_by(cs, symbol) %>%
    dplyr::filter(distance == min(distance)) %>%
    dplyr::group_by(cs) %>%
    dplyr::transmute(cs, symbol,
                     inv_distance = dplyr::case_when(distance == 0 ~ 1,
                                                     TRUE ~ 1 / distance),
                     inv_distance_rank = 1 / rank(distance, ties.method = "min")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  cxg$inv_distance <- distance %>%
    dplyr::transmute(cs, symbol,
                     value = inv_distance)

  cxg$inv_distance_rank <- distance %>%
    dplyr::transmute(cs, symbol,
                     value = inv_distance_rank)

  # prefix names
  names(cxg) <- paste0("cxg_", names(cxg))
  return(cxg)
}
