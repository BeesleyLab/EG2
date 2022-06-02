get_v_level_annotations <- function(variants,
                                    H3K27ac,
                                    annotations,
                                    vxt_master,
                                    DHSs){

  v <- list()

  # intersect with DHS binnings
  v <- v %>%
    intersect_H3K27ac(query = variants,
                      DHSs,
                      H3K27ac = annotations$H3K27ac,
                      variant)

  # annotate open variants (in DHSs)
  v$DHSs <- DHSs %>%
    bed_intersect_left(
      variants, .,
      keepBcoords = F,
      keepBmetadata = F) %>%
    dplyr::distinct(variant)

  # calculate n genes near each variant
  v <- vxt_master %>%
    dplyr::group_by(variant) %>%
    dplyr::summarise(inv_n_genes = 1/dplyr::n_distinct(symbol),
                     inv_n_transcripts = 1/dplyr::n_distinct(enst)) %>%
    tidyr::pivot_longer(c(inv_n_genes, inv_n_transcripts),
                        names_to = "split",
                        values_to = "value") %>%
    split(as.factor(.$split)) %>%
    purrr::map(~ dplyr::select(., -split)) %>%
    c(., v)

  # return
  names(v) <- paste0("v_", names(v))
  return(v)
}
