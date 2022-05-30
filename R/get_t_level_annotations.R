get_t_level_annotations <- function(TSSs,
                                    DHSs,
                                    enriched){

  t <- list()

  # intersect with H3K27ac
  t <- t %>%
    intersect_H3K27ac(query = TSSs %>% dplyr::select(chrom:end, enst),
                      DHSs,
                      H3K27ac = enriched$H3K27ac,
                      enst)

  # annotate open variants (in DHSs)
  t$DHSs <- DHSs %>%
    bed_intersect_left(
      TSSs, .,
      keepBcoords = F,
      keepBmetadata = F) %>%
    dplyr::distinct(enst)

  # return
  names(t) <- paste0("t_", names(t))
  return(t)
}
