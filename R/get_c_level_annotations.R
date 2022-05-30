get_c_level_annotations <- function(variants) {

  c <- list()

  # n variants
  c$n_variants <- variants %>%
    dplyr::group_by(cs) %>%
    dplyr::summarise(value = dplyr::n_distinct(variant))

  # inverse n variants
  c$inv_n_variants <- c$n_variants %>%
    dplyr::transmute(cs,
                     value = 1/value)

  # return
  names(c) <- paste0("c_", names(c))
  return(c)
}
