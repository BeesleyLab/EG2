get_g_level_annotations <- function(vxt_master,
                                    annotations){

  g <- list()

  # gene TPM signal/spec bins
  g <- annotations$expression %>%
    purrr::map( ~ .x %>% dplyr::inner_join(vxt_master %>% dplyr::distinct(ensg), by = "ensg"))
  names(g) <- paste0("expression_", names(g))

  # gene expressed
  g$expressed <- annotations$expressed %>%
    dplyr::inner_join(vxt_master %>% dplyr::distinct(ensg), by = "ensg")

  # return
  names(g) <- paste0("g_", names(g))
  return(g)
}
