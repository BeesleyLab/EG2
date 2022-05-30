get_testable <- function(df, max_n_known_genes_per_CS = Inf){

  testable <- df  %>%
      dplyr::group_by(cs) %>%
      dplyr::filter(known_gene) %>%
      dplyr::distinct(cs, symbol) %>%
      dplyr::filter(dplyr::n_distinct(symbol) > 0,
                    dplyr::n_distinct(symbol) <= max_n_known_genes_per_CS)
  df %>%
    dplyr::filter(cs %in% testable$cs)

}
