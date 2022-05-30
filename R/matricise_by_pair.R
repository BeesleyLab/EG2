matricise_by_pair <- function(df,
                              vxt_master){
  # add value = 1 if there is no value (numeric) column
  value_cols <- df %>% dplyr::ungroup() %>% dplyr::select(where(is.numeric)) %>% names
  if (length(value_cols) == 0) { value_cols <- "value" ; df$value <- 1}

  mat <- vxt_master %>%
    # arrange rows in same order as vxt_master
    dplyr::left_join(df, by = intersect(names(vxt_master), names(df))) %>%
    dplyr::select(dplyr::all_of(value_cols)) %>%
    as.matrix
  mat[is.na(mat)] <- 0

  return(mat)
}
