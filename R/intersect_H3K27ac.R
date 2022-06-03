#' Intersect with the DHS annotations
#'
#' Intersect a query BED with the DHS annotations.
#'
#' @param query A query bed
#' @param H3K27ac A BED tibble (the H3K27ac-in-DHSs annotations) to be intersected by the query bed
#' @param ... query columns to be retained in the output
#'
#' @return A tibble of intersected columns, with .query and .annotation suffixes
intersect_H3K27ac <- function(l,
                              query,
                              DHSs,
                              H3K27ac,
                              ...){

  # intersect with query bed
  intersected_DHSs <- bed_intersect_left(
    bedA = query,
    bedB = DHSs,
    keepBcoords = F
  )
  l <- H3K27ac %>%
    lapply(function(x){
      score_cols <- setdiff(colnames(x), "DHS")
      intersected_DHSs %>%
        dplyr::left_join(x, by = "DHS") %>%
        # group by query column(s) (unit(s) of the annotation)
        dplyr::group_by(...) %>%
        # for overlapping - get max bin to avoid duplicates
        dplyr::summarise(dplyr::across(score_cols, max)) %>%
        dplyr::ungroup()
    })
  names(l) <- paste0("H3K27ac_", names(l))
  return(l)

}

