get_cxt_level_annotations <- function(annotations,
                                      vxt,
                                      variants) {

  multiHiChIP_indep <- annotations$HiChIP %>% names %>%
    sapply(function(ct){
      # Intersect with the HiChIP data
      annotations$HiChIP[[ct]] %>%
        intersect_BEDPE(
          # ! For mutually exclusive intersection with ranges, make variant intervals 1bp long, equal to the end position
          SNPend = variants %>% dplyr::mutate(start = end),
          TSSend = TSSs,
          bedpe = .) %>%
        dplyr::transmute(cs, variant, enst,
                         InteractionID,
                         value = score,
                         celltype = ct)
    }, simplify = F, USE.NAMES = T) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(celltype, cs, enst) %>%
    dplyr::distinct(InteractionID, value) %>%
    dplyr::summarise(n = dplyr::n_distinct(InteractionID),
                     sum = sum(value),
                     binary = as.numeric(n > 1)) %>%
    dplyr::filter(binary == 1) %>%
    tidyr::pivot_longer(cols = c(n, sum, binary),
                        names_to = "annotation",
                        values_to = "value")
  multiHiChIP_indep <- multiHiChIP_indep %>%
    split(multiHiChIP_indep$annotation) %>%
    purrr::map(~ .x %>% tidyr::pivot_wider(id_cols = c(cs, enst),
                                           names_from = celltype,
                                           values_from = value))
  if(length(multiHiChIP_indep)>0){
    names(multiHiChIP_indep) <- names(multiHiChIP_indep) %>% paste0("_multiHiChIP_indep")
  }

  # multiHiChIP statistics within each gene-x-cs-x-experiment combination
  multiHiChIP <- vxt[grepl("HiChIP", names(vxt)) & !grepl("binary", names(vxt))] %>%
    purrr::map(~ dplyr::group_by(., cs, enst))

  # count number of loops between gene and CS
  n_multi <- multiHiChIP %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(!is.na(.x)))))
  names(n_multi) <- names(n_multi) %>% gsub("vxt_", "n_multi", .)

  # binary - multiHiChIP y/n
  binary_multi <- n_multi %>%
    purrr::map(~ dplyr::mutate(., dplyr::across(where(is.numeric), ~ as.numeric(.x > 1))))
  names(binary_multi) <- names(binary_multi) %>% gsub("multiHiChIP", "multiHiChIP_binary", .)

  # sum of loop values between gene and CS
  sum_multi <- multiHiChIP %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(.x, na.rm = T))))
  names(sum_multi) <- names(sum_multi) %>% gsub("vxt_", "sum_multi", .)

  multiHiChIP <- c(n_multi,
                   binary_multi,
                   sum_multi)

  # cxt
  cxt <- c(multiHiChIP_indep,
           multiHiChIP)

  # return
  names(cxt) <- paste0("cxt_", names(cxt))
  return(cxt)
}
