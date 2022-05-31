get_vxt_master <- function(variants, 
                           TSSs,
                           max_variant_to_gene_distance){
  # The transcript-x-variant universe (masterlist of all possible transcript x variant pairs < max_variant_to_gene_distance apart)
  vxt_master <- variants %>%
    # all genes within 2Mb
    valr::bed_slop(both = max_variant_to_gene_distance,
                   genome = ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., TSSs,
                        suffix = c("", ".TSS")) %>%
    # restore coords
    dplyr::select(-c(start, end)) %>%
    dplyr::left_join(variants, by = c("chrom", "variant", "cs")) %>%
    dplyr::transmute(chrom,
                     start.variant = start,
                     end.variant = end,
                     start.TSS, end.TSS,
                     distance = abs(end - end.TSS),
                     variant,
                     pair = paste0(cs, ":", variant, "|", enst.TSS),
                     cs,
                     enst = enst.TSS,
                     symbol = symbol.TSS,
                     ensg = ensg.TSS) %>%
    dplyr::distinct()
  return(vxt_master)
}