get_tissue_enrichment <- function(variants,
                                  DHSs,
                                  H3K27ac_specificity_ranked,
                                  metadata,
                                  out,
                                  ratio_cutoff = 1,
                                  p_value_cutoff = 0.05){
  # for testing: # ratio_cutoff = 1 ; p_value_cutoff = 0.05
  
  cat("Performing enrichment analysis to find enriched celltype(s).\n")
  # Fisher enrichment test
  # Intersect SNPs with DHS sites and compute the mean H3K27ac specificity rank of the intersected DHS.
  # Assuming SNPs overlap sites randomly, compute significance of mean rank of sites where SNPs overlap based on deviation from uniform distribution.
  
  # intersect DHSs
  intersected_DHSs <- bed_intersect_left(variants, DHSs, keepBcoords = F) %>%
  {dplyr::filter(H3K27ac_specificity_ranked, DHS %in% .$DHS)}
  enrichment <- intersected_DHSs %>%
    # mean specificity rank per celltype
    dplyr::summarise(dplyr::across(.cols = where(is.numeric),
                                   .fns = ~ mean(.x))) %>%
    tidyr::pivot_longer(everything(), names_to = "celltype", values_to = "obs_mean_rank") %>%
    dplyr::left_join(metadata %>% dplyr::distinct(celltype, tissue), by = "celltype") %>%
    dplyr::transmute(
      celltype,
      tissue,
      # uniform distribution parameters
      N = dplyr::n_distinct(DHSs$DHS),
      n = dplyr::n_distinct(intersected_DHSs$DHS),
      obs_mean_rank,
      unif_mean_rank = (N + 1)/2,
      unif_variance = sqrt((N^2 - 1)/(12 * n)),
      # deviation from uniform distribution = enrichment
      ratio = obs_mean_rank / unif_mean_rank,
      p_value = pnorm(obs_mean_rank, unif_mean_rank, unif_variance, lower.tail = F),
      p_value_adjust = p_value %>% p.adjust,
      pass = ((p_value < p_value_cutoff) & (ratio > ratio_cutoff))) %>%
    dplyr::arrange(p_value)
  
  enriched <- enrichment %>%
    # Filter to celltypes that pass filters
    dplyr::filter(pass)  %>%
    # Get enriched tissue metadata
    {dplyr::filter(metadata, tissue %in% .$tissue)}
  
  # enrichment message
  if(nrow(enriched) > 0){
    cat("Enriched tissue(s): ") ; enriched$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    cat("Celltype(s) in enriched tissue(s): ") ; enriched$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
  } else {
    message("No enriched cell types found!\n")
  }
  
  write_tibble(enrichment, out$tissue_enrichments)
  cat("Enrichment analysis saved to", out$tissue_enrichments, ".\n")
  
  return(enriched)
}

