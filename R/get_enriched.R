get_enriched <- function(variants,
                         DHSs,
                         H3K27ac_specificity_ranked,
                         H3K27ac,
                         expression,
                         expressed,
                         HiChIP,
                         TADs,
                         metadata,
                         out,
                         celltype_of_interest,
                         tissue_of_interest,
                         celltypes,
                         ratio_cutoff = 1,
                         p_value_cutoff = 0.05){

  
  enriched <- list()

  if (!is.null(tissue_of_interest)) {

    cat("Treating user-provided tissue(s) of interest ", paste(tissue_of_interest, collapse = " + "), " as the enriched tissue(s).\n")
    # user-provided tissue
    enriched[["celltypes"]] <- metadata %>%
      dplyr::filter(tissue %in% tissue_of_interest)

  } else if (!is.null(celltype_of_interest)) {

    cat("Treating user-provided celltype(s) of interest ", paste(celltype_of_interest, collapse = " +"), " as the enriched celltype(s).\n")
    # user-provided celltype
    enriched[["celltypes"]] <- metadata %>%
      dplyr::filter(celltype%in% celltype_of_interest)

  } else if (celltypes %in% c("every_tissue", "all_celltypes")) {
      
    cat("Applying annotations from all available cell types (` celltypes =", celltypes, "`).\n")
    # every_tissue / all_celltypes
    enriched[["celltypes"]] <- metadata
    cat("All tissues: ") ; enriched$celltypes$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    cat("All celltypes: ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
      
  } else {
      
    # enriched_tissues / enriched_celltypes
    
    
    
  }

  # subset annotations
  ## H3K27ac
  enriched$H3K27ac <- H3K27ac %>%
    purrr::map( ~ dplyr::select(., DHS, dplyr::any_of(enriched$celltypes$celltype)))
  
  ## HiChIP
  enriched$HiChIP <- HiChIP[names(HiChIP)[names(HiChIP) %in% enriched$celltypes$celltype]]
  ## expression
  enriched$expression <- expression %>%
    purrr::map( ~ dplyr::select(., ensg, dplyr::any_of(enriched$celltypes$celltype)))
  ## expressed
  enriched$expressed <- expressed %>%
    dplyr::select(ensg, dplyr::any_of(enriched$celltypes$celltype))
  ## TADs
  enriched$TADs <- TADs[names(TADs) %in% enriched$celltypes$celltype]

  return(enriched)
}

