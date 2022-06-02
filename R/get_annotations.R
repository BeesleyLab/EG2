get_annotations <- function(celltypes,
                            celltype_of_interest,
                            tissue_of_interest,
                            metadata,
                            enriched,
                            H3K27ac,
                            expression,
                            expressed,
                            HiChIP,
                            TADs){
  
  if (!is.null(tissue_of_interest)) {
    
    cat("Generating predictions for user-provided tissue(s) of interest ", paste(tissue_of_interest, collapse = " + "), ".\n")
    selected_celltypes <- metadata %>% dplyr::filter(tissue %in% tissue_of_interest)
    
  } else if (!is.null(celltype_of_interest)) {
    
    cat("Generating predictions for user-provided celltype(s) of interest ", paste(celltype_of_interest, collapse = " +"), ".\n")
    selected_celltypes <- metadata %>% dplyr::filter(celltype %in% celltype_of_interest)
    
  } else if (celltypes == "all_celltypes") {
    
    cat("Generating predictions for all available tissues (` celltypes =", celltypes, "`).\n")
    selected_celltypes <- metadata
    cat("All tissues: ") ; enriched$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    cat("All celltypes: ") ; enriched$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    
  } else if (celltypes == "enriched_tissues") {
    
    cat("Generating predictions for enriched tissues (` celltypes =", celltypes, "`).\n")
    selected_celltypes <- enriched
    
  }
    
  # subset annotations accordingly
  annotations <- list(
    selected_celltypes = selected_celltypes,
    H3K27ac = H3K27ac %>% purrr::map(~ dplyr::select(.x, DHS, dplyr::any_of(selected_celltypes$celltype))),
    expression = expression %>% purrr::map(~ dplyr::select(.x, ensg, dplyr::any_of(selected_celltypes$celltype))),
    expressed = expressed %>% dplyr::select(ensg, dplyr::any_of(selected_celltypes$celltype)),
    HiChIP = HiChIP[names(HiChIP)[names(HiChIP) %in% selected_celltypes$celltype]],
    TADs = TADs[names(TADs) %in% selected_celltypes$celltype]
  )
  return(annotations)
  
}
