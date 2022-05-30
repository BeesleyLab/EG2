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

  # for testing: # ratio_cutoff = 1 ; p_value_cutoff = 0.05

  enriched <- list()

  if (!is.null(tissue_of_interest)) {

    cat("Treating user-provided tissue of interest '", tissue_of_interest, "' as the enriched tissue.\n")
    # user-provided tissue
    enriched[["celltypes"]] <- metadata %>%
      dplyr::filter(tissue == tissue_of_interest)

  } else if (!is.null(celltype_of_interest)) {

    cat("Treating user-provided celltype of interest '", celltype_of_interest, "' as the enriched celltype.\n")
    # user-provided celltype
    enriched[["celltypes"]] <- metadata %>%
      dplyr::filter(celltype == celltype_of_interest)

  } else if (celltypes == "all_celltypes") {

    cat("Applying annotations from all available cell types (`celltypes` = 'all_celltypes').\n")
    # Include all celltypes
    enriched[["celltypes"]] <- metadata
    cat("All tissues: ") ; enriched$celltypes$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    cat("All celltypes: ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')

  } else {

    cat("Performing enrichment analysis to find enriched celltype(s).\n")
    # No user-provided tissue, determine via enrichment analysis
    # Fisher enrichment test
    # Intersect SNPs with DHS sites and compute the mean H3K27ac specificity rank of the intersected DHS.
    # Assuming SNPs overlap sites randomly, compute significance of mean rank of sites where SNPs overlap based on deivation from uniform distribution.

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
        unif_mean_rank = (N + 1)/2,
        unif_variance = sqrt((N^2 - 1)/(12 * n)),
        # deviation from uniform distribution = enrichment
        ratio = obs_mean_rank / unif_mean_rank,
        p_value = pnorm(obs_mean_rank, unif_mean_rank, unif_variance, lower.tail = F),
        p_value_adjust = p_value %>% p.adjust,
        pass = ((p_value < p_value_cutoff) & (ratio > ratio_cutoff))) %>%
      dplyr::arrange(p_value)
    write_tibble(enrichment, out$tissue_enrichments)

    # Enriched celltypes/tissues
    enriched[["celltypes"]] <- enrichment %>%
      # Filter to celltypes that pass filters
      dplyr::filter(pass)  %>%
      # Get enriched celltype metadata
      {dplyr::filter(metadata, celltype %in% .$celltype)} %>%
      # Get all samples in the enriched tissue, or only the enriched celltype
      {dplyr::filter(metadata,
                      (celltypes == "enriched_tissues" & tissue %in% .$tissue ) |
                      (celltypes == "enriched_celltypes" & celltype %in% .$celltype))}

    # Error message if no cell types were enriched
    if(nrow(enriched$celltypes)==0){
      stop("No enriched cell types found! Enrichment analysis saved to ", out$tissue_enrichments)}

    # Enriched celltype(s)/tissue(s)
    if(celltypes == "enriched_tissues"){
      cat("Enriched tissue(s): ") ; enriched$celltypes$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
      cat("Celltype(s) in enriched tissue(s): ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    } else {
      cat("Enriched celltype(s): ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')

      full_panel <- metadata$object
      enriched_panel <- dplyr::filter(metadata, celltype %in% enriched$celltypes$celltype)$object
      panel_gaps <- setdiff(full_panel, enriched_panel)
      if(length(panel_gaps) > 0){
      stop("\nOption `celltypes = enriched_celltypes` was given, but the enriched celltype(s) do not constitute a full reference panel.",
           "\nEnriched celltype(s): ", enriched$celltypes$celltype  %>% unique %>% paste(collapse = ", "),
           "\nMissing annotations: ", paste(panel_gaps, collapse = ", "),
           "\nMake sure the enriched celltype(s) have coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table.",
           "\nEither run again with `celltypes = enriched_tissues` or lower the enrichment threshold `p_value_cutoff = ",p_value_cutoff,"`")
      }
    }

  }

  # Subset annotations
  ## H3K27ac
  enriched$H3K27ac <- H3K27ac %>%
    purrr::map( ~ dplyr::select(., DHS,
                                dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "H3K27ac"])))
  ## HiChIP
  enriched$HiChIP <- HiChIP[names(HiChIP)[names(HiChIP) %in% enriched$celltypes$celltype[enriched$celltypes$object == "HiChIP"]]]
  ## expression
  enriched$expression <- expression %>%
    purrr::map( ~ dplyr::select(., ensg,
                                dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "expression"])))
  ## expressed
  enriched$expressed <- expressed %>%
    dplyr::select(ensg, dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "expression"]))
  ## TADs
  enriched$TADs <- TADs[names(TADs) %in% enriched$celltypes$celltype]

  return(enriched)
}

