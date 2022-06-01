check_arguments <- function(metadata,
                            variants_file,
                            known_genes_file,
                            reference_panels_dir,
                            celltype_of_interest,
                            tissue_of_interest,
                            celltypes,
                            do_performance,
                            do_XGBoost){

  if (is.null(variants_file)) {
    stop("Must provide the path to a file of trait variants (`variants_file = path/to/variants_file!") }
  if (is.null(known_genes_file) & do_performance) {
    stop("`do_performance = TRUE` but no known_genes_file provided! Must provide the path to a file of known genes for performance analysis (`known_genes_file = path/to/known_genes_file`).") }
  if (is.null(reference_panels_dir)) {
    stop("Must provide the path to the accompanying reference panels directory (`reference_panels_dir = path/to/EG2_data/`)! Reference data can be downloaded from OSF. Visit https://github.com/BeesleyLab/EG2 for instructions.") }
  if (!is.null(celltype_of_interest)) {
    coi_annotations <- metadata %>% dplyr::filter(celltype == celltype_of_interest) %>% dplyr::pull(object) %>% unique
    if (celltype_of_interest %ni% metadata$celltype) {
      stop("Provided celltype_of_interest '", celltype_of_interest, "' is not represented in the available data. Must be one of...\n",
           paste(unique(metadata$celltype), collapse = ", "))
    } else if (length(setdiff(metadata$object, coi_annotations)) > 0) {
      stop("Provided celltype_of_interest '", celltype_of_interest, "' does not have a full panel of annotations available.\nMissing annotation(s) = ",
           paste(setdiff(metadata$object, coi_annotations), collapse = ", "),
           "\nEither re-run at tissue-level (tissue_of_interest = ", unique(metadata[metadata$celltype == celltype_of_interest,]$tissue),
           ") or generate celltype-level annotations for '", celltype_of_interest, "' in the reference data.") } }
  if (!is.null(tissue_of_interest)) {
    if (tissue_of_interest %ni% metadata$tissue) {
      stop("Provided tissue_of_interest '", tissue_of_interest, "' is not represented in the available data. Must be one of...\n", paste(unique(metadata$tissue), collapse = ", ")) } }
  allowed_celltypes <- c("all_celltypes", "enriched_tissues", "enriched_celltypes", "every_tissue")
  if (celltypes %ni% allowed_celltypes) {
    stop("Provided celltypes argument '", celltypes, "' is not acceptable. Must be one of...\n", paste(allowed_celltypes, collapse = ", ")) }

}
