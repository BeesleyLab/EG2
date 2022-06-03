#' check that all symbols are in the GENCODE data
#'
#' @param known_genes vector of known gene symbols
#' @param known_genes_file name of the known_genes file
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
check_known_genes <- function(known_genes, known_genes_file){
  known_genes_NotFound <- known_genes[known_genes %ni% TSSs$symbol]
  known_genes_NonCoding <- dplyr::filter(TSSs,
                                     symbol %in% known_genes[known_genes %in% TSSs$symbol],
                                     ensg %ni% pcENSGs
                                     )$symbol
  if(length(c(known_genes_NotFound, known_genes_NonCoding)) > 0){
    message(dplyr::n_distinct(known_genes_NotFound),
            " provided known gene(s) from '", basename(known_genes_file), "' are not valid GENCODE gene symbols and therefore cannot be used.",
            "\nUnknown genes: ", paste(known_genes_NotFound, sep="", collapse=", "))
    message(dplyr::n_distinct(known_genes_NonCoding),
            " provided known gene(s) from '", basename(known_genes_file), "' are GENCODE genes but are not protein-coding and therefore cannot be used.",
            "\nNon-coding genes: ", paste(known_genes_NonCoding, sep="", collapse=", "))
    message(length(known_genes[known_genes %ni% c(known_genes_NotFound, known_genes_NonCoding)]),
            " known genes left after removing unknown/non-coding.")
  }

  TSSs %>%
    dplyr::filter(symbol %in% known_genes, ensg %in% pcENSGs)
}
