# Hadley: "All functions in a package should be related to a single problem
# (or a set of closely related problems). Any functions not related to that
# purpose should not be exported. For example, most of my packages have a
# utils.R file that contains many small functions that are useful for me, but
# aren’t part of the core purpose of those packages. I never export these
# functions."
# These simple, internal functions do not need roxygen documentation and will
# not be user-facing.

# Not %in% function
`%ni%` <- Negate('%in%')

# Read in data as tibbles for better printing
#' @importFrom utils read.table
read_tibble <- function(file, header = FALSE) {
  utils::read.delim(file,
                    header = header,
                    stringsAsFactors = FALSE,
                    quote = "") %>%
    dplyr::as_tibble()
}

# Write tibbles to TSV
#' @importFrom utils write.table
write_tibble <- function(tibble, filename) {
  write.table(tibble, filename,
              quote = F, row.names = F, sep = "\t")
}

# Find matches to any of a vector of character strings (for filtering)
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}

# Count the number of distinct genes that meet a set of conditions
condition_n_genes <- function(df, ...){
  df %>% dplyr::filter(...) %>% dplyr::pull(symbol) %>% dplyr::n_distinct()
}

# Count the number of distinct gene-x-CS pairs that meet a set of conditions
condition_n_gene_x_cs_pairs <- function(df, ...){
  df %>% dplyr::ungroup() %>% dplyr::filter(...) %>% dplyr::distinct(symbol, cs) %>% nrow
}
