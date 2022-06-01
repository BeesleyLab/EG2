# default annotation weights
library(EG2)
default_weights_file <- system.file("example_data", "default_weights.tsv", package = "EG2")
# default_weights_file <- "/working/lab_jonathb/alexandT/EG2/inst/example_data/default_weights.tsv"
default_weights <- read.delim(default_weights_file, header = T) %>%
  dplyr::select(annotation, weight) %>%
  tibble::column_to_rownames("annotation") %>%
  as.matrix
usethis::use_data(default_weights, overwrite = T)
