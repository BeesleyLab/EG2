get_weights <- function(weights_file,
                        master){
  
  if(is.null(weights_file)){
    weights_name <- "`default_weights`"
    weights <- default_weights
  } else {
    weights_name <- weights_file
    cat("  > Importing custom weights from ", weights_file, "...\n")
    weights <- read_tibble(weights_file, header = T) %>%
      dplyr::select(annotation, weight) %>%
      tibble::column_to_rownames("annotation") %>%
      as.matrix
  }
  missing_weights <- setdiff(names(master), rownames(weights))
  if(length(missing_weights) > 0){
    message("Annotation(s)\n > ", paste(missing_weights, collapse = "\n > "),
            "\ndo not have a weight in ", weights_name, ". Weighting as 0.")
    zero_weights <- matrix(rep(0, length(missing_weights)))
    rownames(zero_weights) <- missing_weights
    weights <- rbind(weights, zero_weights)
  }
  missing_annotations <- setdiff(rownames(weights), names(master))
  if(length(missing_annotations) > 0){
    message("Weight(s)\n > ", paste(missing_annotations, collapse = "\n > "),
            "\ndo not match the names of any annotations and will not be used. ")
  }
  return(weights)
  
}