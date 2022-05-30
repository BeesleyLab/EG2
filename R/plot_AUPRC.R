plot_AUPRC <- function(performance, ...){
  performance$PR %>%
    dplyr::select(prediction_method, prediction_type, PR_AUC) %>%
    dplyr::filter(prediction_type == "max") %>%
    dplyr::distinct() %>%
    dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .)) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, PR_AUC),
                                 y = PR_AUC,
                                 ...)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Predictor",
                  y = "PR AUC") +
    ggsci::scale_fill_igv() +
    ggplot2::coord_flip()
}
