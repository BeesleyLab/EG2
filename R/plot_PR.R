plot_PR <- function(performance, ...){
  performance$PR %>%
    dplyr::filter(.threshold %ni% c(0, Inf, -Inf),
                  recall != 1,
                  precision != 0) %>%
    dplyr::group_by(prediction_method, prediction_type) %>%
    dplyr::mutate(line = grepl("score", prediction_type),
                  point = !grepl("score", prediction_type)) %>%
    ggplot2::ggplot(ggplot2::aes(x = recall,
                                 y = precision,
                                 ...,
                                 label = prediction_method)) +
    ggplot2::geom_line(data = . %>% dplyr::filter(line),
                       ggplot2::aes(linetype = prediction_type)) +
    ggplot2::geom_point(data = . %>% dplyr::filter(point),
                        ggplot2::aes(shape = prediction_type)) +
    ggsci::scale_color_igv() +
    ggplot2::xlim(0,1) +
    ggplot2::ylim(0,1) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = paste0("recall (n = ", unique(performance$summary$True), ")"))
}
