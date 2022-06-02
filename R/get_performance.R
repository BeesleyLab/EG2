get_performance <- function(annotations, vxt_master, known_genes, pcENSGs, max_n_known_genes_per_CS){

  performance <- list()

  # get all testable CS-gene pairs
  testable <- vxt_master %>%
    # only test protein-coding target predictions against known_genes (assumes all known_genes are protein-coding)
    dplyr::filter(ensg %in% pcENSGs) %>%
    # add known_genes
    dplyr::mutate(known_gene = symbol %in% known_genes$symbol) %>%
    # only test predictions in CSs with a known gene within max prediction distance for performance evaluation
    get_testable(max_n_known_genes_per_CS) %>%
    dplyr::distinct(cs, symbol, known_gene)

  # score, max
  pred <- annotations %>%
    dplyr::select(-dplyr::any_of(c("chrom", "start", "end"))) %>%
    # get testable pairs
    dplyr::right_join(testable, by = c("cs", "symbol")) %>%
    # get maximum score per cxg-x-method
    dplyr::group_by(cs, symbol, known_gene) %>%
    # gather prediction methods
    tidyr::pivot_longer(
      where(is.numeric),
      names_to = "prediction_method",
      values_to = "score"
    ) %>%
    # get maximum score per cxg-x-method
    dplyr::group_by(prediction_method, cs, symbol, known_gene) %>%
    dplyr::summarise(score = max(score)) %>%
    # get maximum score per c-x-method
    dplyr::group_by(prediction_method, cs) %>%
    dplyr::mutate(max = as.numeric(score == max(score) & score > 0)) %>%
    # gather prediction types
    tidyr::pivot_longer(
      c(score, max),
      names_to = "prediction_type",
      values_to = "prediction"
    ) %>%
    dplyr::ungroup()

  # format for PR function input
  PR_in <- pred %>%
    dplyr::select(prediction_type, prediction_method, prediction, known_gene) %>%
    dplyr::group_by(prediction_method, prediction_type) %>%
    # refactor known gene predictions for PR function
    dplyr::mutate(known_gene = ifelse(known_gene, "positive", "negative") %>% factor(c("positive", "negative")))

  # calculate PR curve
  performance$PR <- PR_in %>%
    yardstick::pr_curve(known_gene, prediction) %>%
    dplyr::ungroup()

  # calculate AUPRC
  performance$PR_AUC <- PR_in %>%
    dplyr::filter(prediction_type == "score") %>%
    yardstick::pr_auc(known_gene, prediction) %>%
    dplyr::select(prediction_method,
                  prediction_type,
                  score_PR_AUC = .estimate) %>%
    dplyr::ungroup()

  # get summary statistics (various performance metrics)
  performance$summary <- pred %>%
    # max only (score is not binary, cannot be summarised)
    dplyr::filter(prediction_type == "max") %>%
    dplyr::mutate(prediction = as.logical(prediction)) %>%
    dplyr::group_by(prediction_method, prediction_type) %>%
    # max - confusion matrix
    dplyr::group_modify(
      ~ data.frame(True = .x %>% condition_n_gene_x_cs_pairs(known_gene),
                   Positive = .x %>% condition_n_gene_x_cs_pairs(prediction),
                   TP = .x %>% condition_n_gene_x_cs_pairs(prediction & known_gene),
                   FP = .x %>% condition_n_gene_x_cs_pairs(prediction & !known_gene),
                   TN = .x %>% condition_n_gene_x_cs_pairs(!prediction & !known_gene),
                   FN = .x %>% condition_n_gene_x_cs_pairs(!prediction & known_gene),
                   n_known_genes = .x %>% condition_n_genes(known_gene))
    ) %>%
    dplyr::rowwise() %>%
    # max - performance metrics
    dplyr::mutate(p_value = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FP),
                  Recall = TP / (TP + FN),
                  Sensitivity = TP / (TP + FN),
                  Specificity = TN / (TN + FP),
                  F_score = 2 * ((Precision * Recall) / (Precision + Recall))) %>%
    # score - PR_AUC
    dplyr::left_join(performance$PR_AUC %>%
                       dplyr::select(prediction_method, score_PR_AUC),
                     by = "prediction_method") %>%
    dplyr::ungroup()

  # arrange methods by F score
  performance <- performance %>%
    purrr::map(
      ~ dplyr::left_join(
        performance$summary %>%
          dplyr::arrange(desc(F_score)) %>%
          dplyr::select(prediction_method),
        .x)
    )

  return(performance)

}





