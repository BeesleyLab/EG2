#' Predict target genes of fine-mapped variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param trait Optional. The name of the trait of interest.
#' @param out_dir Optional. The output directory in which to save the predictions. Default is "./out/{trait}/{celltypes}/".
#' @param sub_dir Optional. The output subdirectory to be appended to the default path ("./out/{trait}/{celltypes}/{sub_dir}/").
#' @param variants_file A BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#' @param known_genes_file Optional. The file containing a list of trait known gene symbols. If do_performance is TRUE, must provide a known_genes_file.
#' @param reference_panels_dir The directory containing the external, accompanying reference panels data.
#' @param weights_file A file of alternative weights for annotations. Must contain `annotation` and `weight` columns. Default is data/weights.tsv. If an annotation is missing from the weights_file, it will be weighted 0.
#' @param celltype_of_interest Optional. The celltype(s) of interest for the trait. Only annotations in these celltypes will be used to make predictions. Argument(s) must match the names of celltypes in the metadata. Make sure the celltype of interest has coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table.
#' @param tissue_of_interest Optional. The tissue(s) of interest for the trait. Only annotations in these tissues will be used to make predictions.  Argument(s) must match the names of tissues in the metadata.
#' @param max_variant_to_gene_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Measured as the distance between the variant and the gene's TSS. Default is 2Mb. The HiChIP data is also already filtered to 2Mb.
#' @param max_n_known_genes_per_CS In performance analysis, the maximum number of known genes within max_variant_to_gene_distance of the credible set.
#' @param celltypes Dictates which celltypes' annotations are used. Must be one of c("enriched_celltypes", "enriched_tissues", "all_celltypes"). If "enriched_celltypes", annotations from only the enriched celltype(s) will be used. The enriched celltype(s) must have coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table for this to work. If "enriched_tissues", all annotations from the tissue of the enriched celltype(s) will be used. If "all_celltypes", the enrichment analysis is skipped and annotations from all available celltypes will be used. Default is "enriched_tissues".
#' @param do_performance If TRUE, runs the performance chunk of the script, which measures the performance of the score and each of its constituent annotations in predicting known genes as the targets of nearby variants. Default is FALSE.
#' @param do_XGBoost If TRUE, runs the XGBoost chunk of the script, which generates a model to predict the targets of variants from all available annotations and rates the importance of each annotation. Default is FALSE.
#' @param do_timestamp If TRUE, will save output into a subdirectory timestamped with the data/time of the run.
#' @param HiChIP If you are repeatedly running predict_target_genes, you can load the HiChIP object from the reference_panels_dir into the global environment and pass it to the function to prevent redundant re-loading each call to predict_target_genes.
#' @param H3K27ac If you are repeatedly running predict_target_genes, you can load the H3K27ac object from the reference_panels_dir into the global environment and pass it to the function to prevent redundant re-loading with each call to predict_target_genes.
#' @return The `annotations` tibble with one row per variant-x-transcript pair, one column per annotation and the resulting weighted score for each pair.
#' @export
predict_target_genes <- function(trait = NULL,
                                 out_dir = NULL,
                                 sub_dir = NULL,
                                 variants_file = NULL,
                                 known_genes_file = NULL,
                                 reference_panels_dir = NULL,
                                 weights_file = NULL,
                                 celltype_of_interest = NULL,
                                 tissue_of_interest = NULL,
                                 celltypes = "enriched_tissues",
                                 max_variant_to_gene_distance = 2e6,
                                 max_n_known_genes_per_CS = Inf,
                                 do_performance = T,
                                 do_XGBoost = F,
                                 do_timestamp = F,
                                 HiChIP = NULL,
                                 H3K27ac = NULL) {

  # capture function arguments (do not run when testing internally)
  args <- as.list(environment())[names(as.list(environment())) %ni% c("HiChIP", "H3K27ac")]
  args["HiChIP"] <- list(NULL)
  args["H3K27ac"] <- list(NULL)

  # for testing internally:
  # setwd("/working/lab_jonathb/alexandT/EG2") ; trait="BC_Michailidou2017_FM" ; celltypes = "enriched_celltypes" ; variants_file=paste0("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/",trait,"/variants.bed") ; known_genes_file = paste0("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/",trait,"/known_genes.txt") ; reference_panels_dir = "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/output/" ; weights_file = NULL ; max_variant_to_gene_distance = 2e6 ; max_n_known_genes_per_CS = Inf ; HiChIP = NULL ; H3K27ac = NULL ; celltype_of_interest = NULL ; tissue_of_interest = NULL ; out_dir = NULL ; sub_dir = NULL ; do_scoring = T ; do_performance = T ; do_XGBoost = T ; do_timestamp = F  ; library(devtools) ; load_all()
  # for internally restoring a previous run environment:
  # args <- dget("out/SC3_Chahal2016andSarin2020_LD/SKIN_tissue/arguments_for_predict_target_genes.R") ; list2env(args, envir=.GlobalEnv) ; library(devtools) ; load_all()

  # SETUP ======================================================================================================

  # metadata for all annotations
  metadata <- read_tibble(paste0(reference_panels_dir, "metadata.tsv"), header = T)

  # check arguments
  check_arguments(metadata,
                  variants_file,
                  known_genes_file,
                  reference_panels_dir,
                  celltype_of_interest,
                  tissue_of_interest,
                  celltypes,
                  do_performance,
                  do_XGBoost)

  # define the output files
  celltypes_dir <- {
    if (!is.null(celltype_of_interest)) paste0(paste(celltype_of_interest, collapse = "+"), "_celltype")
    else if (!is.null(tissue_of_interest)) paste0(paste(tissue_of_interest, collapse = "+"), "_tissue")
    else paste0(celltypes)
  }
  out <- list(
    base = "",
    tissue_enrichments = "tissue_enrichments.tsv",
    annotations.rds = "annotations.rds",
    predictions_full = "target_gene_predictions_full.rds",
    predictions_max = "target_gene_predictions_max.rds",
    predictions_max_pc = "target_gene_predictions_max_pc.rds",
    performance = "performance.tsv",
    performance_plot = "performance_plots.pdf",
    args = "arguments_for_predict_target_genes.R",
    XGBoost = "XGBoost_feature_importance.tsv",
    XGBoost_plot = "XGBoost_feature_importance_plot.pdf"
  ) ; if (is.null(out_dir)) {
    out <- "out/" %>%
      paste0(trait, "/", celltypes_dir, "/") %>%
      { if(do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . } %>%
      { if(!is.null(sub_dir)) paste0(., sub_dir, "/") else . } %>%
      { purrr::map(out, function(x) paste0(., x)) }
  } else { out <- out %>% purrr::map(function(x) paste0(out_dir, "/", x)) }
  dir.create(out$base, showWarnings = F, recursive = T)
  cat("Output will be saved to", out$base, "\n")

  # write run arguments to output
  dput(args, file = out$args)

  # ggplot theme settings
  ggplot2::theme_set(
    ggplot2::theme_classic() +
      ggplot2::theme(
        title = ggplot2::element_text(size = 20),
        text  = ggplot2::element_text(size = 18)))

  # import the user-provided variants
  cat("  > Importing variants...\n")
  variants <- import_BED(
    variants_file,
    metadata_cols = c("variant", "cs")) %>%
    dplyr::mutate(variant = as.character(variant),
                  cs = as.character(cs))

  # import the HiChIP data
  if (is.null(HiChIP)) {
    cat("  > Importing HiChIP data...\n")
    HiChIP <- readRDS(paste0(reference_panels_dir, "HiChIP.rds"))
  }

  # import the H3K27ac data
  if (is.null(H3K27ac)) {
    cat("  > Importing H3K27ac data...\n")
    H3K27ac <- readRDS(paste0(reference_panels_dir, "H3K27ac.rds"))
  }
  H3K27ac_specificity_ranked <- readRDS(paste0(reference_panels_dir, "H3K27ac_specificity_rank.rds"))

  # import the expression data
  cat("  > Importing RNA-seq expression data...\n")
  expression <- readRDS(paste0(reference_panels_dir, "expression.rds"))
  expressed <- readRDS(paste0(reference_panels_dir, "expressed.rds"))

  # import DHSs master
  cat("  > Importing DHS data...\n")
  DHSs <- readRDS(paste0(reference_panels_dir, "DHSs.rds"))

  # import the TADs data
  cat("  > Importing TAD data...\n")
  TADs <- readRDS(paste0(reference_panels_dir, "TADs.rds"))

  # 1) CELL TYPE ENRICHMENT ======================================================================================================
  cat("1) Performing cell type enrichment...\n")
  enriched <- get_tissue_enrichment(variants,
                                      DHSs,
                                      H3K27ac_specificity_ranked,
                                      metadata,
                                      ratio_cutoff = 1,
                                      p_value_cutoff = 0.05)
  
  # get annotations 
  annotations <- get_annotations(celltypes,
                                 celltype_of_interest,
                                 tissue_of_interest,
                                 metadata,
                                 enriched,
                                 H3K27ac,
                                 expression,
                                 expressed,
                                 HiChIP,
                                 TADs)

  # 2) VARIANTS ======================================================================================================
  # get variant-to-gene universe
  cat("2) Finding all genes near variants...\n")

  # The transcript-x-variant universe (masterlist of all possible transcript x variant pairs < max_variant_to_gene_distance apart)
  vxt_master <- get_vxt_master(variants,
                               TSSs,
                               max_variant_to_gene_distance)
  vxt_master_basic <- vxt_master %>% dplyr::select(cs, variant, symbol, ensg, enst, protein_coding)

  # 3) ANNOTATING ======================================================================================================
  cat("3) Annotating variant-transcript pairs at every level...\n")

  cat("  > G\tAnnotating genes...\n")
  g <- get_g_level_annotations(
    vxt_master,
    annotations)
  
  cat("  > VxT\tAnnotating variant x transcript pairs...\n")
  vxt <- get_vxt_level_annotations(
    variants,
    DHSs,
    vxt_master,
    annotations)

  cat("  > VxG\tAnnotating variant x gene pairs...\n")
  vxg <- get_vxg_level_annotations(
    variants,
    vxt_master,
    annotations)
 
  # 4) ALL ANNOTATIONS ======================================================================================================
  # Master list of variant-x-transcript annotations
  # -> wide-format (one row per vxt pair, one column per celltype if ct-specific)
  # -> only variant-transcript combinations within 2Mb are included
  # -> rows match vxt_master
  cat("4) Generating master table of transcript x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(vxt, vxg) %>% purrr::map(~ matricise_by_pair(., vxt_master))
  
  # 5) SCORING ======================================================================================================
  cat("5) Scoring variant-gene pairs...\n")

  # get weights
  weights <- get_weights(weights_file, master)
  
  # split annotations into tissues
  tissue_annotations <- annotations$selected_celltypes$tissue %>% unique %>%
    sapply(function(tissue){
      # gather tissue annotations
      annots <- master %>% names %>%
        sapply(function(annotation){
          cols <- colnames(master[[annotation]])
          master[[annotation]][,cols[greplany(c(tissue, "value"), cols)], drop = F] %>% rowMeans
        }, USE.NAMES = T, simplify = T)
      # weight and score
      score <- annots %>% {rowMeans(. * weights[,1][match(colnames(.), rownames(weights))][col(.)])}
      cbind(vxt_master_basic, score, annots)
    }, USE.NAMES = T, simplify = F)
  
  # calculate pair scores
  scores <- tissue_annotations %>%
    purrr::map(
      .f = function(x){
        vxt <- x %>%
          dplyr::as_tibble() %>%
          dplyr::select(dplyr::all_of(colnames(vxt_master_basic)), score, vxt_exon_or_inv_distance) %>%
          dplyr::arrange(-score) 
        # max score per vxg
        vxg <- vxt %>%
        dplyr::group_by(cs, variant, symbol, ensg, protein_coding) %>%
        dplyr::summarise(dplyr::across(where(is.numeric), max))
        # max score per CS
        vxg_max <- vxg %>%
          dplyr::group_by(cs) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), 
                                      ~ as.numeric(.x == max(.x) & .x > 0)))
        # max pc score per CS
        vxg_max_pc <- vxg %>%
          dplyr::group_by(cs, protein_coding) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), 
                                      ~ as.numeric(.x == max(.x) & .x > 0 & protein_coding)))
        return(list(vxt = vxt, 
                    vxg = vxg, 
                    vxg_max = vxg_max, 
                    vxg_max_pc = vxg_max_pc))
      }
    )   
  
  # save
  saveRDS(tissue_annotations, paste0(out$base, "annotations.rds"))
  saveRDS(scores, paste0(out$base, "scores.rds"))

  # 6) PERFORMANCE ======================================================================================================
  if(do_performance){
  cat("6) Measuring EG2 performance...\n")

  # import user-provided known genes and check that all symbols are in the GENCODE database
  cat("  > Importing known genes...\n")
  known_genes <- read_tibble(known_genes_file)$V1 %>%
    check_known_genes(known_genes_file)

  # run performance tests on all annotations
  performance <- tissue_annotations %>%
    purrr::imap(function(x, y){
      performance <- get_performance(x, vxt_master, known_genes, pcENSGs, max_n_known_genes_per_CS)
      
      title_plot <- list(ggplot2::labs(
        title = paste0(
          "\nTrait = ", trait,
          "\nmax n known genes per CS = ", max_n_known_genes_per_CS,
          "; max distance = ", max_variant_to_gene_distance,
          "\nCell types = ", celltypes_dir
        ),
        subtitle = paste0(
          "Tissue(s) = ", y,
          "\n",
          paste(strwrap(
            paste0(
              "Cell type(s) = ", if(celltypes == "all_celltypes"){"all cell types"}
              else{dplyr::filter(metadata, tissue == y)$celltype %>% unique %>% paste(collapse = ", ")}
            ),
            width = 70
          ), collapse = "\n")
        )
      ))
      
      
      
      # plot extras
      weight_facets <- dplyr::tibble(prediction_method = unique(performance$summary$prediction_method)) %>%
        dplyr::full_join(weights %>% dplyr::as_tibble(rownames = "prediction_method"),
                         by = "prediction_method") %>%
        dplyr::mutate(
          weight = factor(dplyr::case_when(
            grepl("^score", prediction_method) ~ "score",
            TRUE ~ as.character(weight)),
            levels = c("score",
                       weights[,"weight"] %>% unique %>% sort(T) %>% as.character)))
      
      {pdf(out$performance_plot, height = 10, width = 15, onefile = T)
        
        # PR score + max
        print(
          performance %>% plot_PR(, colour = prediction_method) + title_plot
        )
        
        # PR score + max facets
        print(
          performance %>% 
                  plot_PR() +
                  ggplot2::facet_wrap(. ~ prediction_method) +
                  ggplot2::geom_text(data = performance$summary %>%
                                       dplyr::transmute(prediction_method,
                                                        label = round(score_PR_AUC, 3)),
                                     mapping = ggplot2::aes(x = -Inf, y = Inf, label = label),
                                     hjust = 0, vjust = 1, size = 5) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                                 axis.text = ggplot2::element_blank()) + 
            title_plot
        )
        
        # PR max
        print(
          performance %>% purrr::map(dplyr::filter, prediction_type == "max") %>%
            plot_PR(colour = prediction_method) +
                ggrepel::geom_text_repel(min.segment.length = 0, max.overlaps = Inf) +
                title_plot
        )
        
        # F score max
        print(
          performance$summary %>%
              dplyr::mutate(F_score = F_score %>% tidyr::replace_na(0)) %>%
              dplyr::distinct() %>%
              ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, F_score),
                                           y = F_score)) +
              ggplot2::geom_col() +
              ggplot2::labs(x = "Predictor",
                            y = "F score") +
              ggsci::scale_fill_igv() +
              ggplot2::coord_flip() +
              title_plot
        )
        
        # performance metrics
        print(
          performance$summary %>%
            dplyr::mutate(dplyr::across(where(is.numeric), tidyr::replace_na, 0),
                          fsc = F_score) %>%
            tidyr::pivot_longer(cols = c(F_score, score_PR_AUC, Precision, Recall),
                                names_to = "metric",
                                values_to = "performance") %>%
            ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, fsc),
                                         y = performance,
                                         fill = prediction_method)) +
            ggplot2::geom_col() +
            ggplot2::facet_grid(. ~ metric,
                                scales = "free", space = "free_y") +
            ggplot2::coord_flip() +
            ggplot2::theme(axis.title = ggplot2::element_blank(),
                           legend.position = "none") +
            ggsci::scale_fill_igv() +
            title_plot
        )
        
        dev.off()}
      
      return(performance)
      
    })
  # write table
  saveRDS(performance, filename = paste0(out$base, "performance.rds"))
 
  }

  # 7) XGBoost MODEL TRAINING ======================================================================================================
  if(do_XGBoost){
  cat("7) Training an XGBoost model...\n")

  # format training set
  full <- annotations %>%
    dplyr::mutate(label = (symbol %in% known_genes$symbol) %>% as.numeric) %>%
    dplyr::group_by(cs) %>%
    dplyr::filter(any(label == 1)) %>%
    dplyr::ungroup()
  train <- list(data = full %>% dplyr::select(names(master)) %>% as.matrix,
                label = full$label)
  dtrain <- xgboost::xgb.DMatrix(data = train$data,
                                 label = train$label)
  # model training
  xgb1 <- xgboost::xgboost(data = dtrain,
                           max.depth = 5,
                           eta = 1,
                           nthread = 1,
                           nrounds = 100,
                           objective = "binary:logistic",
                           verbose = 1)
  # view feature importance plot
  xgb1_feature_importance_mat <- xgboost::xgb.importance(feature_names = colnames(dtrain), model = xgb1)
  xgb1_feature_importance_plot <- xgboost::xgb.ggplot.importance(importance_matrix = xgb1_feature_importance_mat)
  # save plot
  {
    pdf(out$XGBoost_plot, width = 15, height = 15)
    print(xgb1_feature_importance_plot + title_plot)
    dev.off()
  }
  # write table
  write_tibble(xgb1_feature_importance_mat, filename = out$XGBoost)
  }

  cat("\nDONE!\n")
  # 8) SAVE ===
  return(scores)
}




# culled levels: #
# cat("  > C\tAnnotating credible sets...\n")
# c <- get_c_level_annotations(
#   variants)
# cat("  > V\tAnnotating variants...\n")
# levels$v <- get_v_level_annotations(
#   variants,
#   H3K27ac,
#   enriched,
#   vxt_master,
#   DHSs)
# cat("  > T\tAnnotating transcripts...\n")
# levels$t <- get_t_level_annotations(
#   TSSs,
#   DHSs,
#   enriched)
# cat("  > G\tAnnotating genes...\n")
# levels$g <- get_g_level_annotations(
#   vxt_master,
#   enriched)
# cat("  > CxT\tAnnotating credible set x transcript pairs...\n")
# levels$cxt <- get_cxt_level_annotations(
#   enriched,
#   vxt,
#   variants)
# 
# cat("  > CxG\tAnnotating credible set x gene pairs...\n")
# levels$cxg <- get_cxg_level_annotations(
#   vxt_master)


