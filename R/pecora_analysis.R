#' Performs the PeCorA analysis described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/.
#' Take a look at the linked paper to get an idea of the procedure.
#'
#' Takes an msnset, a list of proteins, and a treatment variable. The data is standardized depending on the
#' chosen parameters of 'sample_standardize' and 'peptide_standardize'. Then the PeCorA function from the
#' PeCorA package is used to find discordant peptides. Plots of the significant peptides are made and
#' saved to the given folder. Finally, the output contains a table with the phosphosites mapping to the given
#' proteins, along with their PeCorA analysis p-values.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param treatment_string the name of the column in pData containing treatment information. The data can be numerical of categorical.
#' @param proteins character vector indicating which proteins to analyze.
#' @param sample_standardize logical for whether exprs data should be standardized so that sample wise we have mean zero and variance one.
#' @param peptide_standardize logical for whether exprs data should be normalized so that logratio data is relative to the control group.
#' @param folder folder in which to save the plots. For example folder = "directory/".
#' @param save_plots logical for whether to save boxplots/scatterplots of the significant peptides.
#'
#' @return PeCorA results table for the supplied proteins. Saves plots of the significant peptides to the given folder.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom PeCorA PeCorA
#' @importFrom grid grobTree textGrob gpar
#'
#' @export pecora_analysis
#'
pecora_analysis <- function(m, treatment_string, control_group = NULL, proteins,
                            sample_standardize = TRUE, peptide_standardize = TRUE,
                            folder = "", save_plots = TRUE, median_mod = FALSE){
  mat <- exprs(m)
  peptide_mapping <- fData(m)

  if (!("Protein" %in% colnames(peptide_mapping))){
    stop("fData should have 'Protein' and 'Peptide' columns.\n")
  }

  peptide_mapping <- peptide_mapping %>%
    mutate(Protein = as.character(Protein)) %>%
    select(Protein) %>%
    filter(Protein %in% proteins) %>%
    mutate(Peptide = rownames(.))

  feature_names <- rownames(peptide_mapping)

  if (sample_standardize) {
    cat("Standardizing by sample.\n")
    sample_means <- apply(mat, 2, mean, na.rm = T)
    sample_sds <- apply(mat, 2, sd, na.rm = T)
    mat <- sweep(mat, 2, sample_means, FUN = '-')
    mat <- sweep(mat, 2, sample_sds, FUN = '/')
  }

  metadata <- pData(m) %>%
    select(sym(treatment_string))
  metadata$Sample <- rownames(metadata)

  mat <- mat[feature_names, ]

  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }

  if (peptide_standardize & is.factor(metadata[[1]])) {
    cat("Standardizing to control group.\n")
    if (is.null(control_group)) {
      stop("Must define 'control_group' in order to normalize relative to control\n")
    }

    samples_control <- metadata %>%
      filter(!!sym(treatment_string) == control_group) %>%
      rownames()

    if (length(samples_control) == 0){
      stop("No samples in control group. Please check value of 'control_group' and pData\n")
    }

    ## In case samples_control consists of a single sample we need as.matrix() in the apply call.
    peptide_means <- apply(as.matrix(mat[, samples_control]), 1, mean, na.rm = T)
    mat <- sweep(mat, 1, peptide_means, FUN = '-')

  } else if (peptide_standardize) {

    cat("Mean centering each peptide.\n")
    peptide_means <- apply(mat, 1, mean, na.rm = T)
    mat <- sweep(mat, 1, peptide_means, FUN = '-')

  }

  PeCorA_input <- mat %>%
    as.data.frame() %>%
    mutate(Peptide = rownames(.)) %>%
    pivot_longer(cols = -Peptide, names_to = "Sample", values_to = "LogRatio") %>%
    mutate(modpep_z = Peptide,
           ms1adj = LogRatio) %>%
    merge(peptide_mapping, by = "Peptide") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment_string))

  PeCorA_result <- PeCorA_mod(PeCorA_input, median_mod) %>%
    dplyr::rename(Peptide = peptide,
                  Protein = protein) %>%
    group_by(Protein) %>%
    ## Adjust raw pvalues using only pvalues from the same protein (group by protein)
    mutate(adj_pval2 = p.adjust(pvalue, method = "BH"))

  to_plot <- PeCorA_result %>%
    filter(adj_pval2 < 0.05)

  if (save_plots){
    lapply(1:nrow(to_plot), function(i){
      if (median_mod){
        plot_path <- paste("PeCorA_median_mod", treatment_string, to_plot$Protein[[i]], to_plot$Peptide[[i]], sep = "_") %>%
          paste0(folder, ., ".png")
      } else {
        plot_path <- paste("PeCorA", treatment_string, to_plot$Protein[[i]], to_plot$Peptide[[i]], sep = "_") %>%
          paste0(folder, ., ".png")
      }
      pecora_plot(PeCorA_result, PeCorA_input, to_plot$Protein[[i]], to_plot$Peptide[[i]], plot_path,
                  label = treatment_string)
    })
  }

  return(list("Result" = PeCorA_result, "Input" = PeCorA_input))
}


#' Plots a particular peptide using the PeCorA input and output tables.
#'
#' @param PeCorA_result PeCorA results table from pecora_analysis.
#' @param PeCorA_input PeCorA input table from pecora_analysis
#' @param chosen_protein The protein of interest.
#' @param chosen_peptide The peptide of interest. Must map to the supplied protein.
#' @param plot_path String indicating the path of the plot when saved. Include extension, eg "plot.png". If NULL the plot isn't saved
#' @param label String to include as a custom label in the plots. Should be the treatment variable supplied to pecora_analysis.
#'
#' @return Plot of the specified peptide + protein.
#'
#' @importFrom grid grobTree textGrob gpar
#'
#' @export pecora_analysis
#'

pecora_plot <- function(PeCorA_result, PeCorA_input, chosen_protein, chosen_peptide,
                        plot_path = NULL, label = NULL) {

  results_sig <- PeCorA_result %>%
    filter(Protein == chosen_protein) %>%
    filter(Peptide == chosen_peptide)

  if (nrow(results_sig) == 0){
    stop("Protein/Peptide not found.")
  } else if (nrow(results_sig) > 1){
    ## This probably shouldn't happen, but just in case.
    stop("Given Protein + Peptide match more than one combination.")
  }

  p.value <- results_sig %>%
    pull(adj_pval2) %>% round(12)

  ## Boxplot if condition is a factor. Otherwise, scatterplots when using numerical data.
  if (is.factor(PeCorA_input$Condition)) {

    plot.df <- PeCorA_input %>%
      select(Sample, Protein, Peptide, ms1adj, Condition) %>%
      filter(Protein == chosen_protein) %>%
      mutate(group = case_when(Peptide == chosen_peptide ~ chosen_peptide,
                               TRUE ~ "All other peptides")) %>%
      mutate(group = factor(group, levels = c("All other peptides", chosen_peptide)))

    grob <- grobTree(textGrob(paste0("Adj_pval = ", p.value),
                              x = 0.1,  y = 0.97, hjust = 0,
                              gp = gpar(col = "black", fontsize = 11)))

    p <- ggplot(plot.df, aes(x = Condition, y = ms1adj, fill = group)) +
      geom_boxplot(notch = TRUE, outlier.shape = NA) +
      ggtitle(chosen_peptide) + annotation_custom(grob) +
      xlab(label) + ylab("Log Intensity") +
      theme(plot.title = element_text(hjust = 0.5))

  } else {

    plot.df <- PeCorA_input %>%
      select(Sample, Protein, Peptide, ms1adj, Condition) %>%
      filter(Protein == chosen_protein) %>%
      mutate(group = case_when(Peptide == chosen_peptide ~ chosen_peptide,
                               TRUE ~ "All other peptides"),
             alpha = case_when(group == "All other peptides" ~ 0.005,
                               TRUE ~ 0.25)) %>%
      mutate(group = factor(group, levels = c("All other peptides", chosen_peptide)))

    lm.df <- plot.df %>%
      filter(group == "All other peptides")
    allothers.lm <- lm(lm.df$ms1adj ~ lm.df$Condition)

    lm.df <- plot.df %>%
      filter(group == chosen_peptide)
    chosen.lm <- lm(lm.df$ms1adj ~ lm.df$Condition)

    lm.peptide.plot <- lm(plot.df$ms1adj ~ plot.df$Condition)

    grob <- grobTree(textGrob(paste0("Adj_pval = ", p.value),
                              x = 0.1,  y = 0.97, hjust = 0,
                              gp = gpar(col = "black", fontsize = 11)))

    p <- ggplot(plot.df, aes(x = Condition, y = ms1adj, color = group, alpha = alpha)) + geom_point() +
      ggtitle(chosen_peptide) + annotation_custom(grob) +
      ylab("Log Intensity") + guides(alpha = "none") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_abline(intercept = allothers.lm$coefficients[[1]],
                  slope = allothers.lm$coefficients[[2]],
                  color = "red", size = 1) +
      geom_abline(intercept = chosen.lm$coefficients[[1]],
                  slope = chosen.lm$coefficients[[2]],
                  color = scales::hue_pal()(2)[[2]], size = 1)
  }

  if (!is.null(plot_path)){
    ggsave(plot_path, plot = p, width = 6.5, height = 5)
  }

  return(p)
}













