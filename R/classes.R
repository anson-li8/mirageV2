#' Construct a mirage_result object
#'
#' @param delta_est Data frame of delta estimates.
#' @param eta_est Data frame of eta estimates.
#' @param bf_pp_gene Data frame of gene-level Bayes Factors and Posterior Probabilities.
#' @param bf_all List of variant-level Bayes Factors per gene.
#' @param original_categories The original category labels before internal re-indexing.
#' @param call The matched call.
#' @return An object of class `mirage_result`.
#' @keywords internal
new_mirage_result <- function(delta_est, eta_est, bf_pp_gene, bf_all,
                              original_categories, call) {
  structure(
    list(
      delta.est = delta_est,
      eta.est = eta_est,
      BF.PP.gene = bf_pp_gene,
      BF.all = bf_all,
      original_categories = original_categories,
      call = call
    ),
    class = "mirage_result"
  )
}

#' Construct a mirage_vs_result object
#'
#' @param eta_est Data frame of eta estimates.
#' @param full_info Data frame of variant-level information and Bayes Factors.
#' @param post_prob Data frame of variant-level Posterior Probabilities.
#' @param original_categories The original category labels.
#' @param call The matched call.
#' @return An object of class `mirage_vs_result`.
#' @keywords internal
new_mirage_vs_result <- function(eta_est, full_info, post_prob,
                                 original_categories, call) {
  structure(
    list(
      eta.est = eta_est,
      full.info = full_info,
      post.prob = post_prob,
      original_categories = original_categories,
      call = call
    ),
    class = "mirage_vs_result"
  )
}

#' @export
print.mirage_result <- function(x, ...) {
  cli::cli_h1("MIRAGE Gene-Level Analysis Result")
  cli::cli_text("Number of genes analyzed: {nrow(x$BF.PP.gene)}")
  cli::cli_text("Number of variant categories: {length(x$original_categories)}")
  cli::cli_text("Estimated proportion of risk genes (delta): {round(x$delta.est$delta.est, 4)}")
  cli::cli_h2("Top 5 Genes by Posterior Probability")
  top_genes <- head(x$BF.PP.gene[order(x$BF.PP.gene$post.prob, decreasing = TRUE), ], 5)
  print(top_genes, row.names = FALSE)
  invisible(x)
}

#' @export
print.mirage_vs_result <- function(x, ...) {
  cli::cli_h1("MIRAGE Variant-Set Analysis Result")
  cli::cli_text("Number of variants analyzed: {nrow(x$full.info)}")
  cli::cli_text("Number of variant categories: {length(x$original_categories)}")
  cli::cli_h2("Category-specific risk proportions (eta)")
  print(x$eta.est, row.names = FALSE)
  invisible(x)
}
