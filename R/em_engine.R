#' MIRAGE EM Algorithm
#'
#' @param type Either "gene" or "vs".
#' @param log_var_bf Numeric vector of log Bayes factors for each variant.
#' @param category Integer vector of category indices (1:K) for each variant.
#' @param gene_index Integer vector of gene indices (1:G) for each variant.
#' @param n_genes Number of genes. Required when type = "gene".
#' @param n_categories Number of variant categories (K).
#' @param delta_init Initial value for delta.
#' @param eta_init Numeric vector of initial eta values (length K).
#' @param estimate_delta Logical.
#' @param estimate_eta Logical.
#' @param fixed_eta Numeric vector of fixed eta values.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance.
#' @param verbose Logical.
#' @return A list containing EM results.
#' @keywords internal
em_mirage <- function(type = c("gene", "vs"),
                      log_var_bf,
                      category,
                      gene_index = NULL,
                      n_genes = NULL,
                      n_categories,
                      delta_init = 0.1,
                      eta_init = rep(0.1, n_categories),
                      estimate_delta = TRUE,
                      estimate_eta = TRUE,
                      fixed_eta = NULL,
                      max_iter = 10000L,
                      tol = 1e-5,
                      verbose = TRUE) {
  type <- match.arg(type)
  max_iter <- as.integer(max_iter)

  if (!estimate_eta && is.null(fixed_eta)) {
    cli::cli_abort("{.arg fixed_eta} must be provided when {.arg estimate_eta} is FALSE.")
  }
  if (!estimate_eta) {
    eta_init <- fixed_eta
  }

  if (type == "gene") {
    # C++ expects 0-indexed
    result <- em_mirage_cpp(
      category = as.integer(category - 1L),
      gene_index = as.integer(gene_index - 1L),
      log_var_bf = as.numeric(log_var_bf),
      n_genes = as.integer(n_genes),
      n_categories = as.integer(n_categories),
      delta_init = as.numeric(delta_init),
      eta_init = as.numeric(eta_init),
      estimate_delta = estimate_delta,
      estimate_eta = estimate_eta,
      fixed_eta = if (!estimate_eta) as.numeric(fixed_eta) else rep(0.0, n_categories),
      max_iter = max_iter,
      tol = as.numeric(tol)
    )
  } else {
    # Variant-level: treat each variant as its own "gene"
    n_variants <- length(log_var_bf)
    result <- em_mirage_cpp(
      category = as.integer(category - 1L),
      gene_index = as.integer(seq_len(n_variants) - 1L),
      log_var_bf = as.numeric(log_var_bf),
      n_genes = as.integer(n_variants),
      n_categories = as.integer(n_categories),
      delta_init = 1.0,  # delta=1 means all "genes" are risk genes
      eta_init = as.numeric(eta_init),
      estimate_delta = FALSE,
      estimate_eta = estimate_eta,
      fixed_eta = if (!estimate_eta) as.numeric(fixed_eta) else rep(0.0, n_categories),
      max_iter = max_iter,
      tol = as.numeric(tol)
    )
    # For vs-level, post.prob is just the variant-level posterior
    # EUi is always 1 when delta=1, so post_prob = P_j directly
    # The C++ calculates this correctly since EUi=1 when delta=1
  }

  return(result)
}
