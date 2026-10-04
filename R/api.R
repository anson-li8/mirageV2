#' mirage: MIxture model based Rare variant Analysis on GEnes
#'
#' Gene-level rare variant association test using a Bayesian mixture model.
#'
#' @param data Data frame with columns: ID, Gene, No.case, No.contr, category.
#'   A 4-column input (without ID) is accepted with a warning.
#' @param n1 Sample size in cases.
#' @param n2 Sample size in controls.
#' @param gamma Hyper prior shape parameter(s). Scalar or length-K vector.
#' @param sigma Hyper prior scale parameter(s). Scalar or length-K vector.
#' @param eta.init Initial value for proportion of risk variants.
#' @param delta.init Initial value for proportion of risk genes.
#' @param estimate.delta Logical. Whether to estimate delta.
#' @param estimate.eta Logical. Whether to estimate eta.
#' @param fixed_eta Numeric vector of fixed eta values (required when
#'   estimate.eta = FALSE). Must have length K.
#' @param max.iter Maximum EM iterations.
#' @param tol Convergence tolerance.
#' @param verbose Logical. Print progress messages.
#' @return An object of class \code{mirage_result}.
#' @export
#' @examples
#' # See vignette for full examples
mirage <- function(data, n1, n2, gamma = 3, sigma = 2,
                   eta.init = 0.1, delta.init = 0.1,
                   estimate.delta = TRUE, estimate.eta = TRUE,
                   fixed_eta = NULL, max.iter = 10000L,
                   tol = 1e-5, verbose = TRUE) {

  cl <- match.call()

  # Validate
  dat <- validate_mirage_data(data, type = "gene")
  original_categories <- attr(dat, "original_categories")
  n_categories <- length(original_categories)

  validate_mirage_params(n1, n2, gamma, sigma, eta.init, delta.init,
                         estimate.delta, estimate.eta, fixed_eta,
                         max.iter, tol, n_categories, type = "gene")

  # Variant-level log BFs (category-specific gamma/sigma)
  gamma_cat <- rep_len(gamma, n_categories)
  sigma_cat <- rep_len(sigma, n_categories)

  if (verbose) cli::cli_progress_step("Computing variant-level Bayes factors")

  log_var_bf <- calc_log_bf(
    dat$No.case, dat$No.contr,
    gamma_cat[dat$category], sigma_cat[dat$category],
    n1, n2
  )

  # Gene indexing
  unique_genes <- unique(dat$Gene)
  gene_index <- match(dat$Gene, unique_genes)
  n_genes <- length(unique_genes)

  # EM
  if (verbose) cli::cli_progress_step("Running EM algorithm ({n_genes} genes)")

  em <- em_mirage(
    type = "gene", log_var_bf = log_var_bf, category = dat$category,
    gene_index = gene_index, n_genes = n_genes, n_categories = n_categories,
    delta_init = delta.init, eta_init = rep_len(eta.init, n_categories),
    estimate_delta = estimate.delta, estimate_eta = estimate.eta,
    fixed_eta = fixed_eta, max_iter = max.iter, tol = tol, verbose = verbose
  )

  if (verbose) {
    if (em$converged) {
      cli::cli_alert_success("EM converged in {em$n.iter} iterations")
    } else {
      cli::cli_alert_warning("EM did not converge after {em$n.iter} iterations")
    }
  }

  # Gene-level BFs and posterior probabilities
  delta <- em$delta.est
  eta <- em$eta.est
  gene_bf <- exp(em$gene.log.bf)
  post_prob_gene <- (delta * gene_bf) / (delta * gene_bf + 1 - delta)

  # LRT p-values
  if (verbose) cli::cli_progress_step("Computing LRT statistics and p-values")

  # Per-gene LRT (df = 2: delta + one pooled eta)
  log_lkhd <- vapply(em$gene.log.bf, function(lb) {
    log_sum_exp(log(1 - delta), log(delta) + lb)
  }, numeric(1))
  gene_pvalue <- pchisq(2 * log_lkhd, df = 2, lower.tail = FALSE)

  # Overall delta p-value (df = 1 + K)
  total_pvalue <- pchisq(2 * sum(log_lkhd), df = 1 + n_categories,
                         lower.tail = FALSE)

  # Category-level LRT (df = 1 each)
  cate_pvalue <- vapply(seq_len(n_categories), function(g) {
    gene_log_bf_g <- vapply(seq_len(n_genes), function(i) {
      idx <- which(gene_index == i & dat$category == g)
      if (length(idx) == 0) return(0)
      sum(log_mixture_bf(log_var_bf[idx], eta[g]))
    }, numeric(1))
    log_lkhd_g <- vapply(gene_log_bf_g, function(lb) {
      log_sum_exp(log(1 - delta), log(delta) + lb)
    }, numeric(1))
    stat <- 2 * sum(log_lkhd_g)
    pchisq(stat, df = 1, lower.tail = FALSE)
  }, numeric(1))

  # Per-gene variant tables
  bf_all <- vector("list", n_genes)
  for (i in seq_len(n_genes)) {
    idx <- which(gene_index == i)
    bf_all[[i]] <- data.frame(
      ID = dat$ID[idx],
      Gene = dat$Gene[idx],
      No.case = dat$No.case[idx],
      No.contr = dat$No.contr[idx],
      category = original_categories[dat$category[idx]],
      var.BF = exp(log_var_bf[idx]),
      log_var.BF = log_var_bf[idx],
      stringsAsFactors = FALSE
    )
  }
  names(bf_all) <- unique_genes

  # Assemble output
  delta_est_df <- data.frame(delta.est = delta, delta.pvalue = total_pvalue)

  eta_est_df <- data.frame(eta.est = eta, eta.pvalue = cate_pvalue)
  rownames(eta_est_df) <- original_categories

  bf_pp_gene <- data.frame(Gene = unique_genes, BF = gene_bf,
                           post.prob = post_prob_gene)

  new_mirage_result(delta_est_df, eta_est_df, bf_pp_gene, bf_all,
                    original_categories, cl)
}


#' mirage_vs: Variant-set level analysis
#'
#' Variant-level rare variant association test using a Bayesian mixture model.
#'
#' @inheritParams mirage
#' @param data Data frame with columns: ID, No.case, No.contr, category.
#'   A 3-column input (without ID) is accepted with a warning.
#' @return An object of class \code{mirage_vs_result}.
#' @export
mirage_vs <- function(data, n1, n2, gamma = 3, sigma = 2,
                      eta.init = 0.1, estimate.eta = TRUE,
                      fixed_eta = NULL, max.iter = 10000L,
                      tol = 1e-5, verbose = TRUE) {

  cl <- match.call()

  # Validate
  dat <- validate_mirage_data(data, type = "vs")
  original_categories <- attr(dat, "original_categories")
  n_categories <- length(original_categories)

  validate_mirage_params(n1, n2, gamma, sigma, eta.init,
                         delta.init = 0.1,  # unused for vs
                         estimate.delta = FALSE,
                         estimate.eta, fixed_eta,
                         max.iter, tol, n_categories, type = "vs")

  # Variant-level log BFs
  gamma_cat <- rep_len(gamma, n_categories)
  sigma_cat <- rep_len(sigma, n_categories)

  if (verbose) cli::cli_progress_step("Computing variant-level Bayes factors")

  log_var_bf <- calc_log_bf(
    dat$No.case, dat$No.contr,
    gamma_cat[dat$category], sigma_cat[dat$category],
    n1, n2
  )

  # EM
  if (verbose) cli::cli_progress_step("Running EM algorithm ({nrow(dat)} variants)")

  em <- em_mirage(
    type = "vs", log_var_bf = log_var_bf, category = dat$category,
    n_categories = n_categories,
    eta_init = rep_len(eta.init, n_categories),
    estimate_delta = FALSE, estimate_eta = estimate.eta,
    fixed_eta = fixed_eta, max_iter = max.iter, tol = tol, verbose = verbose
  )

  if (verbose) {
    if (em$converged) {
      cli::cli_alert_success("EM converged in {em$n.iter} iterations")
    } else {
      cli::cli_alert_warning("EM did not converge after {em$n.iter} iterations")
    }
  }

  # LRT p-values
  if (verbose) cli::cli_progress_step("Computing LRT statistics and p-values")

  eta <- em$eta.est

  # Category-level LRT
  cate_pvalue <- vapply(seq_len(n_categories), function(g) {
    idx_g <- which(dat$category == g)
    if (length(idx_g) == 0) return(1)
    log_lkhd_g <- sum(log_mixture_bf(log_var_bf[idx_g], eta[g]))
    pchisq(2 * log_lkhd_g, df = 1, lower.tail = FALSE)
  }, numeric(1))

  # Assemble output
  eta_est_df <- data.frame(eta.est = eta, eta.pvalue = cate_pvalue)
  rownames(eta_est_df) <- original_categories

  full_info <- data.frame(
    ID = dat$ID,
    No.case = dat$No.case,
    No.contr = dat$No.contr,
    category = original_categories[dat$category],
    var.BF = exp(log_var_bf),
    log_var.BF = log_var_bf,
    stringsAsFactors = FALSE
  )

  post_prob_df <- data.frame(
    variant = dat$ID,
    BF = exp(log_var_bf),
    post.prob = em$post.prob
  )

  new_mirage_vs_result(eta_est_df, full_info, post_prob_df,
                       original_categories, cl)
}
