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
  n_variants <- length(log_var_bf)
  max_iter <- as.integer(max_iter)

  # Initialize parameters
  eta <- eta_init
  if (!estimate_eta) eta <- fixed_eta
  delta <- delta_init

  # Precompute factor for rowsum
  category_f <- factor(category, levels = seq_len(n_categories))

  if (type == "gene") {
    gene_index_f <- factor(gene_index, levels = seq_len(n_genes))
  }

  # Convergence tracking
  eta_history <- matrix(NA_real_, nrow = max_iter, ncol = n_categories)
  eta_history[1, ] <- eta
  if (type == "gene") {
    delta_history <- numeric(max_iter)
    delta_history[1] <- delta
  }

  converged <- FALSE
  final_iter <- 1L

  for (iter in 2:max_iter) {
    eta_prev <- eta
    delta_prev <- delta

    if (type == "gene") {
      # Gene-Level EM (vectorized for performance)

      # Gene-level log BFs: sum of log((1-eta_k) + eta_k * BF_j) per gene
      log_mix_bf <- log_mixture_bf(log_var_bf, eta_prev[category])
      gene_log_bf <- as.numeric(rowsum(log_mix_bf, gene_index_f))

      # E-step: EUi = delta*B_i / (delta*B_i + 1-delta)
      log_numer <- log(delta_prev) + gene_log_bf
      log_denom <- log_sum_exp_vec2(log_numer, rep(log1p(-delta_prev), n_genes))
      EUi <- exp(log_numer - log_denom)

      # Variant-level posterior: P_j = eta_k*BF_j / (eta_k*BF_j + 1-eta_k)
      log_eta_cat <- log(eta_prev[category])
      log_one_minus_eta_cat <- log1p(-eta_prev[category])
      a <- log_eta_cat + log_var_bf
      b <- log_one_minus_eta_cat
      log_P_j <- a - log_sum_exp_vec2(a, b)
      P_j <- exp(log_P_j)

      # UiZij = EUi[gene] * P_j
      UiZij <- EUi[gene_index] * P_j

      # M-step: delta
      if (estimate_delta) {
        delta <- mean(EUi)
        delta <- max(delta, .Machine$double.xmin)
        delta <- min(delta, 1 - .Machine$double.eps)
      }

      # M-step: eta (vectorized with table)
      if (estimate_eta) {
        numerator <- as.numeric(rowsum(UiZij, category_f))
        # Denominator: sum over genes of (count_cat_in_gene * EUi[gene])
        counts <- table(gene_index_f, category_f)
        denominator <- as.numeric(colSums(counts * EUi))
        eta <- ifelse(denominator > 0, numerator / denominator, 0)
        eta <- pmax(eta, .Machine$double.xmin)
        eta <- pmin(eta, 1 - .Machine$double.eps)
      }

    } else {
      # Variant-level EM (vectorized for performance)

      # E-step: EZj = eta_k*BF_j / (eta_k*BF_j + 1-eta_k)
      log_eta_cat <- log(eta_prev[category])
      log_one_minus_eta_cat <- log1p(-eta_prev[category])
      a <- log_eta_cat + log_var_bf
      b <- log_one_minus_eta_cat
      EZj <- exp(a - log_sum_exp_vec2(a, b))

      # M-step: eta
      if (estimate_eta) {
        numerator <- as.numeric(rowsum(EZj, category_f))
        denominator <- as.numeric(table(category_f))
        eta <- ifelse(denominator > 0, numerator / denominator, 0)
        eta <- pmax(eta, .Machine$double.xmin)
        eta <- pmin(eta, 1 - .Machine$double.eps)
      }
    }

    # Store history
    eta_history[iter, ] <- eta
    if (type == "gene") delta_history[iter] <- delta

    # Convergence check
    diff <- sum(abs(eta - eta_prev))
    if (type == "gene" && estimate_delta) diff <- diff + abs(delta - delta_prev)
    final_iter <- iter
    if (diff < tol) {
      converged <- TRUE
      break
    }
  }

  # Trim history
  eta_history <- eta_history[seq_len(final_iter), , drop = FALSE]
  if (type == "gene") delta_history <- delta_history[seq_len(final_iter)]

  # Final posteriors at convergence
  if (type == "gene") {
    log_mix_bf_final <- log_mixture_bf(log_var_bf, eta[category])
    gene_log_bf_final <- as.numeric(rowsum(log_mix_bf_final, gene_index_f))

    log_numer_final <- log(delta) + gene_log_bf_final
    log_denom_final <- log_sum_exp_vec2(log_numer_final, rep(log1p(-delta), n_genes))
    EUi_final <- exp(log_numer_final - log_denom_final)

    log_eta_cat_f <- log(eta[category])
    a_f <- log_eta_cat_f + log_var_bf
    b_f <- log1p(-eta[category])
    log_P_j_final <- a_f - log_sum_exp_vec2(a_f, b_f)
    post_prob <- EUi_final[gene_index] * exp(log_P_j_final)

    result <- list(
      delta.est = delta,
      eta.est = eta,
      gene.log.bf = gene_log_bf_final,
      post.prob = post_prob,
      EUi = EUi_final,
      eta.history = eta_history,
      delta.history = delta_history,
      converged = converged,
      n.iter = final_iter
    )
  } else {
    log_eta_cat_f <- log(eta[category])
    a_f <- log_eta_cat_f + log_var_bf
    b_f <- log1p(-eta[category])
    post_prob <- exp(a_f - log_sum_exp_vec2(a_f, b_f))

    result <- list(
      eta.est = eta,
      post.prob = post_prob,
      eta.history = eta_history,
      converged = converged,
      n.iter = final_iter
    )
  }

  return(result)
}
