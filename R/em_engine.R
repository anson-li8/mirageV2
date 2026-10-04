#' MIRAGE EM Algorithm (Log-Space)
#'
#' Unified Expectation-Maximization engine for both gene-level and
#' variant-level mixture models. Runs entirely in log-space to
#' prevent numerical overflow/underflow.
#'
#' @param type Either "gene" or "vs".
#' @param log_var_bf Numeric vector of log Bayes factors for each variant.
#' @param category Integer vector of category indices (1:K) for each variant.
#' @param gene_index Integer vector of gene indices (1:G) for each variant.
#'   Required when type = "gene". Ignored when type = "vs".
#' @param n_genes Number of genes. Required when type = "gene".
#' @param n_categories Number of variant categories (K).
#' @param delta_init Initial value for delta (proportion of risk genes).
#' @param eta_init Numeric vector of initial eta values (length K).
#' @param estimate_delta Logical. Whether to estimate delta.
#' @param estimate_eta Logical. Whether to estimate eta.
#' @param fixed_eta Numeric vector of fixed eta values (used when estimate_eta = FALSE).
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance.
#' @param verbose Logical. Print progress.
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
  if (!estimate_eta) {
    eta <- fixed_eta
  }
  delta <- delta_init

  # Storage for tracking convergence
  eta_history <- matrix(NA_real_, nrow = max_iter, ncol = n_categories)
  eta_history[1, ] <- eta

  if (type == "gene") {
    delta_history <- numeric(max_iter)
    delta_history[1] <- delta
    gene_log_bf <- numeric(n_genes)
  }

  # Gene-level: precompute gene membership
  if (type == "gene") {
    gene_split <- split(seq_len(n_variants), gene_index)
  }

  # EM Loop
  converged <- FALSE
  final_iter <- 1L

  for (iter in 2:max_iter) {
    eta_prev <- eta
    delta_prev <- delta

    if (type == "gene") {
      # GENE-LEVEL EM

      # Compute gene-level log BFs
      # log_B_i = sum_j log((1-eta_kj) + eta_kj * BF_j)
      gene_log_bf <- vapply(seq_len(n_genes), function(i) {
        idx <- gene_split[[i]]
        if (length(idx) == 0) return(0)
        sum(log_mixture_bf(log_var_bf[idx], eta[category[idx]]))
      }, numeric(1))

      # E-step
      # EUi[i] = P(U_i=1 | data) = delta*B_i / (delta*B_i + 1-delta)
      # In log-space: log_EUi = log(delta) + log_B_i - log_sum_exp(log(delta)+log_B_i, log(1-delta))
      log_EUi <- vapply(gene_log_bf, function(lb) {
        log_numer <- log(delta_prev) + lb
        log_denom <- log_sum_exp(log_numer, log(1 - delta_prev))
        log_numer - log_denom
      }, numeric(1))
      EUi <- exp(log_EUi)

      # UiZij: posterior prob variant j is risk = EUi[i] * P_j
      # P_j = eta_k * BF_j / (eta_k * BF_j + (1 - eta_k))
      # log_P_j = log(eta_k) + log_BF_j - log_sum_exp(log(eta_k)+log_BF_j, log(1-eta_k))
      log_P_j <- log(eta_prev[category]) + log_var_bf -
        vapply(seq_len(n_variants), function(j) {
          log_sum_exp(log(eta_prev[category[j]]) + log_var_bf[j],
                      log(1 - eta_prev[category[j]]))
        }, numeric(1))
      P_j <- exp(log_P_j)

      # UiZij = EUi[gene_index] * P_j
      UiZij <- EUi[gene_index] * P_j

      # M-step
      if (estimate_delta) {
        delta <- mean(EUi)
        delta <- max(delta, .Machine$double.xmin)
        delta <- min(delta, 1 - .Machine$double.eps)
      }

      if (estimate_eta) {
        for (g in seq_len(n_categories)) {
          cat_variants <- which(category == g)
          if (length(cat_variants) == 0) next

          # total.Zij: sum of posterior probs of risk variants in cat g
          numerator <- sum(UiZij[cat_variants])

          # total.Ui: expected number of risk variant slots in cat g
          # = sum over genes of (count of cat g variants in gene) * EUi[gene]
          denominator <- 0
          for (i in seq_len(n_genes)) {
            idx <- gene_split[[i]]
            n_cat_in_gene <- sum(category[idx] == g)
            denominator <- denominator + n_cat_in_gene * EUi[i]
          }

          if (denominator > 0) {
            eta[g] <- numerator / denominator
          } else {
            eta[g] <- 0
          }
        }
        eta <- pmax(eta, .Machine$double.xmin)
        eta <- pmin(eta, 1 - .Machine$double.eps)
      }

    } else {
      # VARIANT-LEVEL EM

      # E-step
      # EZj[j] = eta_k * BF_j / (eta_k * BF_j + (1 - eta_k))
      log_EZj <- log(eta_prev[category]) + log_var_bf -
        vapply(seq_len(n_variants), function(j) {
          log_sum_exp(log(eta_prev[category[j]]) + log_var_bf[j],
                      log(1 - eta_prev[category[j]]))
        }, numeric(1))
      EZj <- exp(log_EZj)

      # M-step
      if (estimate_eta) {
        for (g in seq_len(n_categories)) {
          cat_variants <- which(category == g)
          if (length(cat_variants) == 0) {
            eta[g] <- 0
          } else {
            eta[g] <- sum(EZj[cat_variants]) / length(cat_variants)
          }
        }
        eta <- pmax(eta, .Machine$double.xmin)
        eta <- pmin(eta, 1 - .Machine$double.eps)
      }
    }

    # Store history
    eta_history[iter, ] <- eta
    if (type == "gene") {
      delta_history[iter] <- delta
    }

    # Check convergence
    diff <- sum(abs(eta - eta_prev))
    if (type == "gene" && estimate_delta) {
      diff <- diff + abs(delta - delta_prev)
    }

    final_iter <- iter
    if (diff < tol) {
      converged <- TRUE
      break
    }
  }

  # Trim history
  eta_history <- eta_history[seq_len(final_iter), , drop = FALSE]
  if (type == "gene") {
    delta_history <- delta_history[seq_len(final_iter)]
  }

  # Final E-step to get posteriors at convergence
  if (type == "gene") {
    gene_log_bf_final <- vapply(seq_len(n_genes), function(i) {
      idx <- gene_split[[i]]
      if (length(idx) == 0) return(0)
      sum(log_mixture_bf(log_var_bf[idx], eta[category[idx]]))
    }, numeric(1))

    log_EUi_final <- vapply(gene_log_bf_final, function(lb) {
      log_numer <- log(delta) + lb
      log_denom <- log_sum_exp(log_numer, log(1 - delta))
      log_numer - log_denom
    }, numeric(1))

    # Variant-level posteriors
    log_P_j_final <- log(eta[category]) + log_var_bf -
      vapply(seq_len(n_variants), function(j) {
        log_sum_exp(log(eta[category[j]]) + log_var_bf[j],
                    log(1 - eta[category[j]]))
      }, numeric(1))
    post_prob <- exp(log_EUi_final[gene_index] + log_P_j_final)

    result <- list(
      delta.est = delta,
      eta.est = eta,
      gene.log.bf = gene_log_bf_final,
      post.prob = post_prob,
      EUi = exp(log_EUi_final),
      eta.history = eta_history,
      delta.history = delta_history,
      converged = converged,
      n.iter = final_iter
    )
  } else {
    log_EZj_final <- log(eta[category]) + log_var_bf -
      vapply(seq_len(n_variants), function(j) {
        log_sum_exp(log(eta[category[j]]) + log_var_bf[j],
                    log(1 - eta[category[j]]))
      }, numeric(1))
    post_prob <- exp(log_EZj_final)

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
