#' Calculate Log Bayes Factors for Variants
#'
#' Vectorized C++ helper for computing variant-level log Bayes Factors.
#' Uses adaptive Simpson's rule in log-space to prevent numerical over/underflow.
#'
#' @param var_case Integer vector of variant counts in cases.
#' @param var_contr Integer vector of variant counts in controls.
#' @param gamma Numeric vector of shape hyperparameters.
#' @param sigma Numeric vector of scale hyperparameters.
#' @param N1 Total sample size in cases.
#' @param N0 Total sample size in controls.
#' @return Numeric vector of log Bayes Factors.
#' @keywords internal
calc_log_bf <- function(var_case, var_contr, gamma, sigma, N1, N0) {
  n <- length(var_case)
  gamma_vec <- rep_len(gamma, n)
  sigma_vec <- rep_len(sigma, n)

  calc_log_bf_var_cpp(
    var_case = as.integer(var_case),
    var_contr = as.integer(var_contr),
    gamma_vec = as.numeric(gamma_vec),
    sigma_vec = as.numeric(sigma_vec),
    N1 = as.numeric(N1),
    N0 = as.numeric(N0)
  )
}
