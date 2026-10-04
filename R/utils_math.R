#' Log-Sum-Exp for Two Values
#'
#' Calculates log(exp(a) + exp(b)) in a numerically stable way.
#'
#' @param a Numeric scalar.
#' @param b Numeric scalar.
#' @return log(exp(a) + exp(b))
#' @keywords internal
log_sum_exp <- function(a, b) {
  if (is.infinite(a) && a < 0) return(b)
  if (is.infinite(b) && b < 0) return(a)
  m <- max(a, b)
  m + log(exp(a - m) + exp(b - m))
}

#' Log-Sum-Exp for a Vector
#'
#' Calculates log(sum(exp(v))) in a numerically stable way.
#'
#' @param v Numeric vector.
#' @return log(sum(exp(v)))
#' @keywords internal
log_sum_exp_vec <- function(v) {
  v <- v[is.finite(v) | (is.infinite(v) & v > 0)]
  if (length(v) == 0) return(-Inf)
  m <- max(v)
  if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(v - m)))
}

#' Stable Log of (1 - eta) + eta * exp(log_bf)
#'
#' Calculates log((1-eta) + eta * BF) given log(BF), stably.
#' Equivalent to log_sum_exp(log(1-eta), log(eta) + log_bf).
#'
#' @param log_bf Numeric vector of log Bayes factors.
#' @param eta Numeric scalar in (0, 1).
#' @return Numeric vector of log((1-eta) + eta * exp(log_bf)).
#' @keywords internal
log_mixture_bf <- function(log_bf, eta) {
  if (eta <= 0) return(rep(0, length(log_bf)))
  if (eta >= 1) return(log_bf)
  log_one_minus_eta <- log1p(-eta)
  log_eta <- log(eta)
  vapply(log_bf, function(lb) {
    log_sum_exp(log_one_minus_eta, log_eta + lb)
  }, numeric(1))
}
