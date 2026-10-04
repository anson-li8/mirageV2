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
#' Both log_bf and eta can be vectors of the same length,
#' or eta can be a scalar (recycled).
#'
#' @param log_bf Numeric vector of log Bayes factors.
#' @param eta Numeric scalar or vector in \eqn{}\eqn{[0, 1]}.
#' @return Numeric vector of log((1-eta) + eta * exp(log_bf)).
#' @keywords internal
log_mixture_bf <- function(log_bf, eta) {
  if (length(eta) == 1) eta <- rep(eta, length(log_bf))
  vapply(seq_along(log_bf), function(j) {
    if (eta[j] <= 0) return(0)
    if (eta[j] >= 1) return(log_bf[j])
    log_sum_exp(log1p(-eta[j]), log(eta[j]) + log_bf[j])
  }, numeric(1))
}
