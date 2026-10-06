#' Log-Sum-Exp for Two Values
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

#' Vectorized Log-Sum-Exp for Two Vectors
#'
#' Computes log(exp(a) + exp(b)) element-wise for vectors.
#'
#' @param a Numeric vector.
#' @param b Numeric vector (same length as a).
#' @return Numeric vector of log(exp(a) + exp(b)).
#' @keywords internal
log_sum_exp_vec2 <- function(a, b) {
  m <- pmax(a, b)
  # Handle -Inf cases: if both are -Inf, result is -Inf
  result <- m + log(exp(a - m) + exp(b - m))
  result[is.na(result)] <- -Inf
  result
}

#' Log-Sum-Exp for a Vector
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

#' Stable Log of (1 - eta) + eta * exp(log_bf) — Vectorized
#'
#' @param log_bf Numeric vector of log Bayes factors.
#' @param eta Numeric scalar or vector in \eqn{}\eqn{[0, 1]}.
#' @return Numeric vector of log((1-eta) + eta * exp(log_bf)).
#' @keywords internal
log_mixture_bf <- function(log_bf, eta) {
  if (length(eta) == 1) eta <- rep(eta, length(log_bf))
  result <- numeric(length(log_bf))
  valid <- eta > 0 & eta < 1
  if (any(valid)) {
    a <- log1p(-eta[valid])
    b <- log(eta[valid]) + log_bf[valid]
    m <- pmax(a, b)
    result[valid] <- m + log(exp(a - m) + exp(b - m))
  }
  result[eta <= 0] <- 0
  result[eta >= 1] <- log_bf[eta >= 1]
  result
}
