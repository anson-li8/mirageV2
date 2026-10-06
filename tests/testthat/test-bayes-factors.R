test_that("calc_log_bf matches R integrate() for known cases", {
  # Original implementation using R's integrate()
  ref_log_bf <- function(var.case, var.contr, bar.gamma, sig, N1, N0) {
    integrand <- function(aa) {
      exp(dbinom(var.case, var.case + var.contr,
                 aa * N1 / (aa * N1 + N0), log = TRUE) +
            dgamma(aa, bar.gamma * sig, sig, log = TRUE))
    }
    log_marg0 <- dbinom(var.case, var.case + var.contr, N1 / (N1 + N0), log = TRUE)
    log_marg1 <- log(integrate(integrand, lower = 0, upper = 100,
                               stop.on.error = FALSE)$value)
    log_marg1 - log_marg0
  }

  # Test cases
  cases <- c(1L, 5L, 0L, 10L)
  contrs <- c(0L, 2L, 4L, 3L)
  gamma <- 3
  sigma <- 2
  N1 <- 4315
  N0 <- 4315

  v2_result <- calc_log_bf(cases, contrs, gamma, sigma, N1, N0)
  ref_result <- vapply(seq_along(cases), function(i) {
    ref_log_bf(cases[i], contrs[i], gamma, sigma, N1, N0)
  }, numeric(1))

  expect_equal(v2_result, ref_result, tolerance = 1e-4)
})

test_that("calc_log_bf returns 0 for T=0 (no observations)", {
  result <- calc_log_bf(0L, 0L, 3, 2, 1000, 1000)
  expect_equal(result, 0)
})

test_that("calc_log_bf handles category-specific gamma/sigma", {
  # gamma/sigma vectors of length > 1 should be recycled per variant
  result <- calc_log_bf(c(2L, 3L), c(1L, 0L), c(3, 5), c(2, 1), 1000, 1000)
  expect_length(result, 2)
  expect_true(all(is.finite(result)))
})

test_that("calc_log_bf produces positive BF for case-enriched variants", {
  # Variant seen 10x in cases, 0x in controls should have BF > 1 (log > 0)
  result <- calc_log_bf(10L, 0L, 3, 2, 5000, 5000)
  expect_gt(result, 0)
})

test_that("calc_log_bf produces negative BF for control-enriched variants", {
  # Variant seen 0x in cases, 10x in controls should have BF < 1 (log < 0)
  result <- calc_log_bf(0L, 10L, 3, 2, 5000, 5000)
  expect_lt(result, 0)
})
