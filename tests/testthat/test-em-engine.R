test_that("log_sum_exp is numerically stable", {
  # Large values that would overflow unstable exp()
  expect_equal(log_sum_exp(1000, 1000), 1000 + log(2), tolerance = 1e-10)
  # -Inf handling
  expect_equal(log_sum_exp(-Inf, 5), 5)
  expect_equal(log_sum_exp(5, -Inf), 5)
})

test_that("log_sum_exp_vec handles edge cases", {
  expect_equal(log_sum_exp_vec(c(-Inf, -Inf)), -Inf)
  expect_equal(log_sum_exp_vec(c(0, 0)), log(2), tolerance = 1e-10)
  expect_equal(log_sum_exp_vec(numeric(0)), -Inf)
})

test_that("log_mixture_bf handles vector eta correctly", {
  log_bf <- c(0.5, -0.3, 1.2)
  eta <- c(0.1, 0.5, 0.9)
  result <- log_mixture_bf(log_bf, eta)
  expect_length(result, 3)
  # Each result should be log((1-eta) + eta*exp(log_bf))
  expected <- log((1 - eta) + eta * exp(log_bf))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("log_mixture_bf edge cases: eta=0 returns 0, eta=1 returns log_bf", {
  log_bf <- c(1.5, -2.0)
  expect_equal(log_mixture_bf(log_bf, 0), c(0, 0))
  expect_equal(log_mixture_bf(log_bf, 1), log_bf)
})

test_that("em_mirage gene-level converges with valid parameters", {
  dat <- make_tiny_gene_data()
  dat_validated <- validate_mirage_data(dat, type = "gene")
  log_bf <- calc_log_bf(dat_validated$No.case, dat_validated$No.contr,
                        3, 2, 500, 500)
  gene_index <- match(dat_validated$Gene, unique(dat_validated$Gene))

  em <- em_mirage(
    type = "gene", log_var_bf = log_bf, category = dat_validated$category,
    gene_index = gene_index, n_genes = 3L, n_categories = 2L,
    max_iter = 5000L, tol = 1e-5, verbose = FALSE
  )

  expect_true(em$converged)
  expect_gt(em$delta.est, 0)
  expect_lte(em$delta.est, 1)
  expect_true(all(em$eta.est >= 0))
  expect_true(all(em$eta.est <= 1))
  expect_length(em$post.prob, nrow(dat))
  expect_true(all(em$post.prob >= 0 & em$post.prob <= 1))
})

test_that("em_mirage vs-level converges with valid parameters", {
  dat <- make_tiny_vs_data()
  dat_validated <- validate_mirage_data(dat, type = "vs")
  log_bf <- calc_log_bf(dat_validated$No.case, dat_validated$No.contr,
                        3, 2, 500, 500)

  em <- em_mirage(
    type = "vs", log_var_bf = log_bf, category = dat_validated$category,
    n_categories = 2L, max_iter = 5000L, tol = 1e-5, verbose = FALSE
  )

  expect_true(em$converged)
  expect_length(em$eta.est, 2)
  expect_true(all(em$eta.est >= 0 & em$eta.est <= 1))
})

test_that("em_mirage respects fixed eta when estimate_eta=FALSE", {
  dat <- make_tiny_vs_data()
  dat_validated <- validate_mirage_data(dat, type = "vs")
  log_bf <- calc_log_bf(dat_validated$No.case, dat_validated$No.contr,
                        3, 2, 500, 500)
  fixed_eta <- c(0.3, 0.7)

  em <- em_mirage(
    type = "vs", log_var_bf = log_bf, category = dat_validated$category,
    n_categories = 2L, estimate_eta = FALSE, fixed_eta = fixed_eta,
    max_iter = 100L, verbose = FALSE
  )

  expect_equal(em$eta.est, fixed_eta)
})
