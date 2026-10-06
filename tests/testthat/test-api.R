test_that("mirage() runs on toy data and returns valid structure", {
  dat <- make_tiny_gene_data()
  res <- mirage(dat, n1 = 500, n2 = 500, verbose = FALSE)

  expect_s3_class(res, "mirage_result")
  expect_true("delta.est" %in% names(res))
  expect_true("eta.est" %in% names(res))
  expect_true("BF.PP.gene" %in% names(res))
  expect_true("BF.all" %in% names(res))
  expect_equal(nrow(res$BF.PP.gene), 3)
  expect_true(all(res$BF.PP.gene$post.prob >= 0))
  expect_true(all(res$BF.PP.gene$post.prob <= 1))
})

test_that("mirage_vs() runs on toy data and returns valid structure", {
  dat <- make_tiny_vs_data()
  res <- mirage_vs(dat, n1 = 500, n2 = 500, verbose = FALSE)

  expect_s3_class(res, "mirage_vs_result")
  expect_true("eta.est" %in% names(res))
  expect_true("full.info" %in% names(res))
  expect_true("post.prob" %in% names(res))
  expect_equal(nrow(res$post.prob), 6)
})

test_that("mirage() accepts group.index column name", {
  dat <- make_tiny_gene_data()
  names(dat)[names(dat) == "category"] <- "group.index"
  res <- mirage(dat, n1 = 500, n2 = 500, verbose = FALSE)
  expect_s3_class(res, "mirage_result")
})

test_that("mirage() works with estimate.eta=FALSE", {
  dat <- make_tiny_gene_data()
  res <- mirage(dat, n1 = 500, n2 = 500, estimate.eta = FALSE,
                fixed_eta = c(0.2, 0.8), verbose = FALSE)
  expect_s3_class(res, "mirage_result")
})

test_that("mirage() handles single-gene input", {
  dat <- data.frame(
    ID = c("v1", "v2"),
    Gene = c("A", "A"),
    No.case = c(5L, 2L),
    No.contr = c(0L, 1L),
    category = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  res <- mirage(dat, n1 = 1000, n2 = 1000, verbose = FALSE)
  expect_s3_class(res, "mirage_result")
  expect_equal(nrow(res$BF.PP.gene), 1)
})

test_that("mirage() handles single-variant-per-gene input", {
  dat <- data.frame(
    ID = c("v1", "v2"),
    Gene = c("A", "B"),
    No.case = c(3L, 1L),
    No.contr = c(0L, 2L),
    category = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  res <- mirage(dat, n1 = 500, n2 = 500, verbose = FALSE)
  expect_s3_class(res, "mirage_result")
  expect_equal(nrow(res$BF.PP.gene), 2)
})

test_that("mirage() works on full mirage_toy dataset", {
  skip_if_not(file.exists(testthat::test_path("..", "..", "inst", "extdata", "mirage_toy.csv")),
              message = "mirage_toy.csv not found")
  dat <- load_toy_gene()
  res <- mirage(dat, n1 = 4315, n2 = 4315, verbose = FALSE)
  expect_s3_class(res, "mirage_result")
  expect_gt(nrow(res$BF.PP.gene), 100)
  expect_true(all(is.finite(res$BF.PP.gene$BF)))
})
