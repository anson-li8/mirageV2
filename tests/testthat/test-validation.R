test_that("validate_mirage_data rejects non-data.frame input", {
  expect_error(validate_mirage_data(list(a = 1), type = "gene"),
               "must be a data frame")
})

test_that("validate_mirage_data rejects missing columns", {
  bad_df <- data.frame(ID = "v1", Gene = "A", No.case = 1)
  expect_error(validate_mirage_data(bad_df, type = "gene"),
               "must contain columns")
})

test_that("validate_mirage_data rejects negative counts", {
  bad_df <- data.frame(
    ID = "v1", Gene = "A", No.case = -1, No.contr = 0, category = 1
  )
  expect_error(validate_mirage_data(bad_df, type = "gene"),
               "non-negative")
})

test_that("validate_mirage_data rejects NA values", {
  bad_df <- data.frame(
    ID = "v1", Gene = "A", No.case = NA, No.contr = 0, category = 1
  )
  expect_error(validate_mirage_data(bad_df, type = "gene"),
               "missing values")
})

test_that("validate_mirage_params rejects invalid n1/n2", {
  expect_error(
    validate_mirage_params(-1, 100, 3, 2, 0.1, 0.1, TRUE, TRUE, NULL, 100, 1e-5, 2),
    "positive"
  )
})

test_that("validate_mirage_params rejects gamma/sigma length mismatch", {
  expect_error(
    validate_mirage_params(100, 100, c(3, 2, 1), 2, 0.1, 0.1, TRUE, TRUE, NULL, 100, 1e-5, 2),
    "length 1 or match"
  )
})

test_that("validate_mirage_params requires fixed.eta when estimate.eta=FALSE", {
  expect_error(
    validate_mirage_params(100, 100, 3, 2, 0.1, 0.1, TRUE, FALSE, NULL, 100, 1e-5, 2),
    "must be provided"
  )
})
