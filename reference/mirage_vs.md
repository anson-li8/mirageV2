# mirage_vs: Variant-set level analysis

Variant-level rare variant association test using a Bayesian mixture
model.

## Usage

``` r
mirage_vs(
  data,
  n1,
  n2,
  gamma = 3,
  sigma = 2,
  eta.init = 0.1,
  estimate.eta = TRUE,
  fixed_eta = NULL,
  max.iter = 10000L,
  tol = 1e-05,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with columns: ID, No.case, No.contr, category. A 3-column
  input (without ID) is accepted with a warning.

- n1:

  Sample size in cases.

- n2:

  Sample size in controls.

- gamma:

  Hyper prior shape parameter(s). Scalar or length-K vector.

- sigma:

  Hyper prior scale parameter(s). Scalar or length-K vector.

- eta.init:

  Initial value for proportion of risk variants.

- estimate.eta:

  Logical. Whether to estimate eta.

- fixed_eta:

  Numeric vector of fixed eta values (required when estimate.eta =
  FALSE). Must have length K.

- max.iter:

  Maximum EM iterations.

- tol:

  Convergence tolerance.

- verbose:

  Logical. Print progress messages.

## Value

An object of class `mirage_vs_result`.
