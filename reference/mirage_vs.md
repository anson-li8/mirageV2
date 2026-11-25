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

  Prior mean of the relative-risk parameter. Scalar or length-K vector.

- sigma:

  Prior rate of the relative-risk parameter. Scalar or length-K vector.
  The Gamma prior has shape gamma \* sigma and scale 1 / sigma.

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

## Examples

``` r
toy <- data.frame(
  ID = c("v1", "v2", "v3", "v4"),
  No.case = c(3L, 1L, 0L, 2L),
  No.contr = c(1L, 0L, 2L, 0L),
  category = c(1L, 2L, 1L, 2L)
)
mirage_vs(toy, n1 = 500, n2 = 500, verbose = FALSE)
#> 
#> ── MIRAGE Variant-Set Analysis Result ──────────────────────────────────────────
#> Number of variants analyzed: 4
#> Number of variant categories: 2
#> 
#> ── Category-specific risk proportions (eta) ──
#> 
#>  category      eta.est eta.pvalue
#>         1 0.0001999959  1.0000000
#>         2 1.0000000000  0.1323021
```
